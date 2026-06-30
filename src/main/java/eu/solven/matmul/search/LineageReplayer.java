package eu.solven.matmul.search;

import eu.solven.matmul.recombination.BlockSplitSearch;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.concurrent.ConcurrentHashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.Compose;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.recombination.Recombination;
import eu.solven.matmul.catalog.SchemeIO;
import lombok.extern.slf4j.Slf4j;
import tools.jackson.databind.JsonNode;

/**
 * Walks a {@link Lineage.Node} tree and rebuilds the corresponding
 * {@link NonCubicBilinearAlgorithm} by replaying every constructor step
 * via the existing {@link Compose} / {@link Recombination} primitives.
 *
 * <p>This is the inverse of {@code RecursiveMaterialiser}: the materialiser
 * picks a strategy and emits a lineage; the replayer takes a lineage and
 * reconstructs the algorithm without re-searching. It exists primarily to
 * support stub-format scheme files (see {@code MigrateToStubs}) — large
 * derived-recursive schemes are stored as lineage-only stubs, and any
 * consumer that needs the actual tensor materialises it on demand here.</p>
 *
 * <p>Supported ops: {@link Lineage.Atom}, {@link Lineage.KronProduct},
 * {@link Lineage.KronChain}, {@link Lineage.ConcatCols},
 * {@link Lineage.ConcatRows}, {@link Lineage.SumInner}, {@link Lineage.RecombinationN},
 * {@link Lineage.Transpose}. Unsupported (rare in practice):
 * {@link Lineage.RecombinationWithPairN}, {@link Lineage.AugmentSquareDiscard},
 * {@link Lineage.Dce} — these throw {@link UnsupportedOperationException}
 * with a clear message so callers can tally them separately.</p>
 *
 * <p>Caching: two layers, deliberately decoupled to keep memory bounded:</p>
 * <ul>
 *   <li><strong>Per-call memo</strong> ({@link IdentityHashMap} created
 *       fresh on each {@link #replay} entry): shared subtrees within
 *       <em>one</em> lineage (e.g. from {@code @ref} dedup) are computed
 *       once. The memo dies with the call, so retention is bounded by the
 *       size of the largest single tree.</li>
 *   <li><strong>Cross-call shape cache</strong> ({@link #STUB_CACHE},
 *       bounded LRU): when a leaf resolves to a stub file, the
 *       materialised algorithm is cached by canonical file path so
 *       independent callers (e.g. parallel verify across many stubs that
 *       all reference the same ⟨16,16,16⟩) hit it. Bounded so a sweep
 *       across 60+ heavy stubs doesn't accumulate gigabytes.</li>
 * </ul>
 *
 * <p>The replayer instance is otherwise stateless and thread-safe; share
 * one across all worker threads instead of {@code ThreadLocal} per
 * thread.</p>
 */
@Slf4j
public final class LineageReplayer {

	// Accept the shape ref in either marker position: a bare "AxBxC", the suffix
	// forms "AxBxC-direct"/"AxBxC-naive", OR the prefix forms
	// "naive-AxBxC"/"direct-AxBxC" emitted by RecursiveMaterialiser.trivialOneAxis.
	private static final Pattern SHAPE_REF =
			Pattern.compile("(?:(naive|direct)-)?(\\d+)x(\\d+)x(\\d+)(?:-(direct|naive))?");
	/** Loose shape extractors for the resolveLeaf fallback: {@code NxMxP} anywhere
	 *  (source-prefixed / canonical-key refs) and {@code <N,M,P>} (named refs). */
	private static final Pattern SHAPE_X = Pattern.compile("(\\d+)x(\\d+)x(\\d+)");
	private static final Pattern SHAPE_ANGLE =
			Pattern.compile("<\\s*(\\d+)\\s*,\\s*(\\d+)\\s*,\\s*(\\d+)\\s*>");

	private static int[] extractShape(String ref) {
		Matcher mx = SHAPE_X.matcher(ref);
		if (mx.find()) {
			return new int[] {
					Integer.parseInt(mx.group(1)), Integer.parseInt(mx.group(2)), Integer.parseInt(mx.group(3)) };
		}
		Matcher ma = SHAPE_ANGLE.matcher(ref);
		if (ma.find()) {
			return new int[] {
					Integer.parseInt(ma.group(1)), Integer.parseInt(ma.group(2)), Integer.parseInt(ma.group(3)) };
		}
		return null;
	}
	/** Pan/Islam trilinear-aggregation cubic ⟨n,n,n⟩ (DIS09 Appendix Lemma 4),
	 *  emitted by {@code MaterializeDIS09TA} as the lineage of the dis09 cubes. */
	private static final Pattern DIS09_LEMMA4 = Pattern.compile("DIS09Lemma4\\(n=(\\d+)\\)");
	/** KGP 2026 LITA cubic ⟨N,N,N⟩ parametric leaf — replays via LitaTaConstruction.build. */
	private static final Pattern TA_LITA = Pattern.compile("TA_lita\\(n=(\\d+)\\)");

	/**
	 * Cross-call cache of canonical-path → materialised algorithm for stub
	 * replays, bounded by total <strong>weight</strong> (dense {@code double}-cell
	 * count ≈ memory footprint), NOT entry count. A flat entry cap is wrong here
	 * because instances differ in size by orders of magnitude: a small ⟨2,2,2⟩ alg
	 * is a few hundred cells while a ⟨30,32,32⟩ alg is ~40M cells (~330&nbsp;MB), so
	 * 64 of the latter is ~21&nbsp;GB → OOM. The weigher counts
	 * {@code r·(nm+mp+np)} cells; the budget is ~¼ of the heap, so it auto-scales
	 * with {@code -Xmx} and bounds peak retention regardless of shape mix.
	 * (The real fix is a sparse backing for {@link NonCubicBilinearAlgorithm} —
	 * these tensors are ~98% zero, ternary — which would cut this ~50×.)
	 */
	private static final long STUB_CACHE_CELL_BUDGET = Math.max(8_000_000L,
			Runtime.getRuntime().maxMemory() / 8 / 4);  // doubles are 8 bytes; cap at ¼ heap
	private static final com.google.common.cache.Cache<Path, NonCubicBilinearAlgorithm> STUB_CACHE =
			com.google.common.cache.CacheBuilder.newBuilder()
					.maximumWeight(STUB_CACHE_CELL_BUDGET)
					.weigher((Path k, NonCubicBilinearAlgorithm a) -> (int) Math.min(Integer.MAX_VALUE,
							(long) a.r * ((long) a.n * a.m + (long) a.m * a.p + (long) a.n * a.p)))
					.concurrencyLevel(Math.max(1, Runtime.getRuntime().availableProcessors()))
					.build();

	/**
	 * Per-thread set of stub files currently being replayed on this descent, for
	 * cycle detection in {@link #replayFromFile}. A derived stub whose lineage
	 * transitively references itself (a corrupt projection/recombination chain)
	 * would otherwise recurse forever and blow the stack with an UNCATCHABLE
	 * {@link StackOverflowError}. {@code LinkedHashSet} preserves the descent order
	 * so the thrown error names the actual cycle path.
	 */
	private static final ThreadLocal<java.util.LinkedHashSet<String>> ACTIVE_REPLAYS =
			ThreadLocal.withInitial(java.util.LinkedHashSet::new);

	private final FieldAwareLookup lookup;
	private final Map<String, NonCubicBilinearAlgorithm> namedBases;

	public LineageReplayer(FieldAwareLookup lookup,
			Map<String, NonCubicBilinearAlgorithm> namedBases) {
		this.lookup = lookup;
		this.namedBases = namedBases;
	}

	/**
	 * Convenience: build a replayer with the standard named-base pool
	 * (Strassen ⟨2,2,2⟩=7, Laderman ⟨3,3,3⟩=23, …) loaded from disk.
	 */
	public static LineageReplayer withDefaultPool(FieldAwareLookup lookup) {
		Map<String, NonCubicBilinearAlgorithm> pool = new ConcurrentHashMap<>();
		for (BlockSplitSearch.NamedBase nb : BlockSplitSearch.defaultPool()) {
			pool.put(nb.label(), nb.base());
		}
		return new LineageReplayer(lookup, pool);
	}

	/**
	 * Replay a lineage tree into a fully-materialised algorithm. The
	 * call-local {@link IdentityHashMap} memo dies on return, so memory
	 * retention is bounded by the largest single tree.
	 */
	public NonCubicBilinearAlgorithm replay(Lineage.Node node) {
		return replayInternal(node, new IdentityHashMap<>());
	}

	private NonCubicBilinearAlgorithm replayInternal(Lineage.Node node,
			IdentityHashMap<Lineage.Node, NonCubicBilinearAlgorithm> memo) {
		NonCubicBilinearAlgorithm cached = memo.get(node);
		if (cached != null) return cached;
		NonCubicBilinearAlgorithm out = switch (node) {
			case Lineage.Atom l -> resolveLeaf(l.ref());
			case Lineage.KronProduct k -> Compose.kroneckerGeneral(
					replayInternal(k.outer(), memo), replayInternal(k.inner(), memo));
			case Lineage.KronChain c -> kronChain(c.factors(), memo);
			case Lineage.ConcatCols r -> concatCols(
					replayInternal(r.left(), memo), replayInternal(r.right(), memo));
			case Lineage.ConcatRows b -> concatRows(
					replayInternal(b.top(), memo), replayInternal(b.bottom(), memo));
			case Lineage.SumInner s -> concatInner(
					replayInternal(s.left(), memo), replayInternal(s.right(), memo));
			case Lineage.RecombinationN r -> applyRecomb(r, memo);
			case Lineage.Transpose t -> applyTranspose(t, memo);
			case Lineage.OrientAs o -> {
				NonCubicBilinearAlgorithm oc = replayInternal(o.child(), memo);
				// Explicit axis-map → deterministic orientByPerm (no shape search, dim-repeat-safe);
				// legacy null → orientAs inference (first shape match).
				yield o.axisMap() != null
						? oc.orientByPerm(o.axisMap())
						: oc.orientAs(o.n(), o.m(), o.p())
								.orElseThrow(() -> new IllegalStateException("LineageReplayer: OrientAs to ⟨"
										+ o.n() + "," + o.m() + "," + o.p() + "⟩ failed"));
			}
			case Lineage.AxisFlip af -> applyAxisFlip(af, memo);
			case Lineage.AxisPermute ap -> applyAxisPermute(ap, memo);
			case Lineage.RecombinationWithPairN ignored -> throw new UnsupportedOperationException(
					"RecombinationWithPairN replay not implemented (no live usage in catalog as of 2026-05)");
			case Lineage.AugmentSquareDiscard ignored -> throw new UnsupportedOperationException(
					"AugmentSquareDiscard replay not implemented");
			case Lineage.Dce ignored -> throw new UnsupportedOperationException(
					"Dce replay not implemented");
			case Lineage.DisjointSum ignored -> throw new UnsupportedOperationException(
					"DisjointSum materialisation not yet implemented (rank-prediction-only node)");
			// Reconstruct under the CHEAPEST bud ordering (min over ALL_ORDERINGS), matching
			// SerendipitousSearch.bestFor — NOT the DEFAULT_ORDER of productViaBuds, which
			// could build at a higher rank than predicted and get the win discarded.
			case Lineage.SerendipitousProduct sp -> eu.solven.matmul.catalog.SerendipitousBudProduct
					.productViaBudsBest(replayInternal(sp.base(), memo), lookup, sp.n2(), sp.m2(), sp.p2());
			case Lineage.Project pr -> eu.solven.matmul.catalog.Compose.project(
					replayInternal(pr.child(), memo), pr.keepN(), pr.keepM(), pr.keepP());
			// Rectangular Pan-TA peel: re-run the deterministic constructor on the
			// replayed cube + corner (the TA cross-fusion is fixed by (n,s)).
			case Lineage.PeeledViaTa t -> eu.solven.matmul.papers.pan1978.RectangularTrilinearAggregation
					.buildPeeledViaTa(t.n(), t.s(), replayInternal(t.cube(), memo), replayInternal(t.corner(), memo));
			// Generic TA-fused recombination: replay the naïve-grid base + the NC-pinned
			// unpaired leaves, then re-run constructWithTaFusion (which deterministically
			// re-derives the fused cyclic pairs). Leaves resolved by sorted-shape key.
			case Lineage.RecombinationTaN r -> replayRecombinationTa(r, memo);
		};
		memo.put(node, out);
		return out;
	}

	/** Replay a {@link Lineage.RecombinationTaN}: rebuild via {@link
	 *  eu.solven.matmul.recombination.Recombination#constructWithTaFusion} on the replayed
	 *  naïve-grid base, resolving each unpaired leaf from the stored NC-pinned nodes by
	 *  sorted-shape key (re-oriented), and re-deriving the fused pairs deterministically. */
	private NonCubicBilinearAlgorithm replayRecombinationTa(Lineage.RecombinationTaN r,
			IdentityHashMap<Lineage.Node, NonCubicBilinearAlgorithm> memo) {
		NonCubicBilinearAlgorithm baseAlg = replayInternal(r.base(), memo);
		Map<String, NonCubicBilinearAlgorithm> leafByShape = new java.util.HashMap<>();
		for (Lineage.Node lf : r.leaves()) {
			NonCubicBilinearAlgorithm la = replayInternal(lf, memo);
			int[] s = { la.n, la.m, la.p };
			java.util.Arrays.sort(s);
			leafByShape.put(s[0] + "x" + s[1] + "x" + s[2], la);
		}
		eu.solven.matmul.recombination.Recombination.SubResolver resolveSub = (sz) -> {
			int[] k = { sz[0], sz[1], sz[2] };
			java.util.Arrays.sort(k);
			NonCubicBilinearAlgorithm la = leafByShape.get(k[0] + "x" + k[1] + "x" + k[2]);
			if (la == null) {
				return null;
			}
			return la.orientAs(sz[0], sz[1], sz[2]).orElse(la);
		};
		eu.solven.matmul.recombination.Recombination.SotaResolver sota = (a, b, c) -> lookup.findRank(a, b, c);
		return eu.solven.matmul.recombination.Recombination.constructWithTaFusion(
				baseAlg, resolveSub, sota, r.allocA(), r.allocB(), r.allocC()).alg();
	}

	private NonCubicBilinearAlgorithm kronChain(List<Lineage.Node> factors,
			IdentityHashMap<Lineage.Node, NonCubicBilinearAlgorithm> memo) {
		NonCubicBilinearAlgorithm acc = replayInternal(factors.get(0), memo);
		for (int i = 1; i < factors.size(); i++) {
			acc = Compose.kroneckerGeneral(acc, replayInternal(factors.get(i), memo));
		}
		return acc;
	}

	// ── Robust concat replay ───────────────────────────────────────────────
	// Legacy stub lineages (written by an older RecursiveMaterialiser that named
	// direct leaves from the SORTED shape — canon()) encode concat operands in a
	// PERMUTED axis frame: e.g. ConcatCols(2x12x15, 12x15x15) for a ⟨12,15,17⟩
	// result, where the left operand is really ⟨12,15,2⟩ written rotated. The
	// recorded op (Cols/Rows/Inner) is likewise frame-relative and unreliable.
	//
	// Both Compose.concat* and orientAs are correctness-PRESERVING (axis
	// permutation of a matmul scheme is a matmul scheme of the permuted shape),
	// so re-orienting operands to satisfy a concat precondition can only ever
	// yield a valid matmul — never a silently-wrong one. We therefore: try the
	// exact composition first (correct-frame, current files); on a precondition
	// failure, search orientations of both operands for one that concatenates,
	// inferring the axis from the unique differing dimension. The final result is
	// oriented to the file's declared shape and verified by the caller, so a
	// wrong frame surfaces as an honest error, not corruption.

	private NonCubicBilinearAlgorithm concatCols(NonCubicBilinearAlgorithm a, NonCubicBilinearAlgorithm b) {
		try { return Compose.concatRight(a, b); } catch (IllegalArgumentException e) { /* permuted frame */ }
		return robustConcat(a, b, Axis.P);
	}

	private NonCubicBilinearAlgorithm concatRows(NonCubicBilinearAlgorithm a, NonCubicBilinearAlgorithm b) {
		try { return Compose.concatBelow(a, b); } catch (IllegalArgumentException e) { /* permuted frame */ }
		return robustConcat(a, b, Axis.N);
	}

	private NonCubicBilinearAlgorithm concatInner(NonCubicBilinearAlgorithm a, NonCubicBilinearAlgorithm b) {
		try { return Compose.concatInner(a, b); } catch (IllegalArgumentException e) { /* permuted frame */ }
		return robustConcat(a, b, Axis.M);
	}

	/** Which axis a concat adds along: P = ConcatCols, N = ConcatRows, M = SumInner. */
	private enum Axis { N, M, P }

	/**
	 * Orient {@code a} and {@code b} so they share the two non-{@code add} axes
	 * (concatenating along {@code add}), then compose. Tries every orientation
	 * pair; picks the first that satisfies the precondition. Throws if no
	 * orientation pair works (genuinely incompatible operands).
	 */
	private NonCubicBilinearAlgorithm robustConcat(NonCubicBilinearAlgorithm a,
			NonCubicBilinearAlgorithm b, Axis add) {
		for (NonCubicBilinearAlgorithm oa : orientations(a)) {
			for (NonCubicBilinearAlgorithm ob : orientations(b)) {
				switch (add) {
					case P -> { if (oa.n == ob.n && oa.m == ob.m) return Compose.concatRight(oa, ob); }
					case N -> { if (oa.m == ob.m && oa.p == ob.p) return Compose.concatBelow(oa, ob); }
					case M -> { if (oa.n == ob.n && oa.p == ob.p) return Compose.concatInner(oa, ob); }
				}
			}
		}
		throw new IllegalArgumentException("robustConcat(" + add + "): no compatible orientation for ⟨"
				+ a.n + "," + a.m + "," + a.p + "⟩ and ⟨" + b.n + "," + b.m + "," + b.p + "⟩");
	}

	/** All distinct axis-orientations of {@code a} (≤6), as concrete algorithms. */
	private static List<NonCubicBilinearAlgorithm> orientations(NonCubicBilinearAlgorithm a) {
		int[][] perms = { {a.n, a.m, a.p}, {a.n, a.p, a.m}, {a.m, a.n, a.p},
				{a.m, a.p, a.n}, {a.p, a.n, a.m}, {a.p, a.m, a.n} };
		java.util.LinkedHashMap<String, NonCubicBilinearAlgorithm> byShape = new java.util.LinkedHashMap<>();
		for (int[] s : perms) {
			String key = s[0] + "x" + s[1] + "x" + s[2];
			if (byShape.containsKey(key)) continue;
			a.orientAs(s[0], s[1], s[2]).ifPresent(o -> byShape.put(key, o));
		}
		return new java.util.ArrayList<>(byShape.values());
	}

	private NonCubicBilinearAlgorithm resolveLeaf(String ref) {
		// Hash-ref: "{n}x{m}x{p}@{contentHash}" — resolve the EXACT scheme by
		// content hash (precise; avoids the shape-ref picking a bud-poor sibling
		// of the same shape/rank/adds). Falls back to the shape part if the hash
		// no longer resolves (catalog changed).
		int at = ref.indexOf('@');
		if (at > 0) {
			String shapePart = ref.substring(0, at);
			String hash = ref.substring(at + 1);
			Matcher hm = SHAPE_REF.matcher(shapePart);
			if (hm.matches()) {
				int hn = Integer.parseInt(hm.group(2));
				int hmm = Integer.parseInt(hm.group(3));
				int hp = Integer.parseInt(hm.group(4));
				Optional<FieldAwareLookup.WithSource> byHash = lookup.findByHash(hn, hmm, hp, hash);
				if (byHash.isPresent()) {
					return byHash.get().alg();
				}
				// findByHash only hashes materialised candidates — a ref pinned to a
				// STUB's hash never matched it (audit 2026-06-10 blind spot). Check the
				// stubs' stamped "hash" field and replay the matching stub.
				for (java.nio.file.Path cand : lookup.findFiles(hn, hmm, hp)) {
					try {
						tools.jackson.databind.JsonNode candRoot = SchemeIO.parseJson(cand.toFile());
						String stamped = SchemeIO.readHash(candRoot);
						if (stamped != null && stamped.startsWith(hash)
								&& SchemeIO.readLineage(candRoot).isPresent()
								&& !candRoot.has("u") && !candRoot.has("u_sparse")) {
							return replayFromFile(cand.toFile());
						}
					} catch (Exception e) {
						// unreadable candidate — keep scanning
					}
				}
				// Hash gone → fall through to shape-ref. WARN: the replay will use
				// DIFFERENT content than the pin recorded (same shape, possibly same
				// rank) — the phantom-replay mechanism behind task #91's ⟨17,22,29⟩
				// 6129→6138. Callers needing bit-exactness must treat this as failure;
				// RepinDanglingLineageRefs is the repair pass for such pins.
				log.warn("resolveLeaf: pinned ref {} no longer resolves — falling back to "
						+ "shape-best at ⟨{},{},{}⟩ (content will differ from the pin)",
						ref, hn, hmm, hp);
				ref = shapePart;
			}
		}

		// Named base (Strassen<…>=…, Laderman<…>=…, etc.)
		NonCubicBilinearAlgorithm named = namedBases.get(ref);
		if (named != null) return named;

		// Parametric constructor ref, e.g. "DIS09Lemma4(n=26)": reconstruct via
		// the formula-driven constructor it was materialised from, so a stub whose
		// lineage is just this ref replays faithfully (the dis09 even cubes that
		// projection wants as parents are stored this way).
		NonCubicBilinearAlgorithm parametric = resolveParametric(ref);
		if (parametric != null) return parametric;

		// Shape ref: "AxBxC", "AxBxC-direct/-naive", or prefix "naive-/direct-AxBxC".
		Matcher m = SHAPE_REF.matcher(ref);
		if (!m.matches()) {
			// Robust fallback: pull a shape out of a source-prefixed / canonical-key /
			// named ref (e.g. "alphatensor_Z-2x3x3_m15_a58", "perminov_Z-2x5x15_…",
			// "Strassen<2,2,2>=7") and resolve it by shape. This loses the exact
			// scheme (resolves the catalog-best at that shape), which is fine for
			// projection parents and rank-only leaves. Only throw if no shape is
			// recoverable at all.
			int[] sh = extractShape(ref);
			if (sh != null) {
				if (sh[0] == 1 || sh[1] == 1 || sh[2] == 1) {
					return naiveScheme(sh[0], sh[1], sh[2]);
				}
				return resolveShape(sh[0], sh[1], sh[2]);
			}
			throw new IllegalStateException("LineageReplayer: cannot resolve leaf ref '" + ref + "'");
		}
		String prefixMarker = m.group(1);   // "naive"/"direct" or null
		int n = Integer.parseInt(m.group(2));
		int mm = Integer.parseInt(m.group(3));
		int p = Integer.parseInt(m.group(4));
		String suffixMarker = m.group(5);   // "direct"/"naive" or null
		boolean naive = prefixMarker != null || suffixMarker != null;
		// A naive/trivial leaf — or any width-1 axis — is just the elementary
		// n·m·p-product scheme; build it directly (it has no catalog file).
		// NOTE (2026-06-15): a "-direct" leaf (best-at-shape, unpinned) is a CITED
		// BOUND, NOT a buildable explicit scheme — building an explicit scheme over it
		// is invalid (the exact scheme matters for downstream DCE/projection; "take the
		// best is bad"). The principled fix is a separate change (fail on "-direct" in
		// the build path; resolve best-bound only in a cited-bound path) — see the
		// VerifyScheme DAG design. Left as-is here pending that decision.
		if (naive || n == 1 || mm == 1 || p == 1) {
			return naiveScheme(n, mm, p);
		}
		return resolveShape(n, mm, p);
	}

	/**
	 * The elementary "naive" bilinear scheme for {@code ⟨n,m,p⟩}: one product per
	 * (i,j,l) triple, rank {@code n·m·p}. Mirrors
	 * {@code RecursiveMaterialiser.trivialOneAxis} so a {@code naive-NxMxP} leaf
	 * replays identically to how it was emitted.
	 */
	private static NonCubicBilinearAlgorithm naiveScheme(int n, int mm, int p) {
		// Delegate to the single sparse factory so construct and replay agree exactly
		// and neither allocates a dense double[][] for width-1 / trivial blocks.
		return NonCubicBilinearAlgorithm.naive(n, mm, p);
	}

	/**
	 * Reconstruct a parametric constructor leaf ref by invoking the formula it
	 * was materialised from. Returns {@code null} when {@code ref} is not a
	 * recognised parametric form (the caller then tries the shape-ref path).
	 *
	 * <p>Currently handles {@code DIS09Lemma4(n=N)} → {@link
	 * PanTrilinearAggregation#build(int)} and {@code TA_lita(n=N)} → {@link
	 * eu.solven.matmul.papers.khoruzhii2026.LitaTaConstruction#build(int)} (both
	 * non-commutative cubic ⟨N,N,N⟩). Add further families here as their stubs
	 * appear (e.g. RosowskiTheorem2(…), commutative).</p>
	 */
	private NonCubicBilinearAlgorithm resolveParametric(String ref) {
		Matcher d = DIS09_LEMMA4.matcher(ref);
		if (d.matches()) {
			return eu.solven.matmul.papers.dis2009.PanTrilinearAggregation
					.build(Integer.parseInt(d.group(1)));
		}
		Matcher lita = TA_LITA.matcher(ref);
		if (lita.matches()) {
			return eu.solven.matmul.papers.khoruzhii2026.LitaTaConstruction
					.build(Integer.parseInt(lita.group(1)));
		}
		return null;
	}

	/**
	 * Resolve a shape to an algorithm: first via the disk lookup, falling
	 * back to recursive stub-replay if the catalog file is a stub. Throws
	 * if neither path yields an algorithm.
	 *
	 * <p>Cross-call cache: stub replays are cached in {@link #STUB_CACHE}
	 * by canonical file path so a parallel sweep over many shapes that
	 * all reference the same ⟨16,16,16⟩ leaf doesn't re-materialise it
	 * per call. The cache is a bounded LRU.</p>
	 */
	private NonCubicBilinearAlgorithm resolveShape(int n, int mm, int p) {
		Optional<FieldAwareLookup.WithSource> hit = lookup.findWithSource(n, mm, p);
		Optional<Path> bestPath = lookup.findFile(n, mm, p);
		// findWithSource returns the best MATERIALISED (non-stub) entry; findFile
		// returns the best file of ANY kind (rank-ascending, stubs included). When
		// the global best is a stub whose rank beats every materialised sibling —
		// e.g. a dis09 DIS09Lemma4(n) cube, now stored lineage-only — replaying it is
		// correct; returning the worse full file is the bug #68 used to hit. A path
		// inequality means findFile's best is one findWithSource skipped, i.e. a stub.
		boolean betterStubExists = bestPath.isPresent()
				&& (hit.isEmpty() || !bestPath.get().equals(hit.get().path()));
		if (betterStubExists) {
			try {
				return resolveStubFile(bestPath.get(), n, mm, p);
			} catch (RuntimeException e) {
				log.debug("resolveShape ⟨{},{},{}⟩: best stub {} not replayable ({}); "
						+ "falling back to best materialised", n, mm, p,
						bestPath.get().getFileName(), e.toString());
				// fall through to the materialised hit below
			}
		}
		if (hit.isPresent()) return hit.get().alg();
		// No materialised entry at all — the only candidate(s) are stub(s). Replay
		// the best (parseCached returns null on stubs, so findWithSource was empty).
		if (bestPath.isEmpty()) {
			throw new IllegalStateException("LineageReplayer: no catalog entry for ⟨"
					+ n + "," + mm + "," + p + "⟩");
		}
		return resolveStubFile(bestPath.get(), n, mm, p);
	}

	/**
	 * Replay (or STUB_CACHE-hit) a lineage-only file and orient it to the requested
	 * ⟨n,mm,p⟩ (the stub's natural shape is some canonical ordering of the axes).
	 */
	private NonCubicBilinearAlgorithm resolveStubFile(Path path, int n, int mm, int p) {
		NonCubicBilinearAlgorithm cached = STUB_CACHE.getIfPresent(path);
		NonCubicBilinearAlgorithm alg = (cached != null) ? cached : replayFromFile(path.toFile());
		if (cached == null) STUB_CACHE.put(path, alg);
		if (alg.n == n && alg.m == mm && alg.p == p) return alg;
		return alg.orientAs(n, mm, p).orElseThrow(() -> new IllegalStateException(
				"LineageReplayer: cannot orient ⟨" + alg.n + "," + alg.m + "," + alg.p
				+ "⟩ as ⟨" + n + "," + mm + "," + p + "⟩"));
	}

	/**
	 * Read a stub file and replay its lineage. Useful for the recursive
	 * descent through nested stubs (a stub may reference another stub's
	 * shape via a Leaf node).
	 */
	public NonCubicBilinearAlgorithm replayFromFile(File f) {
		// Cycle guard: re-entering a file already on this thread's replay descent
		// means the lineage references itself transitively — it would recurse forever
		// (resolveLeaf → resolveShape → replayFromFile → …) and throw an UNCATCHABLE
		// StackOverflowError that escapes every catch(Exception). Detect it here and
		// throw a catchable IllegalStateException naming the cycle path, so callers
		// (e.g. the projection per-parent catch) skip the corrupt stub and carry on.
		// Path-scoped (added on enter, removed in finally) so legitimate reuse of a
		// shared leaf across sibling branches is NOT flagged as a cycle.
		String key = f.getAbsolutePath();
		java.util.LinkedHashSet<String> active = ACTIVE_REPLAYS.get();
		if (!active.add(key)) {
			throw new IllegalStateException("LineageReplayer: lineage cycle detected: "
					+ active.stream().map(s -> new File(s).getName())
							.collect(java.util.stream.Collectors.joining(" → "))
					+ " → " + f.getName() + " (a derived stub transitively references itself)");
		}
		try {
			JsonNode root = SchemeIO.parseJson(f);
			if (!SchemeIO.isStub(root)) {
				// Not a stub: read directly as bilinear.
				return SchemeIO.readBilinear(f);
			}
			Lineage.Node ln = SchemeIO.readLineage(root)
					.orElseThrow(() -> new IllegalStateException(
							"stub " + f.getName() + " missing lineage"));
			NonCubicBilinearAlgorithm alg = replay(ln);
			// Safety net for robust (permuted-frame) replay: orient to the file's
			// DECLARED shape. A correct-frame replay already matches and this is a
			// no-op; a robust-frame replay lands on a permutation of the declared
			// shape, which orientAs corrects. If the result is NOT even a permutation
			// of the declared shape, the lineage was genuinely unrecoverable — throw
			// rather than hand back a scheme for the wrong shape.
			int[] want = declaredShape(root);
			if (want != null && (alg.n != want[0] || alg.m != want[1] || alg.p != want[2])) {
				alg = alg.orientAs(want[0], want[1], want[2]).orElseThrow(() ->
						new IllegalStateException("LineageReplayer: " + f.getName() + " replayed to ⟨"
								+ rshape(ln) + "⟩ incompatible with declared ⟨"
								+ want[0] + "," + want[1] + "," + want[2] + "⟩"));
			}
			return alg;
		} catch (IOException e) {
			throw new RuntimeException("LineageReplayer: failed to read " + f, e);
		} finally {
			active.remove(key);
		}
	}

	/** The {@code "n":[n,m,p]} shape declared in the file header, or null. */
	private static int[] declaredShape(JsonNode root) {
		JsonNode n = root.get("n");
		if (n != null && n.isArray() && n.size() == 3) {
			return new int[] { n.get(0).asInt(), n.get(1).asInt(), n.get(2).asInt() };
		}
		return null;
	}

	/** Shape of an already-built node, for error messages (cheap: re-reads memo-free). */
	private String rshape(Lineage.Node ln) {
		try { NonCubicBilinearAlgorithm a = replay(ln); return a.n + "," + a.m + "," + a.p; }
		catch (RuntimeException e) { return "?"; }
	}

	private NonCubicBilinearAlgorithm applyRecomb(Lineage.RecombinationN r,
			IdentityHashMap<Lineage.Node, NonCubicBilinearAlgorithm> memo) {
		NonCubicBilinearAlgorithm base = replayInternal(r.base(), memo);
		// The allocation arrays carry the base's axis ORIENTATION (allocA over the base's
		// n-axis, etc.). If the base node replayed to a different axis order than the
		// alloc was recorded for (a task #91-style orientation drift in an older stub —
		// e.g. a ⟨2,2,3⟩ base whose ref now resolves oriented as ⟨2,3,2⟩), re-orient the
		// base to the alloc dims when they are a permutation of the base's. Without this,
		// constructWithAllocation throws "allocation length must equal base dim" and the
		// (otherwise valid) parent/projection is wrongly discarded.
		int la = r.allocA().length, lb = r.allocB().length, lc = r.allocC().length;
		if (base.n != la || base.m != lb || base.p != lc) {
			Optional<NonCubicBilinearAlgorithm> oriented = base.orientAs(la, lb, lc);
			if (oriented.isPresent()) {
				base = oriented.get();
			}
		}
		Recombination.AlgorithmLookup recombLookup = (n, mm, p) -> {
			try { return Optional.of(resolveShape(n, mm, p)); }
			catch (Exception e) { return Optional.empty(); }
		};
		return Recombination.constructWithAllocation(base, recombLookup,
				r.allocA(), r.allocB(), r.allocC());
	}

	private NonCubicBilinearAlgorithm applyTranspose(Lineage.Transpose t,
			IdentityHashMap<Lineage.Node, NonCubicBilinearAlgorithm> memo) {
		NonCubicBilinearAlgorithm child = replayInternal(t.child(), memo);
		// Perm is "SRC->XYZ" where SRC is NMP (current) or ABC (legacy); each output
		// letter names a source axis (A/N = child.n, B/M = child.m, C/P = child.p).
		String perm = t.perm();
		int arrow = perm.indexOf("->");
		if (arrow < 0 || perm.length() - arrow - 2 < 3) {
			return child;  // unparseable legacy perm — let the final orient-to-shape compensate
		}
		// Bit-exact axis-relabel: compose transpose()/cyclicShift() for the EXACT
		// permutation rather than orientAs-by-shape, which is AMBIGUOUS when the base
		// has two equal-sized axes (⟨3,4,4⟩→⟨4,4,3⟩ has two valid orientations giving
		// different U/V/W). The shape-based orientAs picked the wrong one → the recorded
		// recombination base differed from the scored base → predict/build divergence
		// (project_recomb_base_orientation_not_pinned).
		NonCubicBilinearAlgorithm exact = eu.solven.matmul.SymmetryTransforms.permuteAxes(child, perm);
		if (exact != null) {
			return exact;
		}
		// Legacy emitters sometimes wrote a non-bijective perm (e.g. "ABC->ACA")
		// when the source had repeated dims; permuteAxes returns null there. Fall back
		// to the shape-based orient (its previous behaviour) — robust concat + the final
		// orient-to-declared-shape in replayFromFile recover the correct framing.
		Map<Character, Integer> src = Map.of(
				'N', child.n, 'M', child.m, 'P', child.p,
				'A', child.n, 'B', child.m, 'C', child.p);
		String tgt = perm.substring(arrow + 2);
		Integer newN = src.get(tgt.charAt(0)), newM = src.get(tgt.charAt(1)), newP = src.get(tgt.charAt(2));
		if (newN == null || newM == null || newP == null) {
			return child;  // unknown axis letters — compensate downstream
		}
		return child.orientAs(newN, newM, newP).orElse(child);
	}

	/**
	 * Replay an {@link Lineage.AxisFlip} node: reconstruct the variant
	 * from the child by applying the specified axis-reverse mask. The
	 * variant is one of the 8 elements of
	 * {@link eu.solven.matmul.SymmetryTransforms#axisFlipOrbit};
	 * indexing matches that method's mask iteration ({@code bit 0 =
	 * swapA, bit 1 = swapB, bit 2 = swapC}).
	 */
	private NonCubicBilinearAlgorithm applyAxisFlip(Lineage.AxisFlip af,
			IdentityHashMap<Lineage.Node, NonCubicBilinearAlgorithm> memo) {
		NonCubicBilinearAlgorithm child = replayInternal(af.child(), memo);
		List<NonCubicBilinearAlgorithm> orbit =
				eu.solven.matmul.SymmetryTransforms.axisFlipOrbit(child);
		int idx = af.mask();
		if (idx < 0 || idx >= orbit.size()) {
			throw new IllegalStateException("LineageReplayer: axis-flip mask "
					+ idx + " out of orbit range [0," + orbit.size() + ")");
		}
		return orbit.get(idx);
	}

	/**
	 * Replay an {@link Lineage.AxisPermute} node: apply per-axis
	 * permutations to the child's factor matrices via the same logic
	 * used internally by
	 * {@link eu.solven.matmul.SymmetryTransforms#permutationOrbit}.
	 */
	private NonCubicBilinearAlgorithm applyAxisPermute(Lineage.AxisPermute ap,
			IdentityHashMap<Lineage.Node, NonCubicBilinearAlgorithm> memo) {
		NonCubicBilinearAlgorithm child = replayInternal(ap.child(), memo);
		// Reuse the full enumeration to dispatch via the same code path
		// (cheap — permA/B/C are small) and locate the matching entry.
		// Simpler: re-implement the same forward transform here.
		int n = child.n, m = child.m, p = child.p, r = child.r;
		double[][] srcU = child.denseU();
		double[][] srcV = child.denseV();
		double[][] srcW = child.denseW();
		double[][] U2 = new double[n * m][r];
		double[][] V2 = new double[m * p][r];
		double[][] W2 = new double[n * p][r];
		for (int i = 0; i < n; i++) for (int j = 0; j < m; j++)
			for (int k = 0; k < r; k++)
				U2[ap.permA()[i] * m + ap.permB()[j]][k] = srcU[i * m + j][k];
		for (int j = 0; j < m; j++) for (int l = 0; l < p; l++)
			for (int k = 0; k < r; k++)
				V2[ap.permB()[j] * p + ap.permC()[l]][k] = srcV[j * p + l][k];
		for (int i = 0; i < n; i++) for (int l = 0; l < p; l++)
			for (int k = 0; k < r; k++)
				W2[ap.permA()[i] * p + ap.permC()[l]][k] = srcW[i * p + l][k];
		return new NonCubicBilinearAlgorithm(n, m, p, U2, V2, W2);
	}
}
