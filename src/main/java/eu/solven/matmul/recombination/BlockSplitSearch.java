package eu.solven.matmul.recombination;

import eu.solven.matmul.search.MethodCatalog;
import eu.solven.matmul.search.SearchHeuristics;

import eu.solven.matmul.search.AssignmentOptimizer;
import eu.solven.matmul.search.SearchBudget;

import eu.solven.matmul.search.ConcatSplitSearch;
import eu.solven.matmul.search.ConstructiveMethod;
import eu.solven.matmul.search.KnownTauIdentities;
import eu.solven.matmul.search.KroneckerSplitSearch;
import eu.solven.matmul.search.PairFusedRecombination;
import eu.solven.matmul.search.RecombinationPoolConfig;

import lombok.extern.slf4j.Slf4j;

import eu.solven.matmul.catalog.CatalogLimits;

import eu.solven.matmul.isotropy.PairedSubProducts;

import eu.solven.matmul.catalog.Compose;

import eu.solven.matmul.catalog.SchemeIO;

import eu.solven.matmul.recombination.Recombination;

import eu.solven.matmul.papers.dis2009.PanTrilinearAggregation;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.SymmetryTransforms;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.Function;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Phase 1 of the block-split / Sedoglavic procedure:
 * <strong>cheap, formula-only search</strong> over splits {@code n = u + v}
 * (with {@code u > v ≥ 1}) for cubic targets {@code ⟨n,n,n⟩}, scoring each
 * candidate by Sedoglavic's closed-form identity:
 *
 * <pre>R(⟨u+v,u+v,u+v⟩) ≤ R(⟨u,u,u⟩) + 3·R(⟨u,u,v⟩) + 3·R(⟨u,v,v⟩)</pre>
 *
 * <p>Phase 2 ({@link Compose#blockSplitCubic}) takes the chosen {@code (u, v)}
 * and actually constructs the algorithm; see
 * <a href="../../../../../../../TRILINEAR_AGGREGATION.md">TRILINEAR_AGGREGATION.md</a>
 * §5 for the recipe table.
 */
@Slf4j
public class BlockSplitSearch {

	/** Result of evaluating one split candidate. */
	public record SplitCandidate(
			int n,
			int u,
			int v,
			int formulaRank,
			int rUuu,
			int rUuv,
			int rUvv) {

		public String breakdown() {
			return rUuu + " + 3·" + rUuv + " + 3·" + rUvv + " = " + formulaRank;
		}
	}

	/** Non-cubic split candidate: per-axis split + Strassen-via-recombine rank. */
	public record NonCubicSplitCandidate(
			int n, int m, int p,
			int[] allocA, int[] allocB, int[] allocC,
			long rank) {

		public String breakdown() {
			return "Strassen + allocs " + java.util.Arrays.toString(allocA)
					+ "/" + java.util.Arrays.toString(allocB)
					+ "/" + java.util.Arrays.toString(allocC)
					+ " = " + rank;
		}
	}

	/**
	 * Multi-base search candidate: which base + allocation won.
	 */
	public record MultiBaseSplitCandidate(
			int n, int m, int p,
			NonCubicBilinearAlgorithm base,
			String baseLabel,
			int[] allocA, int[] allocB, int[] allocC,
			long rank,
			eu.solven.matmul.catalog.Lineage.Node baseOriginLineage) {

		/** Back-compat: when the search isn't pool-aware (no NamedBase available). */
		public MultiBaseSplitCandidate(int n, int m, int p,
				NonCubicBilinearAlgorithm base, String baseLabel,
				int[] allocA, int[] allocB, int[] allocC, long rank) {
			this(n, m, p, base, baseLabel, allocA, allocB, allocC, rank, null);
		}

		public String breakdown() {
			return baseLabel + " outer + allocs "
					+ java.util.Arrays.toString(allocA) + "/"
					+ java.util.Arrays.toString(allocB) + "/"
					+ java.util.Arrays.toString(allocC)
					+ " → " + rank + " mults";
		}
	}

	/**
	 * Default outer-base pool for {@link #findBestStrategy}: {@link #rootPool}
	 * with its axis-flip orbit expanded (each canonical scheme contributes
	 * up to 8 axis-flip variants).
	 *
	 * <p>Previously this returned {@code rootPool()} verbatim with a
	 * pre-axflipped Winograd variant hand-curated into the pool. After
	 * the axis-flip dedup (task #110), {@code rootPool()} contains
	 * canonical-only entries, and orbit expansion happens here so that
	 * callers retain the ⟨17,17,17⟩=2930 coverage and equivalents at
	 * other unbalanced cubic targets.</p>
	 *
	 * <p>Misses are tolerated: if a scheme file isn't found locally the
	 * pool just omits that base. Callers can always pass their own pool.</p>
	 */
	public static List<NamedBase> defaultPool() {
		List<NamedBase> raw = rootPool();
		List<NamedBase> out = new ArrayList<>();
		java.util.Set<String> sigSeen = new java.util.HashSet<>();
		for (NamedBase nb : raw) {
			NonCubicBilinearAlgorithm a = nb.base();
			eu.solven.matmul.catalog.Lineage.Node canonicalOrigin = nb.originLineage();
			// Record each flip variant as a precise AxisPermute (flips are reversal perms),
			// so replay reconstructs the EXACT base orientation rather than the canonical
			// (task #91) — the AxisFlip orbit-INDEX it used before is ambiguous across the
			// deduped orbit and silently lost the orientation.
			List<SymmetryTransforms.PermutedVariant> variants =
					SymmetryTransforms.internalOrbitWithPerms(a, SymmetryTransforms.InternalOrbitMode.AXIS_FLIP, 8);
			for (int idx = 0; idx < variants.size(); idx++) {
				SymmetryTransforms.PermutedVariant pv = variants.get(idx);
				NonCubicBilinearAlgorithm variant = pv.alg();
				String sig = sigOf(variant);
				if (!sigSeen.add(sig)) continue;
				eu.solven.matmul.catalog.Lineage.Node origin = canonicalOrigin;
				if (origin != null && idx > 0) {
					origin = new eu.solven.matmul.catalog.Lineage.AxisPermute(
							canonicalOrigin, pv.permX(), pv.permY(), pv.permZ());
				}
				String label = idx == 0 ? nb.label() : nb.label() + " :: axflip" + idx;
				out.add(new NamedBase(label, variant, origin));
			}
		}
		return out;
	}

	/**
	 * Hand-curated list of <strong>historical root</strong> templates:
	 * schemes whose factor matrices were tabulated or discovered as the
	 * first non-trivial bound for their shape, NOT derived from any
	 * composition of smaller schemes in this catalog.
	 *
	 * <p>Each root contributes its full S₃-orbit (up to 6 distinct
	 * shape-permuted variants via {@link SymmetryTransforms#s3Orbit}),
	 * so e.g. {@code ⟨2,2,3⟩=11} also enables outer splits along the
	 * {@code ⟨3,2,2⟩} and {@code ⟨2,3,2⟩} block layouts.</p>
	 *
	 * <p>Bases here are <strong>NC, Z-arithmetic</strong> so the result
	 * lifts to recursive matmul over arbitrary rings. Commutative-only
	 * roots (Waksman, Makarov, Rosowski) and field-restricted ones
	 * (AT/F2, AE/C-only) are excluded.</p>
	 */
	/**
	 * Single config-driven entry point. Replaces ad-hoc rootPool/extendedPool
	 * calls in drivers. Filters the catalog by
	 * {@link RecombinationPoolConfig#cubicOnly}, {@link RecombinationPoolConfig#rootsOnly},
	 * {@link RecombinationPoolConfig#maxBaseDim}, then expands each retained scheme
	 * by the configured {@link RecombinationPoolConfig#orbitMode}.
	 *
	 * <p>The two pre-existing helpers {@link #rootPool()} and
	 * {@link #extendedPool(int)} remain for back-compat; both are now
	 * special cases of {@code buildPool}.</p>
	 *
	 * @see RecombinationPoolConfig#simple
	 */
	public static List<NamedBase> buildPool(RecombinationPoolConfig config) {
		return buildPool(config, null);
	}

	/**
	 * Field-aware pool: the extended outer bases are filtered to those VALID over
	 * {@code fieldTag} (i.e. {@code fieldTag ∈ fields[]} — the catalog encodes every
	 * field a scheme is valid over: an integer scheme carries {@code [F2,F3,Z,Q,R,C]},
	 * an F₂-native scheme carries {@code [F2]}). So an {@code F2} sweep gets F₂-native
	 * bases (AlphaTensor ⟨4,4,4⟩=47) AND the integer roots that reduce mod 2, while a
	 * {@code Q} sweep gets only char-0-valid bases. {@code fieldTag == null} keeps the
	 * legacy char-0 (Z/Q/R) filter. The hand-curated {@code rootPool} templates are
	 * integer (Z) so they are valid over every field and are always included.
	 */
	public static List<NamedBase> buildPool(RecombinationPoolConfig config, String fieldTag) {
		// Start from rootPool's hand-curated list of historical
		// templates (Strassen + AT⟨2,2,3⟩ + AT⟨2,3,3⟩ + Laderman) —
		// these are always Leaf-lineage. Always applied first so the
		// well-known templates lead the enumeration order.
		List<NamedBase> raw = rootPool();
		if (!config.rootsOnly()) {
			// Append every Leaf-lineage NC catalog scheme valid over the field, up to
			// maxBaseDim. extendedPool already inlines rootPool, so we dedup-merge.
			// Dedup by CONTENT, not just (shape, rank): two DISTINCT schemes at the same
			// ⟨n,m,p⟩=r recombine DIFFERENTLY (different product supports → different
			// per-block effective dims), so dropping one as a "duplicate" silently loses the
			// better recombination. ⟨5,20,26⟩ via the alphatensor_Z ⟨2,4,4⟩=26 base reaches
			// 1702, but the hk71 ⟨2,4,4⟩=26 only 1716 — keying on "2x4x4:r=26" kept whichever
			// came first and lost master's 1702 (the residual large-unbalanced regressions —
			// the missing base, NOT the optimizer's cost model, which is exact).
			java.util.Set<String> have = new java.util.HashSet<>();
			for (NamedBase nb : raw) have.add(poolContentKey(nb.base()));
			for (NamedBase nb : extendedPool(config.maxBaseDim(), fieldTag)) {
				if (have.add(poolContentKey(nb.base()))) raw.add(nb);
			}
		}
		List<NamedBase> out = new ArrayList<>();
		java.util.Set<String> sigSeen = new java.util.HashSet<>();
		for (NamedBase nb : raw) {
			NonCubicBilinearAlgorithm a = nb.base();
			boolean cubicBase = (a.n == a.m && a.m == a.p);
			// Two-threshold cap (user 2026-06-26): non-cubic bases ≤ maxBaseDim; cubic bases
			// may reach maxCubicBaseDim (≥ maxBaseDim) so e.g. ⟨7,7,7⟩ stays in a pool whose
			// rectangular bases stop at 5. cubicOnly still drops every non-cubic base.
			int dimCap = cubicBase
					? Math.max(config.maxBaseDim(), config.maxCubicBaseDim())
					: config.maxBaseDim();
			if (Math.max(a.n, Math.max(a.m, a.p)) > dimCap) continue;
			if (config.cubicOnly() && !cubicBase) continue;
			eu.solven.matmul.catalog.Lineage.Node canonicalOrigin = nb.originLineage();
			// Each variant carries the EXACT per-axis permutation that reconstructs it from
			// the canonical scheme (flips are reversal permutations) — so the origin lineage
			// is a precise, replayable AxisPermute, not a lossy AxisFlip orbit-INDEX or an
			// identity placeholder. This is the task #91 fix: previously a recombination that
			// won on a permuted/flipped base recorded a lineage that replayed to the CANONICAL
			// base, yielding a valid-but-different-rank scheme the write-guard then discarded.
			List<SymmetryTransforms.PermutedVariant> variants =
					SymmetryTransforms.internalOrbitWithPerms(a, config.orbitMode(), config.permutationOrbitCap());
			for (int idx = 0; idx < variants.size(); idx++) {
				SymmetryTransforms.PermutedVariant pv = variants.get(idx);
				NonCubicBilinearAlgorithm variant = pv.alg();
				String sig = sigOf(variant);
				if (!sigSeen.add(sig)) continue;
				eu.solven.matmul.catalog.Lineage.Node origin = canonicalOrigin;
				if (origin != null && idx > 0) {
					origin = new eu.solven.matmul.catalog.Lineage.AxisPermute(
							canonicalOrigin, pv.permX(), pv.permY(), pv.permZ());
				}
				out.add(new NamedBase(nb.label() + " :: " + config.orbitMode().name(),
						variant, origin));
			}
		}
		return out;
	}

	private static int[] identity(int n) {
		int[] p = new int[n];
		for (int i = 0; i < n; i++) p[i] = i;
		return p;
	}

	/** Pool dedup key = (shape, rank) + the CANONICAL MULTISET of per-product SUPPORTS.
	 *  The recombination cost {@code Σ R(per-product effective dims)} depends ONLY on which
	 *  blocks each product touches on each axis — NOT on the coefficient values. So two schemes
	 *  with the same support multiset tile EVERY target identically (dedup), while two with
	 *  different supports recombine differently and must BOTH be kept (the hk71 vs alphatensor_Z
	 *  ⟨2,4,4⟩=26 case: 1716 vs 1700 on ⟨5,20,26⟩). This is coarser than a content hash (which
	 *  over-splits on gauge/coefficient differences that don't affect tiling) and far finer than
	 *  the old {@code shape:r} key (which merged recombination-distinct bases — the bug). The
	 *  per-product tuple is sorted within each axis-support and the products are sorted, so the
	 *  key is invariant to product reordering. [[project_optimizer_pairing_blind]] */
	private static String poolContentKey(NonCubicBilinearAlgorithm b) {
		AnalyticalMaskSearch.SchemeSupports sup = AnalyticalMaskSearch.SchemeSupports.extract(b);
		java.util.List<String> products = new java.util.ArrayList<>(b.r);
		for (int k = 0; k < b.r; k++) {
			products.add(sortedSig(sup.uRowSupport[k]) + "/" + sortedSig(sup.uColSupport[k]) + "/"
					+ sortedSig(sup.vRowSupport[k]) + "/" + sortedSig(sup.vColSupport[k]) + "/"
					+ sortedSig(sup.wRowSupport[k]) + "/" + sortedSig(sup.wColSupport[k]));
		}
		java.util.Collections.sort(products);
		return b.n + "x" + b.m + "x" + b.p + ":r=" + b.r + ":" + String.join(";", products);
	}

	private static String sortedSig(int[] support) {
		int[] s = support.clone();
		java.util.Arrays.sort(s);
		return java.util.Arrays.toString(s);
	}

	private static String sigOf(NonCubicBilinearAlgorithm a) {
		StringBuilder sb = new StringBuilder();
		sb.append(a.n).append('x').append(a.m).append('x').append(a.p).append('|');
		for (double[] row : a.denseU()) for (double v : row) sb.append(v).append(',');
		sb.append('|');
		for (double[] row : a.denseV()) for (double v : row) sb.append(v).append(',');
		sb.append('|');
		for (double[] row : a.denseW()) for (double v : row) sb.append(v).append(',');
		return sb.toString();
	}

	public static List<NamedBase> rootPool() {
		List<NamedBase> pool = new ArrayList<>();
		// ── Trivial axis-splitters (DIS09 mul211/mul121/mul112) ──
		// Each is a rank-2 ⟨1,2,1⟩ / ⟨2,1,1⟩ / ⟨1,1,2⟩ scheme that lets the
		// recombination framework express a pure single-axis split (no
		// Strassen-style mixing). DIS09 includes these in its pool because
		// they unify axis-splits with the rest of the search; the S₃ orbit
		// generator below produces the other 2 axis variants automatically.
		pool.add(new NamedBase("AxisSplit<2,1,1>=2", eu.solven.matmul.AxisSplitBases.mul211()));
		pool.add(new NamedBase("AxisSplit<1,2,1>=2", eu.solven.matmul.AxisSplitBases.mul121()));
		pool.add(new NamedBase("AxisSplit<1,1,2>=2", eu.solven.matmul.AxisSplitBases.mul112()));
		// ── Rank-4 naïve grids ⟨1,2,2⟩ / ⟨2,1,2⟩ / ⟨2,2,1⟩ ──
		// Single-block 2×2 grids with one axis unsplit. Unlike the rank-2 AxisSplits,
		// these carry a disjoint cyclic-rotation product PAIR (e.g. ⟨n,b,c⟩ & its rot²
		// ⟨n,c,b⟩ under a symmetric m=p split) that Pan trilinear aggregation fuses at
		// abc+ab+bc+ca < 2·R — the in-recombination TA path (Recombination.tryBuildTaFusion,
		// gated on isNaiveGrid). Without them the search cannot reach the naive-1x2x2
		// RecombinationTa schemes master held (⟨26,29,29⟩=11693, ⟨28,31,31⟩=14043): TA is a
		// saving WITHIN a recombination's multiplications, not a separate strategy, so the
		// grid must be IN the pool for the saving to be weighed. All three axis orientations
		// are added so the unsplit axis can align with any target's heaviest axis.
		pool.add(new NamedBase("Naive<1,2,2>=4", NonCubicBilinearAlgorithm.naive(1, 2, 2)));
		pool.add(new NamedBase("Naive<2,1,2>=4", NonCubicBilinearAlgorithm.naive(2, 1, 2)));
		pool.add(new NamedBase("Naive<2,2,1>=4", NonCubicBilinearAlgorithm.naive(2, 2, 1)));
		// ── ⟨2,*,*⟩ ──
		// Strassen / Winograd 2×2×2 have the full Burichenko-order-36
		// stabilizer including S₃ axis-permutation → cubicSymmetric.
		tryAddOrbit(pool, "Strassen<2,2,2>=7",
				"src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json",
				true);
		// Winograd 1971 — canonical form of the rank-7 ⟨2,2,2⟩ algorithm
		// that reaches ⟨17,17,17⟩=2930 (sub-shape distribution
		// 1×⟨9,9,9⟩+4×⟨8,9,9⟩-cyc+1×⟨8,8,9⟩+1×⟨8,8,8⟩ = 486+4·430+388+336).
		// NO base-matrix flip is required: for a 2×2×2 base the axis-flip is
		// the per-axis allocation ORDER, so canonical Winograd at the (8,9)³
		// ordering already gives 2930 (the mirrored (9,8)³ ordering gives
		// 2954). Verified by TestStrassenCousinHunt (canonical @ (8,9)³=2930,
		// @ (9,8)³=2954) and TestProbe17ConfigCompare (RecombinationPoolConfig.simple()
		// CANONICAL and auditAxisFlip BOTH return 2930 from base
		// "Winograd :: CANONICAL/AXIS_FLIP" at alloc [8,9]³). The axis-flip
		// orbit machinery (#106 AnalyticalMaskSearch, #108 axisFlipCanonical,
		// #133 lifting the 2×2×2 mask-sweep gate, now native in
		// findBestMultiBaseSplit) is what made this (8,9)³ route reachable
		// from the canonical pool entry — it is not used to flip the matrix
		// at the winning allocation. Earlier this entry was the pre-axflipped
		// `solven-winograd-cousin-axflip1`; #110 replaced it with canonical
		// Winograd-1971 since the orbit/allocation machinery now covers the
		// same ground without pool duplication. Discovery: TestStrassenCousinHunt.
		tryAddOrbit(pool, "Winograd<2,2,2>=7",
				"src/main/resources/schemes/known/section2/winograd_1971-2x2x2_m7_a24.json",
				true);
		tryAddOrbit(pool, "AT<2,2,3>=11",
				"src/main/resources/schemes/known/section3/alphatensor_Z-2x2x3_m11_a25.json");
		tryAddOrbit(pool, "AT<2,2,4>=14",
				"src/main/resources/schemes/known/section4/alphatensor_Z-2x2x4_m14_a48.json");
		// ── ⟨2,3,*⟩, ⟨3,3,*⟩ ──
		tryAddOrbit(pool, "AT<2,3,3>=15",
				"src/main/resources/schemes/known/section3/alphatensor_Z-2x3x3_m15_a58.json");
		tryAddOrbit(pool, "Laderman<3,3,3>=23",
				"src/main/resources/schemes/known/section3/laderman_1976-3x3x3_m23_a98.json");
		tryAddOrbit(pool, "AT<3,3,4>=29",
				"src/main/resources/schemes/known/section4/alphatensor_Z-3x3x4_m29_a148.json");
		// ── ⟨3,3,6⟩=40: Smirnov 2013 (was imported via fmm-lille; re-attributed
		//     by the 2026-06 rename — hint updated to the current filename) ──
		tryAddOrbit(pool, "Smirnov<3,3,6>=40",
				"src/main/resources/schemes/known/section6/3x3x6-r40-smirnov13-d246191.json");
		// ── cubic ⟨4,4,4⟩=49 (AT-Z) — critical for the `simple` config
		//     where cubicOnly filters everything else above out. ──
		tryAddOrbit(pool, "AT<4,4,4>=49",
				"src/main/resources/schemes/known/section4/alphatensor_Z-4x4x4_m49_a468.json");
		// ── ⟨5,5,5⟩ at Z: best is AlphaEvolve 93 (cubic) ──
		tryAddOrbit(pool, "AE<5,5,5>=93",
				"src/main/resources/schemes/known/section5/alphaevolve-5x5x5_m93_a846.json");
		// ── ⟨7,7,7⟩=250 (Sedoglavic via Perminov-ZT mirror; cubic) ──
		tryAddOrbit(pool, "Sedoglavic<7,7,7>=250",
				"src/main/resources/schemes/known/section7/perminov_ZT-7x7x7_m250_a2417.json");
		return pool;
	}

	/**
	 * Exhaustive pool: every catalog scheme whose lineage is a
	 * {@code Leaf} (or has no lineage — older root entries), filtered to
	 * NC Z/Q/ZT-arithmetic. Used for thorough audits where
	 * {@link #rootPool} doesn't surface the closing template.
	 *
	 * <p><strong>Slow.</strong> Enumeration cost in
	 * {@link #findBestMultiBaseSplit} scales with C(n-1, bn-1)·
	 * C(m-1, bm-1)·C(p-1, bp-1) per base — for n, m, p ≤ 32 and bases
	 * up to ⟨5,5,5⟩ this can explode without the
	 * {@code maxCombosPerBase} auto-fallback to balanced enumeration.</p>
	 *
	 * @param maxBaseDim cap on {@code max(bn, bm, bp)} of each base —
	 *                   keeps the pool focused on actually-useful outer
	 *                   templates (recommend 5–8). Anything above is
	 *                   better used as a Kronecker factor / SOTA leaf.
	 */
	public static List<NamedBase> extendedPool(int maxBaseDim) {
		return extendedPool(maxBaseDim, null);
	}

	/** Field-aware {@link #extendedPool}: keep only leaf bases VALID over {@code fieldTag}
	 *  ({@code fieldTag ∈ fields[]}); {@code null} = legacy char-0 (Z/Q/R) filter. */
	public static List<NamedBase> extendedPool(int maxBaseDim, String fieldTag) {
		List<NamedBase> pool = new ArrayList<>();
		java.util.Set<String> seenShape = new java.util.HashSet<>();
		// Add the rootPool first so the well-known templates appear at
		// the front of enumeration.
		for (NamedBase nb : rootPool()) {
			NonCubicBilinearAlgorithm a = nb.base();
			if (Math.max(a.n, Math.max(a.m, a.p)) > maxBaseDim) continue;
			seenShape.add(poolContentKey(a));
			pool.add(nb);
		}
		// Walk the WHOLE schemes tree, not just top-level sectionN: the 2026-06-09
		// known/derived/curated split moved every sectionN dir under those subfolders,
		// so the old listFiles("section*") on the root silently loaded NOTHING (the
		// extended pool collapsed to just the seed templates). Recurse so the pool
		// sees every catalog leaf again, regardless of which subtree it lives in.
		java.nio.file.Path schemesDir = java.nio.file.Path.of("src/main/resources/schemes");
		List<java.io.File> jsonFiles = new java.util.ArrayList<>();
		try (var walk = java.nio.file.Files.walk(schemesDir)) {
			walk.filter(p -> p.getFileName().toString().endsWith(".json"))
					.sorted()
					.forEach(p -> jsonFiles.add(p.toFile()));
		} catch (java.io.IOException e) {
			return pool;
		}
		for (java.io.File f : jsonFiles) {
			try {
				log.debug("Loading {}", f);
				if (!isFieldValidLeafNC(f, maxBaseDim, fieldTag)) continue;
				NonCubicBilinearAlgorithm alg = SchemeIO.read(f);
				if (Math.max(alg.n, Math.max(alg.m, alg.p)) > maxBaseDim) continue;
				// Same orbit-expansion semantics as rootPool: shape orbit only,
				// one representative per distinct shape. Algorithmic variants
				// (axis-flip, permutation) are opt-in.
				// Pin the canonical by CONTENT HASH (resolves to the EXACT loaded scheme,
				// not the shape-ambiguous filename atom) and record each oriented variant
				// with its EXACT axis-relabel. Without a pinned origin the variant fell back
				// to an Atom(label) whose replay re-orients by shape — AMBIGUOUS for equal
				// axes (⟨3,4,4⟩→⟨4,4,3⟩) → the recombination base replayed to a DIFFERENT
				// orientation than the one scored → predict/build divergence
				// (project_recomb_base_orientation_not_pinned).
				eu.solven.matmul.catalog.Lineage.Node canonicalLeaf =
						new eu.solven.matmul.catalog.Lineage.Atom(alg.n + "x" + alg.m + "x" + alg.p
								+ "@" + SchemeIO.contentHash(alg));
				List<SymmetryTransforms.S3Variant> orbit = SymmetryTransforms.s3OrbitWithPerms(alg);
				java.util.Map<String, SymmetryTransforms.S3Variant> byShape =
						new java.util.LinkedHashMap<>();
				for (SymmetryTransforms.S3Variant v : orbit) {
					byShape.putIfAbsent(v.alg().n + "x" + v.alg().m + "x" + v.alg().p, v);
				}
				for (SymmetryTransforms.S3Variant v : byShape.values()) {
					NonCubicBilinearAlgorithm a = v.alg();
					if (Math.max(a.n, Math.max(a.m, a.p)) > maxBaseDim) continue;
					String key = poolContentKey(a);
					if (seenShape.add(key)) {
						String suffix = byShape.size() == 1
								? ""
								: " / <" + a.n + "," + a.m + "," + a.p + ">=" + a.r;
						eu.solven.matmul.catalog.Lineage.Node origin = "ABC->ABC".equals(v.perm())
								? canonicalLeaf
								: new eu.solven.matmul.catalog.Lineage.Transpose(canonicalLeaf, v.perm());
						pool.add(new NamedBase("ext[" + f.getName() + "]" + suffix, a, origin));
					}
				}
			} catch (Exception ignored) {
				// silent: skip un-readable
			}
		}
		return pool;
	}

	/**
	 * Quick text check: scheme file is Leaf-lineage (or no lineage), NC,
	 * and Z/Q/ZT field. We grep the file header rather than fully parsing
	 * to skip the heavy expansion on schemes we'll drop anyway.
	 */
	static boolean isFieldValidLeafNC(java.io.File f, int maxBaseDim, String fieldTag) {
		try {
			// Read the header (everything before the u/v/w payload, which the
			// canonical formatter writes AFTER n/m/hash/fields/source). 4 KB covers
			// even a long lineage_str.
			byte[] buf = new byte[4096];
			try (java.io.FileInputStream in = new java.io.FileInputStream(f)) {
				int read = in.read(buf);
				if (read <= 0) return false;
			}
			String header = new String(buf, java.nio.charset.StandardCharsets.UTF_8);
			// PAYLOAD-FIRST files (imports whose key order puts u/v/w before the
			// metadata, e.g. the flips_mod2 family) keep fields[]/commutative at the
			// TAIL — outside this window. Treating them as "legacy, no fields[]"
			// silently admitted F₂/C-only bases into char-0 pools (verify-gate caught
			// the bogus results, but the pool was polluted). If the window shows no
			// fields[] and the file is bigger than the window, scan the WHOLE file.
			if (!header.contains("\"fields\"") && f.length() > buf.length) {
				header = java.nio.file.Files.readString(f.toPath());
			}
			// Reject commutative (cannot be an NC recombination base); reject composed
			// lineage (an outer base must be a Leaf-lineage atom).
			if (header.contains("\"commutative\": true")) return false;
			if (header.contains("\"op\":\"RecombinationN\"")
					|| header.contains("\"op\":\"RecombinationTa")
					|| header.contains("\"op\":\"RecombinationWithPair")
					|| header.contains("\"op\":\"KronProduct\"")
					|| header.contains("\"op\":\"KronChain\"")
					|| header.contains("\"op\":\"ConcatCols\"")
					|| header.contains("\"op\":\"ConcatRows\"")
					|| header.contains("\"op\":\"SumInner\"")
					|| header.contains("\"op\":\"ConcatRight\"")  // legacy alias
					|| header.contains("\"op\":\"ConcatBelow\"")  // legacy alias
					|| header.contains("\"op\":\"Transpose\"")) {
				return false;
			}
			// Field filter. The catalog uses a `fields[]` ARRAY (the old singular `"field"`
			// key is gone — checking it silently let EVERY F₂-only / C-only base leak into
			// the char-0 pool, e.g. AlphaTensor ⟨4,4,4⟩=47/F₂ as a Q-sweep base → a phantom
			// 9316). A base is valid for a sweep over `fieldTag` iff `fieldTag ∈ fields[]`
			// (the catalog records every field a scheme is valid over: an integer scheme
			// carries all of [F2,F3,Z,Q,R,C]). `fieldTag == null` ⇒ legacy char-0 (Z/Q/R).
			int fi = header.indexOf("\"fields\"");
			if (fi >= 0) {
				int lb = header.indexOf('[', fi);
				int rb = lb >= 0 ? header.indexOf(']', lb) : -1;
				if (lb >= 0 && rb > lb) {
					String arr = header.substring(lb, rb);
					if (fieldTag != null) {
						// INCLUSION-aware (Z⊂Q⊂R⊂C): a scheme valid over a SUB-field is valid
						// over fieldTag. A bare exact `fieldTag ∈ fields[]` starved the pool — a
						// Z/Q base stamped ["Z"] was rejected from an R sweep, collapsing the R
						// extended pool from 103 to 5 (Strassen-only). Mirror FieldAwareLookup's
						// field.fallbackChain() so the pool matches the lookup.
						for (eu.solven.matmul.algebra.Field allowed
								: eu.solven.matmul.algebra.Field.fromTag(fieldTag).fallbackChain()) {
							if (arr.contains("\"" + allowed.tag() + "\"")) return true;
						}
						return false;
					}
					return arr.contains("\"Z\"") || arr.contains("\"Q\"") || arr.contains("\"R\"");
				}
			}
			// Legacy / no fields[]: such files are the catalog's hand-curated INTEGER (Z)
			// atoms (e.g. the un-stamped 2×2×2 shims), valid over EVERY field — accept
			// unless an explicit singular C/F₂ marker says otherwise.
			String normalised = header.replaceAll("\"field\"\\s*:\\s*", "\"field\":");
			if (normalised.contains("\"field\":\"C\"")
					|| normalised.contains("\"field\":\"F2\"")
					|| normalised.contains("\"field\":\"F_2\"")) {
				return false;
			}
			return true;
		} catch (Exception e) {
			return false;
		}
	}

	private static void tryAdd(List<NamedBase> pool, String label, String path) {
		try {
			java.io.File f = new java.io.File(path);
			if (f.exists()) pool.add(new NamedBase(label, SchemeIO.read(f)));
		} catch (Exception ignored) {
			// silent: the pool is a best-effort default
		}
	}

	/**
	 * Like {@link #tryAdd} but adds every distinct shape in the base's
	 * S₃-orbit (up to 6 entries). The label is suffixed by the
	 * orbit member's shape to disambiguate in reports.
	 */
	private static void tryAddOrbit(List<NamedBase> pool, String baseLabel, String path) {
		tryAddOrbit(pool, baseLabel, path, false);
	}

	/** The stable part of a scheme filename: everything up to and including the
	 *  {@code _m{rank}} token. The trailing tokens ({@code _a{adds}},
	 *  {@code _b{border}}, field suffix) drift with filename-convention changes
	 *  (#173, the {@code _b0} border token, additions re-optimisation), so the
	 *  hand-curated {@link #rootPool} entries must NOT hard-code them. */
	private static final java.util.regex.Pattern SCHEME_STEM =
			java.util.regex.Pattern.compile("^.*_m\\d+");

	/**
	 * Resolve a hand-curated pool file by its stable
	 * {@code {source}-{shape}_m{rank}} stem, tolerant of drifting trailing tokens.
	 * Returns the exact path if it exists; otherwise globs the directory for a file
	 * sharing the stem (token boundary respected, so {@code _m7} never matches
	 * {@code _m70}) and picks the shortest name — the one with the fewest extra
	 * tokens. Falls back to the literal path so {@link #tryAddOrbit}'s
	 * {@code exists()} check still skips cleanly when nothing matches.
	 */
	static java.io.File resolvePoolFile(String hintPath) {
		// Content base-resolver (shape from the catalog index, source token disambiguates):
		// tolerant of the known/derived/curated folder split, the 2026 filename rename, and
		// historical _a/_b token drift. Replaces the old parent-dir stem-glob, which the
		// rename broke (new names no longer start with "{source}-{shape}_m{rank}").
		return eu.solven.matmul.catalog.SchemeResolver.byHint(hintPath);
	}

	/**
	 * @param cubicSymmetric mark every NamedBase added by this call with
	 *     {@code cubicSymmetric=true}. Only set for bases known to have
	 *     a stabilizer that includes the full S₃ axis-permutation subgroup
	 *     (Strassen / Winograd 2×2×2 — Burichenko-order-36 stabilizer).
	 */
	private static void tryAddOrbit(List<NamedBase> pool, String baseLabel, String path,
			boolean cubicSymmetric) {
		try {
			java.io.File f = resolvePoolFile(path);
			if (!f.exists()) return;
			NonCubicBilinearAlgorithm root = SchemeIO.read(f);
			// Pin the base by CONTENT HASH, not the shape-ambiguous filename atom. The
			// filename atom resolves to the rank-BEST scheme of this shape, which need NOT
			// be the exact scheme this file holds (e.g. the winograd ⟨2,2,2⟩ file hashes to
			// one variant, while the best ⟨2,2,2⟩ is another). A recombination on this base
			// then replayed against a DIFFERENT base → a valid-but-wrong-rank scheme that the
			// write-guard discards (task #91). The hash-ref resolves to the exact loaded base.
			eu.solven.matmul.catalog.Lineage.Node canonicalLeaf =
					new eu.solven.matmul.catalog.Lineage.Atom(root.n + "x" + root.m + "x" + root.p
							+ "@" + eu.solven.matmul.catalog.SchemeIO.contentHash(root));
			List<SymmetryTransforms.S3Variant> orbit = SymmetryTransforms.s3OrbitWithPerms(root);
			java.util.Map<String, SymmetryTransforms.S3Variant> byShape =
					new java.util.LinkedHashMap<>();
			for (SymmetryTransforms.S3Variant v : orbit) {
				byShape.putIfAbsent(v.alg().n + "x" + v.alg().m + "x" + v.alg().p, v);
			}
			if (byShape.size() == 1) {
				NamedBase nb = new NamedBase(baseLabel, byShape.values().iterator().next().alg(), canonicalLeaf);
				if (cubicSymmetric) nb = nb.withCubicSymmetric();
				pool.add(nb);
			} else {
				for (SymmetryTransforms.S3Variant v : byShape.values()) {
					NonCubicBilinearAlgorithm a = v.alg();
					String suffix = "<" + a.n + "," + a.m + "," + a.p + ">=" + a.r;
					// EXACT axis-relabel (not size-based axisLetter, ambiguous for equal axes)
					// so replay reproduces this precise orientation.
					eu.solven.matmul.catalog.Lineage.Node origin = "ABC->ABC".equals(v.perm())
							? canonicalLeaf
							: new eu.solven.matmul.catalog.Lineage.Transpose(canonicalLeaf, v.perm());
					NamedBase nb = new NamedBase(baseLabel + " / " + suffix, a, origin);
					if (cubicSymmetric) nb = nb.withCubicSymmetric();
					pool.add(nb);
				}
			}
		} catch (Exception ignored) {
			// silent: best-effort
		}
	}

	/**
	 * Strategy-agnostic wrapper: either a multi-base recombination
	 * candidate, or an additive concat split — whichever has the
	 * smaller predicted rank. Returned by {@link #findBestStrategy}.
	 */
	public record NonCubicStrategy(
			String label,
			long rank,
			MultiBaseSplitCandidate recombination,
			ConcatSplitSearch.ConcatSplit concat,
			KroneckerSplitSearch.KroneckerSplit kronecker,
			PairFusedRecombination.Prediction pairFused,
			KnownTauIdentities.Identity tauIdentity) {
		public static NonCubicStrategy ofRecombination(MultiBaseSplitCandidate c) {
			return new NonCubicStrategy(
					"recombination[" + c.baseLabel() + "]", c.rank(), c, null, null, null, null);
		}
		public static NonCubicStrategy ofConcat(ConcatSplitSearch.ConcatSplit c) {
			String axisName = c.axis() == 2 ? "p" : "n";
			return new NonCubicStrategy(
					"concat-" + axisName + "[" + c.leftSize() + "," + c.rightSize() + "]",
					c.totalRank(), null, c, null, null, null);
		}
		public static NonCubicStrategy ofKronecker(KroneckerSplitSearch.KroneckerSplit k) {
			return new NonCubicStrategy(
					"kronecker[⟨" + k.n1() + "," + k.m1() + "," + k.p1() + "⟩=" + k.r1()
							+ " ⊗ ⟨" + k.n2() + "," + k.m2() + "," + k.p2() + "⟩=" + k.r2() + "]",
					k.totalRank(), null, null, k, null, null);
		}
		public static NonCubicStrategy ofPairFused(PairFusedRecombination.Prediction p, String baseLabel) {
			return new NonCubicStrategy(
					"pair-fused[" + baseLabel + "; " + p.pairs() + "×Pan + " + p.solos() + "×solo ⟨"
							+ p.subK() + "⟩³]",
					p.totalRank(), null, null, null, p, null);
		}
		public static NonCubicStrategy ofTauIdentity(KnownTauIdentities.Identity id) {
			return new NonCubicStrategy(
					"τ-identity[" + id.attribution() + "]",
					id.predictedRank(), null, null, null, null, id);
		}
		/**
		 * Generic ConstructiveMethod wrapper. After the
		 * {@link MethodCatalog} refactor, this is the preferred factory:
		 * every {@link ConstructiveMethod.Prediction} becomes a candidate
		 * with its method-supplied label and lineage.
		 */
		public static NonCubicStrategy ofMethodPrediction(ConstructiveMethod.Prediction p) {
			return new NonCubicStrategy(
					p.label(), p.predictedRank(), null, null, null, null, null);
		}
	}

	/**
	 * Combined search: try {@link #findBestMultiBaseSplit} AND
	 * {@link ConcatSplitSearch#findBest}, return whichever predicts the
	 * smaller rank.
	 *
	 * <p>For "narrow" shapes (any axis = 2) concat usually wins because
	 * the additive decomposition tracks the Hopcroft-Kerr formula
	 * closely. For "thick" shapes (all axes ≥ 3) Strassen-recombination
	 * usually wins because it exploits sub-cubic recursion.</p>
	 */
	public static Optional<NonCubicStrategy> findBestStrategy(
			int n, int m, int p,
			List<NamedBase> basePool, Recombination.SotaResolver sota, boolean balancedOnly) {
		return findBestStrategy(n, m, p, basePool, sota, balancedOnly, Integer.MAX_VALUE);
	}

	/**
	 * Variant of {@link #findBestStrategy(int, int, int, List, Recombination.SotaResolver, boolean)}
	 * that caps per-axis allocation imbalance to {@code maxImbalance}.
	 * See {@link #findBestMultiBaseSplit(int, int, int, List, Recombination.SotaResolver, boolean, long, int)}.
	 */
	public static Optional<NonCubicStrategy> findBestStrategy(
			int n, int m, int p,
			List<NamedBase> basePool, Recombination.SotaResolver sota, boolean balancedOnly,
			int maxImbalance) {
		return findBestStrategy(n, m, p, basePool, sota, balancedOnly, maxImbalance, Integer.MAX_VALUE);
	}

	public static Optional<NonCubicStrategy> findBestStrategy(
			int n, int m, int p,
			List<NamedBase> basePool, Recombination.SotaResolver sota, boolean balancedOnly,
			int maxImbalance, int maxCombinations) {
		return findBestStrategy(n, m, p, basePool, sota, balancedOnly, maxImbalance, maxCombinations, 0);
	}

	public static Optional<NonCubicStrategy> findBestStrategy(
			int n, int m, int p,
			List<NamedBase> basePool, Recombination.SotaResolver sota, boolean balancedOnly,
			int maxImbalance, int maxCombinations, int maxPadding) {
		return findBestStrategy(n, m, p, basePool, sota, balancedOnly,
				maxImbalance, maxCombinations, maxPadding, Long.MAX_VALUE);
	}

	/**
	 * As {@link #findBestStrategy(int, int, int, List, Recombination.SotaResolver, boolean, int, int, int)},
	 * but seeds the recombination B&amp;B with an external {@code upperBoundHint} — e.g.
	 * the catalog INCUMBENT rank in improve mode. The allocation optimiser then prunes
	 * every branch {@code ≥} this bound up front, so the common "no allocation beats what
	 * we already have" case terminates fast instead of grinding {@code ALLOC_STAGNATION}
	 * nodes per base. EXACT: a recombination {@code ≥} the incumbent is never kept, so
	 * this can never drop a real improvement. Pass {@code Long.MAX_VALUE} for no hint.
	 */
	public static Optional<NonCubicStrategy> findBestStrategy(
			int n, int m, int p,
			List<NamedBase> basePool, Recombination.SotaResolver sota, boolean balancedOnly,
			int maxImbalance, int maxCombinations, int maxPadding, long upperBoundHint) {
		// Exclude same-ORDERED-shape bases: a base ⟨n,m,p⟩ can only reach the target ⟨n,m,p⟩
		// via a no-op all-1s tiling — a degenerate SELF-recombination the write-guard refuses
		// (SELF-SHAPE). Left in the pool it WINS on rank (it just re-emits the import) and
		// shadows the best STRICTLY-smaller-base derivation, so the refusal leaves the shape
		// with no derived witness at all (the ⟨4,4,4⟩=49 Strassen² / ⟨5,5,5⟩=93 gap). Drop it
		// here so the search keeps the genuine recombination. Orientation/transpose bases (same
		// multiset, DIFFERENT order — e.g. ⟨4,4,2⟩ for ⟨2,4,4⟩) are legitimate and kept.
		{
			List<NamedBase> trimmed = new java.util.ArrayList<>(basePool.size());
			for (NamedBase nb : basePool) {
				NonCubicBilinearAlgorithm b = nb.base();
				if (b.n == n && b.m == m && b.p == p) continue;
				trimmed.add(nb);
			}
			basePool = trimmed;
		}
		// Evaluate cheap candidates FIRST so the heavy recombination sweep
		// has a tight upperBound to prune against. Concat + Kronecker are
		// O(factor_pairs) — milliseconds vs the multi-minute recombination
		// grind on heavy bases (AE 5×5×5 at ⟨12,12,12⟩ explored 500k+ allocs
		// without ever beating Kronecker's ⟨2,4,4⟩=26 ⊗ ⟨6,3,3⟩=40 = 1040).
		Optional<ConcatSplitSearch.ConcatSplit> concat =
				ConcatSplitSearch.findBest(n, m, p, sota);
		Optional<KroneckerSplitSearch.KroneckerSplit> kron =
				KroneckerSplitSearch.findBest(n, m, p, sota);
		// Seed with the incumbent hint so the recombination B&B never explores branches
		// that can't beat what the catalog already has (improve mode).
		long cheapUpperBound = upperBoundHint;
		if (concat.isPresent()) cheapUpperBound = Math.min(cheapUpperBound, concat.get().totalRank());
		if (kron.isPresent()) cheapUpperBound = Math.min(cheapUpperBound, kron.get().totalRank());

		// Recombination. For the no-peel, unbalanced case use AllocationOptimizer
		// (exact B&B) — it reports the rank-minimising allocation directly, so it
		// REPLACES the legacy enumeration (no double work) and finds the unbalanced
		// optima a capped/balanced enumeration would miss (2026-06-04 balance work).
		// The legacy findBestMultiBaseSplit is kept for the two cases the optimizer
		// doesn't model: output-zero peel (maxPadding > 0) and the balancedOnly
		// restriction.
		Optional<MultiBaseSplitCandidate> recomb;
		if (!balancedOnly && maxPadding <= 0) {
			Optional<MultiBaseSplitCandidate> opt = USE_ASSIGNMENT_OPTIMIZER
					? bestViaAssignmentOptimizer(n, m, p, basePool, sota, cheapUpperBound)
					: bestViaAllocationOptimizer(n, m, p, basePool, sota, cheapUpperBound);
			// The allocation/assignment optimizer minimises the UN-PAIRED leaf-sum, so
			// it is blind to Pan trilinear-aggregation savings: a disjoint cyclic-rotation
			// cross-pair ⟨a,b,c⟩+⟨b,c,a⟩ (e.g. the two off-diagonal blocks of a symmetric
			// ⟨1,2,2⟩ peel) costs abc+ab+bc+ca instead of 2·R(⟨a,b,c⟩). TA is not a
			// separate strategy — it is a saving WITHIN decomposition that must be weighed.
			// A disjoint cyclic-rotation cross-pair (the structure TA fuses) only arises
			// from a NAÏVE-GRID base — every product a single block (the ⟨1,2,2⟩ peel, FMM's
			// ⟨2,3,3⟩ grids, …); a block-combining base like Strassen/Laderman cannot carry
			// one (cubic same-shape pairs are covered by the balanced PairFusedRecombination
			// path below). So run the pairing-aware mask sweep ONLY over the naïve-grid bases —
			// zero cost when the pool has none (the common Strassen/Laderman case) — pruned
			// against the same incumbent, keeping whichever recombination is cheaper. Bounded
			// by PAIRING_SWEEP_COMBO_BUDGET per base so a grid degrades to balanced enumeration
			// rather than exploding.
			recomb = opt;
			List<NamedBase> grids = basePool.stream()
					.filter(nb -> eu.solven.matmul.recombination.Recombination.isNaiveGrid(nb.base()))
					.toList();
			if (!grids.isEmpty()) {
				Optional<MultiBaseSplitCandidate> paired = findBestMultiBaseSplit(
						n, m, p, grids, sota, false, PAIRING_SWEEP_COMBO_BUDGET,
						maxImbalance, maxCombinations, 0, cheapUpperBound);
				if (paired.isPresent() && (recomb.isEmpty() || paired.get().rank() < recomb.get().rank())) {
					recomb = paired;
				}
			}
		} else {
			recomb = findBestMultiBaseSplit(n, m, p, basePool, sota, balancedOnly,
					Long.MAX_VALUE, maxImbalance, maxCombinations, maxPadding, cheapUpperBound);
		}

		// Collect candidates and pick min.
		List<NonCubicStrategy> candidates = new java.util.ArrayList<>();
		recomb.ifPresent(r -> candidates.add(NonCubicStrategy.ofRecombination(r)));
		concat.ifPresent(c -> candidates.add(NonCubicStrategy.ofConcat(c)));
		kron.ifPresent(k -> candidates.add(NonCubicStrategy.ofKronecker(k)));

		// Pair-fused recombination (Pan TA): cubic target ⟨n,n,n⟩ only,
		// each cubic base in the pool, balanced [k,k]³ allocation. Profitable
		// when 2·R(⟨k,k,k⟩) > k³+3k² — see MaterializeViaPanPair table.
		if (n == m && m == p) {
			for (NamedBase nb : basePool) {
				PairFusedRecombination.predictBalancedCubic(n, nb.base(), sota)
						.ifPresent(pred -> candidates.add(NonCubicStrategy.ofPairFused(pred, nb.label())));
			}
		}

		// Hand-extracted τ-identities (Sedoglavic 2017 Prop 1, FMM-Lille
		// recipes, etc.). For each identity targeting this shape (under
		// axis-permutation orbit) whose arithmetic verifies against the
		// current sota resolver, add as a candidate.
		for (KnownTauIdentities.Identity id : KnownTauIdentities.applicableTo(n, m, p, sota)) {
			candidates.add(NonCubicStrategy.ofTauIdentity(id));
		}

		// ConstructiveMethod registry (#161). Iterates every known
		// formula-driven constructor: HK71 and Pan TA today; Waksman,
		// Rosowski, etc. as adopters land. (Sedoglavic Prop 1 is NOT a
		// registered method — its bound is the Strassen ⟨2,2,2⟩ recombination
		// multiset at a [u,v] split, already covered by the recomb +
		// PairFused candidates above; see MethodCatalog.)
		// Each method's prediction (when applicable to this shape) becomes
		// a search candidate competing with Kron / Concat / recombination.
		for (ConstructiveMethod.Prediction pred : MethodCatalog.predictAll(n, m, p, sota)) {
			// Only ELECT predictions the method can actually construct. Unverified
			// (theoretical / formula-only) bounds — e.g. the HK1971 ⟨2,m,n⟩ closed
			// form for parities we have not ported — must not drive the search:
			// for small matrices we require a realisable scheme, not an unproven
			// extrapolation. (Such bounds may still be DISPLAYED as cited bounds.)
			if (!pred.verified()) continue;
			candidates.add(NonCubicStrategy.ofMethodPrediction(pred));
		}

		if (TRACE_SHAPE != null && TRACE_SHAPE.equals(n + "x" + m + "x" + p)) {
			String pick = candidates.stream().min((a, b) -> Long.compare(a.rank(), b.rank()))
					.map(s -> s.label() + "=" + s.rank()).orElse("none");
			log.info("[bss-trace {}] sota={} upperHint={} cheapUB={} | concat={} kron={} recomb={} |"
					+ " candidates=[{}] PICK={}",
					n + "x" + m + "x" + p, sota.getClass().getSimpleName(), upperBoundHint, cheapUpperBound,
					concat.map(ConcatSplitSearch.ConcatSplit::totalRank).orElse(-1L),
					kron.map(KroneckerSplitSearch.KroneckerSplit::totalRank).orElse(-1L),
					recomb.map(MultiBaseSplitCandidate::rank).orElse(-1L),
					candidates.stream().map(s -> s.label() + "=" + s.rank())
							.collect(java.util.stream.Collectors.joining(", ")),
					pick);
		}
		if (candidates.isEmpty()) return Optional.empty();
		return candidates.stream().min((a, b) -> Long.compare(a.rank(), b.rank()));
	}

	/** Per-shape diagnostic trace gate: {@code -Dbss.trace=12x13x13} logs every
	 *  {@link #findBestStrategy} candidate's rank for that shape (and the sota it ran
	 *  under), to pin predict-vs-build / sota-sensitivity divergences. Off by default. */
	private static final String TRACE_SHAPE = System.getProperty("bss.trace");

	/**
	 * Exact single-base, no-peel, non-recursive allocation optimum via
	 * {@link AllocationOptimizer} (branch-and-bound). Iterates the pool; each
	 * base's incumbent is seeded with {@code upperBound} and then with the
	 * running best across bases, so a base whose root lower bound already meets
	 * the incumbent is dropped outright. Returns the best
	 * {@code (base, optimal-allocation)} as a {@link MultiBaseSplitCandidate}
	 * (so it flows through the same recombination materialisation path), or
	 * empty if no base beats {@code upperBound}.
	 *
	 * <p>This is the exact replacement for the no-peel portion of
	 * {@link #findBestMultiBaseSplit}: it reports the rank-minimising allocation
	 * directly, so no separate allocation search is needed downstream.</p>
	 */
	/** Optional anytime budgets for the default (flat {@link AllocationOptimizer})
	 *  recombination path. {@code MAX_VALUE} = unbounded exact search (the historic
	 *  behaviour). {@code ALLOC_STAGNATION} stops a base once that many nodes pass
	 *  with no incumbent improvement — the cheap way to make large-base targets
	 *  (band ≥11) terminate fast at a near-optimal (anytime) rank instead of paying
	 *  the multi-hour exact proof. Set via {@code SchemeSweep --maxNodes/--stagnation}. */
	public static volatile long ALLOC_MAX_NODES = Long.MAX_VALUE;
	/** Default stagnation cap: stop a base after this many allocations with no
	 *  incumbent improvement. 100k is far past where the balance-first incumbent
	 *  is found (empirically &lt; 5k nodes), so it keeps the proven-good answer while
	 *  making even large-base targets terminate in seconds instead of hours.
	 *  Override with {@code SchemeSweep --stagnation} (set to a huge value for an
	 *  exact, unbounded proof). */
	public static final long DEFAULT_ALLOC_STAGNATION = 100_000L;
	public static volatile long ALLOC_STAGNATION = DEFAULT_ALLOC_STAGNATION;

	public static Optional<MultiBaseSplitCandidate> bestViaAllocationOptimizer(
			int n, int m, int p, List<NamedBase> basePool,
			Recombination.SotaResolver sota, long upperBound) {
		MultiBaseSplitCandidate best = null;
		long bestRank = upperBound;
		for (NamedBase nb : basePool) {
			NonCubicBilinearAlgorithm b = nb.base();
			if (b.n > n || b.m > m || b.p > p) continue;
			AllocationOptimizer.Result r = AllocationOptimizer.optimize(
					b, sota, n, m, p, new SearchBudget(bestRank, ALLOC_MAX_NODES, ALLOC_STAGNATION), null);
			if (r.rank() < bestRank) {
				bestRank = r.rank();
				best = new MultiBaseSplitCandidate(n, m, p, b, nb.label(),
						r.allocA(), r.allocB(), r.allocC(), r.rank(), nb.originLineage());
			}
		}
		return Optional.ofNullable(best);
	}

	/** Opt-in: route recombination through the partition+assignment exact B&B
	 *  ({@link AssignmentOptimizer}) instead of the flat {@link AllocationOptimizer}.
	 *  Set via {@code SchemeSweep --optimizer=assignment}. Default off because the
	 *  exact assignment search is slow for large bases (e.g. AlphaEvolve ⟨5,5,5⟩ at
	 *  target ≥12) even with a cross-base bound; the node budget below caps it so a
	 *  sweep can't hang (it then returns best-found, anytime). */
	public static volatile boolean USE_ASSIGNMENT_OPTIMIZER = false;

	/** Per-base allocation-combination budget for the pairing-aware mask sweep that
	 *  {@link #findBestStrategy} runs alongside the (pairing-blind) optimizer to catch
	 *  Pan trilinear-aggregation savings. A base whose full unbalanced enumeration
	 *  exceeds this degrades to balanced-only allocations there (its un-paired optimum
	 *  is still covered by the optimizer), bounding the broad sweep's added cost. */
	public static final long PAIRING_SWEEP_COMBO_BUDGET = 200_000L;
	/** Per-base node budget for the assignment optimizer (safety against the slow
	 *  large-base case). Beyond this it returns best-found rather than the proven
	 *  optimum. */
	public static volatile long ASSIGNMENT_MAX_NODES = SearchHeuristics.DEFAULT_ASSIGNMENT_MAX_NODES;

	/** Mirror of {@link #bestViaAllocationOptimizer} using {@link AssignmentOptimizer}
	 *  (partition + exact arrangement B&B), pruned by the cross-base incumbent and
	 *  capped by {@link #ASSIGNMENT_MAX_NODES}. */
	public static Optional<MultiBaseSplitCandidate> bestViaAssignmentOptimizer(
			int n, int m, int p, List<NamedBase> basePool,
			Recombination.SotaResolver sota, long upperBound) {
		MultiBaseSplitCandidate best = null;
		long bestRank = upperBound;
		for (NamedBase nb : basePool) {
			NonCubicBilinearAlgorithm b = nb.base();
			if (b.n > n || b.m > m || b.p > p) continue;
			AssignmentOptimizer.Result r =
					AssignmentOptimizer.optimize(b, sota, n, m, p, new SearchBudget(bestRank, ASSIGNMENT_MAX_NODES, Long.MAX_VALUE));
			if (r.rank() < bestRank) {
				bestRank = r.rank();
				best = new MultiBaseSplitCandidate(n, m, p, b, nb.label(),
						r.allocA(), r.allocB(), r.allocC(), r.rank(), nb.originLineage());
			}
		}
		return Optional.ofNullable(best);
	}

	/**
	 * Multi-base block-split search. For target {@code ⟨n,m,p⟩}, iterate
	 * over a pool of outer bases (Strassen, Laderman, Hopcroft-Kerr,
	 * trivial axis-splitters, etc.) and for each, iterate all
	 * non-degenerate allocations summing to the target. Pick the
	 * (base, allocation) tuple with the smallest total rank.
	 *
	 * <p>For larger targets the allocation enumeration grows
	 * combinatorially. Pass {@code balancedOnly = true} to restrict to
	 * "balanced" splits per DIS09's heuristic for n ≥ 7 — distribute
	 * the budget into nearly-equal parts.</p>
	 */
	public static Optional<MultiBaseSplitCandidate> findBestMultiBaseSplit(int n, int m, int p,
			List<NamedBase> basePool, Recombination.SotaResolver sota, boolean balancedOnly) {
		return findBestMultiBaseSplit(n, m, p, basePool, sota, balancedOnly, Long.MAX_VALUE);
	}

	/**
	 * Variant of {@link #findBestMultiBaseSplit} that auto-switches to
	 * balanced enumeration for any base whose unbalanced combo count
	 * exceeds {@code maxCombosPerBase}. Lets you enable full enumeration
	 * for cheap bases (Strassen, axis-split: 2 blocks/axis) while keeping
	 * Laderman (3 blocks/axis, C(n-1,2)³ combos at large n) on the
	 * balanced heuristic so the search stays feasible.
	 */
	public static Optional<MultiBaseSplitCandidate> findBestMultiBaseSplit(int n, int m, int p,
			List<NamedBase> basePool, Recombination.SotaResolver sota,
			boolean balancedOnly, long maxCombosPerBase) {
		return findBestMultiBaseSplit(n, m, p, basePool, sota, balancedOnly, maxCombosPerBase, Integer.MAX_VALUE);
	}

	/**
	 * Variant that additionally accepts {@code maxImbalance} — drop any
	 * per-axis allocation whose {@code max(parts) − min(parts) > maxImbalance}.
	 * Set to {@code Integer.MAX_VALUE} for "no cap" (fully unbalanced).
	 * Set to {@code 1} for the historical "balanced" behaviour.
	 * Set to {@code 0} for "all parts equal" (only legal when budget % blocks == 0).
	 *
	 * <p>Composes with {@code balancedOnly} — if {@code balancedOnly=true},
	 * {@code maxImbalance} is irrelevant (we already only enumerate the
	 * balanced multiset's permutations).</p>
	 */
	public static Optional<MultiBaseSplitCandidate> findBestMultiBaseSplit(int n, int m, int p,
			List<NamedBase> basePool, Recombination.SotaResolver sota,
			boolean balancedOnly, long maxCombosPerBase, int maxImbalance) {
		return findBestMultiBaseSplit(n, m, p, basePool, sota, balancedOnly, maxCombosPerBase,
				maxImbalance, Integer.MAX_VALUE);
	}

	/**
	 * Variant that caps the total Cartesian-product size per base to
	 * {@code maxCombinations}. When the cap fires, allocations are
	 * sorted by per-axis imbalance ascending and the top-K from each
	 * axis are taken so the resulting product is bounded; the chosen
	 * tuples are themselves the {@code maxCombinations} smallest by total
	 * imbalance. Set to {@code Integer.MAX_VALUE} for "no cap".
	 *
	 * <p>This is the right knob for keeping unbalanced search tractable
	 * on heavy bases (Pan ⟨3,4,7⟩=63 produces 538M combos at ⟨17,17,17⟩
	 * with unbalanced enumeration — totally infeasible). Setting
	 * {@code maxCombinations=16} probes just the 16 most-balanced
	 * tuples per base, finishing in milliseconds.</p>
	 */
	public static Optional<MultiBaseSplitCandidate> findBestMultiBaseSplit(int n, int m, int p,
			List<NamedBase> basePool, Recombination.SotaResolver sota,
			boolean balancedOnly, long maxCombosPerBase, int maxImbalance,
			int maxCombinations) {
		return findBestMultiBaseSplit(n, m, p, basePool, sota, balancedOnly, maxCombosPerBase,
				maxImbalance, maxCombinations, 0);
	}

	/**
	 * Variant that additionally enumerates over-allocations with output-side
	 * peel ({@code maxPadding} ≥ 0). For each axis we additionally try
	 * over-allocations summing to {@code target + 1 .. target + maxPadding};
	 * the excess is peeled off the LAST block via the
	 * {@link Recombination#recombineWithAllocation(NonCubicBilinearAlgorithm,
	 * Recombination.SotaResolver, int[], int[], int[], int[], int[], int[])}
	 * peel-aware overload.
	 *
	 * <p>This captures DIS09's γ5 recipe: at ⟨3,3,3⟩ via Strassen ⟨2,2,2⟩=7,
	 * over-allocation (2,2)³ + peel (0,1)³ collapses 4 of the 7 sub-products
	 * to smaller shapes, dropping the total from 49 to 25 — the exact bound
	 * cited in DIS09 §6. Similarly at ⟨17,17,17⟩ → padded ⟨18,18,18⟩.</p>
	 *
	 * <p>Search-space impact: multiplies enumeration by up to
	 * {@code (maxPadding+1)^3} per base. Setting {@code maxPadding=base.n−1}
	 * is usually enough to capture the dominant peel pattern.</p>
	 */
	public static Optional<MultiBaseSplitCandidate> findBestMultiBaseSplit(int n, int m, int p,
			List<NamedBase> basePool, Recombination.SotaResolver sota,
			boolean balancedOnly, long maxCombosPerBase, int maxImbalance,
			int maxCombinations, int maxPadding) {
		return findBestMultiBaseSplit(n, m, p, basePool, sota,
				balancedOnly, maxCombosPerBase, maxImbalance, maxCombinations, maxPadding,
				Long.MAX_VALUE);
	}

	/**
	 * Bounded variant: seed {@code bestRank} with an externally-known
	 * {@code upperBound} (e.g.\ from cheap Kronecker/Concat strategies)
	 * and prune bases whose trivial lower bound ({@code base.r}, since
	 * every sub-product costs ≥1 and pairing never reduces below this
	 * for non-pair-fusion-amenable shapes) already meets/exceeds the
	 * upper bound. Returns {@code Optional.empty()} when no recombination
	 * candidate strictly beats {@code upperBound}.
	 *
	 * <p>Sound: every prune below is provably correct under the strict
	 * SAT-style discipline — no heuristic slack. Pairing-savings analysis
	 * would let us prune inside the sub-product loop too, but only the
	 * trivial LB is unconditionally safe across pairing/non-pairing
	 * regimes.</p>
	 */
	public static Optional<MultiBaseSplitCandidate> findBestMultiBaseSplit(int n, int m, int p,
			List<NamedBase> basePool, Recombination.SotaResolver sota,
			boolean balancedOnly, long maxCombosPerBase, int maxImbalance,
			int maxCombinations, int maxPadding, long upperBound) {
		long bestRank = upperBound;
		MultiBaseSplitCandidate best = null;
		// Intra-search progress (5-second cadence). Gated on max(n,m,p)≥12
		// so cheap shapes stay silent.
		final boolean progressEnabled = Math.max(n, Math.max(m, p)) >= 12;
		final long searchStartNs = System.nanoTime();
		long lastEmitNs = searchStartNs;
		long allocsExplored = 0;
		long allocsExploredSinceEmit = 0;
		final long emitEveryNs = 5_000_000_000L;
		for (int bi = 0; bi < basePool.size(); bi++) {
			NamedBase nb = basePool.get(bi);
			NonCubicBilinearAlgorithm base = nb.base;
			// Sound per-base lower-bound prune. Every base of rank r
			// produces r sub-products, each costing ≥1 (sota ≥ 1 for any
			// non-degenerate shape). applyPairing can only pair WHEN
			// savings > 0, and pair_cost(a,b,c) = abc + ab + bc + ca ≥ 4
			// for ⟨1,1,1⟩, so trivially-shaped sub-products are never
			// reduced below the unpaired sum. Hence total ≥ r is a sound
			// lower bound under the strict (SAT-style) discipline. Skip
			// any base whose r already meets/exceeds the running best.
			if (base.r >= bestRank) {
				if (progressEnabled) {
					long elapsedMs = (System.nanoTime() - searchStartNs) / 1_000_000L;
					System.err.printf("    [search ⟨%d,%d,%d⟩ +%dms] base %d/%d = %s "
									+ "SKIPPED (r=%d ≥ bestRank=%d)%n",
							n, m, p, elapsedMs, bi + 1, basePool.size(), nb.label(),
							base.r, bestRank);
				}
				continue;
			}
			if (progressEnabled) {
				long elapsedMs = (System.nanoTime() - searchStartNs) / 1_000_000L;
				System.err.printf("    [search ⟨%d,%d,%d⟩ +%dms] base %d/%d = %s "
								+ "(allocs explored: %,d, best: %s)%n",
						n, m, p, elapsedMs, bi + 1, basePool.size(), nb.label(),
						allocsExplored,
						bestRank == Long.MAX_VALUE ? "—" : String.valueOf(bestRank));
				lastEmitNs = System.nanoTime();
			}
			boolean useBalanced = balancedOnly;
			if (!useBalanced) {
				long combos = (long) compositions(base.n, n)
						* compositions(base.m, m)
						* compositions(base.p, p);
				if (combos > maxCombosPerBase) useBalanced = true;
			}
			// Per-base cached block-supports — reused across every
			// per-allocation recombine call below via
			// recombineWithAllocationFast. Drops per-call cost from
			// O(r·dim²) factor-matrix scan to O(r·dim) bitmask lookup.
			AnalyticalMaskSearch.SchemeSupports baseSupports =
					AnalyticalMaskSearch.SchemeSupports.extract(base);
			// Per-base shape-multiset signature dedup. Two allocations
			// that produce the same multiset of sub-shapes give the same
			// cost; the second one wastes recombine + applyPairing work.
			// We compute a 64-bit FNV-1a hash over the sorted multiset
			// (cheap: r long packs + sort + mix) and skip pair-fusion
			// when a duplicate is detected. Cap the set to avoid OOM.
			final int SIGNATURE_CAP = 5_000_000;
			java.util.Set<Long> seenSignatures = new java.util.HashSet<>();
			long[] sigBuf = new long[base.r];
			long dedupSkipped = 0;
			// Outer loop: over-allocation per axis. dN=dM=dP=0 is the original
			// no-peel allocation (back-compat); positive d enables DIS09 γ5.
			for (int dN = 0; dN <= maxPadding; dN++) {
				for (int dM = 0; dM <= maxPadding; dM++) {
					for (int dP = 0; dP <= maxPadding; dP++) {
						int padN = n + dN, padM = m + dM, padP = p + dP;
						if (padN < base.n || padM < base.m || padP < base.p) continue;
						List<int[]> allocsA = useBalanced
								? balancedSplits(base.n, padN)
								: filterByImbalance(nonDegenerateFillings(base.n, padN), maxImbalance);
						List<int[]> allocsB = useBalanced
								? balancedSplits(base.m, padM)
								: filterByImbalance(nonDegenerateFillings(base.m, padM), maxImbalance);
						List<int[]> allocsC = useBalanced
								? balancedSplits(base.p, padP)
								: filterByImbalance(nonDegenerateFillings(base.p, padP), maxImbalance);
						if (maxCombinations != Integer.MAX_VALUE
								&& (long) allocsA.size() * allocsB.size() * allocsC.size() > maxCombinations) {
							int perAxisK = (int) Math.ceil(Math.cbrt(maxCombinations)) + 1;
							allocsA = topKByImbalance(allocsA, perAxisK);
							allocsB = topKByImbalance(allocsB, perAxisK);
							allocsC = topKByImbalance(allocsC, perAxisK);
						}
						if (allocsA.isEmpty() || allocsB.isEmpty() || allocsC.isEmpty()) continue;
						// Per-axis reverse canonicalisation: skip allocations
						// that are lex-greater than their reverse. Safe iff
						// the axis-flip orbit is covered elsewhere — by the
						// mask sweep below, which runs for every base when
						// peel is null (task #133 lifted the 2×2×2 gate).
						final boolean canonicaliseAxisFlip =
								(dN == 0 && dM == 0 && dP == 0);
						// S₃ axis-permutation canonicalisation. Safe iff the
						// base is cubic-symmetric (T[a,a,a] invariant under
						// axis permutation) AND the target is cubic AND
						// peel is null. Up to 6× reduction on top of axis-flip
						// for cubic-symmetric pool entries (Strassen, Winograd).
						final boolean canonicaliseAxisPerm =
								(dN == 0 && dM == 0 && dP == 0)
										&& nb.cubicSymmetric()
										&& base.n == base.m && base.m == base.p
										&& n == m && m == p;
						for (int[] aA : allocsA) {
							// Peel = excess concentrated in LAST block. Skip when peel meets
							// OR exceeds that block (#98): aA[last]==dN leaves an EFFECTIVE
							// block of size 0 (entirely padding → degenerate sub-product).
							if (dN > 0 ? aA[aA.length - 1] <= dN : aA[aA.length - 1] < dN) continue;
							int[] peelA = (dN == 0) ? null : tailPeel(aA.length, dN);
							if (canonicaliseAxisFlip && lexGreaterThanReverse(aA)) continue;
							for (int[] aB : allocsB) {
								if (dM > 0 ? aB[aB.length - 1] <= dM : aB[aB.length - 1] < dM) continue;
								int[] peelB = (dM == 0) ? null : tailPeel(aB.length, dM);
								if (canonicaliseAxisFlip && lexGreaterThanReverse(aB)) continue;
								for (int[] aC : allocsC) {
									if (dP > 0 ? aC[aC.length - 1] <= dP : aC[aC.length - 1] < dP) continue;
									int[] peelC = (dP == 0) ? null : tailPeel(aC.length, dP);
									if (canonicaliseAxisFlip && lexGreaterThanReverse(aC)) continue;
									if (canonicaliseAxisPerm && !isS3Canonical(aA, aB, aC)) continue;
									Recombination.Result r = Recombination.recombineWithAllocationFast(
											base, baseSupports, sota, aA, aB, aC,
											peelA, peelB, peelC);
									// Skip duplicates by shape-signature (peel=null only).
									if (peelA == null && peelB == null && peelC == null
											&& seenSignatures.size() < SIGNATURE_CAP) {
										long sig = signatureHash(r.smallMatrixSizes, sigBuf);
										if (!seenSignatures.add(sig)) {
											dedupSkipped++;
											continue;
										}
									}
									long paired = PairedSubProducts.applyPairing(r.smallMatrixSizes, sota);
									long rk = Math.min(r.totalRank, paired);
									if (rk < bestRank) {
										bestRank = rk;
										best = new MultiBaseSplitCandidate(n, m, p, base, nb.label,
												aA.clone(), aB.clone(), aC.clone(), rk, nb.originLineage());
									}
									allocsExplored++;
									allocsExploredSinceEmit++;
									if (progressEnabled) {
										long nowNs = System.nanoTime();
										if (nowNs - lastEmitNs >= emitEveryNs) {
											long elapsedMs = (nowNs - searchStartNs) / 1_000_000L;
											long rate = allocsExploredSinceEmit * 1_000_000_000L
													/ Math.max(1, nowNs - lastEmitNs);
											System.err.printf("    [search ⟨%d,%d,%d⟩ +%dms] base %d/%d=%s "
															+ "alloc=%s×%s×%s "
															+ "(allocs: %,d, %,d sig-skipped, %,d/s, best: %s)%n",
													n, m, p, elapsedMs, bi + 1, basePool.size(), nb.label(),
													java.util.Arrays.toString(aA),
													java.util.Arrays.toString(aB),
													java.util.Arrays.toString(aC),
													allocsExplored, dedupSkipped, rate,
													bestRank == Long.MAX_VALUE ? "—" : String.valueOf(bestRank));
											lastEmitNs = nowNs;
											allocsExploredSinceEmit = 0;
										}
									}
									// Analytical axis-flip mask sweep — restores orbit
									// coverage when the pool wasn't expanded (e.g.\ RecombinationPoolConfig
									// with CANONICAL orbit mode). Gated on peel == null because
									// AnalyticalMaskSearch doesn't model peel; for peel > 0 the
									// pool-expansion path is the only coverage. Task #133:
									// previously gated on base.n==base.m==base.p==2; lifted so
									// axis-flip canonicalisation (above) is safe for all bases.
									// Uses cached baseSupports — per-call overhead is O(8·r·dim).
									if (peelA == null && peelB == null && peelC == null) {
										List<AnalyticalMaskSearch.MaskCandidate> maskTop =
												AnalyticalMaskSearch.topKMasks(baseSupports, aA, aB, aC, sota, 8);
										for (AnalyticalMaskSearch.MaskCandidate mc : maskTop) {
											if (mc.mask == 0) continue; // already evaluated above
											if (mc.cost < bestRank) {
												NonCubicBilinearAlgorithm masked =
														SymmetryTransforms.applyAxisFlipMask(base, mc.mask);
												eu.solven.matmul.catalog.Lineage.Node origin =
														nb.originLineage();
												if (origin != null) {
													origin = new eu.solven.matmul.catalog.Lineage.AxisFlip(
															origin, mc.mask);
												}
												bestRank = mc.cost;
												best = new MultiBaseSplitCandidate(n, m, p, masked,
														nb.label + " :: mask" + mc.mask,
														aA.clone(), aB.clone(), aC.clone(),
														mc.cost, origin);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		return Optional.ofNullable(best);
	}

	/** Peel pattern: 0 for all blocks except the last, which gets {@code d}. */
	/**
	 * True iff {@code alloc} is lexicographically strictly greater than its
	 * reverse. Used to skip the reverse-image of each allocation when the
	 * axis-flip orbit of the base is covered elsewhere (e.g. by the
	 * 2×2×2 mask sweep). Palindromes return false (kept as canonical).
	 */
	private static boolean lexGreaterThanReverse(int[] alloc) {
		int n = alloc.length;
		for (int i = 0; i < n / 2; i++) {
			int j = n - 1 - i;
			if (alloc[i] != alloc[j]) return alloc[i] > alloc[j];
		}
		return false;
	}

	/**
	 * True iff {@code (aA, aB, aC)} is the lexicographically smallest of its
	 * 6 S₃ axis-permutations. Used to canonicalise the cubic axis-permutation
	 * orbit when the base is cubic-symmetric and the target is cubic
	 * (DIS09 §3: {@code T[a,a,a]} is invariant under axis-permutation).
	 */
	private static boolean isS3Canonical(int[] aA, int[] aB, int[] aC) {
		if (compareTripleLex(aA, aB, aC, aA, aC, aB) > 0) return false;
		if (compareTripleLex(aA, aB, aC, aB, aA, aC) > 0) return false;
		if (compareTripleLex(aA, aB, aC, aB, aC, aA) > 0) return false;
		if (compareTripleLex(aA, aB, aC, aC, aA, aB) > 0) return false;
		if (compareTripleLex(aA, aB, aC, aC, aB, aA) > 0) return false;
		return true;
	}

	private static int compareTripleLex(int[] a1, int[] a2, int[] a3,
			int[] b1, int[] b2, int[] b3) {
		int c = compareIntArr(a1, b1);
		if (c != 0) return c;
		c = compareIntArr(a2, b2);
		if (c != 0) return c;
		return compareIntArr(a3, b3);
	}

	private static int compareIntArr(int[] a, int[] b) {
		int n = Math.min(a.length, b.length);
		for (int i = 0; i < n; i++) {
			if (a[i] != b[i]) return Integer.compare(a[i], b[i]);
		}
		return Integer.compare(a.length, b.length);
	}

	/**
	 * Order-invariant 64-bit hash of the multiset of sub-shape triples.
	 * Sum of Murmur3-64-finalised per-element packs — commutative,
	 * associative, multiplicity-aware. O(r), no allocation, no sort.
	 * The {@code buf} parameter is kept for source compatibility but no
	 * longer used.
	 */
	private static long signatureHash(int[][] shapes, long[] buf) {
		long h = 0;
		int n = shapes.length;
		for (int i = 0; i < n; i++) {
			int[] x = shapes[i];
			long packed = ((long) x[0] << 32) | ((long) x[1] << 16) | (long) x[2];
			h += murmur3Finalizer(packed);
		}
		return h;
	}

	private static long murmur3Finalizer(long x) {
		x ^= x >>> 33;
		x *= 0xff51afd7ed558ccdL;
		x ^= x >>> 33;
		x *= 0xc4ceb9fe1a85ec53L;
		x ^= x >>> 33;
		return x;
	}

	private static int[] tailPeel(int blocks, int d) {
		int[] p = new int[blocks];
		p[blocks - 1] = d;
		return p;
	}

	/** Per-allocation imbalance = max(parts) - min(parts). */
	private static int imbalance(int[] alloc) {
		int min = alloc[0], max = alloc[0];
		for (int v : alloc) { if (v < min) min = v; if (v > max) max = v; }
		return max - min;
	}

	/**
	 * Returns the K allocations with smallest imbalance (most balanced
	 * first). Stable wrt input order on ties.
	 */
	private static List<int[]> topKByImbalance(List<int[]> allocs, int k) {
		if (allocs.size() <= k) return allocs;
		List<int[]> copy = new ArrayList<>(allocs);
		copy.sort(java.util.Comparator.comparingInt(BlockSplitSearch::imbalance));
		return new ArrayList<>(copy.subList(0, k));
	}

	/**
	 * Drop any allocation whose {@code max(parts) − min(parts) > maxImbalance}.
	 * No-op when {@code maxImbalance == Integer.MAX_VALUE}.
	 */
	private static List<int[]> filterByImbalance(List<int[]> allocs, int maxImbalance) {
		if (maxImbalance == Integer.MAX_VALUE) return allocs;
		List<int[]> out = new ArrayList<>();
		for (int[] a : allocs) {
			int min = a[0], max = a[0];
			for (int v : a) { if (v < min) min = v; if (v > max) max = v; }
			if (max - min <= maxImbalance) out.add(a);
		}
		return out;
	}

	/** Number of compositions of {@code budget} into {@code blocks} positive parts: C(budget-1, blocks-1). */
	private static long compositions(int blocks, int budget) {
		if (blocks <= 0 || budget < blocks) return 0;
		// C(budget-1, blocks-1)
		int n = budget - 1, k = blocks - 1;
		if (k < 0 || k > n) return 0;
		if (k > n - k) k = n - k;
		long c = 1;
		for (int i = 0; i < k; i++) {
			c = c * (n - i) / (i + 1);
		}
		return c;
	}

	/** Named base for {@link #findBestMultiBaseSplit} reports. */
	/**
	 * A pool entry: the scheme to use as outer template plus a human label.
	 * The optional {@code originLineage} tells the materialiser how this
	 * variant was derived from the canonical leaf — typically a
	 * {@link eu.solven.matmul.catalog.Lineage.Atom} for canonical
	 * entries, wrapped in {@link eu.solven.matmul.catalog.Lineage.AxisFlip}
	 * or {@link eu.solven.matmul.catalog.Lineage.AxisPermute} when the
	 * pool was built with a non-{@code CANONICAL} orbit mode. Null is
	 * allowed and means "no origin recorded" — back-compat for callers
	 * built with the 2-arg constructor.
	 */
	public record NamedBase(String label, NonCubicBilinearAlgorithm base,
			eu.solven.matmul.catalog.Lineage.Node originLineage,
			boolean cubicSymmetric) {
		public NamedBase(String label, NonCubicBilinearAlgorithm base) {
			this(label, base, null, false);
		}
		public NamedBase(String label, NonCubicBilinearAlgorithm base,
				eu.solven.matmul.catalog.Lineage.Node originLineage) {
			this(label, base, originLineage, false);
		}
		/** Same base with {@code cubicSymmetric=true}. */
		public NamedBase withCubicSymmetric() {
			return new NamedBase(label, base, originLineage, true);
		}
	}

	/** All non-degenerate compositions of {@code budget} into {@code blocks} positive parts. */
	private static List<int[]> nonDegenerateFillings(int blocks, int budget) {
		List<int[]> out = new ArrayList<>();
		for (int[] alloc : Recombination.blockFillings(blocks, budget)) {
			boolean ok = true;
			for (int x : alloc) if (x == 0) { ok = false; break; }
			if (ok) out.add(alloc);
		}
		return out;
	}

	/**
	 * DIS09's balanced-split heuristic: distribute {@code budget} into
	 * {@code blocks} parts as evenly as possible, then return all distinct
	 * permutations of that composition (and only those).
	 *
	 * <p>E.g. budget=10, blocks=3 → base composition is (4,3,3); returns
	 * the 3 distinct permutations. For budget=15, blocks=3 → (5,5,5) only.</p>
	 */
	private static List<int[]> balancedSplits(int blocks, int budget) {
		if (blocks <= 0 || budget < blocks) return List.of(); // can't fill all with positive
		int small = budget / blocks;
		int remainder = budget - small * blocks;  // remainder large = small+1, rest = small
		int largeCount = remainder;
		// Compose: largeCount copies of (small+1) followed by (blocks-largeCount) copies of small
		int[] template = new int[blocks];
		for (int i = 0; i < blocks; i++) template[i] = (i < largeCount) ? (small + 1) : small;
		// All distinct permutations
		List<int[]> out = new ArrayList<>();
		permute(template, 0, out, new java.util.HashSet<>());
		return out;
	}

	private static void permute(int[] arr, int start, List<int[]> out, java.util.Set<String> seen) {
		if (start == arr.length - 1) {
			String key = java.util.Arrays.toString(arr);
			if (seen.add(key)) out.add(arr.clone());
			return;
		}
		java.util.Set<Integer> swappedValues = new java.util.HashSet<>();
		for (int i = start; i < arr.length; i++) {
			if (!swappedValues.add(arr[i])) continue;
			int tmp = arr[start]; arr[start] = arr[i]; arr[i] = tmp;
			permute(arr, start + 1, out, seen);
			tmp = arr[start]; arr[start] = arr[i]; arr[i] = tmp;
		}
	}

	/**
	 * Non-cubic generalisation: for target {@code ⟨n,m,p⟩}, enumerate all
	 * non-degenerate per-axis splits {@code (u_n+v_n, u_m+v_m, u_p+v_p)}
	 * and pick the min total rank under
	 * {@link Recombination#recombineWithAllocation} with a Strassen
	 * {@code ⟨2,2,2⟩} outer base.
	 *
	 * @param strassen the canonical Strassen ⟨2,2,2⟩=7 algorithm (caller-pinned
	 *                 because choice of specific rank-7 file matters; see
	 *                 TRILINEAR_AGGREGATION.md §3ter)
	 * @param sota     rank lookup for sub-formats — should respect field
	 *                 discipline (e.g. R-pure or F₂-pure)
	 */
	public static Optional<NonCubicSplitCandidate> findBestSplitNonCubic(int n, int m, int p,
			NonCubicBilinearAlgorithm strassen, Recombination.SotaResolver sota) {
		long bestRank = Long.MAX_VALUE;
		int[] bestAllocA = null, bestAllocB = null, bestAllocC = null;
		for (int vN = 1; vN < n; vN++) {
			int[] allocA = { n - vN, vN };
			for (int vM = 1; vM < m; vM++) {
				int[] allocB = { m - vM, vM };
				for (int vP = 1; vP < p; vP++) {
					int[] allocC = { p - vP, vP };
					Recombination.Result r = Recombination.recombineWithAllocation(strassen, sota,
							allocA, allocB, allocC);
					if (r.totalRank < bestRank) {
						bestRank = r.totalRank;
						bestAllocA = allocA.clone();
						bestAllocB = allocB.clone();
						bestAllocC = allocC.clone();
					}
				}
			}
		}
		if (bestAllocA == null) return Optional.empty();
		return Optional.of(new NonCubicSplitCandidate(n, m, p,
				bestAllocA, bestAllocB, bestAllocC, bestRank));
	}

	/**
	 * For target {@code ⟨n,n,n⟩}, try every split {@code n = u + v} with
	 * {@code u > v ≥ 1} and return the one with the smallest formula rank,
	 * or empty if no split has all three sub-formats available.
	 */
	public static Optional<SplitCandidate> findBestSplit(int n,
			Function<int[], Optional<Integer>> rankLookup) {
		SplitCandidate best = null;
		for (int v = 1; v < n; v++) {
			int u = n - v;
			if (u <= v) break;
			Optional<Integer> rUuu = rankLookup.apply(canonical(u, u, u));
			Optional<Integer> rUuv = rankLookup.apply(canonical(u, u, v));
			Optional<Integer> rUvv = rankLookup.apply(canonical(u, v, v));
			if (rUuu.isEmpty() || rUuv.isEmpty() || rUvv.isEmpty()) continue;
			int total = rUuu.get() + 3 * rUuv.get() + 3 * rUvv.get();
			if (best == null || total < best.formulaRank) {
				best = new SplitCandidate(n, u, v, total, rUuu.get(), rUuv.get(), rUvv.get());
			}
		}
		return Optional.ofNullable(best);
	}

	/**
	 * Walks the entire scheme catalog and builds a best-rank map keyed by
	 * canonical (sorted) format, mixing schemes from ALL fields. Includes
	 * trivial entries for {@code ⟨1,n,m⟩} = {@code n·m} (matrix-vector
	 * product). The result is a <strong>cross-field upper bound</strong> —
	 * not realisable as a single-field algorithm. For field-pure bounds,
	 * use {@link #loadCatalogBestRanksForField}.
	 */
	public static Map<String, Integer> loadCatalogBestRanks() {
		return loadCatalogBestRanksFiltered(name -> true);
	}

	/**
	 * Field discipline: load best-rank map restricted to a single field-class.
	 * <p>Field-classes (per {@code TRILINEAR_AGGREGATION.md} §3bis):</p>
	 * <ul>
	 * <li>{@code "R"}: includes Q, Z, R schemes (the single classical cluster) —
	 *     excludes F₂/Z₂-tagged files and excludes C-only schemes.</li>
	 * <li>{@code "C"}: R-class ∪ C-only schemes (since {@code R ⊂ C}, any
	 *     R-algorithm trivially works over C).</li>
	 * <li>{@code "F2"}: only F₂/Z₂-tagged schemes.</li>
	 * </ul>
	 */
	/** Type-safe entry point. */
	public static Map<String, Integer> loadCatalogBestRanksForField(eu.solven.matmul.algebra.Field field) {
		return switch (field) {
			// Z / Q / R: same filter for now — most schemes are Z-coefficient
			// and trivially valid over Q and R. Once we have schemes that are
			// strictly Q-only (e.g. with `1/3` coefficients) and not Z-realizable,
			// we'll split these cases.
			case Z, Q, R -> loadCatalogBestRanksFiltered(e ->
					isNonCommutative(e) && (hasField(e, "Z") || hasField(e, "Q") || hasField(e, "R")));
			case C -> loadCatalogBestRanksFiltered(e -> isNonCommutative(e) && hasField(e, "C"));
			case F2 -> loadCatalogBestRanksFiltered(e -> hasField(e, "F2"));
			case F3 -> throw new IllegalArgumentException("F3 catalog filter not yet implemented");
		};
	}

	/** @deprecated use {@link #loadCatalogBestRanksForField(eu.solven.matmul.algebra.Field)}. */
	@Deprecated
	public static Map<String, Integer> loadCatalogBestRanksForField(String fieldClass) {
		return loadCatalogBestRanksForField(eu.solven.matmul.algebra.Field.fromTag(fieldClass));
	}

	/** Shared Jackson mapper for reading the generated catalog manifest. */
	private static final tools.jackson.databind.json.JsonMapper MAPPER =
			tools.jackson.databind.json.JsonMapper.builder().build();

	/** {@code true} if the catalog entry's {@code fields[]} names the given algebra tag. */
	private static boolean hasField(tools.jackson.databind.JsonNode entry, String tag) {
		for (tools.jackson.databind.JsonNode f : entry.path("fields")) {
			if (tag.equals(f.asString())) {
				return true;
			}
		}
		return false;
	}

	/**
	 * {@code true} unless the entry is commutative-only (Waksman, Rosowski,
	 * Makarov 1986) — those do not lift to recursive matmul over a
	 * non-commutative ring, so they must be excluded from the R/C bound maps.
	 */
	private static boolean isNonCommutative(tools.jackson.databind.JsonNode entry) {
		return !entry.path("commutative").asBoolean(false);
	}

	/**
	 * Content-driven best-rank map. Reads {@code docs/catalog.json} (the
	 * Pareto-deduped manifest — already best-per-shape), keeping entries that
	 * pass {@code entryFilter} (a field/commutativity predicate over the JSON
	 * {@code fields[]}/{@code commutative} content, never the filename). Keyed
	 * by canonical (sorted) format.
	 */
	private static Map<String, Integer> loadCatalogBestRanksFiltered(
			java.util.function.Predicate<tools.jackson.databind.JsonNode> entryFilter) {
		Map<String, Integer> best = new HashMap<>();
		Path catalog = Path.of("docs/catalog.json");
		tools.jackson.databind.JsonNode root;
		try (java.io.Reader r = Files.newBufferedReader(catalog)) {
			root = MAPPER.readTree(r);
		} catch (IOException e) {
			throw new RuntimeException("reading " + catalog, e);
		}
		for (tools.jackson.databind.JsonNode entry : root.path("schemes")) {
			if (!entryFilter.test(entry)) {
				continue;
			}
			tools.jackson.databind.JsonNode fmt = entry.path("format");
			if (!fmt.isArray() || fmt.size() != 3) {
				continue;
			}
			int rank = entry.path("rank").asInt(-1);
			if (rank < 0) {
				continue;
			}
			String key = canonicalKey(fmt.get(0).asInt(), fmt.get(1).asInt(), fmt.get(2).asInt());
			best.merge(key, rank, Integer::min);
		}
		// Trivial: ⟨1, n, m⟩ has rank n·m (matrix-vector product); same in any field.
		for (int n = 1; n <= 64; n++) {
			for (int m = n; m <= 64; m++) {
				best.putIfAbsent(canonicalKey(1, n, m), n * m);
			}
		}
		// Pan trilinear aggregation: closed-form ⟨n,n,n⟩ upper bound, applies in any field
		// (genuinely non-commutative construction). Takes the MIN against existing catalog
		// entries so we never override a stricter bound from explicit schemes.
		for (int n = 2; n <= 64; n++) {
			long pan = PanTrilinearAggregation.cubicBound(n);
			if (pan <= Integer.MAX_VALUE) {
				best.merge(canonicalKey(n, n, n), (int) pan, Integer::min);
			}
		}
		return best;
	}

	/** Builds a rank-lookup function from a pre-computed catalog map. */
	public static Function<int[], Optional<Integer>> rankLookupFromMap(Map<String, Integer> best) {
		return key -> Optional.ofNullable(best.get(canonicalKey(key[0], key[1], key[2])));
	}

	private static int[] canonical(int n, int m, int p) {
		int[] a = { n, m, p };
		Arrays.sort(a);
		return a;
	}

	private static String canonicalKey(int n, int m, int p) {
		int[] s = canonical(n, m, p);
		return s[0] + "x" + s[1] + "x" + s[2];
	}

	/**
	 * CLI: for each cubic target {@code n ∈ [4, 32]}, print the best split
	 * and the formula breakdown.
	 */
	public static void main(String[] args) {
		Map<String, Integer> ranks = loadCatalogBestRanks();
		Function<int[], Optional<Integer>> lookup = rankLookupFromMap(ranks);

		log.info(String.format("%10s | %8s | %10s | %30s | %8s%n",
				"target", "(u,v)", "formula", "breakdown", "direct"));
		log.info("-".repeat(80));
		for (int n = 4; n <= CatalogLimits.MAX_DIM; n++) {
			Optional<SplitCandidate> best = findBestSplit(n, lookup);
			Optional<Integer> direct = lookup.apply(new int[] { n, n, n });
			String tgt = "⟨" + n + "," + n + "," + n + "⟩";
			if (best.isEmpty()) {
				log.info(String.format("%10s | %8s | %10s | %30s | %8s%n",
						tgt, "—", "—", "", direct.map(String::valueOf).orElse("—")));
			} else {
				SplitCandidate c = best.get();
				log.info(String.format("%10s | %8s | %10d | %30s | %8s%n",
						tgt, c.u + "+" + c.v, c.formulaRank, c.breakdown(),
						direct.map(String::valueOf).orElse("—")));
			}
		}
	}
}
