package eu.solven.matmul.search;

import eu.solven.matmul.recombination.BlockSplitSearch;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.Compose;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.recombination.Recombination;
import eu.solven.matmul.catalog.SchemeIO;
import lombok.extern.slf4j.Slf4j;

/**
 * Materialise scheme files on disk for any {@code ⟨n,m,p⟩} reachable via
 * {@link BlockSplitSearch#findBestStrategy}. Closes the gap between
 * "the search returns a rank prediction" and "the scheme exists on
 * disk".
 *
 * <p>Pipeline for a target shape:</p>
 * <ol>
 *   <li>Direct catalog hit ({@link FieldAwareLookup#find}) → return as
 *       {@link Lineage.Atom}.</li>
 *   <li>Otherwise consult {@link BlockSplitSearch#findBestStrategy} to
 *       pick the best composition. Recursively materialise every
 *       sub-shape the strategy references (via an internal
 *       {@link RecursiveLookup} overlay that adds derived schemes to
 *       {@link Recombination.AlgorithmLookup#find}'s search space).</li>
 *   <li>Glue the recursively-materialised pieces via the appropriate
 *       {@link Compose} / {@link Recombination} primitive, attach a
 *       structural {@link Lineage.Node}.</li>
 *   <li>Verify via {@link Verifier#passesRandomMatmulSpotCheck},
 *       optionally write to disk via
 *       {@link SchemeIO#write(NonCubicBilinearAlgorithm, File, boolean, Lineage.Node)}
 *       (the write notifies {@link FieldAwareLookup#onSchemeWritten}
 *       so the index stays in sync).</li>
 * </ol>
 *
 * <p>This subsumes {@link MaterialiseGaps} (which only knew about the
 * Strassen-recombine strategy). The recursive materialiser handles
 * Kronecker, concat (both axes), and recombination; pair-fused
 * (Pan TA fusion) is left as a follow-up because it's the most
 * specialised case and rarely wins over the others.</p>
 */
@Slf4j
public final class RecursiveMaterialiser {

	private final FieldAwareLookup diskLookup;
	private final List<BlockSplitSearch.NamedBase> pool;
	private final Recombination.SotaResolver sota;
	private final Path writeRoot;
	private final boolean writeNewSchemes;
	private final boolean balancedOnly;
	/** Improve-mode: when a disk scheme already exists, still run the search and
	 *  replace it if a Verifier-passing composition is STRICTLY better. Off by
	 *  default (fill-only — the historic behaviour). */
	private final boolean improveExisting;
	/** Leaf resolution mode (user 2026-06-26). DEFAULT {@code false} = LOAD: a sub-scheme
	 *  (recombination/Kron/concat leaf) is resolved from the DISK catalog as-is, NEVER
	 *  re-derived — so a sweep is NOT implicitly recursive. The catalog is populated by
	 *  re-deriving whole BANDS bottom-up (smallest max-axis first), so a leaf is already
	 *  on disk by the time a larger shape needs it. {@code true} = OPTIMIZE: leaves are
	 *  recursively materialised (the old implicit behaviour) — an EXPLICIT opt-in, intended
	 *  for introducing a NEW base whose sub-products aren't on disk yet. A constructor
	 *  parameter (immutable per instance); {@code SchemeSweep --recursive-derive} sets it. */
	private final boolean recursiveDerive;
	/** Derive-best mode: ALWAYS materialise the best DERIVED scheme per shape and
	 *  persist it when it ties-or-beats the catalog incumbent — even if a hand-crafted
	 *  / imported direct entry already exists at the same rank. The derived (replayable)
	 *  scheme is a first-class artifact: we want the best reproducible construction for
	 *  every shape, not just where composition strictly beats an import. Unlike
	 *  improve-mode this does NOT bound the search by the incumbent (so it can find the
	 *  tie) and persists on {@code <=}, not {@code <}. */
	private final boolean deriveBest;

	/** Compose-phase strategy tokens (default: all). Gate via {@link #setStrategies}. */
	public static final String STRAT_RECOMBINATION = "recombination";
	public static final String STRAT_SERENDIPITOUS = "serendipitous";
	public static final String STRAT_PROJECTION = "projection";
	/** Plain Kronecker products — no block-split search, no bud fusion. Elects only the
	 *  {@link BlockSplitSearch.CandidateKind#KRONECKER} kind of the upward search.
	 *  Typical use is a restricted {@code --strategies=kron} sweep (e.g. A/B-ing how
	 *  much serendipitous bud fusion buys over the plain product). */
	public static final String STRAT_KRONECKER = "kronecker";
	/** Additive axis-concat splits — the {@link BlockSplitSearch.CandidateKind#CONCAT}
	 *  kind of the upward search, selectable on its own. */
	public static final String STRAT_CONCAT = "concat";
	/** Default = every operator. Each token now maps 1:1 to its candidate kind(s) —
	 *  {@code recombination} no longer implicitly elects Kronecker/concat picks, so
	 *  Kron and concat must be (and are) first-class members of the default set. */
	private java.util.Set<String> strategies = java.util.Set.of(
			STRAT_RECOMBINATION, STRAT_KRONECKER, STRAT_CONCAT,
			STRAT_SERENDIPITOUS, STRAT_PROJECTION);

	/** Restrict the compose phase to the given strategy tokens (default: all three). */
	public void setStrategies(java.util.Set<String> strategies) {
		this.strategies = java.util.Set.copyOf(strategies);
	}

	/**
	 * The {@link BlockSplitSearch.CandidateKind}s electable under a strategy-token
	 * set — the 1:1 mapping that keeps {@code --strategies} honest. The τ-identity /
	 * MethodCatalog recipes travel with the {@code recombination} token until they
	 * earn a token of their own (Pan-TA pair fusion is not a kind at all — it is a
	 * saving WITHIN recombination). Shared by {@link #materialise} and SchemeSweep's
	 * closure mode so the mapping cannot drift.
	 */
	public static java.util.EnumSet<BlockSplitSearch.CandidateKind> kindsFor(
			java.util.Set<String> strategies) {
		java.util.EnumSet<BlockSplitSearch.CandidateKind> kinds =
				java.util.EnumSet.noneOf(BlockSplitSearch.CandidateKind.class);
		if (strategies.contains(STRAT_RECOMBINATION)) {
			kinds.add(BlockSplitSearch.CandidateKind.RECOMBINATION);
			kinds.add(BlockSplitSearch.CandidateKind.METHOD);
		}
		if (strategies.contains(STRAT_KRONECKER)) {
			kinds.add(BlockSplitSearch.CandidateKind.KRONECKER);
		}
		if (strategies.contains(STRAT_CONCAT)) {
			kinds.add(BlockSplitSearch.CandidateKind.CONCAT);
		}
		return kinds;
	}

	/**
	 * Whether commutative-only schemes (Waksman/Rosowski/Makarov) may be used as
	 * sub-algorithms / projection parents. Default {@code false} = non-commutative:
	 * a commutative scheme does NOT lift to recursive matmul over a non-commutative
	 * ring, so it is invalid as a recombination/projection ingredient for NC matmul.
	 * This is the commutativity axis (orthogonal to the field axis) — transported from
	 * {@code SchemeSweep --commutative}, defaulted to NC. (A commutative-matmul sweep,
	 * e.g. scalar / Rosowski / Waksman work, sets it true.)
	 */
	private boolean allowCommutative = false;

	/** Set the commutativity mode (default false = non-commutative). */
	public void setAllowCommutative(boolean allow) {
		this.allowCommutative = allow;
	}

	/** Per-base allocation-enumeration caps for the recombination search. Default
	 *  {@code Integer.MAX_VALUE} = UNBOUNDED (optimal-within-rule-set). Bounding is
	 *  OPT-IN (--maxImbalance / --maxCombinations) — silently capping would turn the
	 *  search into a bounded heuristic, violating the don't-silently-prune rule. The
	 *  one place this matters: a ⟨2,5,5⟩-class base splits two axes into 5 parts each,
	 *  so its Cartesian allocation product is ~10⁹ tuples (the ~40 min/shape stall).
	 *  {@code --maxCombinations} caps it to the top-K most-balanced tuples — where the
	 *  near-balanced winners (e.g. [6,6,6,5]) live — making the rich pool tractable,
	 *  at the honest cost of dropping to optimal-WITHIN-SCOPE for that base. */
	private int maxImbalance = Integer.MAX_VALUE;
	private int maxCombinations = Integer.MAX_VALUE;

	/** Opt-in bound on per-axis allocation imbalance (default unbounded). */
	public void setMaxImbalance(int v) {
		this.maxImbalance = v;
	}

	/** Opt-in bound on the per-base allocation Cartesian product (default unbounded). */
	public void setMaxCombinations(int v) {
		this.maxCombinations = v;
	}

	/** Disk-only leaf resolution: the best on-disk scheme for ⟨n,m,p⟩ wrapped as a
	 *  {@link Result} (fromDisk), or empty if absent — NO search, NO recursion. This is the
	 *  default leaf resolver (LOAD mode); the catalog is filled by re-deriving bands
	 *  bottom-up so the leaf is already present. Mirrors the direct-hit branch of
	 *  {@link #materialise} without the compose fall-through. */
	private Optional<Result> diskBest(int n, int m, int p) {
		if (n == 1 || m == 1 || p == 1) {
			return Optional.of(trivialOneAxis(n, m, p));
		}
		Optional<FieldAwareLookup.WithSource> disk = diskLookup.findWithSource(n, m, p);
		if (disk.isEmpty()) {
			return Optional.empty();
		}
		Lineage.Node leaf;
		try {
			Optional<Lineage.Node> deep = SchemeIO.readLineage(disk.get().path().toFile());
			leaf = deep.isPresent() ? deep.get()
					: Lineage.atomFromFilename(disk.get().path().getFileName().toString());
		} catch (java.io.IOException e) {
			leaf = Lineage.atomFromFilename(disk.get().path().getFileName().toString());
		}
		int[] nat = shapeFromName(disk.get().path().getFileName().toString());
		if (nat != null && (nat[0] != n || nat[1] != m || nat[2] != p)) {
			leaf = Lineage.orientAs(leaf, nat[0], nat[1], nat[2], n, m, p);
		}
		return Optional.of(new Result(disk.get().alg(), leaf, true));
	}

	/** Whether the downward projection-closure phase should run for this config. */
	public boolean runsProjection() {
		return strategies.contains(STRAT_PROJECTION);
	}

	/** Schemes materialised in this run, keyed by canonical "nxmxp". */
	private final Map<String, NonCubicBilinearAlgorithm> derived = new HashMap<>();
	/** Lineage trees for materialised schemes. */
	private final Map<String, Lineage.Node> derivedLineage = new HashMap<>();
	/** Cycle guard. */
	private final java.util.Set<String> inFlight = new java.util.HashSet<>();
	/** Lazily-built replayer used to resolve stub (maxDim&gt;16) projection parents
	 *  whose actual matrices are not on disk — see {@link #resolveParent}. */
	private LineageReplayer replayer;

	/** LRU cache of resolved projection parents, keyed by ORIENTED shape "NxMxP" (the
	 *  exact value {@link #resolveParent} returns). A single sweep projects MANY children
	 *  that share the same slightly-larger parent — e.g. ⟨24,27,29⟩ was replayed 7× in
	 *  one 2-minute window — and each stub replay (lineage → dense expand) costs seconds.
	 *  Memoising is correctness-neutral (same parent, no re-replay); invalidated for a
	 *  shape's orbit the moment {@link #persist} writes a better scheme there, so a
	 *  cascade win in the same pass is never masked. PER-INSTANCE (not static): the
	 *  resolved parent depends on this materialiser's {@code schemesRoot}, so a shared
	 *  cache would leak parents across roots. Access-ordered + size-capped LRU. */
	private static final int PARENT_CACHE_MAX = 48;
	/** A resolved projection/recombination parent paired with the on-disk file it
	 *  came from. The {@code path} is what lets a building-block leaf be pinned by the
	 *  SOURCE file's (always-resolvable) stamped hash rather than the in-memory
	 *  oriented/replayed hash — see {@link #durableLeafRef}. */
	private record ParentHit(NonCubicBilinearAlgorithm alg, Path path) {}

	private final Map<String, ParentHit> parentCache =
			java.util.Collections.synchronizedMap(
					new java.util.LinkedHashMap<String, ParentHit>(64, 0.75f, true) {
						@Override
						protected boolean removeEldestEntry(
								Map.Entry<String, ParentHit> eldest) {
							return size() > PARENT_CACHE_MAX;
						}
					});

	public RecursiveMaterialiser(FieldAwareLookup diskLookup,
			List<BlockSplitSearch.NamedBase> pool,
			Recombination.SotaResolver sota,
			Path writeRoot, boolean writeNewSchemes, boolean balancedOnly) {
		this(diskLookup, pool, sota, writeRoot, writeNewSchemes, balancedOnly, false);
	}

	public RecursiveMaterialiser(FieldAwareLookup diskLookup,
			List<BlockSplitSearch.NamedBase> pool,
			Recombination.SotaResolver sota,
			Path writeRoot, boolean writeNewSchemes, boolean balancedOnly, boolean improveExisting) {
		this(diskLookup, pool, sota, writeRoot, writeNewSchemes, balancedOnly, improveExisting, false);
	}

	public RecursiveMaterialiser(FieldAwareLookup diskLookup,
			List<BlockSplitSearch.NamedBase> pool,
			Recombination.SotaResolver sota,
			Path writeRoot, boolean writeNewSchemes, boolean balancedOnly,
			boolean improveExisting, boolean deriveBest) {
		this(diskLookup, pool, sota, writeRoot, writeNewSchemes, balancedOnly,
				improveExisting, deriveBest, false);
	}

	/** Full constructor. {@code recursiveDerive} (last) selects LOAD (false, default) vs
	 *  OPTIMIZE (true) leaf resolution — see {@link #recursiveDerive}. */
	public RecursiveMaterialiser(FieldAwareLookup diskLookup,
			List<BlockSplitSearch.NamedBase> pool,
			Recombination.SotaResolver sota,
			Path writeRoot, boolean writeNewSchemes, boolean balancedOnly,
			boolean improveExisting, boolean deriveBest, boolean recursiveDerive) {
		this.diskLookup = diskLookup;
		this.pool = pool;
		this.sota = sota;
		this.writeRoot = writeRoot;
		this.writeNewSchemes = writeNewSchemes;
		this.balancedOnly = balancedOnly;
		this.improveExisting = improveExisting;
		this.deriveBest = deriveBest;
		this.recursiveDerive = recursiveDerive;
	}

	/** When true (with {@link #deriveBest}), persist the best DERIVED scheme for a
	 *  shape EVEN WHEN an imported/hand-crafted atom beats it — so every shape gets a
	 *  reproducible derivation on disk for derived-vs-atom comparison (e.g. ⟨4,4,4⟩=49
	 *  Strassen⊗Strassen alongside the imported 47). Off by default: the normal
	 *  derive-best gate only keeps ties-or-better to avoid dominated clutter. */
	private boolean registerDerivedAnyway = false;

	/**
	 * Field-aware exactness spot-check. For an F₂/F₃ sweep a scheme must be verified
	 * over GF(p) — the char-0 (random-real) check would WRONGLY reject a valid F_p-only
	 * scheme (e.g. AlphaTensor ⟨4,4,4⟩=47, which computes matmul only mod 2). For Z/Q/R/C
	 * the char-0 spot-check is correct (and an integer scheme that happens to be F_p-valid
	 * also passes char-0). Returns false if the scheme isn't even representable mod p.
	 */
	private boolean verifies(NonCubicBilinearAlgorithm alg) {
		eu.solven.matmul.algebra.Field f = diskLookup.field();
		if (f == eu.solven.matmul.algebra.Field.F2) return Verifier.passesRandomMatmulSpotCheckFp(alg, 2);
		if (f == eu.solven.matmul.algebra.Field.F3) return Verifier.passesRandomMatmulSpotCheckFp(alg, 3);
		return Verifier.passesRandomMatmulSpotCheck(alg);
	}

	public RecursiveMaterialiser registerDerivedAnyway(boolean v) {
		this.registerDerivedAnyway = v;
		return this;
	}

	/** True iff the shape already has a scheme WE produced (derived / constructed /
	 *  curated) on disk — i.e. it is NOT atom-only. Used to gate the
	 *  {@link #registerDerivedAnyway} fill so we never write a derivation that is
	 *  redundant with an existing derived scheme. */
	private boolean hasOursScheme(int n, int m, int p) {
		for (Path path : diskLookup.findFiles(n, m, p)) {
			String s = path.toString();
			if (s.contains("/derived/") || s.contains("/constructed/") || s.contains("/curated/")) {
				return true;
			}
		}
		return false;
	}

	/** True iff a reproducible <em>recombination</em>-derived scheme already exists for this
	 *  shape (the {@code derived/} folder only). Unlike {@link #hasOursScheme}, a
	 *  {@code constructed/} entry (e.g. a COMMUTATIVE Waksman scheme that does NOT lift to
	 *  non-commutative matmul) does NOT count — so {@code --derive-all} still writes a genuine
	 *  NC derived witness for every shape that lacks one ("a derived for any shape"). */
	private boolean hasDerivedScheme(int n, int m, int p) {
		for (Path path : diskLookup.findFiles(n, m, p)) {
			if (path.toString().contains("/derived/")) return true;
		}
		return false;
	}

	public record Result(NonCubicBilinearAlgorithm alg, Lineage.Node lineage, boolean fromDisk) {}

	public Optional<Result> materialise(int n, int m, int p) {
		String key = canon(n, m, p);

		// (0) Trivial axis: ⟨1,m,p⟩ is just m·p naïve scalar multiplications.
		// Same for any-axis-1; build the naïve scheme directly so concat strategies
		// that split off a width-1 axis can compose.
		if (n == 1 || m == 1 || p == 1) {
			return Optional.of(trivialOneAxis(n, m, p));
		}

		// (1) Direct catalog hit. Reach for the source file too so we can
		// pick up the file's own `lineage` field if present (avoids the
		// shallow-leaf bug where composed leaves lose their deep history).
		Optional<FieldAwareLookup.WithSource> disk = diskLookup.findWithSource(n, m, p);
		Result diskResult = null;
		if (disk.isPresent()) {
			Lineage.Node leaf;
			try {
				Optional<Lineage.Node> deepLineage =
						SchemeIO.readLineage(disk.get().path().toFile());
				if (deepLineage.isPresent()) {
					leaf = deepLineage.get();
				} else {
					// File has no lineage field — use its canonical key
					// {NxMxP}_m{R}_a{A} as the Atom ref (drops the source
					// prefix so the lineage stays stable when the catalog
					// later re-tags the file's attribution). See
					// Lineage.canonicalKey.
					leaf = Lineage.atomFromFilename(disk.get().path().getFileName().toString());
				}
			} catch (java.io.IOException e) {
				// Fall back to filename-as-leaf if read fails.
				leaf = Lineage.atomFromFilename(disk.get().path().getFileName().toString());
			}
			// findWithSource returns the alg oriented to ⟨n,m,p⟩, but `leaf` (the
			// file's own lineage / canonical-key atom) replays to the file's NATIVE
			// orientation. Record the reorientation explicitly so replay reproduces
			// the exact ⟨n,m,p⟩ alg — the replayable-lineage invariant.
			int[] nat = shapeFromName(disk.get().path().getFileName().toString());
			if (nat != null && (nat[0] != n || nat[1] != m || nat[2] != p)) {
				leaf = Lineage.orientAs(leaf, nat[0], nat[1], nat[2], n, m, p);
			}
			diskResult = new Result(disk.get().alg(), leaf, true);
			// Fill-mode (default): the disk scheme is the answer. Derive-best still
			// composes — we want the best DERIVED scheme even when a direct entry exists.
			if (!improveExisting && !deriveBest) return Optional.of(diskResult);
			// Improve / derive-best: fall through to compose.
		}

		// (2) Already-derived this run (memoised).
		NonCubicBilinearAlgorithm cached = derived.get(key);
		if (cached != null) {
			Optional<NonCubicBilinearAlgorithm> oriented = cached.orientAs(n, m, p);
			if (oriented.isPresent()) {
				Lineage.Node ln = derivedLineage.get(key);
				if (cached.n != n || cached.m != m || cached.p != p) {
					ln = Lineage.orientAs(ln, cached.n, cached.m, cached.p, n, m, p);
				}
				return Optional.of(new Result(oriented.get(), ln, false));
			}
		}

		// (3) Compose via the search (Kronecker / concat / recombination).
		Optional<Result> composed = compose(n, m, p);

		// Best rank ALREADY known in the committed catalog for this shape —
		// including lineage-only STUBS, which findWithSource (above) skips because
		// they carry no matrices, yet whose rank lives in the index. This is the
		// authoritative "do we already have something at least this good?" gate:
		// gating on diskResult.alg.r alone is WRONG when an expandable-but-worse
		// sibling (e.g. dis09 ⟨18,18,18⟩=3306) shadows a better stub (=3200) — the
		// sweep would re-derive 3200 as a phantom "win" over 3306. (2026-06-06
		// serendipity-band incident.) findRank reads only the committed index, so
		// in-run staging writes stay invisible until promoted — by design.
		int bestKnownRank = (improveExisting || deriveBest)
				? diskLookup.findRank(n, m, p)
				: Recombination.SotaResolver.UNKNOWN_RANK;

		// Derive-best: persist the best DERIVED scheme whenever it TIES-OR-BEATS the
		// catalog incumbent (<=, not <) — a hand-crafted/imported equal is no reason to
		// skip; we want a reproducible derivation on disk for every shape. If the
		// composition is strictly worse (or empty), keep the existing disk entry and
		// don't write dominated clutter.
		if (deriveBest) {
			// Persist when the derivation ties-or-beats the incumbent (the normal
			// derive-best gate) OR — under registerDerivedAnyway — when the shape is
			// ATOM-ONLY (no derived/constructed/curated scheme yet), so every shape
			// gets a reproducible derived witness for derived-vs-atom comparison
			// (e.g. ⟨4,4,4⟩=49 next to the imported 47). We deliberately do NOT write a
			// derivation that is dominated by an EXISTING derived scheme — that is pure
			// clutter (e.g. a ⟨16,16,16⟩=2304 recursion when a derived 2209 is on disk).
			boolean tieOrBetter = composed.isPresent() && composed.get().alg.r <= bestKnownRank;
			// Fill a reproducible NC derived witness for any shape that lacks one in derived/
			// — even if a (possibly commutative) constructed/imported scheme already covers it
			// ("a derived for any shape", user 2026-06-28). Uses hasDerivedScheme, not
			// hasOursScheme, so a commutative Waksman constructed/ entry doesn't suppress it.
			boolean fillAtomOnly = registerDerivedAnyway && composed.isPresent() && !hasDerivedScheme(n, m, p);
			if (tieOrBetter || fillAtomOnly) {
				persist(n, m, p, composed.get());
				return composed;
			}
			return diskResult != null ? Optional.of(diskResult) : composed;
		}

		if (diskResult != null) {
			// Improve-mode: replace only if the composition is strictly better than
			// EVERYTHING already catalogued (stubs included) AND verified (compose()
			// already spot-checked it).
			if (composed.isPresent() && composed.get().alg.r < bestKnownRank) {
				persist(n, m, p, composed.get());
				return composed;
			}
			return Optional.of(diskResult);
		}
		if (composed.isEmpty()) return Optional.empty();
		// No expandable disk hit — but a stub may still cover this shape (its rank
		// is in bestKnownRank). Skip re-deriving a tie or worse. (Fill-mode keeps
		// MAX so it falls through and materialises the stub into explicit matrices.)
		if (bestKnownRank <= composed.get().alg.r) {
			return Optional.empty();
		}
		persist(n, m, p, composed.get());
		return composed;
	}

	/** Run the search at {@code ⟨n,m,p⟩} and build+verify the best composition.
	 *  Returns the built scheme (NOT written/memoised) or empty. Sub-shapes are
	 *  resolved through {@link #materialise}, so they pick up disk atoms (and, in
	 *  improve-mode, their own improvements). Cycle-guarded. */
	private Optional<Result> compose(int n, int m, int p) {
		String key = canon(n, m, p);
		if (!inFlight.add(key)) {
			log.warn("cycle detected at ⟨{}⟩ — bailing", key);
			return Optional.empty();
		}
		try {
			Result built = null;
			// Upward composite search (findBestStrategy), electing ONLY the candidate
			// kinds mapped 1:1 from the selected strategy tokens — `recombination` no
			// longer implicitly returns Kronecker/concat picks (their cheap bounds are
			// still computed inside as B&B pruning seeds).
			java.util.EnumSet<BlockSplitSearch.CandidateKind> kinds = kindsFor(strategies);
			if (!kinds.isEmpty()) {
				// In improve mode, bound the (expensive) recombination B&B by the catalog
				// incumbent: we only keep strict improvements, so any allocation ≥ the
				// incumbent is wasted exploration. findRank sees stubs; sentinel/0 → no
				// bound (fill mode, or a shape we don't have yet). EXACT prune.
				long incumbentBound = Long.MAX_VALUE;
				if (improveExisting && !deriveBest) {
					int known = diskLookup.findRank(n, m, p);
					if (known > 0 && known < Recombination.SotaResolver.UNKNOWN_RANK) incumbentBound = known;
				}
				// derive-best: do NOT bound by the incumbent — we want to FIND the best
				// derivation even when it only ties an import (a strict-improvement bound
				// would prune the tie and we'd never materialise it).
				//
				// Retry over kinds: findBestStrategy commits to a SINGLE min-rank pick;
				// if that pick fails to build or verify, drop ITS kind and re-elect among
				// the rest. This rescues e.g. the Kron fallback when a recombination pick
				// fails its spot-check, and ends the old dead-end where an unbuildable
				// pair-fused / method pick nulled the whole upward phase.
				java.util.EnumSet<BlockSplitSearch.CandidateKind> remaining = kinds.clone();
				while (!remaining.isEmpty()) {
					Optional<BlockSplitSearch.NonCubicStrategy> picked =
							BlockSplitSearch.findBestStrategy(n, m, p, pool, sota, balancedOnly,
									maxImbalance, maxCombinations, 0, incumbentBound, remaining);
					if (picked.isEmpty()) break;
					BlockSplitSearch.NonCubicStrategy s = picked.get();
					if (s.kronecker() != null) {
						built = buildKronecker(s.kronecker());
					} else if (s.concat() != null) {
						built = buildConcat(s.concat());
					} else if (s.recombination() != null) {
						built = buildRecombination(n, m, p, s.recombination());
					}
					// pair-fused / method picks are not yet buildable here → built stays
					// null and their kind is dropped below (was: whole phase dead-ended).
					if (built != null && !verifies(built.alg)) {
						log.warn("⟨{}⟩ materialised at r={} but FAILED spot-check (strategy={})",
								key, built.alg.r, s.label());
						built = null;
					}
					if (built != null) break;
					remaining.remove(BlockSplitSearch.kindOf(s));
				}
			}

			// Serendipitous product (#159, Phase 1): bud-decompose the catalog-best
			// base of each divisor shape and multiply at best per-block rank. Competes
			// with the findBestStrategy pick; kept only if strictly better and verified.
			if (strategies.contains(STRAT_SERENDIPITOUS)) {
				long upper = built != null ? built.alg.r : Long.MAX_VALUE;
				Result serendip = trySerendipitous(n, m, p, upper);
				if (serendip != null && (built == null || serendip.alg.r < built.alg.r)) {
					built = serendip;
				}
			}

			// Projection (#159 / ROADMAP): restrict a slightly-larger parent scheme
			// down to ⟨n,m,p⟩ + DCE. This is the *downward* operator — it improves a
			// smaller shape from a better-than-divisible larger one (FMM's cube→cube
			// [[1,15],[15]] pattern). It is what makes the closure multipass.
			if (strategies.contains(STRAT_PROJECTION)) {
				long upper2 = built != null ? built.alg.r : Long.MAX_VALUE;
				Result projected = tryProjection(n, m, p, upper2);
				if (projected != null && (built == null || projected.alg.r < built.alg.r)) {
					built = projected;
				}
			}

			// Direct dis09 / Pan-TA formula cube for cubic ⟨n,n,n⟩ above the materialise
			// cap. tryProjection only projects from LARGER parents (identity is skipped), so
			// the SAME-shape buildable cube at cubicBound(n) — a DIS09Lemma4(n) atom that
			// replays exactly — is otherwise unreachable: materialise then settles for a
			// worse Project(⟨n,n,n+1⟩) that ProjectionSearch mis-prices AT cubicBound while it
			// actually replays higher (the Schwartz-Zwecher rank). That split-source result is
			// the ⟨2k,2k,2k+2⟩ phantom/divergence GENERATOR (⟨22,22,22⟩ claimed 5566, replayed
			// 5596). Offer the cube directly and PREFER it on tie (≤): same rank, but genuinely
			// buildable. verifies()+replaysConsistently gate out fields where the Q-strict cube
			// is invalid. [[project_findrank_poisoned_by_nonbuildable_stubs]]
			if (n == m && m == p && n > eu.solven.matmul.catalog.CatalogPolicy.MATERIALISE_MAX_DIM) {
				Result cube = tryDis09Cube(n);
				if (cube != null && (built == null || cube.alg.r <= built.alg.r)) {
					built = cube;
				}
				// KGP-2026 LITA cube (n≥19): strictly beats dis09 + catalog for odd n
				// and large even n.
				Result lita = tryLitaCube(n);
				if (lita != null && (built == null || lita.alg.r < built.alg.r)) {
					built = lita;
				}
			}

			if (built == null) {
				log.info("⟨{}⟩ no materialisable strategy", key);
				return Optional.empty();
			}
			return Optional.of(built);
		} finally {
			inFlight.remove(key);
		}
	}

	/**
	 * Phase-1 serendipitous attempt: for each proper divisor factorisation
	 * {@code ⟨n,m,p⟩ = ⟨n1,m1,p1⟩ ⊗ ⟨n2,m2,p2⟩}, use the catalog-best base of
	 * the first factor; {@link eu.solven.matmul.catalog.SerendipitousSearch}
	 * bud-decomposes it and keeps the cheapest verified product below
	 * {@code upper}. Lineage = {@code SerendipitousProduct(Atom(n1xm1xp1), n2,m2,p2)},
	 * which replays deterministically (shape-ref → catalog best → re-bud).
	 */
	private Result trySerendipitous(int n, int m, int p, long upper) {
		List<NonCubicBilinearAlgorithm> bases = new java.util.ArrayList<>();
		// Track each base's SOURCE FILE (by identity) so the winning base can be pinned
		// durably — by the file's stamped hash + an OrientAs, NOT contentHash(orientedAlg),
		// which is stamped on no file and replays to nothing (the 747 dangling stubs).
		java.util.IdentityHashMap<NonCubicBilinearAlgorithm, Path> basePaths = new java.util.IdentityHashMap<>();
		for (int n1 = 1; n1 <= n; n1++) {
			if (n % n1 != 0) continue;
			for (int m1 = 1; m1 <= m; m1++) {
				if (m % m1 != 0) continue;
				for (int p1 = 1; p1 <= p; p1++) {
					if (p % p1 != 0) continue;
					if (n1 * m1 * p1 == 1) continue;                 // trivial base
					if ((n / n1) * (m / m1) * (p / p1) == 1) continue; // base == target
					diskLookup.findWithSource(n1, m1, p1).ifPresent(ws -> {
						bases.add(ws.alg());
						basePaths.put(ws.alg(), ws.path());
					});
					// Serendipity feeds on BUDS, not minimal rank — so also offer EVERY
					// bud-bearing base at this shape (the rank-best one is usually bud-free,
					// e.g. Strassen). Crucially we offer ALL of them, NOT just the single
					// budScore-MAX one: budScore is a misleading proxy for the actual saving
					// σ, and the σ-paying base is often a LOWER-budScore sibling whose bud is
					// on the axis THIS factorisation enlarges. The ⟨8,9,9⟩=430 case: its
					// size-3 V-bud base has budScore 4 and was masked by a budScore-11 U-bud
					// sibling (σ_V=0 here) → the old picker lost the win at 432. bestFor
					// prices every candidate by σ over all orderings in its cheap PREDICT
					// phase and builds only the winner, so feeding all bud-bearers is ~free.
					for (BudBase b : budBasesAt(n1, m1, p1)) {
						bases.add(b.alg());
						basePaths.put(b.alg(), b.path());
					}
				}
			}
		}
		if (bases.isEmpty()) return null;
		// Stub-capable build resolver: an enlarged fusion target (e.g. ⟨4,4,20⟩=230)
		// often exists only as a lineage stub — resolveParentHit replays it (cached,
		// corrupt-over-claim-guarded), where the default findWithSource resolver
		// would throw and the candidate would be silently dropped.
		eu.solven.matmul.catalog.SerendipitousBudProduct.InnerResolver resolver = (a, b, c) -> {
			ParentHit ph = resolveParentHit(a, b, c);
			return Optional.ofNullable(ph == null ? null : ph.alg());
		};
		var hit = eu.solven.matmul.catalog.SerendipitousSearch.bestFor(
				n, m, p, bases, diskLookup, upper, resolver);
		if (hit.isEmpty()) return null;
		var h = hit.get();
		if (!verifies(h.scheme())) return null;
		// Reference the base by a DURABLE file-pin (hash + OrientAs) so replay resolves
		// the EXACT (bud-rich) base — a bare shape-ref would resolve to the rank-best
		// bud-free sibling and reconstruct a different scheme, while contentHash(oriented)
		// matched no file. If the winning base has no source file we cannot pin it
		// reproducibly → discard rather than write a dangling stub.
		Path basePath = basePaths.get(h.base());
		if (basePath == null) {
			return null;
		}
		Lineage.Node baseNode = durableLeafRef(basePath, h.base().n, h.base().m, h.base().p);
		Lineage.Node tree = new Lineage.SerendipitousProduct(baseNode, h.n2(), h.m2(), h.p2());
		return new Result(h.scheme(), tree, false);
	}

	/** A bud-rich base together with the SOURCE FILE it was loaded+oriented from — the
	 *  file is needed to pin the serendipitous base DURABLY (file-hash + OrientAs), not by
	 *  the oriented alg's content-hash, which matches NO file → a dangling lineage ref. */
	private record BudBase(NonCubicBilinearAlgorithm alg, Path path) {}

	/** Per-shape cache of ALL bud-bearing bases, for serendipitous bases. */
	private final Map<String, List<BudBase>> budBasesCache =
			new java.util.concurrent.ConcurrentHashMap<>();

	/**
	 * EVERY scheme on disk at ⟨n,m,p⟩ (oriented) that carries at least one bud
	 * ({@code budScore > 0}), with its source file. We hand ALL of them to
	 * {@link eu.solven.matmul.catalog.SerendipitousSearch#bestFor} rather than the
	 * single budScore-maximal one, because budScore COUNTS buds while the real
	 * currency σ PRICES them against the chosen inner: the σ-paying base for a given
	 * product is frequently a lower-budScore sibling whose bud sits on the axis that
	 * factorisation enlarges (⟨8,9,9⟩=430's size-3 V-bud, budScore 4, masked by a
	 * budScore-11 U-bud sibling). bestFor prices each candidate over all orderings in
	 * its cheap predict-only phase and builds only the winner, so offering every
	 * bud-bearer costs ~nothing and lets the σ-right base win. See
	 * references/BUD_STRUCTURE_THEORY.md (count vs value).
	 */
	private List<BudBase> budBasesAt(int n, int m, int p) {
		return budBasesCache.computeIfAbsent(n + "x" + m + "x" + p, k -> {
			List<BudBase> out = new java.util.ArrayList<>();
			for (Path path : diskLookup.findFiles(n, m, p)) {
				NonCubicBilinearAlgorithm a;
				try {
					a = SchemeIO.read(path.toFile());
				} catch (Exception e) {
					continue;
				}
				Optional<NonCubicBilinearAlgorithm> oriented = a.orientAs(n, m, p);
				if (oriented.isEmpty()) {
					continue;
				}
				if (budScore(oriented.get()) > 0) {
					out.add(new BudBase(oriented.get(), path));
				}
			}
			return out;
		});
	}

	/** Independent-partition bud score: Σ over U/V/W class sizes ≥ 2. */
	private static int budScore(NonCubicBilinearAlgorithm a) {
		int[][] cs = eu.solven.matmul.catalog.SerendipitousBudProduct.independentClassSizes(a);
		int s = 0;
		for (int[] c : cs) {
			for (int x : c) {
				if (x >= 2) {
					s += x;
				}
			}
		}
		return s;
	}

	/** Max indices dropped per axis when projecting a parent down to the target. */
	private static final int PROJECTION_MAX_DELTA = 1;

	/** Count of projection parents skipped by the build-time margin prune (no replay
	 *  paid). Observability only — a high count vs slow throughput confirms the prune
	 *  is doing its job. Process-wide; logged by callers via {@link #prunedParents()}. */
	private static final java.util.concurrent.atomic.AtomicLong PRUNED_PARENTS =
			new java.util.concurrent.atomic.AtomicLong();

	/** Total disk parents skipped by the projection margin prune since JVM start. */
	public static long prunedParents() { return PRUNED_PARENTS.get(); }

	/**
	 * Downward-only closure step (#159 / ROADMAP, multipass): attempt to improve
	 * {@code ⟨n,m,p⟩} by projecting a slightly-larger catalog parent down to it
	 * (restrict indices + DCE). Unlike {@link #materialise}, this does NOT re-run
	 * the upward search (Kronecker / concat / recombination / serendipitous) — it
	 * is purely the projection operator, so a projection closure pass stays cheap.
	 *
	 * <p>Returns the new scheme (already persisted) iff it is a strict improvement
	 * over the current catalog entry for {@code ⟨n,m,p⟩}; otherwise empty. The
	 * strictness comes from {@link eu.solven.matmul.catalog.ProjectionSearch},
	 * which only keeps projected ranks below the supplied upper bound.</p>
	 */
	public Optional<Result> projectInto(int n, int m, int p) {
		if (n == 1 || m == 1 || p == 1) return Optional.empty(); // trivial axis: nothing to project
		long upper = currentUpper(n, m, p);
		Result projected = tryProjection(n, m, p, upper);
		if (projected == null) return Optional.empty();
		persist(n, m, p, projected); // tryProjection guarantees projected.r < upper → strict win
		return Optional.of(projected);
	}

	/**
	 * Edge-driven projection: project ONE explicitly-named parent {@code ⟨N,M,P⟩}
	 * down to child {@code ⟨n,m,p⟩}, bypassing {@link #tryProjection}'s automatic
	 * ≤7-neighbour discovery. This is the "drive the edge to visit" operator — the
	 * caller picks the DAG edge (e.g. replicate FMM's cascade
	 * ⟨26,27,32⟩→⟨26,26,32⟩→⟨25,26,32⟩) rather than letting the gather choose. The
	 * parent shape resolves to the catalog-best/stub at {@code ⟨N,M,P⟩} (same ref
	 * convention as auto-projection, so the lineage replays identically). The
	 * per-axis drop is whatever the edge requires (not capped at the auto-sweep's
	 * small delta), since the caller chose it deliberately. Persisted only if it
	 * strictly beats the child's current best. Fail-loud: an unresolvable parent
	 * propagates (we no longer swallow).
	 *
	 * @return the persisted projection, or empty if the parent can't reach the
	 *         child / the result isn't a strict improvement / isn't replayable.
	 */
	public Optional<Result> projectEdge(int[] parent, int[] child) {
		int N = parent[0], M = parent[1], P = parent[2];
		int n = child[0], m = child[1], p = child[2];
		if (n == 1 || m == 1 || p == 1) return Optional.empty();
		if (N < n || M < m || P < p) {
			throw new IllegalArgumentException("projectEdge: parent ⟨" + N + "," + M + "," + P
					+ "⟩ is not ≥ child ⟨" + n + "," + m + "," + p + "⟩ on every axis");
		}
		ParentHit parHit = resolveParentHit(N, M, P);
		if (parHit == null) return Optional.empty();
		NonCubicBilinearAlgorithm par = parHit.alg();
		long upper = currentUpper(n, m, p);
		int delta = Math.max(N - n, Math.max(M - m, P - p));
		var hit = eu.solven.matmul.catalog.ProjectionSearch.bestFor(
				n, m, p, java.util.List.of(par), upper, Math.max(delta, PROJECTION_MAX_DELTA));
		if (hit.isEmpty()) return Optional.empty();
		var h = hit.get();
		if (!verifies(h.scheme())) return Optional.empty();
		// Pin the EXACT parent (N×M×P@hash); projection DCE is parent-sensitive, so a bare
		// "@sota" parent would make the projection non-reproducible.
		String ref = preciseParentRef(parHit.path(), N, M, P);
		Lineage.Node tree = new Lineage.Project(new Lineage.Atom(ref), h.keepN(), h.keepM(), h.keepP());
		if (!replaysConsistently(tree, h.scheme())) {
			log.warn("projectEdge ⟨{},{},{}⟩=r{} via ⟨{},{},{}⟩ produced a NON-replayable lineage"
					+ " — discarding.", n, m, p, h.scheme().r, N, M, P);
			return Optional.empty();
		}
		Result r = new Result(h.scheme(), tree, false);
		persist(n, m, p, r);
		return Optional.of(r);
	}

	/**
	 * Strict upper bound for a new ⟨n,m,p⟩ projection: the best rank we ALREADY know,
	 * so only a strict improvement is ever kept. Bounds by BOTH {@code findWithSource}
	 * (real matrices) AND {@code findRank} (which also sees stub/lineage-only entries
	 * {@code findWithSource} skips) — without the latter, projection would happily
	 * write a scheme DOMINATED by an existing stub (e.g. ⟨25,26,32⟩=11530 when =11343
	 * is on disk): pure pollution. {@code Long.MAX_VALUE} when the shape is unknown.
	 */
	private long currentUpper(int n, int m, int p) {
		Optional<FieldAwareLookup.WithSource> disk = diskLookup.findWithSource(n, m, p);
		long upper = disk.map(ws -> (long) ws.alg().r).orElse(Long.MAX_VALUE);
		int known = diskLookup.findRank(n, m, p);
		if (known > 0 && known < Recombination.SotaResolver.UNKNOWN_RANK) upper = Math.min(upper, known);
		return upper;
	}

	/** A resolved projection parent plus the lineage ref ({@code "NxMxP"} for a disk
	 *  scheme, {@code "DIS09Lemma4(n=N)"} for a PanTA cube) that replays it. */
	private record NamedParent(NonCubicBilinearAlgorithm alg, String ref) {}

	/**
	 * Parent-centric (scatter) projection sweep over a set of target {@code children}.
	 * Inverts {@link #projectInto}'s child-centric gather: instead of, for each child,
	 * re-resolving its ≤7 slightly-larger parents (loading the same big stub once per
	 * child), it builds the {@code parentShape → coveredChildren} map and resolves
	 * each parent shape EXACTLY ONCE per pass, fanning it out to every child it can
	 * reach via {@link eu.solven.matmul.catalog.ProjectionSearch#projectToMany} (one
	 * {@code Supports} build per parent). Parent shapes are processed largest-first so
	 * a child improved at a high tier is available — through the next pass — as a
	 * parent for smaller children (cascade). The build-time margin prune skips
	 * resolving (replaying) a parent shape outright when NO covered child could beat
	 * its incumbent. Each child's best projection of the pass is persisted once.
	 *
	 * @return total strict-improvement wins persisted across all passes
	 */
	public int projectScatter(List<int[]> children, int maxPasses) {
		int total = 0;
		for (int pass = 1; pass <= maxPasses; pass++) {
			int wins = scatterPass(children, pass);
			total += wins;
			log.info("=== scatter pass {}: {} win(s) (margin-pruned parent replays so far: {}) ===",
					pass, wins, PRUNED_PARENTS.get());
			if (wins == 0) break;
		}
		log.info("projectScatter done: {} win(s) over {} target shape(s); {} parent replays skipped by margin prune.",
				total, children.size(), PRUNED_PARENTS.get());
		return total;
	}

	private int scatterPass(List<int[]> children, int pass) {
		int nC = children.size();
		// 1. coverage: each non-identity parent shape ⟨N,M,P⟩ (= child + δ, δ∈[0,1]³\0)
		//    → the indices of every child it can project down to.
		// Per-child FMM gap (ours − fmm) when the target row carries it ({n,m,p,ours,fmm});
		// 0 for explicit-shape inputs (no gap → falls back to size ordering).
		int[] gap = new int[nC];
		for (int ci = 0; ci < nC; ci++) {
			int[] c = children.get(ci);
			gap[ci] = (c.length >= 5 && c[3] > 0 && c[4] > 0) ? Math.max(0, c[3] - c[4]) : 0;
		}
		Map<String, List<Integer>> coverage = new HashMap<>();
		Map<String, Integer> parentMaxGap = new HashMap<>();
		for (int ci = 0; ci < nC; ci++) {
			int[] c = children.get(ci);
			if (c[0] == 1 || c[1] == 1 || c[2] == 1) continue; // trivial axis
			for (int dn = 0; dn <= PROJECTION_MAX_DELTA; dn++)
				for (int dm = 0; dm <= PROJECTION_MAX_DELTA; dm++)
					for (int dp = 0; dp <= PROJECTION_MAX_DELTA; dp++) {
						if (dn == 0 && dm == 0 && dp == 0) continue;
						String ps = (c[0] + dn) + "x" + (c[1] + dm) + "x" + (c[2] + dp);
						coverage.computeIfAbsent(ps, k -> new ArrayList<>()).add(ci);
						parentMaxGap.merge(ps, gap[ci], Math::max);
					}
		}
		// 2. parent shapes by LARGEST covered-child gap first (user 2026-06-08: close the
		//    biggest FMM gaps first, so an interrupted sweep has already landed the most
		//    rank). Every parent of a high-gap child has maxGap ≥ that gap, so ALL of a
		//    child's parents fall in its gap tier → it is fully resolved before any
		//    smaller-gap work. Tie-break larger shape first (cascade within a gap tier).
		List<int[]> pshapes = new ArrayList<>();
		for (String ps : coverage.keySet()) pshapes.add(shapeFromName(ps));
		pshapes.sort((a, b) -> {
			int ga = parentMaxGap.getOrDefault(a[0] + "x" + a[1] + "x" + a[2], 0);
			int gb = parentMaxGap.getOrDefault(b[0] + "x" + b[1] + "x" + b[2], 0);
			if (ga != gb) return gb - ga;
			int ma = Math.max(a[0], Math.max(a[1], a[2])), mb = Math.max(b[0], Math.max(b[1], b[2]));
			if (ma != mb) return mb - ma;
			return (b[0] + b[1] + b[2]) - (a[0] + a[1] + a[2]);
		});
		// 3. per-child live upper, updated in place as wins are persisted.
		long[] upper = new long[nC];
		for (int ci = 0; ci < nC; ci++) {
			int[] c = children.get(ci);
			upper[ci] = currentUpper(c[0], c[1], c[2]);
		}
		// Wins are PERSISTED IMMEDIATELY (not staged to end-of-pass): a 30-min sweep that
		// is interrupted must keep the high-gap wins it already found, and an immediate
		// write + parentCache-invalidation lets a freshly-won shape serve as a (fresh)
		// parent for smaller children later in the SAME pass (cascade). A child improved
		// by a later, better parent is simply re-persisted; the dominated earlier file is
		// harmless (findRank takes the min; the manifest dedups by shape).
		java.util.Set<Integer> improved = new java.util.HashSet<>();
		long t0 = System.currentTimeMillis();
		int processed = 0;
		for (int[] P : pshapes) {
			processed++;
			String pkey = P[0] + "x" + P[1] + "x" + P[2];
			List<Integer> kidIdx = coverage.get(pkey);
			// 3a. shape-level margin prune: skip resolving (replaying) this parent shape
			//     iff EVERY covered child is hopeless (R − k·μ ≥ that child's live upper).
			//     A child whose μ is unknown (−1) keeps the parent (never drop on missing
			//     data). One hopeful child → resolve.
			int R = diskLookup.findRank(P[0], P[1], P[2]);
			int mu = diskLookup.projectionMarginUpperBound(P[0], P[1], P[2]);
			boolean canPrune = R > 0 && R < Recombination.SotaResolver.UNKNOWN_RANK && mu >= 0;
			boolean anyHopeful = !canPrune;
			if (canPrune) {
				for (int ci : kidIdx) {
					int[] c = children.get(ci);
					int k = (P[0] > c[0] ? 1 : 0) + (P[1] > c[1] ? 1 : 0) + (P[2] > c[2] ? 1 : 0);
					if ((long) R - (long) k * mu < upper[ci]) { anyHopeful = true; break; }
				}
			}
			if (!anyHopeful) { PRUNED_PARENTS.incrementAndGet(); continue; }
			// 3b. resolve candidate parents: the disk best ONCE (cached replay) + the
			//     structured PanTA cube for a cube shape (high, un-manifested margin).
			List<NamedParent> cands = new ArrayList<>();
			// Pin the EXACT disk parent (P@hash), not the bare pkey: projection DCE is
			// parent-sensitive, so a bare/@sota parent makes the projection non-reproducible.
			ParentHit diskHit = resolveParentHit(P[0], P[1], P[2]);
			if (diskHit != null) {
				cands.add(new NamedParent(diskHit.alg(),
						preciseParentRef(diskHit.path(), P[0], P[1], P[2])));
			}
			if (P[0] == P[1] && P[1] == P[2]) {
				NonCubicBilinearAlgorithm panta = pantaCube(P[0]);
				if (panta != null) cands.add(new NamedParent(panta, "DIS09Lemma4(n=" + P[0] + ")"));
			}
			if (cands.isEmpty()) continue;
			// 3c. project each candidate down to all covered children, Supports once.
			int[][] kidShapes = new int[kidIdx.size()][];
			long[] kidUppers = new long[kidIdx.size()];
			for (int j = 0; j < kidIdx.size(); j++) {
				int ci = kidIdx.get(j);
				int[] c = children.get(ci);
				kidShapes[j] = new int[] { c[0], c[1], c[2] };
				kidUppers[j] = upper[ci];
			}
			for (NamedParent np : cands) {
				List<Optional<eu.solven.matmul.catalog.ProjectionSearch.Hit>> hits =
						eu.solven.matmul.catalog.ProjectionSearch.projectToMany(
								np.alg(), kidShapes, kidUppers, PROJECTION_MAX_DELTA);
				for (int j = 0; j < hits.size(); j++) {
					if (hits.get(j).isEmpty()) continue;
					var h = hits.get(j).get();
					int ci = kidIdx.get(j);
					int[] c = children.get(ci);
					if (h.scheme().r >= upper[ci]) continue; // strict improvement over live best
					Lineage.Node tree = new Lineage.Project(new Lineage.Atom(np.ref()),
							h.keepN(), h.keepM(), h.keepP());
					if (!replaysConsistently(tree, h.scheme())) continue;
					persist(c[0], c[1], c[2], new Result(h.scheme(), tree, false)); // immediate, crash-safe
					upper[ci] = h.scheme().r;   // live bound for the rest of the pass
					kidUppers[j] = h.scheme().r; // tighten for the next candidate parent
					improved.add(ci);
					log.info("  [scatter win p{}] ⟨{},{},{}⟩ → {} (gap {})",
							pass, c[0], c[1], c[2], h.scheme().r, gap[ci]);
				}
			}
			if (processed % 50 == 0) {
				log.info("[progress][scatter pass {}] {}/{} parent-shapes, {} children improved, {} pruned, {}ms",
						pass, processed, pshapes.size(), improved.size(),
						PRUNED_PARENTS.get(), System.currentTimeMillis() - t0);
			}
		}
		return improved.size();
	}

	/**
	 * Projection attempt: gather catalog-best schemes for every shape
	 * {@code ⟨N,M,P⟩} with {@code n≤N≤n+δ}, {@code m≤M≤m+δ}, {@code p≤P≤p+δ}
	 * (excluding the target itself), and let
	 * {@link eu.solven.matmul.catalog.ProjectionSearch} restrict each down to
	 * {@code ⟨n,m,p⟩} (exhaustive over drop positions) + DCE, keeping the
	 * cheapest verified result below {@code upper}. Lineage =
	 * {@code Project(Atom(NxMxP), keepN, keepM, keepP)}, which replays
	 * deterministically (shape-ref → catalog best → restrict + DCE).
	 *
	 * <p>This is the downward operator that makes the closure multipass: a round
	 * that improves ⟨N,M,P⟩ enables a later round to project that gain down to
	 * ⟨n,m,p⟩.</p>
	 */
	private Result tryProjection(int n, int m, int p, long upper) {
		List<NonCubicBilinearAlgorithm> parents = new ArrayList<>();
		// Lineage ref per candidate parent, so the emitted Project node replays the
		// EXACT parent — crucial because the projection-best parent is NOT always the
		// rank-best one (a structured higher-rank cube can project lower: see the
		// projection-margin R−μ argument, paper §projmargin). A bare shape-ref would
		// resolve to the rank-best sibling and reconstruct a different scheme.
		java.util.IdentityHashMap<NonCubicBilinearAlgorithm, Lineage.Node> parentRef =
				new java.util.IdentityHashMap<>();
		for (int dn = 0; dn <= PROJECTION_MAX_DELTA; dn++)
			for (int dm = 0; dm <= PROJECTION_MAX_DELTA; dm++)
				for (int dp = 0; dp <= PROJECTION_MAX_DELTA; dp++) {
					if (dn == 0 && dm == 0 && dp == 0) continue; // identity, not a projection
					int N = n + dn, M = m + dm, P = p + dp;
					// Build-time prune: skip a disk parent that provably cannot beat the
					// incumbent — BEFORE the (expensive, for >16 stubs) replay. A parent of
					// rank R, dropping one index on each of k increased axes, projects to
					// rank ≥ R − k·μ (each dropped index DCEs ≤ μ products; rigorous for the
					// ≤1-row-per-axis drops here). If even that LB ≥ upper, no projection of
					// any same-shape parent wins, so we never pay the replay. μ=-1 (shape
					// has no manifest margin) disables the prune for that shape → replay, so
					// nothing is ever dropped on missing data. Only the disk parent is gated;
					// the structured PanTA cube (high, un-manifested margin) is still offered.
					boolean pruneDisk = false;
					if (dn <= 1 && dm <= 1 && dp <= 1) {
						int R = diskLookup.findRank(N, M, P);
						int mu = diskLookup.projectionMarginUpperBound(N, M, P);
						if (R > 0 && R < Recombination.SotaResolver.UNKNOWN_RANK && mu >= 0) {
							int k = (dn > 0 ? 1 : 0) + (dm > 0 ? 1 : 0) + (dp > 0 ? 1 : 0);
							if ((long) R - (long) k * mu >= upper) {
								pruneDisk = true;
								PRUNED_PARENTS.incrementAndGet();
							}
						}
					}
					ParentHit ph = pruneDisk ? null : resolveParentHit(N, M, P);
					NonCubicBilinearAlgorithm parent = (ph == null) ? null : ph.alg();
					if (parent != null) {
						parents.add(parent);
						// Pin the EXACT parent scheme (N×M×P@hash), NOT a bare best-at-shape
						// ref: projection DCE is sensitive to the precise parent — a different
						// (even better-rank) parent eliminates different products, so a bare
						// "@sota" parent makes the projection non-reproducible.
						// Pin by the file's NATIVE shape@hash with orientation made EXPLICIT in the
						// lineage (an exact-perm Transpose), never an ORIENTED "NxMxP@hash" the replay
						// must re-derive — for an equal-axis parent (⟨19,20,20⟩, two 20-axes) that
						// re-derivation is ambiguous → picks a worse-projecting axis → predict/build
						// divergence (⟨19,19,20⟩ 4154 vs 4237). [[project_projection_parent_orientation_not_pinned]]
						parentRef.put(parent, projectionParentNode(ph.path(), parent, N, M, P));
					}
					// Structured parent: a Pan-trilinear-aggregation cube projects far
					// better than a flat rank-best cube (high projection margin), so
					// offer it even though its rank is higher. Replays via the
					// parametric ref DIS09Lemma4(n=N) (LineageReplayer.resolveParametric).
					if (N == M && M == P) {
						NonCubicBilinearAlgorithm panta = pantaCube(N);
						if (panta != null) {
							parents.add(panta);
							parentRef.put(panta, new Lineage.Atom("DIS09Lemma4(n=" + N + ")"));
						}
					}
				}
		if (parents.isEmpty()) return null;
		var hit = eu.solven.matmul.catalog.ProjectionSearch.bestFor(
				n, m, p, parents, upper, PROJECTION_MAX_DELTA);
		if (hit.isEmpty()) return null;
		var h = hit.get();
		if (!verifies(h.scheme())) return null;
		Lineage.Node baseNode = parentRef.getOrDefault(h.parent(),
				new Lineage.Atom(h.parent().n + "x" + h.parent().m + "x" + h.parent().p));
		Lineage.Node tree = new Lineage.Project(baseNode, h.keepN(), h.keepM(), h.keepP());
		// Enforce replayability before we persist: a written scheme (especially a
		// stub above MATERIALISE_MAX_DIM) is only its lineage on disk, so the
		// emitted Project node MUST reconstruct an equivalent scheme. Reject (and
		// warn) rather than persist a scheme whose lineage we cannot replay.
		if (!replaysConsistently(tree, h.scheme())) {
			log.warn("projection ⟨{},{},{}⟩=r{} via ⟨{},{},{}⟩ produced a NON-replayable lineage"
					+ " — discarding the win to keep the catalog replay-clean.",
					n, m, p, h.scheme().r, h.parent().n, h.parent().m, h.parent().p);
			return null;
		}
		return new Result(h.scheme(), tree, false);
	}

	/**
	 * Replayability guard: replay {@code tree} and confirm it reconstructs a
	 * scheme of the right shape that still passes the matmul spot-check, with a
	 * rank that is <strong>at least as good as the evaluated rank</strong>.
	 * Used to certify a {@link Lineage.Project} node before it is persisted, so
	 * the catalog never carries a stub whose lineage cannot be replayed.
	 *
	 * <p>Rank policy (the asymmetry is load-bearing):</p>
	 * <ul>
	 *   <li>{@code replay.r < expected.r} — the rebuild is BETTER than the
	 *       evaluated rank. Legitimate: projection dead-code elimination can drop
	 *       products a flat evaluation counted. Accepted.</li>
	 *   <li>{@code replay.r == expected.r} — exact agreement. Accepted.</li>
	 *   <li>{@code replay.r > expected.r} — the rebuild is WORSE than evaluated.
	 *       This cannot happen for an honest construction (DCE only ever removes
	 *       products); it means the evaluated rank was a phantom over-claim — a
	 *       predict/build divergence BUG. {@code FAIL LOUD} (throw), do NOT warn +
	 *       silently discard, which let over-claiming stubs survive at their
	 *       inflated rank (see {@code feedback_fail_loud_dont_swallow} and the
	 *       ⟨4,19,20⟩=1002 phantom that replayed to 1012/1016).</li>
	 *   <li>wrong shape / fails spot-check / replay threw — cannot certify the
	 *       rebuild at all; discard (returns {@code false}). Not a rank bug.</li>
	 * </ul>
	 */
	/**
	 * Fail-loud rank policy for a rebuilt lineage: a rebuild may be BETTER than the
	 * evaluated rank (projection dead-code elimination can drop products) but NEVER
	 * worse. A worse rebuild ({@code replayedR > evaluatedR}) is a predict/build
	 * divergence — the evaluated rank was a phantom over-claim — so we throw rather
	 * than warn + silently discard (which let over-claiming stubs survive). Package
	 * -private so the invariant is unit-testable in isolation (no catalog load).
	 */
	static void assertRebuildNotWorse(int n, int m, int p, int evaluatedR, int replayedR) {
		if (replayedR > evaluatedR) {
			throw new IllegalStateException(String.format(
					"rebuild WORSE than evaluated for ⟨%d,%d,%d⟩: evaluated r=%d but lineage"
					+ " replayed to r=%d. Projection DCE can make a rebuild BETTER, never worse —"
					+ " this is a predict/build divergence (the evaluated rank was a phantom"
					+ " over-claim), not a tolerable discard.", n, m, p, evaluatedR, replayedR));
		}
	}

	private boolean replaysConsistently(Lineage.Node tree, NonCubicBilinearAlgorithm expected) {
		NonCubicBilinearAlgorithm rep;
		try {
			if (replayer == null) replayer = LineageReplayer.withDefaultPool(diskLookup);
			rep = replayer.replay(tree);
		} catch (RuntimeException e) {
			// Log the lineage + a stack for non-trivial replay crashes (AIOOBE off-by-one in a
			// replay op silently discards a possibly-valid strategy — we want to fix the op).
			if (e instanceof ArrayIndexOutOfBoundsException || e instanceof IndexOutOfBoundsException) {
				log.warn("[replay-crash] ⟨{},{},{}⟩ replay threw {} on lineage={}",
						expected.n, expected.m, expected.p, e.toString(), Lineage.prettyString(tree), e);
			} else {
				log.warn("lineage replay check failed: {}", e.toString());
			}
			return false;
		}
		boolean shapeOk = rep.n == expected.n && rep.m == expected.m && rep.p == expected.p;
		if (!shapeOk || !verifies(rep)) {
			log.warn("[replay-diag] ⟨{},{},{}⟩=r{} lineage replayed to ⟨{},{},{}⟩=r{} but is NON-replayable"
					+ " (shapeOk={}, spotcheck={}) — discarding to keep the catalog replay-clean.",
					expected.n, expected.m, expected.p, expected.r,
					rep.n, rep.m, rep.p, rep.r, shapeOk, verifies(rep));
			return false;
		}
		if (rep.r > expected.r) {
			log.error("[divergence-diag] ⟨{},{},{}⟩ evaluated r={} replayed r={} (+{}); top-node={}; lineage={}",
					expected.n, expected.m, expected.p, expected.r, rep.r, rep.r - expected.r,
					tree.getClass().getSimpleName(), Lineage.prettyString(tree));
			// Walk the IMMEDIATE children of any node type and show each child's replayed
			// rank vs findRank — localises the divergence even for Kron/Concat/Project/recursion
			// (not just RecombinationN). project_recomb_base_orientation_not_pinned residual.
			for (Lineage.Node child : Lineage.childrenOf(tree)) {
				try {
					NonCubicBilinearAlgorithm la = replayer.replay(child);
					int fr = diskLookup.findRank(la.n, la.m, la.p);
					log.error("    child kind={} ⟨{},{},{}⟩ replayed r={} | findRank={}{}",
							child.getClass().getSimpleName(), la.n, la.m, la.p, la.r, fr,
							la.r > fr ? "  <<< OVER vs findRank" : "");
				} catch (RuntimeException e) {
					log.error("    child kind={} replay failed: {}", child.getClass().getSimpleName(), e.toString());
				}
			}
			if (tree instanceof Lineage.RecombinationN rn) {
				log.error("    allocA={} allocB={} allocC={} (#leaves={})",
						java.util.Arrays.toString(rn.allocA()), java.util.Arrays.toString(rn.allocB()),
						java.util.Arrays.toString(rn.allocC()), rn.leaves().size());
				int idx = 0;
				for (Lineage.Node leaf : rn.leaves()) {
					try {
						NonCubicBilinearAlgorithm la = replayer.replay(leaf);
						int fr = diskLookup.findRank(la.n, la.m, la.p);
						String flag = la.r > fr ? "  <<< OVER-CLAIM" : "";
						log.error("    leaf[{}] kind={} ⟨{},{},{}⟩ replayed r={} | findRank={}{}",
								idx, leaf.getClass().getSimpleName(), la.n, la.m, la.p, la.r, fr, flag);
					} catch (RuntimeException e) {
						log.error("    leaf[{}] kind={} replay failed: {}", idx, leaf.getClass().getSimpleName(), e.toString());
					}
					idx++;
				}
			}
		}
		assertRebuildNotWorse(expected.n, expected.m, expected.p, expected.r, rep.r); // FAIL LOUD on worse
		if (rep.r < expected.r) {
			log.info("[replay-dce] ⟨{},{},{}⟩ evaluated r={} but lineage replays to a BETTER r={}"
					+ " (projection DCE) — accepting.",
					expected.n, expected.m, expected.p, expected.r, rep.r);
		}
		return true;
	}

	/**
	 * Resolve the actual factor matrices for a projection parent {@code ⟨N,M,P⟩}.
	 * Tries the on-disk lookup first (real matrices, maxDim ≤ 16); if that is
	 * empty the shape is likely a lineage-only stub (maxDim &gt; 16), so it falls
	 * back to replaying the stub's lineage via {@link LineageReplayer} — this is
	 * what lets projection fire on the 17–32 band (e.g. ⟨25,25,25⟩ ← ⟨26,26,26⟩).
	 * Returns {@code null} when no catalog entry exists or replay fails.
	 */
	/** Cache of structured Pan-trilinear-aggregation cubes ⟨N,N,N⟩, used as
	 *  HIGH-projection-margin parents: a structured cube projects far better than
	 *  the rank-best flat cube (paper §projmargin). Built lazily; empty if out of
	 *  range or the construction fails. */
	// STATIC: PanTA cubes are immutable + deterministic, and large (~300 MB at
	// n≥30); a per-instance cache would duplicate them across all worker
	// materialisers and blow the heap. Shared across threads (ConcurrentHashMap).
	private static final Map<Integer, Optional<NonCubicBilinearAlgorithm>> pantaCache =
			new java.util.concurrent.ConcurrentHashMap<>();

	private NonCubicBilinearAlgorithm pantaCube(int dim) {
		if (dim < 2 || dim > eu.solven.matmul.catalog.CatalogLimits.MAX_CUBIC_DIM) {
			return null;
		}
		return pantaCache.computeIfAbsent(dim, d -> {
			try {
				return Optional.of(
						eu.solven.matmul.papers.dis2009.PanTrilinearAggregation.build(d));
			} catch (RuntimeException e) {
				log.warn("PanTA build({}) failed — not offered as a projection parent: {}",
						d, e.toString());
				return Optional.empty();
			}
		}).orElse(null);
	}

	/**
	 * The buildable dis09 / Pan-TA Lemma-4 cube for ⟨n,n,n⟩ as a {@link Result}:
	 * rank {@code cubicBound(n)}, lineage the parametric atom
	 * {@code DIS09Lemma4(n=N)} (which {@link LineageReplayer} reconstructs via
	 * {@link eu.solven.matmul.papers.dis2009.PanTrilinearAggregation#build}). This
	 * is the SAME-shape buildable cube that {@link #tryProjection} cannot reach (it
	 * only projects from LARGER parents), so offering it directly is what stops
	 * materialise settling for a worse, mis-priced {@code Project(⟨n,n,n+1⟩)}.
	 * Returns {@code null} when the cube can't be built, doesn't verify over this
	 * field (the Q-strict ÷3 makes it invalid over e.g. F₂/ℤ), or doesn't replay.
	 */
	private Result tryDis09Cube(int n) {
		NonCubicBilinearAlgorithm cube = pantaCube(n);
		if (cube == null || !verifies(cube)) return null;
		Lineage.Node tree = new Lineage.Atom("DIS09Lemma4(n=" + n + ")");
		if (!replaysConsistently(tree, cube)) return null;
		return new Result(cube, tree, false);
	}

	/**
	 * The buildable KGP-2026 LITA cube for ⟨n,n,n⟩ as a {@link Result}: rank
	 * {@code LitaTrilinearAggregation.cubicRank(n)}, lineage the parametric atom
	 * {@code TA_lita(n=N)} (which {@link LineageReplayer} reconstructs via
	 * {@link eu.solven.matmul.papers.khoruzhii2026.LitaTaConstruction#build}). LITA
	 * beats {@link #tryDis09Cube} (and the catalog) for odd {@code n≥19} and large
	 * even n.
	 *
	 * <p>{@code verifies()} is the random spot-check (not the exact term-map), so
	 * the dense even-n cubes (the φ-embedding) are handled here too — consistent with
	 * how the catalog verifiers ({@code VerifyScheme}/{@code VerifyAllSchemes})
	 * spot-check char-0 schemes.</p>
	 */
	private Result tryLitaCube(int n) {
		if (n < eu.solven.matmul.papers.khoruzhii2026.LitaTrilinearAggregation.MIN_N) {
			return null;
		}
		NonCubicBilinearAlgorithm cube;
		try {
			cube = eu.solven.matmul.papers.khoruzhii2026.LitaTaConstruction.build(n);
		} catch (RuntimeException e) {
			log.warn("LITA build({}) failed — not offered as a cube: {}", n, e.toString());
			return null;
		}
		if (!verifies(cube)) return null;
		Lineage.Node tree = new Lineage.Atom("TA_lita(n=" + n + ")");
		if (!replaysConsistently(tree, cube)) return null;
		return new Result(cube, tree, false);
	}


	/** Like {@link #resolveParent} but also returns the SOURCE file the parent was
	 *  loaded/replayed from, so a building-block leaf can be pinned by that file's
	 *  stamped hash (always resolvable) instead of the in-memory oriented hash. */
	private ParentHit resolveParentHit(int N, int M, int P) {
		String key = N + "x" + M + "x" + P;
		ParentHit cached = parentCache.get(key);
		if (cached != null) return cached;
		ParentHit resolved = resolveParentUncached(N, M, P);
		if (resolved != null) parentCache.put(key, resolved);
		return resolved;
	}

	private ParentHit resolveParentUncached(int N, int M, int P) {
		Optional<FieldAwareLookup.WithSource> direct = diskLookup.findWithSource(N, M, P);
		// Stub-blindness guard: findWithSource (like find) SKIPS lineage-only stubs, so
		// it can return a dense Schwartz-Zwecher ⟨2k,2k,2k⟩=5596 while a strictly-cheaper
		// BUILDABLE dis09 cube STUB (5566, a DIS09Lemma4(n) formula atom) exists that
		// findRank sees. Returning the dense 5596 makes the concat operand EXCEED what
		// the upward search priced it at (findRank=5566) → the ⟨2k,2k,2k+2⟩ predict/build
		// divergence (⟨22,22,24⟩ evaluated 6303, built 6333). Only short-circuit on the
		// dense hit when findRank agrees nothing cheaper exists; otherwise fall through to
		// the rank-ordered scan below, which replays the cheaper stub (dis09) first.
		// [[project_findrank_poisoned_by_nonbuildable_stubs]] / find-vs-findRank stub-blindness.
		if (direct.isPresent() && diskLookup.findRank(N, M, P) >= direct.get().alg().r) {
			return new ParentHit(direct.get().alg(), direct.get().path());
		}
		// No usable dense file (or a cheaper stub is known) → try the lineage-only stubs
		// (maxDim>16) in rank order and return the first that REPLAYS TO ≤ ITS CLAIMED
		// RANK — don't give up on the single lowest-rank file when it happens to be a
		// non-replayable stub (e.g. a Schwartz-Zwecher or Strassen-mask cube), since a
		// slightly higher-rank but replayable sibling (e.g. the dis09 cube) may follow.
		// CRITICAL: a candidate whose lineage replays WORSE than it claims is a CORRUPT
		// phantom (the ⟨2k,2k,2k⟩ derived dup: claims cubicBound 5566, its Project lineage
		// replays 5596). It "replays" fine, so the old "first that replays" rule returned
		// IT and the concat operand exceeded what findRank priced → the divergence. Skip
		// any over-claimer so the honest dis09 cube (claims 5566, replays 5566) wins.
		// Honour the commutativity mode (default NC): an NC projection must lift over
		// NC rings, so a commutative scheme (e.g. the Rosowski ⟨N,3,3⟩ Algorithm-1,
		// which is also non-bilinear and so un-replayable by the bilinear replayer) is
		// invalid as a parent. A commutative-matmul sweep (allowCommutative) keeps them.
		java.util.List<Path> files = allowCommutative
				? diskLookup.findFiles(N, M, P)
				: diskLookup.findFilesNonCommutative(N, M, P);
		if (files.isEmpty()) return null; // no entry at all
		if (replayer == null) replayer = LineageReplayer.withDefaultPool(diskLookup);
		boolean anyFailed = false;
		for (Path path : files) {
			try {
				long t0 = System.nanoTime();
				NonCubicBilinearAlgorithm parent = replayer.replayFromFile(path.toFile());
				long ms = (System.nanoTime() - t0) / 1_000_000L;
				if (ms > 1_000) {
					log.info("projection: replayed stub parent ⟨{},{},{}⟩ r={} in {}ms ({})",
							N, M, P, parent.r, ms, path.getFileName());
				}
				int claimed = recordedRankOf(path);
				if (claimed >= 0 && parent.r > claimed) {
					// CORRUPT_RANK phantom: lineage replays worse than the file claims. Using it
					// would silently inflate the parent → predict/build divergence. Skip it.
					log.warn("projection parent {} CLAIMS r={} but lineage REPLAYS WORSE to r={} "
							+ "(corrupt over-claim) — skipping so an honest sibling wins.",
							path.getFileName(), claimed, parent.r);
					anyFailed = true;
					continue;
				}
				return new ParentHit(parent.orientAs(N, M, P).orElse(parent), path);
			} catch (RuntimeException e) {
				anyFailed = true;
				log.warn("projection parent {} is NOT replayable: {} — skipping to next candidate.",
						path.getFileName(), e.toString());
			}
		}
		if (anyFailed) {
			log.warn("projection parent ⟨{},{},{}⟩: ALL {} catalog file(s) failed to replay"
					+ " (or were corrupt over-claimers) — fix the stub lineages so projection can use"
					+ " this shape.", N, M, P, files.size());
		}
		return null;
	}

	/** The rank ({@code m}) a scheme file records on disk, or −1 if unreadable. Used to
	 *  detect CORRUPT_RANK over-claimers (lineage replays worse than the file claims). */
	private static int recordedRankOf(Path path) {
		try {
			tools.jackson.databind.JsonNode d = SchemeIO.parseJson(path.toFile());
			return d.has("m") ? d.get("m").asInt() : -1;
		} catch (Exception e) {
			return -1;
		}
	}

	/** Memoise a freshly-composed scheme and (if configured) write it to disk. */
	private void persist(int n, int m, int p, Result built) {
		String key = canon(n, m, p);
		// WRITE-TIME CYCLE GUARD (user 2026-06-25): never memoise/persist a scheme whose lineage is
		// cyclic or carries an unresolvable "@ref?:Ln" fallback. A cycle can NEVER be uniquely optimal
		// — project-then-concat-back reconstructs the removed slice, so the rank is ≥ the original — so
		// FAILING here loses no genuine win: a real rank is reachable acyclically (the search keeps the
		// clean alternative), a phantom is correctly skipped. Loud, so a "cyclic = claimed-better" case
		// is surfaced for investigation rather than silently written. [[cyclic lineage → silent SOE]]
		if (built.lineage != null) {
			String corruption = lineageCorruption(built.lineage, n, m, p,
					SchemeIO.contentHash(built.alg).substring(0, 7));
			if (corruption != null) {
				log.error("REFUSING ⟨{},{},{}⟩ r={}: cyclic/corrupt lineage [{}] — SKIPPED (a cycle is "
						+ "never uniquely optimal; if this rank was thought best, investigate)",
						n, m, p, built.alg.r, corruption);
				return;
			}
		}
		derived.put(key, built.alg);
		derivedLineage.put(key, built.lineage);
		// A better scheme now exists at ⟨n,m,p⟩, so any cached projection parent for
		// that shape (in ANY orientation) is stale — drop the whole orbit so the next
		// child that projects through ⟨n,m,p⟩ re-resolves the freshly-written winner
		// (cascade wins within the same pass must not be masked by the cache).
		invalidateParentCache(n, m, p);
		if (writeNewSchemes && writeRoot != null) {
			writeToDisk(n, m, p, built);
			// Drop in-memory copy now that the scheme is on disk; the
			// FieldAwareLookup.onSchemeWritten hook has already added it to the
			// index, so subsequent lookups read it from cache.
			derived.remove(key);
			derivedLineage.remove(key);
		}
	}

	/** Drop every orientation of {@code ⟨n,m,p⟩} from the projection-parent cache. */
	private void invalidateParentCache(int n, int m, int p) {
		int[][] perms = { { n, m, p }, { n, p, m }, { m, n, p }, { m, p, n }, { p, n, m }, { p, m, n } };
		for (int[] q : perms) parentCache.remove(q[0] + "x" + q[1] + "x" + q[2]);
	}

	/** Standalone plain-Kronecker probe: the cheapest {@code ⟨n1,m1,p1⟩⊗⟨n2,m2,p2⟩}
	 *  over all divisor factorisations, priced by the sota resolver and materialised
	 *  via {@link #buildKronecker} (durable lineage). {@link #materialise} routes the
	 *  {@link #STRAT_KRONECKER} token through findBestStrategy's KRONECKER kind (same
	 *  underlying {@link KroneckerSplitSearch#findBest}); this direct entry serves
	 *  probes/tests. Returns null unless the verified product lands strictly below
	 *  {@code upper}. */
	Result tryKronecker(int n, int m, int p, long upper) {
		Optional<KroneckerSplitSearch.KroneckerSplit> split =
				KroneckerSplitSearch.findBest(n, m, p, sota);
		if (split.isEmpty() || split.get().totalRank() >= upper) {
			return null;
		}
		Result built = buildKronecker(split.get());
		if (built == null || built.alg.r >= upper || !verifies(built.alg)) {
			return null;
		}
		return built;
	}

	private Result buildKronecker(KroneckerSplitSearch.KroneckerSplit k) {
		// Target = the Kronecker product shape; resolveSubScheme gives the factor
		// schemes' matrices even when they are stubs improve-mode declines to re-derive
		// (the same stub-withholding that broke concat), and guards shape-acyclicity.
		int tn = k.n1() * k.n2(), tm = k.m1() * k.m2(), tp = k.p1() * k.p2();
		Optional<Result> outer = resolveSubScheme(k.n1(), k.m1(), k.p1(), tn, tm, tp);
		Optional<Result> inner = resolveSubScheme(k.n2(), k.m2(), k.p2(), tn, tm, tp);
		if (outer.isEmpty() || inner.isEmpty()) return null;
		NonCubicBilinearAlgorithm composed = Compose.kroneckerGeneral(outer.get().alg, inner.get().alg);
		Lineage.Node tree = new Lineage.KronProduct(outer.get().lineage, inner.get().lineage);
		return new Result(composed, tree, false);
	}

	/**
	 * Resolve a composition SUB-shape to a usable {@link Result} (matrices + a
	 * replayable lineage). Tries {@link #materialise} first (improve-mode may better
	 * it); if that returns empty — which happens when the sub-shape is a lineage-only
	 * STUB that improve-mode declines to re-derive (materialise line ~227 returns
	 * empty rather than expanding the stub) — it falls back to REPLAYING the stub so
	 * the composition can still use its matrices. Without this, concat/Kronecker
	 * silently fail to build on the 17–32 band: e.g. ⟨26,26,32⟩ = ⟨26,26,6⟩ +
	 * ⟨26,26,26⟩ predicts 11165 but the ⟨26,26,26⟩=8658 stub's matrices were withheld,
	 * so the obvious win was reported as "no-improvement". The lineage is the shape-ref
	 * atom, which replays to the same catalog best the matrices came from.
	 */
	private Optional<Result> resolveSubScheme(int n, int m, int p, int tn, int tm, int tp) {
		// LOAD by default (no implicit recursion): a leaf comes from the disk catalog as-is.
		// OPTIMIZE (recursiveDerive) only when explicitly asked — e.g. introducing a new base.
		Optional<Result> mat = recursiveDerive ? materialise(n, m, p) : diskBest(n, m, p);
		// A composed sub-block whose own lineage is CORRUPT (dangling/cyclic) must NOT be
		// handed to a parent: writeToDisk would refuse to persist it, yet the parent would
		// pin it by a file that never lands → the parent itself dangles. Reject it here so
		// the parent falls back to a clean alternative (or another split). This is the
		// production-side twin of the write-time DAG guard.
		if (mat.isPresent() && mat.get().lineage != null && !mat.get().fromDisk()) {
			String h = SchemeIO.contentHash(mat.get().alg).substring(0, 7);
			if (lineageCorruption(mat.get().lineage, n, m, p, h) != null) {
				mat = Optional.empty();
			}
		}
		// Fallback-ref guard — runs even for fromDisk leaves (the !fromDisk gate above skips them):
		// a disk scheme whose STORED lineage is cyclic parses to an unresolvable Atom("@ref?:L0")
		// IN MEMORY (parseLineageNode forward-ref fallback). Embedding it propagates the corruption
		// into the parent (the ⟨17,19,20⟩ self-referential ⟨2,3,3⟩ memo node). Cheap in-memory scan.
		if (mat.isPresent() && mat.get().lineage != null && hasFallbackRef(mat.get().lineage)) {
			mat = Optional.empty();
		}
		// findRank may know a strictly-better REPLAYABLE block than materialise returns:
		// a lower-rank STUB that findWithSource skips (maxDim>16) and compose can't rebuild
		// — e.g. a DIS09 ⟨20,20,20⟩=4340 cube vs the explicit Schwartz-Zwecher 4378. The
		// upward search SCORES blocks with findRank, so building the worse block makes the
		// concat/Kron rebuild EXCEED its evaluated rank → predict/build divergence (AND a
		// lost win: ⟨20,20,22⟩=610+4340=4950 vs the built 610+4378=4988). Prefer the stub
		// when it actually REPLAYS cheaper than materialise (so build == what scoring priced).
		int matR = mat.map(r -> r.alg.r).orElse(Integer.MAX_VALUE);
		if (diskLookup.findRank(n, m, p) < matR && !subSchemeReaches(n, m, p, tn, tm, tp)) {
			ParentHit better = resolveParentHit(n, m, p);
			if (better != null) {
				NonCubicBilinearAlgorithm a = better.alg();
				Optional<NonCubicBilinearAlgorithm> oriented =
						(a.n == n && a.m == m && a.p == p) ? Optional.of(a) : a.orientAs(n, m, p);
				if (oriented.isPresent() && oriented.get().r < matR) {
					return oriented.map(o -> new Result(o, durableLeafRef(better.path(), n, m, p), true));
				}
			}
		}
		if (mat.isPresent()) return mat;
		// Fallback: the sub-shape is a stub improve-mode declined to re-derive — replay
		// it so the composition can use its matrices. But FIRST guard the shape-level
		// acyclicity invariant: a scheme for target T must not transitively reference T.
		// resolveParent's raw replay bypasses compose()'s inFlight guard, so referencing
		// a stub that (via the catalog graph) reaches T would close a cycle AFTER we
		// overwrite T — the ⟨15,26,32⟩→⟨15,26,30⟩→⟨15,26,31⟩→⟨15,26,32⟩ p-band loop. Skip
		// such a sub-block (the composition tries another split, or reports no win).
		if (subSchemeReaches(n, m, p, tn, tm, tp)) return Optional.empty();
		ParentHit hit = resolveParentHit(n, m, p);
		if (hit == null) return Optional.empty();
		NonCubicBilinearAlgorithm alg = hit.alg();
		Optional<NonCubicBilinearAlgorithm> oriented =
				(alg.n == n && alg.m == m && alg.p == p) ? Optional.of(alg) : alg.orientAs(n, m, p);
		// Pin the building block by the SOURCE FILE's stamped hash + an OrientAs node,
		// NOT by contentHash(orientedAlg). The oriented (and, for a replayed STUB, the
		// in-memory-reconstructed) scheme is stamped on NO file, so pinning its hash
		// produced a ref that could never resolve — resolveLeaf then silently fell back
		// to shape-best and the rank drifted, so the write-guard discarded the win (the
		// ~46% refusal rate observed in the non-5 sweep, all on stub-sourced parents).
		// durableLeafRef references the file we actually loaded (always resolvable) and
		// reconstructs the (n,m,p) orientation bit-exactly — honouring "either the ref
		// resolves, or there is no ref". [[lineage replay must be bit-exact]]
		return oriented.map(a -> new Result(a, durableLeafRef(hit.path(), n, m, p), true));
	}

	/**
	 * Build an always-resolvable, bit-exact leaf ref for a building block loaded from
	 * {@code src}: pin the file's NATIVE shape + its stamped content hash (which a
	 * later {@code resolveLeaf} resolves via {@code findByHash} for dense files or the
	 * stamped-hash stub scan), wrapped in an {@link Lineage.OrientAs} when the file's
	 * native shape differs from the requested ⟨n,m,p⟩. Falls back to a bare (still
	 * resolvable) shape ref if the file carries no stamped hash.
	 */
	/** Precise {@code N×M×P@hash} ref-string for a projection parent file; fails rather
	 *  than emit a bare/@sota parent (projection DCE needs the exact parent scheme). */
	private String preciseParentRef(Path src, int N, int M, int P) {
		String hash;
		try {
			hash = SchemeIO.readHash(SchemeIO.parseJson(src.toFile()));
		} catch (Exception e) {
			throw new IllegalStateException("projection parent " + src + " unreadable — refusing to"
					+ " emit a bare/@sota parent ref for ⟨" + N + "," + M + "," + P + "⟩", e);
		}
		if (hash == null || hash.isBlank()) hash = hashFromFilename(src); // older atoms hash-in-name only
		if (hash == null || hash.isBlank()) {
			throw new IllegalStateException("projection parent " + src.getFileName() + " has no content"
					+ " hash (neither JSON nor filename) — refusing a bare/@sota parent ref for ⟨" + N
					+ "," + M + "," + P + "⟩ (the exact parent is load-bearing for projection DCE).");
		}
		return N + "x" + M + "x" + P + "@" + hash;
	}

	/**
	 * Pin a projection parent as a replayable {@link Lineage.Node}: the file's NATIVE
	 * (as-stored) {@code shape@hash} Atom, wrapped in an EXACT-perm {@link Lineage.Transpose}
	 * whenever the search oriented the parent into a different ⟨N,M,P⟩ frame. This keeps the
	 * {@code @hash} a true content ref of the file as stored and makes the orientation EXPLICIT
	 * in the lineage, so replay reconstructs the exact parent the search projected — never an
	 * ambiguous re-{@code orientAs}. For an EQUAL-AXIS parent (e.g. ⟨19,20,20⟩, two 20-axes) that
	 * re-orientation is ambiguous and can pick a worse-projecting axis → predict/build divergence
	 * (⟨19,19,20⟩ evaluated 4154, replayed 4237). The exact perm is found by matching the search's
	 * oriented parent against {@link SymmetryTransforms#s3OrbitWithPerms} (by content hash).
	 * [[project_projection_parent_orientation_not_pinned]]
	 */
	private Lineage.Node projectionParentNode(Path src, NonCubicBilinearAlgorithm orientedParent,
			int N, int M, int P) {
		String hash;
		tools.jackson.databind.JsonNode root;
		try {
			root = SchemeIO.parseJson(src.toFile());
			hash = SchemeIO.readHash(root);
		} catch (Exception e) {
			throw new IllegalStateException("projection parent " + src + " unreadable — refusing a"
					+ " bare parent ref for ⟨" + N + "," + M + "," + P + "⟩", e);
		}
		if (hash == null || hash.isBlank()) hash = hashFromFilename(src);
		if (hash == null || hash.isBlank()) {
			throw new IllegalStateException("projection parent " + src.getFileName() + " has no content"
					+ " hash — refusing a bare parent ref (exact parent is load-bearing for projection DCE).");
		}
		tools.jackson.databind.JsonNode shp = root.get("n");
		if (shp == null || !shp.isArray() || shp.size() != 3) {
			throw new IllegalStateException("projection parent " + src.getFileName() + " has no n[] shape.");
		}
		int nn = shp.get(0).asInt(), nm = shp.get(1).asInt(), np = shp.get(2).asInt();
		Lineage.Node atom = new Lineage.Atom(nn + "x" + nm + "x" + np + "@" + hash);
		if (nn == N && nm == M && np == P) return atom; // already in the projected frame; no Transpose
		// Orientation needed: pin the EXACT axis-perm reproducing the search's oriented parent.
		if (replayer == null) replayer = LineageReplayer.withDefaultPool(diskLookup);
		NonCubicBilinearAlgorithm nativeAlg = replayer.replayFromFile(src.toFile());
		String wanted = SchemeIO.contentHash(orientedParent);
		for (eu.solven.matmul.SymmetryTransforms.S3Variant v
				: eu.solven.matmul.SymmetryTransforms.s3OrbitWithPerms(nativeAlg)) {
			NonCubicBilinearAlgorithm a = v.alg();
			if (a.n == N && a.m == M && a.p == P && SchemeIO.contentHash(a).equals(wanted)) {
				return new Lineage.Transpose(atom, v.perm());
			}
		}
		throw new IllegalStateException("projection parent ⟨" + N + "," + M + "," + P + "⟩: no S3 perm of "
				+ src.getFileName() + " reproduces the search's oriented parent — refusing an ambiguous ref.");
	}

	/** Content hash from a scheme filename's {@code -<hash7>.json} suffix — the canonical
	 *  label encodes it even when the JSON body carries no {@code hash} field (older imports). */
	private static String hashFromFilename(Path src) {
		java.util.regex.Matcher m = java.util.regex.Pattern
				.compile("-([0-9a-f]{4,})\\.json$").matcher(src.getFileName().toString());
		return m.find() ? m.group(1) : null;
	}

	private Lineage.Node durableLeafRef(Path src, int n, int m, int p) {
		int sn, sm, sp;
		String hash;
		try {
			var root = SchemeIO.parseJson(src.toFile());
			var nNode = root.get("n");
			sn = nNode.get(0).asInt(); sm = nNode.get(1).asInt(); sp = nNode.get(2).asInt();
			hash = SchemeIO.readHash(root);
		} catch (Exception e) {
			// FAIL LOUD, do NOT emit a bare leaf. A bare "NxMxP" ref is a @sota / best-at-
			// shape leaf — a CITED BOUND, not a precise scheme — and building an explicit
			// scheme over it is invalid (the exact sub-scheme matters for downstream
			// DCE/projection). "We must not fall back to naive/best by default: fail."
			throw new IllegalStateException("durableLeafRef: cannot read building block " + src
					+ " to pin a precise leaf for ⟨" + n + "," + m + "," + p + "⟩ — refusing to"
					+ " emit a bare/@sota leaf in an explicit build", e);
		}
		if (hash == null || hash.isBlank()) hash = hashFromFilename(src); // older atoms hash-in-name only
		if (hash == null || hash.isBlank()) {
			throw new IllegalStateException("durableLeafRef: building block " + src.getFileName()
					+ " carries no content hash (neither JSON nor filename) — refusing to emit a"
					+ " bare/@sota leaf for ⟨" + n + "," + m + "," + p + "⟩ in an explicit build.");
		}
		String shapeRef = sn + "x" + sm + "x" + sp;
		Lineage.Node atom = new Lineage.Atom(shapeRef + "@" + hash);
		return (sn == n && sm == m && sp == p) ? atom : Lineage.orientAs(atom, sn, sm, sp, n, m, p);
	}

	private static final java.util.regex.Pattern REF_SHAPE =
			java.util.regex.Pattern.compile("(\\d+)x(\\d+)x(\\d+)");

	/** Does ⟨n,m,p⟩'s current catalog construction transitively reference shape
	 *  ⟨tn,tm,tp⟩ (orientation-insensitive)? BFS over the catalog file graph, exactly
	 *  the relation {@code DetectCyclicStubs} flags — used to refuse a composition
	 *  sub-block that would close a shape-level cycle with the target being composed. */
	private boolean subSchemeReaches(int n, int m, int p, int tn, int tm, int tp) {
		int[] tgt = { tn, tm, tp };
		java.util.Arrays.sort(tgt);
		java.util.Set<String> seen = new java.util.HashSet<>();
		java.util.Deque<int[]> stack = new java.util.ArrayDeque<>();
		stack.push(new int[] { n, m, p });
		while (!stack.isEmpty()) {
			int[] s = stack.pop();
			int[] ss = { s[0], s[1], s[2] };
			java.util.Arrays.sort(ss);
			if (!seen.add(ss[0] + "x" + ss[1] + "x" + ss[2]) || seen.size() > 8192) continue;
			Optional<java.nio.file.Path> f = diskLookup.findFile(s[0], s[1], s[2]);
			if (f.isEmpty()) continue;
			for (int[] ref : stubRefShapes(f.get())) {
				int[] rs = { ref[0], ref[1], ref[2] };
				java.util.Arrays.sort(rs);
				if (java.util.Arrays.equals(rs, tgt)) return true;
				stack.push(ref);
			}
		}
		return false;
	}

	/** All ⟨n,m,p⟩ shapes referenced in a stub file's lineage ("ref" string tokens). */
	private List<int[]> stubRefShapes(java.nio.file.Path file) {
		List<int[]> out = new java.util.ArrayList<>();
		try {
			collectRefShapes(SchemeIO.parseJson(file.toFile()).get("lineage"), out);
		} catch (Exception e) {
			// unreadable lineage → treat as no refs
		}
		return out;
	}

	private void collectRefShapes(tools.jackson.databind.JsonNode node, List<int[]> out) {
		if (node == null) return;
		if (node.isObject()) {
			tools.jackson.databind.JsonNode ref = node.get("ref");
			if (ref != null && ref.isString()) {
				java.util.regex.Matcher m = REF_SHAPE.matcher(ref.asString());
				if (m.find()) out.add(new int[] { Integer.parseInt(m.group(1)),
						Integer.parseInt(m.group(2)), Integer.parseInt(m.group(3)) });
			}
			for (var e : node.properties()) collectRefShapes(e.getValue(), out);
		} else if (node.isArray()) {
			for (tools.jackson.databind.JsonNode c : node) collectRefShapes(c, out);
		}
	}

	private Result buildConcat(ConcatSplitSearch.ConcatSplit c) {
		// axis: 0 = n-axis (concatBelow), 1 = m-axis (concatInner), 2 = p-axis (concatRight)
		int leftL = c.leftSize(), rightL = c.rightSize();
		if (c.axis() == 0) {
			Optional<Result> top = resolveSubScheme(leftL, c.m(), c.p(), c.n(), c.m(), c.p());
			Optional<Result> bot = resolveSubScheme(rightL, c.m(), c.p(), c.n(), c.m(), c.p());
			if (top.isEmpty() || bot.isEmpty()) return null;
			NonCubicBilinearAlgorithm composed = Compose.concatBelow(top.get().alg, bot.get().alg);
			Lineage.Node tree = new Lineage.ConcatRows(top.get().lineage, bot.get().lineage);
			return new Result(composed, tree, false);
		}
		if (c.axis() == 1) {
			// m-axis contraction sum: split m, share n,p; C = A1·B1 + A2·B2.
			Optional<Result> left = resolveSubScheme(c.n(), leftL, c.p(), c.n(), c.m(), c.p());
			Optional<Result> right = resolveSubScheme(c.n(), rightL, c.p(), c.n(), c.m(), c.p());
			if (left.isEmpty() || right.isEmpty()) return null;
			NonCubicBilinearAlgorithm composed = Compose.concatInner(left.get().alg, right.get().alg);
			Lineage.Node tree = new Lineage.SumInner(left.get().lineage, right.get().lineage);
			return new Result(composed, tree, false);
		}
		// p-axis
		Optional<Result> left = resolveSubScheme(c.n(), c.m(), leftL, c.n(), c.m(), c.p());
		Optional<Result> right = resolveSubScheme(c.n(), c.m(), rightL, c.n(), c.m(), c.p());
		if (left.isEmpty() || right.isEmpty()) return null;
		NonCubicBilinearAlgorithm composed = Compose.concatRight(left.get().alg, right.get().alg);
		Lineage.Node tree = new Lineage.ConcatCols(left.get().lineage, right.get().lineage);
		return new Result(composed, tree, false);
	}

	private Result buildRecombination(int n, int m, int p,
			BlockSplitSearch.MultiBaseSplitCandidate r) {
		// Generic Pan-TA fusion: when the candidate base is a naïve grid (⟨1,2,2⟩ peel,
		// FMM's ⟨2,3,3⟩ grids, …), its disjoint cyclic-rotation single-block product pairs
		// FUSE via Pan trilinear aggregation (fusedRank=abc+ab+bc+ca per pair), exactly the
		// saving PairedSubProducts already prices in scoring. The generic
		// constructWithAllocation path below glues those pairs INDEPENDENTLY and so builds
		// at the un-fused rank, never realising the scored saving — route through
		// Recombination.constructWithTaFusion so the BUILT scheme matches the SCORED rank
		// and is replayable (Lineage.RecombinationTaN).
		Result taFused = tryBuildTaFusion(n, m, p, r);
		if (taFused != null) {
			return taFused;
		}
		// The recombination uses a base ⟨baseA, baseB, baseC⟩ outer scheme +
		// allocations. Sub-shapes consulted are determined by the base's U/V/W
		// structure and the allocation. We use a RecursiveLookup overlay so
		// any missing sub-shape gets materialised on demand.
		RecursiveLookup overlay = new RecursiveLookup();
		try {
			NonCubicBilinearAlgorithm composed = Recombination.constructWithAllocation(
					r.base(), overlay, r.allocA(), r.allocB(), r.allocC());
			// Collect lineage leaves from the overlay (in order seen).
			List<Lineage.Node> leaves = new ArrayList<>(overlay.lineageOrder);
			// Record the base via its ORIGIN lineage (a resolvable Atom, or an
			// AxisFlip/AxisPermute node) rather than the bare pool label — the label
			// can be a synthetic, unresolvable string like
			// "Winograd<2,2,2>=7 :: AXIS_FLIP" that breaks replay. (Task #91.)
			//
			// Subtlety verified empirically (TestStrassenCousinHunt): for a 2×2×2
			// base the axis-flip IS the per-axis allocation ORDERING — canonical
			// winograd_1971 at (8,9)³ = 2930, at (9,8)³ = 2954. So an axis-flip-orbit
			// WINNER that the search reached is reproduced by the CANONICAL atom at
			// the recorded ordered allocation; baseOriginLineage() is the canonical
			// Atom and the lineage replays exactly. No AxisFlip wrapper is needed for
			// this family — the allocation tuple carries the orientation.
			//
			// CAVEAT (task #91 residual): if a genuine matrix flip (one NOT expressible
			// as an allocation reorder, e.g. for ≥3-part bases) ever wins, the upstream
			// AxisFlip node carries the orbit-list INDEX, not the true flip mask
			// (axisFlipOrbit dedups/compacts the list) — that case would need the true
			// mask threaded through. The stub-replay check below would catch it.
			Lineage.Node baseNode = r.baseOriginLineage();
			if (baseNode == null) {
				// Pool entry without an origin lineage (every stock pool stamps one — this is
				// an ad-hoc caller, e.g. a GL-member base). Best effort: pin by content hash
				// if the exact base resolves in the catalog; otherwise fall back to the label
				// and WARN — a label ref is a best-at-shape cited bound, the stub will not be
				// explicitable, and its replay is hostage to future catalog state.
				String hash = SchemeIO.contentHash(r.base());
				if (diskLookup.findByHash(r.base().n, r.base().m, r.base().p, hash).isPresent()) {
					baseNode = new Lineage.Atom(r.base().n + "x" + r.base().m + "x" + r.base().p + "@" + hash);
				} else {
					baseNode = new Lineage.Atom(r.baseLabel());
					log.warn("recombination base ⟨{},{},{}⟩ '{}' has no origin lineage and no catalog match"
							+ " — writing an UNPINNED label ref (stub will not be explicitable)",
							r.base().n, r.base().m, r.base().p, r.baseLabel());
				}
			}
			// CYCLE PREVENTION: the base's origin lineage can itself carry inherited cycle corruption
			// (an "@ref?:L0" fallback read from a legacy cyclic disk scheme). Embedding it spreads the
			// corruption into the parent recombination (the ⟨17,19,20⟩ ⟨2,3,3⟩ base case). Pin the base
			// by content hash instead — replay loads it; corruption stays contained to the legacy scheme.
			if (hasFallbackRef(baseNode)) {
				String pinned = r.base().n + "x" + r.base().m + "x" + r.base().p + "@" + SchemeIO.contentHash(r.base());
				log.warn("recombination BASE ⟨{},{},{}⟩ had corrupt (cyclic) origin lineage — pinned by hash ({})",
						r.base().n, r.base().m, r.base().p, pinned);
				baseNode = new Lineage.Atom(pinned);
			}
			Lineage.Node tree = new Lineage.RecombinationN(
					baseNode,
					r.allocA().clone(), r.allocB().clone(), r.allocC().clone(),
					leaves);
			return new Result(composed, tree, false);
		} catch (Exception e) {
			log.warn("recombination ⟨{},{},{}⟩ failed: {}", n, m, p, e.getMessage());
			return null;
		}
	}

	/**
	 * Generic Pan-TA fusion over ANY naïve-grid base (⟨1,2,2⟩ peel, FMM's ⟨2,3,3⟩
	 * grids, …): {@link Recombination#constructWithTaFusion} fuses every disjoint
	 * cyclic-rotation single-block product pair at {@code fusedRank} instead of the
	 * two leaves' summed rank. Returns {@code null} (→ caller falls back to the glue
	 * path) when the base is not a naïve grid or no pair fuses.
	 *
	 * <p>Built by REPLAYING the recorded {@link Lineage.RecombinationTaN}: the base is
	 * the naïve grid (Atom {@code "naive-NxMxP"}) and each unpaired leaf is pinned to
	 * the NON-COMMUTATIVE best by {@link #ncLeafNode} (the global rank-best can be a
	 * COMMUTATIVE Waksman/Rosowski scheme — e.g. ⟨26,3,3⟩=159 vs the NC 175 — invalid
	 * as an NC leaf and not what scoring priced). Building via replay guarantees
	 * build ≡ replay, so the registered rank is reproducible.</p>
	 */
	private Result tryBuildTaFusion(int n, int m, int p, BlockSplitSearch.MultiBaseSplitCandidate r) {
		NonCubicBilinearAlgorithm base = r.base();
		if (!Recombination.isNaiveGrid(base)) {
			return null;   // only naïve grids carry TA-fusable disjoint cyclic pairs
		}
		if (replayer == null) {
			replayer = LineageReplayer.withDefaultPool(diskLookup);
		}
		// Base node = the naïve grid; build the base alg BY REPLAY so build ≡ replay.
		Lineage.Node baseNode = new Lineage.Atom("naive-" + base.n + "x" + base.m + "x" + base.p);
		NonCubicBilinearAlgorithm baseAlg = replayer.replay(baseNode);

		// Lazy NC sub-resolver for UNPAIRED leaves: cache (sorted-shape → node, alg). The
		// recorded leaves are exactly the nodes resolved here, so replay re-resolves the
		// same schemes (it keys by sorted shape + orients identically).
		java.util.LinkedHashMap<String, Lineage.Node> leafNodes = new java.util.LinkedHashMap<>();
		java.util.Map<String, NonCubicBilinearAlgorithm> leafAlgs = new java.util.HashMap<>();
		Recombination.SubResolver resolver = (sz) -> {
			int[] s = { sz[0], sz[1], sz[2] };
			java.util.Arrays.sort(s);
			String key = s[0] + "x" + s[1] + "x" + s[2];
			NonCubicBilinearAlgorithm a = leafAlgs.get(key);
			if (a == null) {
				Lineage.Node node = ncLeafNode(s[0], s[1], s[2]);   // sorted orientation, NC-pinned
				if (node == null) return null;
				try {
					a = replayer.replay(node);
				} catch (RuntimeException e) {
					return null;
				}
				leafNodes.put(key, node);
				leafAlgs.put(key, a);
			}
			return a.orientAs(sz[0], sz[1], sz[2]).orElse(a);
		};

		Recombination.TaFusedConstruction tc;
		try {
			tc = Recombination.constructWithTaFusion(baseAlg, resolver, sota,
					r.allocA(), r.allocB(), r.allocC());
		} catch (RuntimeException e) {
			log.info("⟨{},{},{}⟩ TA-fusion build failed: {} — falling back to glue", n, m, p, e.toString());
			return null;
		}
		if (tc.fusedPairs().isEmpty()) {
			return null;   // no cyclic pair fused — the glue path is equivalent
		}
		NonCubicBilinearAlgorithm alg = tc.alg();
		if (alg.n != n || alg.m != m || alg.p != p) {
			log.warn("⟨{},{},{}⟩ TA-fusion produced ⟨{},{},{}⟩ — discarding", n, m, p, alg.n, alg.m, alg.p);
			return null;
		}
		if (!verifies(alg)) {
			log.warn("⟨{},{},{}⟩ TA-fusion built at r={} but FAILED spot-check — discarding", n, m, p, alg.r);
			return null;
		}
		Lineage.Node tree = new Lineage.RecombinationTaN(baseNode,
				r.allocA().clone(), r.allocB().clone(), r.allocC().clone(),
				new java.util.ArrayList<>(leafNodes.values()));
		log.info("⟨{},{},{}⟩ TA-fusion via base ⟨{},{},{}⟩: {}",
				n, m, p, base.n, base.m, base.p, describeFusionSummary(baseAlg, sota, r, tc, alg));
		return new Result(alg, tree, false);
	}

	/**
	 * The Pan-TA breakdown one-liner for logs (and the same string the catalog stamps
	 * as {@code ta_fusion.summary}). Falls back to a bare pair-count if rank-only
	 * re-pricing can't reproduce a naïve grid (should not happen for a fused build).
	 */
	private String describeFusionSummary(NonCubicBilinearAlgorithm baseAlg,
			Recombination.SotaResolver sota, BlockSplitSearch.MultiBaseSplitCandidate r,
			Recombination.TaFusedConstruction tc, NonCubicBilinearAlgorithm alg) {
		Recombination.TaFusionBreakdown bd = Recombination.describeTaFusion(
				baseAlg, sota, r.allocA(), r.allocB(), r.allocC());
		return bd != null ? bd.summary()
				: tc.fusedPairs().size() + " cross-pair(s) fused -> r" + alg.r;
	}

	/**
	 * A replayable lineage leaf pinned to the NON-COMMUTATIVE best scheme of shape
	 * ⟨a,b,c⟩, oriented to exactly ⟨a,b,c⟩. Pins by content hash ({@code shape@hash7})
	 * so replay resolves THIS file — not a commutative rank-best sibling
	 * ({@link FieldAwareLookup#findFilesNonCommutative} excludes commutative; the hash
	 * pin defends against a future lower-rank commutative entry being added). Walks the
	 * NC files rank-ascending and returns the first that actually replays to ⟨a,b,c⟩
	 * (skipping e.g. a non-bilinear Rosowski-Algorithm-1 file that fails to load).
	 * Returns {@code null} when no NC scheme at this shape is replayable.
	 */
	private Lineage.Node ncLeafNode(int a, int b, int c) {
		if (replayer == null) {
			replayer = LineageReplayer.withDefaultPool(diskLookup);
		}
		for (Path path : diskLookup.findFilesNonCommutative(a, b, c)) {
			int[] nat = shapeFromName(path.getFileName().toString());
			if (nat == null) {
				continue;
			}
			// Prefer the file's OWN lineage (the exact reproducible construction — what
			// materialise embeds); else pin by content hash; else the canonical-key atom.
			Lineage.Node leaf;
			try {
				Optional<Lineage.Node> own = SchemeIO.readLineage(path.toFile());
				if (own.isPresent()) {
					leaf = own.get();
				} else {
					String hash = SchemeIO.readHash(SchemeIO.parseJson(path.toFile()));
					leaf = (hash != null && hash.length() >= 7)
							? new Lineage.Atom(nat[0] + "x" + nat[1] + "x" + nat[2] + "@" + hash.substring(0, 7))
							: Lineage.atomFromFilename(path.getFileName().toString());
				}
			} catch (java.io.IOException e) {
				continue;
			}
			if (nat[0] != a || nat[1] != b || nat[2] != c) {
				leaf = Lineage.orientAs(leaf, nat[0], nat[1], nat[2], a, b, c);
			}
			try {
				NonCubicBilinearAlgorithm rep = replayer.replay(leaf);
				if (rep.n == a && rep.m == b && rep.p == c) {
					return leaf;
				}
			} catch (RuntimeException ignore) {
				// non-replayable NC candidate (e.g. non-bilinear file) — try the next
			}
		}
		return null;
	}


	private static final java.util.regex.Pattern SHAPE_REF_PAT =
			java.util.regex.Pattern.compile("(\\d+)x(\\d+)x(\\d+)");

	/** Write-time DAG corruption check for a stub's lineage. Returns a reason string when
	 *  the lineage would verify CORRUPT — a {@code shape@hash} ref that resolves to NO file
	 *  (DANGLING), or a ref that transitively reaches the WRITTEN shape (CYCLE) — else null.
	 *  Bounded transitive walk over {@code @hash} refs (the precise edges); bare/named refs
	 *  get only the cheap self-shape check (they resolve to catalog-best, not a fixed file). */
	// Package-private for white-box testing of the write-time guard (TestRecursiveMaterialiser).
	String lineageCorruption(Lineage.Node root, int n, int m, int p, String selfHash7) {
		// DFS over the @hash-pinned dependency closure. Reports the FIRST corruption:
		//   DANGLING — a shape@hash ref resolves to no file;
		//   CYCLE    — a back-edge (a ref reaches a file already on the DFS path, OR the
		//              new stub's own hash) → persisting it closes an INFINITE replay loop.
		// Detecting ANY closure cycle (not just the self-hash) is what catches a NEW stub
		// that INHERITS corruption by referencing a pre-existing cyclic stub (⟨10,11,20⟩ →
		// ⟨3,3,7⟩@h ← ⟨3,3,8⟩ ← ⟨3,3,7⟩@h). Bare/named refs resolve to catalog-best (not a
		// fixed file), so they are not followed — only @hash edges are precise.
		java.util.List<String> rootRefs = new java.util.ArrayList<>();
		collectAtomRefs(root, rootRefs);
		// SELF-SHAPE guard: a same-ORDERED-shape Atom leaf is a degenerate self-derivation —
		// "build ⟨n,m,p⟩ from a (different-hash) ⟨n,m,p⟩ scheme" via project-back/concat. It
		// TERMINATES (the atom is a concrete, acyclic scheme), so the @hash-cycle DFS below never
		// catches it — yet it can never be uniquely optimal (project-then-concat-back reconstructs
		// the removed slice, so rank ≥ the original). No legitimate op (concat/project/Kron/
		// recombination) ever yields a same-ORDERED-shape leaf; an orientation/transpose leaf
		// (same multiset, different order — e.g. ⟨2,3,2⟩ from ⟨2,2,3⟩) is fine and NOT flagged.
		// [[cyclic lineage → silent SOE]] — this is the terminating-cousin of that family.
		for (String ref : rootRefs) {
			// naive-NxMxP is the TERMINAL ground-truth leaf (n·m·p scalar products,
			// synthesised by trivialOneAxis, no catalog file) — same-shape by
			// construction and trivially acyclic, NOT a self-derivation. Flagging it
			// silently killed every Kron/concat build with a unit-axis factor
			// (⟨3,3,18⟩ = ⟨1,1,3⟩⊗⟨3,3,6⟩ predicted 120 then failed to materialise).
			// Guard: TestKroneckerOnlyStrategy.kron_only_handles_prime_axes_via_unit_factors.
			if (ref != null && ref.startsWith("naive-")) continue;
			int[] s = shapeOfRef(ref);
			if (s != null && s[0] == n && s[1] == m && s[2] == p) {
				return "SELF-SHAPE: lineage references an Atom of its own shape ⟨"
						+ n + "," + m + "," + p + "⟩ (" + ref + ")";
			}
		}
		java.util.Set<Path> onPath = new java.util.HashSet<>();
		java.util.Set<Path> cleanVisited = new java.util.HashSet<>();
		int[] budget = { 50_000 };
		for (String ref : rootRefs) {
			String c = dfsCorruption(ref, selfHash7, onPath, cleanVisited, budget);
			if (c != null) return c;
		}
		return null;
	}

	private String dfsCorruption(String ref, String selfHash7, java.util.Set<Path> onPath,
			java.util.Set<Path> cleanVisited, int[] budget) {
		if (budget[0]-- <= 0) return null; // bounded; treat over-budget as clean
		// An "@ref?:L0" Atom is the parse-time FALLBACK for an unresolvable internal lineage-id
		// (a self-referential / dedup-broken node — the ⟨17,19,20⟩ degenerate-memo cycle). It can
		// never resolve to a scheme, so a sub-block carrying it is corrupt and must be rejected (the
		// composition then uses a clean alternative instead of embedding the broken node).
		if (ref.startsWith("@ref?:")) return "FALLBACK-REF " + ref;
		int at = ref.indexOf('@');
		if (at <= 0) return null; // bare/named/naive — resolves to catalog-best, not a file edge
		String tail = ref.substring(at + 1);
		if (tail.equals("sota") || tail.equals("naive")) return null;
		String shape = ref.substring(0, at);
		String hash7 = tail.substring(0, Math.min(7, tail.length()));
		if (selfHash7 != null && hash7.equals(selfHash7)) {
			return "CYCLE: closure references the stub's own hash @" + hash7;
		}
		Path f = findFileByHash(shape, hash7);
		if (f == null) {
			return "DANGLING: ref " + ref + " resolves to no on-disk file";
		}
		if (onPath.contains(f)) {
			return "CYCLE: closure loops back to " + f.getFileName();
		}
		if (cleanVisited.contains(f)) return null; // sub-closure already proven clean
		onPath.add(f);
		for (String childRef : jsonRefs(f)) {
			String c = dfsCorruption(childRef, selfHash7, onPath, cleanVisited, budget);
			if (c != null) return c;
		}
		onPath.remove(f);
		cleanVisited.add(f);
		return null;
	}

	private static String sortedShapeKey(int n, int m, int p) {
		int[] d = { n, m, p };
		java.util.Arrays.sort(d);
		return d[0] + "x" + d[1] + "x" + d[2];
	}

	private static String sortedShapeKey(String shape) {
		String[] s = shape.split("x");
		return sortedShapeKey(Integer.parseInt(s[0]), Integer.parseInt(s[1]), Integer.parseInt(s[2]));
	}

	/** Collect every Atom ref in an in-memory lineage tree. */
	/** Ordered shape {@code [n,m,p]} parsed from a lineage Atom ref (e.g. {@code "2x2x3@h"},
	 *  {@code "naive-6x14x1"}), or {@code null} if the ref carries no NxMxP token. Reuses the
	 *  class-level {@link #REF_SHAPE} pattern. */
	private static int[] shapeOfRef(String ref) {
		if (ref == null) return null;
		java.util.regex.Matcher m = REF_SHAPE.matcher(ref);
		if (!m.find()) return null;
		return new int[] { Integer.parseInt(m.group(1)), Integer.parseInt(m.group(2)),
				Integer.parseInt(m.group(3)) };
	}

	private static void collectAtomRefs(Lineage.Node node, java.util.Collection<String> out) {
		if (node instanceof Lineage.Atom a) {
			out.add(a.ref());
		} else {
			for (Lineage.Node ch : Lineage.childrenOf(node)) collectAtomRefs(ch, out);
		}
	}

	/** True if any leaf is an unresolvable {@code "@ref?:Ln"} fallback Atom — the parse-time
	 *  placeholder for a cyclic/dedup-broken lineage node (never resolves to a scheme). */
	private static boolean hasFallbackRef(Lineage.Node node) {
		if (node instanceof Lineage.Atom a) {
			return a.ref() != null && a.ref().startsWith("@ref?:");
		}
		for (Lineage.Node ch : Lineage.childrenOf(node)) {
			if (hasFallbackRef(ch)) return true;
		}
		return false;
	}

	/** Every {@code ref} string in an on-disk scheme's JSON lineage. */
	private static java.util.List<String> jsonRefs(Path f) {
		java.util.List<String> out = new java.util.ArrayList<>();
		try {
			collectJsonRefs(SchemeIO.parseJson(f.toFile()).get("lineage"), out);
		} catch (Exception ignore) { /* unreadable → no edges (caught elsewhere) */ }
		return out;
	}

	private static void collectJsonRefs(tools.jackson.databind.JsonNode node, java.util.List<String> out) {
		if (node == null) return;
		if (node.isObject()) {
			tools.jackson.databind.JsonNode ref = node.get("ref");
			if (ref != null && ref.isString()) out.add(ref.asString());
			for (var e : node.properties()) collectJsonRefs(e.getValue(), out);
		} else if (node.isArray()) {
			for (tools.jackson.databind.JsonNode c : node) collectJsonRefs(c, out);
		}
	}

	/** The on-disk file at {@code shape} whose canonical name carries {@code hash7}, or null.
	 *  When several files share the hash (a derived DUPLICATE of a known/constructed/curated
	 *  scheme — same content, same hash), prefer the CANONICAL non-derived one: it is the
	 *  real leaf, so a degenerate {@code derived} identity-stub that references its own hash
	 *  ({@code ⟨2,2,2⟩=Recombination(base=2x2x2@<self>)}) does NOT masquerade as a self-loop. */
	private Path findFileByHash(String shape, String hash7) {
		String[] s = shape.split("x");
		int sn = Integer.parseInt(s[0]), sm = Integer.parseInt(s[1]), sp = Integer.parseInt(s[2]);
		Path derivedMatch = null;
		for (Path f : diskLookup.findFiles(sn, sm, sp)) {
			if (f.getFileName().toString().contains(hash7)) {
				if (!f.toString().contains("/derived/")) return f; // canonical leaf — prefer
				derivedMatch = f;
			}
		}
		if (derivedMatch != null) return derivedMatch;
		// Filename token is cosmetic — old-convention files (fmm-lille_5x7x7_r176_a3315)
		// carry the hash ONLY in the stamped JSON "hash" field, which is also what
		// LineageReplayer resolves pins by. Fall back to content so the audit agrees
		// with the replayer instead of flagging a resolvable pin as DANGLING.
		for (Path f : diskLookup.findFiles(sn, sm, sp)) {
			try {
				String stamped = SchemeIO.readHash(SchemeIO.parseJson(f.toFile()));
				if (stamped != null && stamped.startsWith(hash7)) return f;
			} catch (Exception e) {
				// unreadable candidate — keep scanning
			}
		}
		return null;
	}

	private void writeToDisk(int n, int m, int p, Result r) {
		int maxDim = Math.max(Math.max(n, m), p);
		// Everything this materialiser emits is DERIVED (regenerable) output, so it
		// always lands under the catalog's derived/ subtree — never top-level
		// schemes/sectionN/ (known/derived/curated split, 2026-06-09). writeRoot is
		// the schemes root; we inject derived/ here so no caller can misplace output.
		File dir = writeRoot.resolve("derived").resolve("section" + maxDim).toFile();
		dir.mkdirs();
		// Above the materialise cap, write a lineage-only stub (no factor
		// matrices) per CatalogPolicy — the explicit U/V/W are reproduced on
		// demand by LineageReplayer. Keeps the catalog small for large derived
		// shapes while still recording construction + rank.
		boolean stub = maxDim > eu.solven.matmul.catalog.CatalogPolicy.MATERIALISE_MAX_DIM
				&& r.lineage != null;
		// Canonical filename via the shared helper (2026-06 convention). Addition
		// count and bud-score live in JSON content, never the name. The content hash
		// dedups: identical (shape, rank, U/V/W) → identical name → exists()-skip below.
		// Note "derived" (not "derived_recursive") to match the established canonical
		// convention — the concrete mechanism (concat / Kron / recombination) lives in
		// lineage_compact, and RenameSchemes already folds derived_recursive→derived.
		String fn = SchemeIO.canonicalName(r.alg, "derived");
		File out = new File(dir, fn);
		if (out.exists()) {
			log.info("⟨{},{},{}⟩=r{} already on disk: {}", n, m, p, r.alg.r, out.getName());
			return;
		}
		// A stub IS only its lineage on disk, so an un-replayable one is a phantom catalog
		// entry. Verify replay BEFORE writing — this is the write-boundary guard that would
		// have caught the 2026-05-31 sorted-frame concat bug (operands recorded as the
		// sorted shape with no OrientAs → replay reconstructs a different shape; fixed in
		// materialise on 2026-06-04, but ~150 phantom stubs were already written). The
		// projection paths already discard un-replayable results before persist; this
		// extends the same discipline to the compose path (concat / Kron / recombination).
		if (stub && r.lineage != null && !replaysConsistently(r.lineage, r.alg)) {
			log.warn("REFUSING un-replayable STUB {} — lineage {} does not reconstruct ⟨{},{},{}⟩=r{};"
					+ " discarding rather than write a phantom catalog entry.",
					fn, Lineage.prettyString(r.lineage), n, m, p, r.alg.r);
			return;
		}
		// WRITE-TIME DAG GUARD: replaysConsistently passes even for a lineage whose @hash
		// refs DANGLE (resolveLeaf silently falls back to shape-best) or form a cross-file
		// CYCLE (⟨3,3,7⟩←⟨3,3,8⟩←⟨3,3,7⟩). Both produce a catalog VerifyScheme flags corrupt.
		// Refuse to write either: a stub must have every ref resolve EXACTLY and never reach
		// its own shape transitively. (project_findrank_poisoned… / cyclic-lineage class.)
		if (r.lineage != null) {
			String selfHash7 = SchemeIO.contentHash(r.alg).substring(0, 7);
			String corruption = lineageCorruption(r.lineage, n, m, p, selfHash7);
			if (corruption != null) {
				log.warn("REFUSING corrupt STUB {} — {} ; lineage={}. Discarding rather than"
						+ " write a catalog entry that would verify CORRUPT.",
						fn, corruption, Lineage.prettyString(r.lineage));
				return;
			}
		}
		// Field set = INTERSECTION of the leaf atoms' stamped fields[] (content-only). A
		// derived scheme without fields[] is INVISIBLE to every field-aware query, so BOTH
		// stubs AND explicit schemes must be born-stamped — the explicit branch used to skip
		// this, leaving all ≤16 derived field-invisible (the ⟨16,16,16⟩=2304 "vanished" bug).
		// Empty means a leaf the build resolved can't be re-resolved → FAIL LOUD rather than
		// ship an over-claiming [Z] default or an invisible unstamped scheme.
		java.util.List<String> fields = diskLookup.fieldNamesFromLineage(r.lineage);
		if (fields.isEmpty()) {
			throw new IllegalStateException("field inference returned EMPTY for " + fn
					+ " — a leaf could not be resolved though the build resolved it. lineage="
					+ Lineage.prettyString(r.lineage));
		}
		// COEFFICIENT OVERRIDE: fieldNamesFromLineage intersects the leaves' fields,
		// assuming the composition preserves their ring — coefficient-BLIND. A
		// recombination / ½-polarization step (e.g. an ABC->ACA projection) divides
		// out of that ring, so the materialized matrices can hold a 1/8 that excludes
		// Z/F2 even though every integer leaf claimed them. Narrow the inferred set to
		// what the ACTUAL coefficients support, so we never BORN-STAMP an over-claim
		// the field gate would later reject (project-derived-schemes-overclaim-z).
		java.util.List<String> narrowed = SchemeIO.narrowFieldsToCoefficients(r.alg, fields);
		if (narrowed.size() != fields.size()) {
			log.info("narrowed fields {} → {} for {} (coefficients exclude the dropped tag(s); "
					+ "lineage inference was ring-blind)", fields, narrowed, fn);
			fields = narrowed;
		}
		try {
			if (stub) {
				SchemeIO.writeStub(r.alg, out, r.lineage, fields);
				log.info("wrote STUB {} (r={}, fields={}, lineage={})", out.getName(), r.alg.r,
						fields, Lineage.prettyString(r.lineage));
			} else {
				SchemeIO.write(r.alg, out, r.lineage);
				// Stamp fields[] onto the just-written explicit scheme (write() doesn't).
				SchemeIO.updateFields(out, java.util.Map.of("fields", fields),
						java.util.List.of(), true);
				log.info("wrote {} (r={}, fields={}, lineage={})", out.getName(), r.alg.r,
						fields, Lineage.prettyString(r.lineage));
			}
		} catch (IOException e) {
			log.warn("write failed for ⟨{},{},{}⟩: {}", n, m, p, e.getMessage());
		}
	}

	/**
	 * Naïve scheme for any {@code ⟨n,m,p⟩} where at least one axis equals 1.
	 * Rank is exactly the product of the two non-unit dimensions (e.g.
	 * {@code ⟨1,m,p⟩} = m·p scalar multiplications). Used as a leaf when
	 * a concat-search picks a width-1 split.
	 */
	private Result trivialOneAxis(int n, int m, int p) {
		int r = n * m * p;  // for n=1 or m=1 or p=1 this is the non-unit product
		double[][] U = new double[n * m][r];
		double[][] V = new double[m * p][r];
		double[][] W = new double[n * p][r];
		int k = 0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				for (int l = 0; l < p; l++) {
					U[i * m + j][k] = 1;
					V[j * p + l][k] = 1;
					W[i * p + l][k] = 1;
					k++;
				}
			}
		}
		NonCubicBilinearAlgorithm alg = new NonCubicBilinearAlgorithm(n, m, p, U, V, W);
		Lineage.Node leaf = new Lineage.Atom(String.format("naive-%dx%dx%d", n, m, p));
		return new Result(alg, leaf, false);
	}

	private static String canon(int n, int m, int p) {
		int[] s = { n, m, p };
		Arrays.sort(s);
		return s[0] + "x" + s[1] + "x" + s[2];
	}

	/** Native shape ⟨n,m,p⟩ from a scheme filename's {@code NxMxP} token, or null. */
	private static final java.util.regex.Pattern SHAPE_IN_NAME =
			java.util.regex.Pattern.compile("(\\d+)x(\\d+)x(\\d+)");

	private static int[] shapeFromName(String name) {
		java.util.regex.Matcher mm = SHAPE_IN_NAME.matcher(name);
		if (!mm.find()) return null;
		return new int[] { Integer.parseInt(mm.group(1)), Integer.parseInt(mm.group(2)),
				Integer.parseInt(mm.group(3)) };
	}

	/**
	 * Lookup overlay that wraps {@link #diskLookup} but recurses into
	 * {@link #materialise} on miss. Records lineage nodes in insertion
	 * order so the outer recombination can attach them as leaves.
	 */
	private final class RecursiveLookup implements Recombination.AlgorithmLookup {
		final List<Lineage.Node> lineageOrder = new ArrayList<>();
		final Map<String, Lineage.Node> seen = new HashMap<>();

		@Override
		public Optional<NonCubicBilinearAlgorithm> find(int n, int m, int p) {
			String key = canon(n, m, p);
			// LOAD by default (no implicit recursion) — see resolveSubScheme / recursiveDerive.
			Optional<Result> r = recursiveDerive ? materialise(n, m, p) : diskBest(n, m, p);
			if (r.isEmpty()) return Optional.empty();
			Lineage.Node leaf = r.get().lineage;
			// CYCLE PREVENTION: materialise() can return a catalog-best block whose lineage carries
			// INHERITED cycle corruption (a legacy 565-class stub — an unresolvable "@ref?:L0" fallback
			// from a projection cycle). This overlay path (unlike resolveSubScheme) would EMBED that
			// broken lineage into the parent, spreading the corruption. Instead, pin the block by its
			// content hash — the parent stays clean, replay loads the block by hash, and the corruption
			// is contained to the (separately purged) legacy scheme. [[cyclic lineage → silent SOE]]
			if (leaf != null && hasFallbackRef(leaf)) {
				String pinned = n + "x" + m + "x" + p + "@" + SchemeIO.contentHash(r.get().alg);
				log.warn("recombination leaf ⟨{},{},{}⟩ had corrupt (cyclic) lineage — pinned by hash "
						+ "instead of embedding ({})", n, m, p, pinned);
				leaf = new Lineage.Atom(pinned);
			}
			if (seen.putIfAbsent(key, leaf) == null) {
				lineageOrder.add(leaf);
			}
			return Optional.of(r.get().alg);
		}
	}
}
