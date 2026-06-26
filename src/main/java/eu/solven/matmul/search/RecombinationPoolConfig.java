package eu.solven.matmul.search;

import eu.solven.matmul.SymmetryTransforms.InternalOrbitMode;

/**
 * Configuration for the outer-template base pool used by
 * {@link BlockSplitSearch#buildPool(RecombinationPoolConfig)}. Three independent
 * axes — {@code (cubicOnly, rootsOnly, orbitMode)} — span the 3·3 = 9
 * meaningful tradeoff points between "fast default sweep" and
 * "thorough audit". Defaults match the user-stated "DIS-like simple
 * default": cubic outer templates, genuine roots only, no algorithmic
 * orbit expansion (canonical shape only).
 *
 * <p>For routine use call {@link #simple()}; for audit sweeps see the
 * other factory methods or supply your own combination via the
 * canonical constructor.</p>
 *
 * <h2>Axes</h2>
 * <ul>
 *   <li>{@code cubicOnly} — outer template must be ⟨n,n,n⟩.
 *       Excludes AT⟨2,2,3⟩, AT⟨2,3,3⟩, AT⟨2,3,4⟩, etc. Filters at the
 *       root node only; sub-shape SOTA lookups still use the full
 *       catalog (cost lookup, not recursion).</li>
 *   <li>{@code rootsOnly} — outer template must be Leaf-lineage
 *       (genuine root scheme, not derived via composition). Excludes
 *       {@code derived-recursive_*}, {@code solven-strassen-*},
 *       {@code dis09-Q_*}, etc. Same "filter at root node only"
 *       semantics.</li>
 *   <li>{@code orbitMode} — see {@link InternalOrbitMode}.</li>
 * </ul>
 *
 * <h2>Cost</h2>
 * <p>Pool size grows multiplicatively with each axis loosened:</p>
 * <ul>
 *   <li>simple = O(10) outer templates.</li>
 *   <li>cubicOnly=false → O(30) (add rectangular AT roots).</li>
 *   <li>rootsOnly=false → O(thousands at maxBaseDim=5).</li>
 *   <li>orbitMode=AXIS_FLIP → ×4-8 per template (axis reversal).</li>
 *   <li>orbitMode=PERMUTATION_BOUNDED → ×24-216 per template at
 *       maxBaseDim=3, intractable above.</li>
 * </ul>
 * Recombination search cost is roughly linear in pool size, so a
 * loose config can multiply wall-clock by 100× or more.
 */
public record RecombinationPoolConfig(
		boolean cubicOnly,
		boolean rootsOnly,
		InternalOrbitMode orbitMode,
		int maxBaseDim,
		long permutationOrbitCap,
		boolean balancedOnly,
		int maxImbalance,
		int maxCombinations,
		int maxPadding,
			int maxCubicBaseDim) {

	/** Sentinel for "no imbalance cap"; equivalent to fully-unbalanced search. */
	public static final int UNBOUNDED_IMBALANCE = Integer.MAX_VALUE;

	/** Sentinel for "no per-base combination cap". */
	public static final int UNBOUNDED_COMBINATIONS = Integer.MAX_VALUE;

	/** Default {@code maxPadding} = 0: no over-allocation (no DIS09 γ5 peel pattern tried). */
	public static final int NO_PADDING = 0;

	public RecombinationPoolConfig {
		if (maxBaseDim < 2) throw new IllegalArgumentException("maxBaseDim must be ≥ 2, got " + maxBaseDim);
		if (permutationOrbitCap < 8) throw new IllegalArgumentException(
				"permutationOrbitCap must be ≥ 8 (axis-flip subset size), got " + permutationOrbitCap);
		if (orbitMode == null) throw new IllegalArgumentException("orbitMode required");
		if (maxImbalance < 0) throw new IllegalArgumentException("maxImbalance must be ≥ 0, got " + maxImbalance);
		if (maxCombinations < 1) throw new IllegalArgumentException("maxCombinations must be ≥ 1, got " + maxCombinations);
		if (maxPadding < 0) throw new IllegalArgumentException("maxPadding must be ≥ 0, got " + maxPadding);
	}

	/**
	 * Back-compat 5-arg constructor: defaults the new search-shape knobs
	 * to {@code balancedOnly = false, maxImbalance = UNBOUNDED_IMBALANCE}
	 * — i.e. the broadest legal search. Callers that want the old
	 * "balanced only" behaviour must opt in explicitly.
	 *
	 * <p>Rationale: silent pruning of the search space is a footgun
	 * (see project memory feedback_dont_silently_prune_search_space).
	 * Pruning heuristics must be opt-in, not the default.</p>
	 */
	public RecombinationPoolConfig(boolean cubicOnly, boolean rootsOnly, InternalOrbitMode orbitMode,
			int maxBaseDim, long permutationOrbitCap) {
		this(cubicOnly, rootsOnly, orbitMode, maxBaseDim, permutationOrbitCap,
				false, UNBOUNDED_IMBALANCE, UNBOUNDED_COMBINATIONS, NO_PADDING, maxBaseDim);
	}

	/** 7-arg back-compat (pre-maxCombinations / maxPadding). */
	public RecombinationPoolConfig(boolean cubicOnly, boolean rootsOnly, InternalOrbitMode orbitMode,
			int maxBaseDim, long permutationOrbitCap,
			boolean balancedOnly, int maxImbalance) {
		this(cubicOnly, rootsOnly, orbitMode, maxBaseDim, permutationOrbitCap,
				balancedOnly, maxImbalance, UNBOUNDED_COMBINATIONS, NO_PADDING, maxBaseDim);
	}

	/** Returns a copy of {@code this} with {@code balancedOnly} replaced. */
	public RecombinationPoolConfig withBalancedOnly(boolean v) {
		return new RecombinationPoolConfig(cubicOnly, rootsOnly, orbitMode, maxBaseDim, permutationOrbitCap, v, maxImbalance, maxCombinations, maxPadding, maxCubicBaseDim);
	}

	/** Returns a copy of {@code this} with {@code maxImbalance} replaced. */
	public RecombinationPoolConfig withMaxImbalance(int v) {
		return new RecombinationPoolConfig(cubicOnly, rootsOnly, orbitMode, maxBaseDim, permutationOrbitCap, balancedOnly, v, maxCombinations, maxPadding, maxCubicBaseDim);
	}

	public RecombinationPoolConfig withMaxBaseDim(int v) {
		return new RecombinationPoolConfig(cubicOnly, rootsOnly, orbitMode, v, permutationOrbitCap, balancedOnly, maxImbalance, maxCombinations, maxPadding, maxCubicBaseDim);
	}

	public RecombinationPoolConfig withOrbitMode(InternalOrbitMode v) {
		return new RecombinationPoolConfig(cubicOnly, rootsOnly, v, maxBaseDim, permutationOrbitCap, balancedOnly, maxImbalance, maxCombinations, maxPadding, maxCubicBaseDim);
	}

	public RecombinationPoolConfig withMaxCombinations(int v) {
		return new RecombinationPoolConfig(cubicOnly, rootsOnly, orbitMode, maxBaseDim, permutationOrbitCap, balancedOnly, maxImbalance, v, maxPadding, maxCubicBaseDim);
	}

	public RecombinationPoolConfig withMaxPadding(int v) {
		return new RecombinationPoolConfig(cubicOnly, rootsOnly, orbitMode, maxBaseDim, permutationOrbitCap, balancedOnly, maxImbalance, maxCombinations, v, maxCubicBaseDim);
	}

	/**
	 * "DIS-like simple default": cubic, root-only, canonical orbit,
	 * cap 5 (covers ⟨2,2,2⟩ through ⟨5,5,5⟩). The single recommended
	 * starting point for routine sweeps; closes most catalog gaps with
	 * sub-15-minute wall-clock on ~200 shapes.
	 */
	/** Copy with {@code maxCubicBaseDim} replaced — the dim ceiling for CUBIC bases, which
	 *  may exceed {@code maxBaseDim} so cubic bases (⟨6,6,6⟩…⟨8,8,8⟩, e.g. Sedoglavic
	 *  ⟨7,7,7⟩) stay in a pool whose NON-cubic bases are capped lower. */
	public RecombinationPoolConfig withMaxCubicBaseDim(int v) {
		return new RecombinationPoolConfig(cubicOnly, rootsOnly, orbitMode, maxBaseDim, permutationOrbitCap,
				balancedOnly, maxImbalance, maxCombinations, maxPadding, v);
	}

	/**
	 * The sensible-but-large DEFAULT recombination pool (user 2026-06-26): non-cubic bases
	 * up to {@code maxBaseDim=5}, cubic bases up to {@code maxCubicBaseDim=8} (so ⟨6,6,6⟩…
	 * ⟨7,7,7⟩ Sedoglavic stay in), {@code cubicOnly=false}, extended catalog included,
	 * {@code AXIS_FLIP} orbit. Replaces {@link #simple()} ({@code cubicOnly=true}) as the
	 * materialize default: simple() silently dropped EVERY non-cubic base (the AxisSplits
	 * and the naïve grids), which regressed the fullDerive's non-cubic recombinations
	 * (⟨7,7,30⟩, ⟨13,13,32⟩, and the Pan-TA ⟨26,29,29⟩/⟨28,31,31⟩). Users narrow it
	 * ({@code --cubicOnly}, smaller {@code --maxBaseDim}) or widen it for specific bigger bases.
	 *
	 * <p>Orbit is {@code AXIS_FLIP} (the serious-sweep default): the flip orbit explores
	 * each base's axis ORIENTATIONS, which is exactly what the per-base frontier-multiset
	 * mechanism (Phase 2) will make efficient — iterating the frontier instead of re-deriving
	 * every flipped base. {@code CANONICAL} is an A/B-test edge case (fast, single
	 * orientation), useful to check whether the flip orbit adds anything, but NOT a default
	 * for any serious sweep — use {@link #withOrbitMode}{@code (CANONICAL)} to probe.</p>
	 */
	public static RecombinationPoolConfig defaultLarge() {
		return new RecombinationPoolConfig(false, false, InternalOrbitMode.AXIS_FLIP, 5, 8,
				false, UNBOUNDED_IMBALANCE, UNBOUNDED_COMBINATIONS, NO_PADDING, 8);
	}

	public static RecombinationPoolConfig simple() {
		return new RecombinationPoolConfig(true, true, InternalOrbitMode.CANONICAL, 5, 8);
	}

	/**
	 * Same as simple but with axis-flip orbit expansion. ×4-8 pool size,
	 * proportionally slower. Marks lineage as
	 * {@link eu.solven.matmul.catalog.Lineage.AxisFlip approximate}.
	 */
	public static RecombinationPoolConfig auditAxisFlip() {
		return new RecombinationPoolConfig(true, true, InternalOrbitMode.AXIS_FLIP, 5, 8);
	}

	/** Cubic, root, full permutation orbit bounded at 216 (tractable through ⟨3,3,3⟩). */
	public static RecombinationPoolConfig auditPermutation() {
		return new RecombinationPoolConfig(true, true, InternalOrbitMode.PERMUTATION_BOUNDED, 5, 216);
	}

	/**
	 * Same as {@link #simple()} but with {@link InternalOrbitMode#AXIS_FLIP_ONLY} —
	 * the canonical scheme is EXCLUDED from the pool, leaving only the
	 * 7 axis-flipped variants. Used as an A/B probe: if a search under
	 * this config matches or beats the {@link #simple()} rank, the
	 * canonical scheme is contributing nothing the flipped variants don't
	 * already cover. If it's strictly worse, the canonical scheme is the
	 * essential one and the orbit expansion is decorative.
	 */
	public static RecombinationPoolConfig axisFlipOnly() {
		return new RecombinationPoolConfig(true, true, InternalOrbitMode.AXIS_FLIP_ONLY, 5, 8);
	}

	/** Allow non-cubic roots (rectangular AT⟨2,2,3⟩, etc.). Canonical orbit. */
	public static RecombinationPoolConfig rectangular() {
		return new RecombinationPoolConfig(false, true, InternalOrbitMode.CANONICAL, 5, 8);
	}

	/** Include derived catalog schemes — the full "extendedPool" universe at cap 5. */
	public static RecombinationPoolConfig includeDerived() {
		return new RecombinationPoolConfig(false, false, InternalOrbitMode.CANONICAL, 5, 8);
	}

	/**
	 * Most thorough realistic config: rectangular + derived + axis-flip.
	 * Permutation NOT enabled here because at maxBaseDim=5 the pool
	 * would explode past the cap on every ⟨4,4,4⟩ / ⟨5,5,5⟩ scheme
	 * and silently fall back to axis-flip anyway. For true permutation
	 * sweeps, lower {@code maxBaseDim} to 3.
	 */
	public static RecombinationPoolConfig thorough() {
		return new RecombinationPoolConfig(false, false, InternalOrbitMode.AXIS_FLIP, 5, 8);
	}

	/** A compact handle for logs / report headers. */
	public String shortLabel() {
		String alloc = balancedOnly ? "bal"
				: (maxImbalance == UNBOUNDED_IMBALANCE ? "unbal" : "imb" + maxImbalance);
		return (cubicOnly ? "cub" : "any")
				+ "+" + (rootsOnly ? "root" : "all")
				+ "+" + orbitMode.name().toLowerCase().replace('_', '-')
				+ "+d" + maxBaseDim
				+ "+" + alloc;
	}
}
