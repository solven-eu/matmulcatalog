package eu.solven.matmul.search;

import eu.solven.matmul.SymmetryTransforms.InternalOrbitMode;

/**
 * Configuration for the outer-template base pool used by
 * {@link BlockSplitSearch#buildPool(PoolConfig)}. Three independent
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
public record PoolConfig(
		boolean cubicOnly,
		boolean rootsOnly,
		InternalOrbitMode orbitMode,
		int maxBaseDim,
		long permutationOrbitCap,
		boolean balancedOnly,
		int maxImbalance,
		int maxCombinations,
		int maxPadding) {

	/** Sentinel for "no imbalance cap"; equivalent to fully-unbalanced search. */
	public static final int UNBOUNDED_IMBALANCE = Integer.MAX_VALUE;

	/** Sentinel for "no per-base combination cap". */
	public static final int UNBOUNDED_COMBINATIONS = Integer.MAX_VALUE;

	/** Default {@code maxPadding} = 0: no over-allocation (no DIS09 γ5 peel pattern tried). */
	public static final int NO_PADDING = 0;

	public PoolConfig {
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
	public PoolConfig(boolean cubicOnly, boolean rootsOnly, InternalOrbitMode orbitMode,
			int maxBaseDim, long permutationOrbitCap) {
		this(cubicOnly, rootsOnly, orbitMode, maxBaseDim, permutationOrbitCap,
				false, UNBOUNDED_IMBALANCE, UNBOUNDED_COMBINATIONS, NO_PADDING);
	}

	/** 7-arg back-compat (pre-maxCombinations / maxPadding). */
	public PoolConfig(boolean cubicOnly, boolean rootsOnly, InternalOrbitMode orbitMode,
			int maxBaseDim, long permutationOrbitCap,
			boolean balancedOnly, int maxImbalance) {
		this(cubicOnly, rootsOnly, orbitMode, maxBaseDim, permutationOrbitCap,
				balancedOnly, maxImbalance, UNBOUNDED_COMBINATIONS, NO_PADDING);
	}

	/** Returns a copy of {@code this} with {@code balancedOnly} replaced. */
	public PoolConfig withBalancedOnly(boolean v) {
		return new PoolConfig(cubicOnly, rootsOnly, orbitMode, maxBaseDim, permutationOrbitCap, v, maxImbalance, maxCombinations, maxPadding);
	}

	/** Returns a copy of {@code this} with {@code maxImbalance} replaced. */
	public PoolConfig withMaxImbalance(int v) {
		return new PoolConfig(cubicOnly, rootsOnly, orbitMode, maxBaseDim, permutationOrbitCap, balancedOnly, v, maxCombinations, maxPadding);
	}

	public PoolConfig withMaxBaseDim(int v) {
		return new PoolConfig(cubicOnly, rootsOnly, orbitMode, v, permutationOrbitCap, balancedOnly, maxImbalance, maxCombinations, maxPadding);
	}

	public PoolConfig withOrbitMode(InternalOrbitMode v) {
		return new PoolConfig(cubicOnly, rootsOnly, v, maxBaseDim, permutationOrbitCap, balancedOnly, maxImbalance, maxCombinations, maxPadding);
	}

	public PoolConfig withMaxCombinations(int v) {
		return new PoolConfig(cubicOnly, rootsOnly, orbitMode, maxBaseDim, permutationOrbitCap, balancedOnly, maxImbalance, v, maxPadding);
	}

	public PoolConfig withMaxPadding(int v) {
		return new PoolConfig(cubicOnly, rootsOnly, orbitMode, maxBaseDim, permutationOrbitCap, balancedOnly, maxImbalance, maxCombinations, v);
	}

	/**
	 * "DIS-like simple default": cubic, root-only, canonical orbit,
	 * cap 5 (covers ⟨2,2,2⟩ through ⟨5,5,5⟩). The single recommended
	 * starting point for routine sweeps; closes most catalog gaps with
	 * sub-15-minute wall-clock on ~200 shapes.
	 */
	public static PoolConfig simple() {
		return new PoolConfig(true, true, InternalOrbitMode.CANONICAL, 5, 8);
	}

	/**
	 * Same as simple but with axis-flip orbit expansion. ×4-8 pool size,
	 * proportionally slower. Marks lineage as
	 * {@link eu.solven.matmul.catalog.Lineage.AxisFlip approximate}.
	 */
	public static PoolConfig auditAxisFlip() {
		return new PoolConfig(true, true, InternalOrbitMode.AXIS_FLIP, 5, 8);
	}

	/** Cubic, root, full permutation orbit bounded at 216 (tractable through ⟨3,3,3⟩). */
	public static PoolConfig auditPermutation() {
		return new PoolConfig(true, true, InternalOrbitMode.PERMUTATION_BOUNDED, 5, 216);
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
	public static PoolConfig axisFlipOnly() {
		return new PoolConfig(true, true, InternalOrbitMode.AXIS_FLIP_ONLY, 5, 8);
	}

	/** Allow non-cubic roots (rectangular AT⟨2,2,3⟩, etc.). Canonical orbit. */
	public static PoolConfig rectangular() {
		return new PoolConfig(false, true, InternalOrbitMode.CANONICAL, 5, 8);
	}

	/** Include derived catalog schemes — the full "extendedPool" universe at cap 5. */
	public static PoolConfig includeDerived() {
		return new PoolConfig(false, false, InternalOrbitMode.CANONICAL, 5, 8);
	}

	/**
	 * Most thorough realistic config: rectangular + derived + axis-flip.
	 * Permutation NOT enabled here because at maxBaseDim=5 the pool
	 * would explode past the cap on every ⟨4,4,4⟩ / ⟨5,5,5⟩ scheme
	 * and silently fall back to axis-flip anyway. For true permutation
	 * sweeps, lower {@code maxBaseDim} to 3.
	 */
	public static PoolConfig thorough() {
		return new PoolConfig(false, false, InternalOrbitMode.AXIS_FLIP, 5, 8);
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
