package eu.solven.matmul.search.flip;

import eu.solven.matmul.recombination.Recombination;

import java.util.HashMap;
import java.util.Map;

/**
 * Stock {@link FlipObjective}s plus the native integer re-implementations of
 * the two catalog structure metrics they optimize. The native scores MUST agree
 * with the catalog reference implementations —
 * {@code SerendipitousBudProduct.independentClassSizes} (Σ class sizes ≥2) and
 * {@code ProjectionSearch.projectionMargin} — which {@code TestFlipObjectives}
 * cross-checks on real catalog schemes. Both metrics are HIGHER-is-better in
 * the catalog, so the costs below negate them under a rank-first lexicographic
 * key: {@code cost = rank·SCALE − score}.
 */
public final class FlipObjectives {

	private FlipObjectives() {}

	/**
	 * Lexicographic headroom for the structure score. Both scores are bounded by
	 * {@code 3·rank} (budScore) and {@code rank} (margin), far below this, so a
	 * rank drop always dominates any structure gain.
	 */
	static final long SCALE = 1_000_000L;

	/** Pure rank descent — the classic Kauers–Moosbauer objective. */
	public static FlipObjective minRank() {
		return new FlipObjective() {
			@Override
			public long cost(FlipScheme s) {
				return s.rank();
			}

			@Override
			public String describe() {
				return "rank↓";
			}
		};
	}

	/** Rank first, then maximize bud-richness (serendipitous-product fuel). */
	public static FlipObjective maxBudScore() {
		return new FlipObjective() {
			@Override
			public long cost(FlipScheme s) {
				return s.rank() * SCALE - budScore(s);
			}

			@Override
			public String describe() {
				return "rank↓ then budScore↑";
			}
		};
	}

	/** Rank first, then maximize the projection margin μ (downward-parent strength). */
	public static FlipObjective maxProjectionMargin() {
		return new FlipObjective() {
			@Override
			public long cost(FlipScheme s) {
				return s.rank() * SCALE - projectionMargin(s);
			}

			@Override
			public String describe() {
				return "rank↓ then projectionMargin↑";
			}
		};
	}

	/**
	 * Linear tradeoff cost — rank is NOT lexicographic here: it trades against
	 * structure at user-chosen exchange rates,
	 * {@code cost = rank·wRank − budScore·wBud − margin·wMargin}.
	 * E.g. {@code weighted(10, 1, 0)} accepts rank+1 whenever it buys ≥10 bud
	 * points (the "rank+1 but bud−10-cost" regime); {@code weighted(0, 1, 0)}
	 * ignores rank entirely (pair with a tight {@code maxRankAbove} budget, or
	 * the walk will happily inflate rank for structure). Weights must not be
	 * all zero.
	 */
	public static FlipObjective weighted(long wRank, long wBud, long wMargin) {
		if (wRank == 0 && wBud == 0 && wMargin == 0) {
			throw new IllegalArgumentException("at least one weight must be nonzero");
		}
		return new FlipObjective() {
			@Override
			public long cost(FlipScheme s) {
				long c = (long) s.rank() * wRank;
				if (wBud != 0) {
					c -= (long) budScore(s) * wBud;
				}
				if (wMargin != 0) {
					c -= (long) projectionMargin(s) * wMargin;
				}
				return c;
			}

			@Override
			public String describe() {
				return "weighted(rank·" + wRank + " − bud·" + wBud + " − margin·" + wMargin + ")";
			}
		};
	}

	/**
	 * DIRECT serendipitous-product cost (the real currency, not the budScore
	 * proxy): {@code cost = SerendipitousBudProduct.serendipitousCost(s, lookup,
	 * n2, m2, p2)} — the predicted rank of {@code s ⊗ ⟨n2,m2,p2⟩} under the
	 * catalog rank oracle. The 2026-06-11 probe ({@code ProbeFlipBudHarvest})
	 * showed budScore maximization WORSENS this (an extra term costs a full
	 * R(inner) while a size-2 bud saves only {@code 2·R(inner) − R(enlarged)}),
	 * so harvest walks must optimize this objective directly.
	 *
	 * <p>Each evaluation re-runs the greedy decomposition over all 6 orderings
	 * (O(r·dim) hashing ×18) plus rank-oracle lookups — fine at small shapes,
	 * the incremental version is the known phase-3 perf follow-up. Unknown
	 * enlarged-inner ranks surface as a huge cost, so the walk simply avoids
	 * bud sizes the catalog cannot price.</p>
	 */
	public static FlipObjective serendipitous(
			eu.solven.matmul.catalog.FieldAwareLookup lookup, int n2, int m2, int p2) {
		return new FlipObjective() {
			@Override
			public long cost(FlipScheme s) {
				return eu.solven.matmul.catalog.SerendipitousBudProduct
						.serendipitousCost(s.toAlgorithm(), lookup, n2, m2, p2);
			}

			@Override
			public String describe() {
				return "serendipitousCost(inner ⟨" + n2 + "," + m2 + "," + p2 + "⟩)↓";
			}
		};
	}

	/**
	 * Self-serendipitous cost — the canonical, parameter-free "serendipity
	 * potential" of a base: {@link #serendipitous} with the inner shape equal to
	 * the scheme's OWN shape, i.e. the predicted rank of {@code base ⊗ base}
	 * (⟨n²,m²,p²⟩ — the quantity that drives recursive squaring and ω). This is
	 * the principled form of "bud structure summed over axes, minus rank": by
	 * Smith 2002 eq. (69) as generalised in the paper (§serendipitous),
	 * <pre>
	 *   r_s = Σ_buds R(enlargedᵢ) + trivial·R(inner)
	 *       = R_base·R(inner) − Σ_buds σ(kᵢ),   σ(k) = k·R(inner) − R(k-enlarged)
	 * </pre>
	 * so minimizing {@code r_s} ≡ maximizing (priced bud savings − rank·R(inner)).
	 * Rank carries the exchange rate {@code R(inner)} per unit — for inner
	 * ⟨3,3,3⟩ that is 23 per rank vs σ(2)=2·23−R⟨3,3,6⟩=6 per size-2 bud, which
	 * is exactly why raw budScore maximization loses (ProbeFlipBudHarvest).
	 * Works across formats (MetaFlipWalk): the inner follows the current shape.
	 */
	public static FlipObjective selfSerendipitous(
			eu.solven.matmul.catalog.FieldAwareLookup lookup) {
		return new FlipObjective() {
			@Override
			public long cost(FlipScheme s) {
				return eu.solven.matmul.catalog.SerendipitousBudProduct
						.serendipitousCost(s.toAlgorithm(), lookup, s.n, s.m, s.p);
			}

			@Override
			public String describe() {
				return "selfSerendipitousCost(base ⊗ base)↓";
			}
		};
	}

	/**
	 * DIRECT projected cost against a concrete target ⟨tn,tm,tp⟩ (the projection
	 * analogue of {@link #serendipitous}): {@code cost = exact projected rank} =
	 * {@code rank − μ(best drop combo)} via
	 * {@code ProjectionSearch.projectedRank}. The real currency for "is this base
	 * a good projecting parent for THAT target". Rank's exchange rate against
	 * margin is exactly 1 (each base rank point costs one projected point), so —
	 * unlike serendipity (rate {@code R(inner)}) — there is NO a-priori ceiling
	 * on profitable rank relaxation: rank+δ pays whenever it buys margin {@code
	 * > δ}, a purely realizability-bound question. Walks at {@code maxRankAbove >
	 * 0} are therefore the natural mode here. Schemes that don't reach the
	 * target (shape too small / over-delta) cost a huge constant, steering the
	 * walk back.
	 */
	public static FlipObjective projectedTo(int tn, int tm, int tp, int maxDelta) {
		return new FlipObjective() {
			@Override
			public long cost(FlipScheme s) {
				long r = eu.solven.matmul.catalog.ProjectionSearch
						.projectedRank(s.toAlgorithm(), tn, tm, tp, maxDelta);
				return r < 0 ? Long.MAX_VALUE / 8 : r;
			}

			@Override
			public String describe() {
				return "projectedRank(→⟨" + tn + "," + tm + "," + tp + "⟩)↓";
			}
		};
	}

	/**
	 * The absolute serendipitous SAVING of a base vs the plain Kronecker
	 * product with {@code ⟨n2,m2,p2⟩}: {@code Σ_buds σ(kᵢ) = R_base·R(inner) −
	 * r_s} — the "priced bud structure" term of the potential, for reporting.
	 * 0 means the bud structure buys nothing at this inner (e.g. ANY size-2 bud
	 * with inner ⟨2,2,2⟩: σ(2)=2·7−R⟨2,2,4⟩=14−14=0). A bound, like everything
	 * the greedy decomposition produces.
	 */
	public static long serendipitySaving(FlipScheme s,
			eu.solven.matmul.catalog.FieldAwareLookup lookup, int n2, int m2, int p2) {
		long inner = lookup.findRank(n2, m2, p2);
		long cost = eu.solven.matmul.catalog.SerendipitousBudProduct
				.serendipitousCost(s.toAlgorithm(), lookup, n2, m2, p2);
		if (inner >= Recombination.SotaResolver.UNKNOWN_RANK || cost >= Long.MAX_VALUE / 8) {
			return 0;
		}
		return (long) s.rank() * inner - cost;
	}

	/**
	 * Per-axis serendipity profile {@code {σ_U, σ_V, σ_W}} — the priced bud
	 * savings split by bud type. This is the multi-dimensional form of the
	 * potential (user 2026-06-11): U-buds amplify the inner's {@code p}-axis,
	 * V-buds its {@code n}-axis, W-buds its {@code m}-axis, so WHICH axis's
	 * richness pays off depends on the partner scheme — one very-rich axis and
	 * two moderately-rich axes are different goods, and a scalar total only
	 * ranks bases for one fixed inner. Computed from the default greedy
	 * decomposition (a bound; ordering-sensitive like every greedy bud figure).
	 */
	public static long[] serendipitySavingByAxis(FlipScheme s,
			eu.solven.matmul.catalog.FieldAwareLookup lookup, int n2, int m2, int p2) {
		long inner = lookup.findRank(n2, m2, p2);
		long[] out = new long[3];
		if (inner >= Recombination.SotaResolver.UNKNOWN_RANK) {
			return out;
		}
		eu.solven.matmul.catalog.SerendipitousBudProduct.BudDecomposition dec =
				eu.solven.matmul.catalog.SerendipitousBudProduct.findBuds(s.toAlgorithm());
		for (eu.solven.matmul.catalog.SerendipitousBudProduct.Bud b : dec.buds()) {
			int k = b.terms().length;
			long enlarged = switch (b.type()) {
				case U -> lookup.findRank(n2, m2, k * p2);
				case V -> lookup.findRank(k * n2, m2, p2);
				case W -> lookup.findRank(n2, k * m2, p2);
			};
			if (enlarged >= Recombination.SotaResolver.UNKNOWN_RANK) {
				continue;  // unpriceable bud: contributes no saving
			}
			out[b.type().ordinal()] += k * inner - enlarged;
		}
		return out;
	}

	// ── native metric implementations ───────────────────────────────────────

	/**
	 * Bud score: for each of U/V/W independently, partition the products by
	 * proportional vector direction (gcd+sign normalization — exact for integer
	 * vectors) and sum all class sizes ≥2. Matches
	 * {@code Σ classes≥2 of SerendipitousBudProduct.independentClassSizes}.
	 * Zero vectors are excluded from classes (transient walk states only — a
	 * reduced scheme has none).
	 */
	public static int budScore(FlipScheme s) {
		int total = 0;
		for (FlipScheme.Slot slot : FlipScheme.Slot.values()) {
			Map<FlipScheme.VecKey, Integer> sizes = new HashMap<>();
			for (int l = 0; l < s.rank(); l++) {
				int[] vec = s.vec(slot, l);
				if (FlipScheme.isZero(vec)) {
					continue;
				}
				sizes.merge(new FlipScheme.VecKey(directionNormalize(vec)), 1, Integer::sum);
			}
			for (int size : sizes.values()) {
				if (size >= 2) {
					total += size;
				}
			}
		}
		return total;
	}

	/** Vector divided by gcd of its entries, negated to a positive leading nonzero. */
	private static int[] directionNormalize(int[] vec) {
		int g = 0;
		for (int x : vec) {
			g = gcd(g, Math.abs(x));
		}
		int[] norm = FlipScheme.signNormalize(vec);
		if (g <= 1) {
			return norm;
		}
		int[] out = norm == vec ? vec.clone() : norm;
		for (int i = 0; i < out.length; i++) {
			out[i] /= g;
		}
		return out;
	}

	private static int gcd(int a, int b) {
		while (b != 0) {
			int t = a % b;
			a = b;
			b = t;
		}
		return a;
	}

	/**
	 * Projection margin μ: max over (axis, index) of the number of products
	 * whose support on that axis — in either of the two matrices carrying it
	 * (n→{U,W}, m→{U,V}, p→{V,W}) — collapses to that single index. Mirrors
	 * {@code ProjectionSearch.projectionMargin}/{@code axisMargin} including the
	 * count-once-per-index rule.
	 */
	public static int projectionMargin(FlipScheme s) {
		int best = 0;
		int r = s.rank();
		// n-axis: U row index = idx/m, W row index = idx/p
		best = Math.max(best, axisMargin(s, r, s.n,
				l -> singleAxisIndex(s.u(l), idx -> idx / s.m),
				l -> singleAxisIndex(s.w(l), idx -> idx / s.p)));
		// m-axis: U col index = idx%m, V row index = idx/p
		best = Math.max(best, axisMargin(s, r, s.m,
				l -> singleAxisIndex(s.u(l), idx -> idx % s.m),
				l -> singleAxisIndex(s.v(l), idx -> idx / s.p)));
		// p-axis: V col index = idx%p, W col index = idx%p
		best = Math.max(best, axisMargin(s, r, s.p,
				l -> singleAxisIndex(s.v(l), idx -> idx % s.p),
				l -> singleAxisIndex(s.w(l), idx -> idx % s.p)));
		return best;
	}

	private interface AxisIndexOf {
		/** The product's single axis-index in one carrying matrix, or −1 if its
		 *  support there spans 0 or ≥2 indices. */
		int get(int product);
	}

	private static int axisMargin(FlipScheme s, int r, int dim,
			AxisIndexOf first, AxisIndexOf second) {
		int[] death = new int[dim];
		for (int l = 0; l < r; l++) {
			int i1 = first.get(l);
			int i2 = second.get(l);
			if (i1 >= 0) {
				death[i1]++;
			}
			if (i2 >= 0 && i2 != i1) {
				death[i2]++;
			}
		}
		int max = 0;
		for (int d : death) {
			max = Math.max(max, d);
		}
		return max;
	}

	private interface IndexMap {
		int axisIndex(int flatIndex);
	}

	private static int singleAxisIndex(int[] vec, IndexMap map) {
		int found = -1;
		for (int idx = 0; idx < vec.length; idx++) {
			if (vec[idx] == 0) {
				continue;
			}
			int ai = map.axisIndex(idx);
			if (found < 0) {
				found = ai;
			} else if (found != ai) {
				return -1;
			}
		}
		return found;
	}
}
