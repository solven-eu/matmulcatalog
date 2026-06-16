package eu.solven.matmul.papers.sedoglavic2017;

import eu.solven.matmul.catalog.Recombination;

import java.util.Optional;

import eu.solven.matmul.catalog.Recombination.SotaResolver;

/**
 * <strong>Sedoglavic 2017 Proposition 1</strong> — the closed-form
 * recipe behind FMM-Lille's cubic catalog entries.
 *
 * <p>The proposition (hal-01572046v2, Dec 2017) states:</p>
 *
 * <pre>
 *   ⟨u+v, u+v, u+v⟩ ≤ ⟨u,u,u⟩ + 3·⟨u,u,v⟩ + 3·⟨v,v,u⟩    when u &gt; v
 * </pre>
 *
 * <p>Plug in any (u, v) with u &gt; v ≥ 1 and read off a cubic upper
 * bound from the three sub-shape ranks. The famous (u,v)=(4,3) gives
 * the paper's title bound {@code ⟨7,7,7⟩=250}; (u,v)=(6,5) gives the
 * FMM-Lille bound {@code ⟨11,11,11⟩=873}; (u,v)=(9,8) gives
 * {@code ⟨17,17,17⟩=2940}; (u,v)=(10,9) gives {@code ⟨19,19,19⟩=4044}.</p>
 *
 * <h2>Doubling extension: u = v = k</h2>
 *
 * <p>The published proposition requires u &gt; v. The natural extension
 * to u = v = k uses Pan 1980's TA pair-fusion to absorb the six
 * same-shape cubic copies into three paired computations:</p>
 *
 * <pre>
 *   ⟨2k, 2k, 2k⟩ ≤ ⟨k,k,k⟩ + 3·pair_cost(k,k,k)
 *                where pair_cost(k,k,k) = k³ + 3k²
 * </pre>
 *
 * <p>k=7 gives ⟨14,14,14⟩=1719, k=11 gives ⟨22,22,22⟩=5955 etc. —
 * both confirmed against catalog imports.</p>
 *
 * <h2>How to use</h2>
 *
 * <p>{@link #predict(int, SotaResolver)} for the cubic target ⟨n,n,n⟩
 * enumerates all (u,v) with u+v=n and u≥v, returning the minimum
 * predicted rank with its attribution. {@link Prediction} carries the
 * winning (u,v) tuple and the constructive recipe in
 * {@code lineage_compact} form.</p>
 *
 * <p>Belongs to {@code eu.solven.matmul.papers.sedoglavic2017} per
 * project convention (one package per paper) and is wired into
 * {@code BlockSplitSearch.findBestStrategy} alongside Kronecker,
 * Concat, recombination, pair-fused, and {@code KnownTauIdentities}
 * candidates. Closes the cubic gaps the hand-encoded identities
 * inside {@code KnownTauIdentities} can't cover individually.</p>
 *
 * @see <a href="https://hal.science/hal-01572046v2">Sedoglavic 2017 (hal-01572046v2)</a>
 */
public final class SedoglavicProp1 {

	private SedoglavicProp1() {}

	/**
	 * Result of the (u,v) enumeration at a cubic target. {@link #lineageCompact}
	 * captures the recipe in the same canonical-key notation the rest of
	 * the catalog uses.
	 */
	public record Prediction(
			int n,
			int u, int v,
			long predictedRank,
			boolean usesPanPairCost,
			String lineageCompact) {}

	/**
	 * Pan 1980 pair-fusion cost: {@code a·b·c + a·b + b·c + c·a}. The
	 * formula computes two cyclically-related products jointly with
	 * fewer multiplications than naïve sum. We use it for the
	 * u = v = k doubling extension.
	 */
	public static long pairCost(int a, int b, int c) {
		return (long) a * b * c + (long) a * b + (long) b * c + (long) c * a;
	}

	/**
	 * Predict the rank of {@code ⟨n,n,n⟩} via Sedoglavic Prop 1 by
	 * enumerating all (u,v) splits with u+v=n and u≥v.
	 *
	 * @param n target cubic side length, n ≥ 2
	 * @param sota catalog rank resolver
	 * @return the prediction with the smallest rank, or empty if no
	 *         valid (u,v) split has all sub-ranks resolvable
	 */
	public static Optional<Prediction> predict(int n, SotaResolver sota) {
		if (n < 2) return Optional.empty();
		Prediction best = null;
		// u ≥ v ≥ 1; u + v = n
		for (int v = 1; v <= n / 2; v++) {
			int u = n - v;
			long pred;
			boolean doubling;
			String formula;
			if (u == v) {
				// Doubling extension via Pan TA pair-fusion.
				int rUUU = sota.getRank(u, u, u);
				if (rUUU >= Recombination.SotaResolver.UNKNOWN_RANK) continue;
				long pc = pairCost(u, u, u);
				pred = rUUU + 3L * pc;
				doubling = true;
				formula = String.format(
						"SedoglavicDoubling(k=%d) = %dx%dx%d_m%d + 3·PanPair(%d,%d,%d)",
						u, u, u, u, rUUU, u, u, u);
			} else {
				int rUUU = sota.getRank(u, u, u);
				int rUUV = sota.getRank(u, u, v);
				int rVVU = sota.getRank(v, v, u);
				if (rUUU >= Recombination.SotaResolver.UNKNOWN_RANK
						|| rUUV >= Recombination.SotaResolver.UNKNOWN_RANK
						|| rVVU >= Recombination.SotaResolver.UNKNOWN_RANK) continue;
				pred = (long) rUUU + 3L * rUUV + 3L * rVVU;
				doubling = false;
				formula = String.format(
						"SedoglavicProp1(u=%d,v=%d) = %dx%dx%d_m%d + 3·%dx%dx%d_m%d + 3·%dx%dx%d_m%d",
						u, v,
						u, u, u, rUUU,
						u, u, v, rUUV,
						v, v, u, rVVU);
			}
			if (best == null || pred < best.predictedRank) {
				best = new Prediction(n, u, v, pred, doubling, formula);
			}
		}
		return Optional.ofNullable(best);
	}

	/**
	 * Closed-form rank when both ⟨u,u,v⟩ and ⟨v,v,u⟩ resolve to their
	 * trivial-cubic baselines. Useful for spot-checks; not used in the
	 * live search (which always defers to the SOTA resolver).
	 */
	public static long naiveBound(int u, int v) {
		// Trivial: ⟨u,u,u⟩ ≤ u³, ⟨u,u,v⟩ ≤ u²v, ⟨v,v,u⟩ ≤ v²u.
		return (long) u * u * u + 3L * u * u * v + 3L * v * v * u;
	}
}
