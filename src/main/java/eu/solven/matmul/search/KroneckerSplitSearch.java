package eu.solven.matmul.search;

import eu.solven.matmul.catalog.Compose;

import eu.solven.matmul.recombination.Recombination;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Search over Kronecker-product decompositions of a target
 * {@code ⟨n, m, p⟩}: for every factorisation
 * {@code n = n₁·n₂, m = m₁·m₂, p = p₁·p₂} with both factors ≥ 2,
 * consult {@code sota} for {@code R(⟨n₁,m₁,p₁⟩)} and
 * {@code R(⟨n₂,m₂,p₂⟩)}; the resulting Kronecker product (Sedoglavic
 * Lemma 13 in his 2017 ⟨7,7,7⟩=250 paper, also Compose.kroneckerGeneral)
 * achieves rank {@code r₁·r₂}. Returns the (n₁,m₁,p₁; n₂,m₂,p₂) tuple
 * minimising the product rank.
 *
 * <p>By itself this is just an inventory: it doesn't construct the
 * scheme (use {@link #materialise} for that). But fed into the
 * recursive closure {@link RecursiveClosureSota}, every newly-reached
 * shape unlocks further Kronecker splits as a side-effect — that's
 * how DIS09's table-fill discovers compositions like
 * {@code ⟨18,18,18⟩ = ⟨2,2,2⟩⊗⟨9,9,9⟩} without hardcoding.</p>
 */
public final class KroneckerSplitSearch {

	private KroneckerSplitSearch() {}

	public record KroneckerSplit(
			int targetN, int targetM, int targetP,
			int n1, int m1, int p1, int r1,
			int n2, int m2, int p2, int r2,
			long totalRank) {
		public String breakdown() {
			return String.format("⟨%d,%d,%d⟩=%d ⊗ ⟨%d,%d,%d⟩=%d → %d",
					n1, m1, p1, r1, n2, m2, p2, r2, totalRank);
		}
	}

	/**
	 * Find the minimum-rank Kronecker decomposition. Returns
	 * {@link Optional#empty()} when no non-trivial factorisation exists
	 * (any of {@code n, m, p} prime).
	 */
	public static Optional<KroneckerSplit> findBest(int targetN, int targetM, int targetP,
			Recombination.SotaResolver sota) {
		long bestRank = Long.MAX_VALUE;
		KroneckerSplit best = null;

		// Include UNIT factor pairs (1, n) / (n, 1): they enable the "stack k
		// copies along one axis" decomposition, e.g. ⟨3,3,18⟩ = ⟨1,1,3⟩⊗⟨3,3,6⟩
		// = 3·40 = 120 — which the ≥2-only enumeration missed entirely (and which
		// is the only Kronecker route when a target dim is prime, like n=m=3).
		List<int[]> nFactors = factorPairsWithUnit(targetN);
		List<int[]> mFactors = factorPairsWithUnit(targetM);
		List<int[]> pFactors = factorPairsWithUnit(targetP);

		for (int[] nf : nFactors) {
			for (int[] mf : mFactors) {
				for (int[] pf : pFactors) {
					int n1 = nf[0], n2 = nf[1];
					int m1 = mf[0], m2 = mf[1];
					int p1 = pf[0], p2 = pf[1];
					// Skip the no-op: one factor being ⟨1,1,1⟩ ⇒ the other is the
					// whole target (would just echo R(target) and risk recursion).
					if ((n1 == 1 && m1 == 1 && p1 == 1) || (n2 == 1 && m2 == 1 && p2 == 1)) continue;
					int r1 = rankOf(sota, n1, m1, p1);
					int r2 = rankOf(sota, n2, m2, p2);
					if (r1 <= 0 || r2 <= 0) continue;
					if (r1 >= n1 * m1 * p1 && r2 >= n2 * m2 * p2) continue;  // both naive — skip
					long tot = (long) r1 * r2;
					if (tot < bestRank) {
						bestRank = tot;
						best = new KroneckerSplit(targetN, targetM, targetP,
								n1, m1, p1, r1, n2, m2, p2, r2, tot);
					}
				}
			}
		}
		return Optional.ofNullable(best);
	}

	/**
	 * Materialise the Kronecker scheme implied by a split. Looks up
	 * the actual sub-algorithms via {@code lookup}.
	 */
	public static NonCubicBilinearAlgorithm materialise(KroneckerSplit split,
			Recombination.AlgorithmLookup lookup) {
		NonCubicBilinearAlgorithm outer = factorScheme(lookup, split.n1, split.m1, split.p1);
		NonCubicBilinearAlgorithm inner = factorScheme(lookup, split.n2, split.m2, split.p2);
		return Compose.kroneckerGeneral(outer, inner);
	}

	/** Resolve a Kronecker factor to an actual scheme. Degenerate shapes (any
	 *  axis = 1) have no catalog file but their optimal scheme IS the naive one,
	 *  so synthesise it via {@link eu.solven.matmul.NaiveMatMul#ofNonCubic}. */
	private static NonCubicBilinearAlgorithm factorScheme(
			Recombination.AlgorithmLookup lookup, int n, int m, int p) {
		if (n == 1 || m == 1 || p == 1) {
			return eu.solven.matmul.NaiveMatMul.ofNonCubic(n, m, p);
		}
		return lookup.find(n, m, p).orElseThrow(() -> new IllegalStateException(
				"missing Kronecker factor ⟨" + n + "," + m + "," + p + "⟩"));
	}

	/** Rank of a sub-factor: for a degenerate shape (any axis = 1) the bilinear
	 *  rank is exactly the naive product n·m·p (no Strassen saving possible), so
	 *  compute it directly rather than relying on the catalog resolver to carry
	 *  ⟨1,…⟩ entries. Otherwise defer to {@code sota}. */
	private static int rankOf(Recombination.SotaResolver sota, int a, int b, int c) {
		if (a == 1 || b == 1 || c == 1) return a * b * c;
		return sota.getRank(a, b, c);
	}

	/** All ordered factor pairs (a, b) with {@code a·b = n}, including the UNIT
	 *  pairs (1, n) and (n, 1). Always non-empty (≥ the unit pairs), so a prime
	 *  target dim still admits the trivial ⟨1,…⟩ Kronecker factor. */
	private static List<int[]> factorPairsWithUnit(int n) {
		List<int[]> out = new ArrayList<>();
		out.add(new int[] { 1, n });
		if (n != 1) out.add(new int[] { n, 1 });
		for (int a = 2; a * a <= n; a++) {
			if (n % a == 0) {
				int b = n / a;
				if (b >= 2) {
					out.add(new int[] { a, b });
					if (a != b) out.add(new int[] { b, a });
				}
			}
		}
		return out;
	}
}
