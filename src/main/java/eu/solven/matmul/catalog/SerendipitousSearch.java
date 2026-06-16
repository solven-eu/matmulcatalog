package eu.solven.matmul.catalog;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;

/**
 * Discovery driver for serendipitous products (#159). For a target
 * {@code ⟨n,m,p⟩}, enumerates factorizations {@code ⟨n,m,p⟩ = base ⊗ ⟨n₂,m₂,p₂⟩}
 * over a set of candidate base schemes, decomposes each base by its buds
 * ({@link SerendipitousBudProduct}), predicts the serendipitous rank
 * {@code r_s = Σ Sᵢ·R(⟨Nᵢn₂,Mᵢm₂,Pᵢp₂⟩)}, and — for the cheapest predictions —
 * builds and {@link Verifier#isExactNonCubic verifies} the scheme, returning the
 * best one strictly below {@code upperBound}.
 *
 * <p>The win comes entirely from the base's bud structure: a base with no buds
 * reduces to the naive Kronecker product, so this only improves on shapes where
 * some candidate base carries buds.</p>
 */
public final class SerendipitousSearch {

	private SerendipitousSearch() {}

	/** A verified serendipitous construction for the target. */
	public record Hit(NonCubicBilinearAlgorithm scheme, NonCubicBilinearAlgorithm base,
			int n2, int m2, int p2, long rank) {}

	/**
	 * Best verified serendipitous product for {@code ⟨n,m,p⟩} strictly below
	 * {@code upperBound}, over the given candidate bases.
	 */
	public static Optional<Hit> bestFor(int n, int m, int p,
			List<NonCubicBilinearAlgorithm> bases, FieldAwareLookup lookup, long upperBound) {
		// Phase 1 — PREDICT only (cheap): collect every viable (base, inner)
		// factorization with its predicted serendipitous rank, WITHOUT building or
		// verifying. Building the dense scheme + verifying is the expensive part, so
		// we defer it to the single best candidate (Phase 2). The serendipitous
		// product is a proven identity (Smith 2002 / KMW 2026), so prediction == the
		// built rank whenever the inputs are well-formed; we still confirm per-build.
		List<Candidate> cands = new ArrayList<>();
		for (NonCubicBilinearAlgorithm base : bases) {
			int n1 = base.n, m1 = base.m, p1 = base.p;
			if (n % n1 != 0 || m % m1 != 0 || p % p1 != 0) continue;
			int n2 = n / n1, m2 = m / m1, p2 = p / p1;
			if (n2 * m2 * p2 == 1) continue;       // would be the base itself
			if (n1 * m1 * p1 == 1) continue;       // trivial base
			long naiveKron = (long) base.r * rank(lookup, n2, m2, p2);
			if (naiveKron < 0) continue;           // second scheme unavailable

			// The greedy decomposition is order-sensitive (a term shared between a
			// U-class and a V-class goes to whichever type is processed first), so a
			// larger, cheaper bud of a later type can be masked. Try every type
			// ordering and keep the cheapest — this is what surfaces e.g. the size-3
			// V-bud of ⟨4,3,3⟩ that yields ⟨8,9,9⟩=430 (U-first would only see 434).
			SerendipitousBudProduct.BudDecomposition bestDec = null;
			long predicted = Long.MAX_VALUE;
			for (SerendipitousBudProduct.BudType[] order : SerendipitousBudProduct.ALL_ORDERINGS) {
				SerendipitousBudProduct.BudDecomposition dec = SerendipitousBudProduct.findBuds(base, order);
				if (dec.buds().isEmpty()) continue;  // no buds → no saving over Kronecker
				long pr = predictRank(dec, lookup, n2, m2, p2);
				if (pr >= 0 && pr < predicted) { predicted = pr; bestDec = dec; }
			}
			if (bestDec == null || predicted >= upperBound || predicted >= naiveKron) continue;

			cands.add(new Candidate(base, n2, m2, p2, predicted, bestDec));
		}
		if (cands.isEmpty()) return Optional.empty();

		// Phase 2 — BUILD + VERIFY best-first: only the cheapest-predicted candidate
		// is built and checked; if it passes it is returned immediately (it is the
		// global best, since the list is sorted by predicted rank ascending). A
		// candidate that fails to build at the predicted rank or fails the spot-check
		// falls through to the next-best. Verification is a cheap random matmul probe
		// (catches a wrong index map / bad inner with overwhelming probability);
		// EXACT verification is reserved for promote-time, not the hot search loop.
		cands.sort(Comparator.comparingLong(Candidate::predicted));
		for (Candidate c : cands) {
			NonCubicBilinearAlgorithm built = SerendipitousBudProduct.productFromDecomposition(
					c.base(), c.dec(), lookup, c.n2(), c.m2(), c.p2(),
					java.util.EnumSet.allOf(SerendipitousBudProduct.BudType.class));
			if (built.r != c.predicted()) continue;                 // sanity
			if (!Verifier.passesRandomMatmulSpotCheck(built)) continue;
			return Optional.of(new Hit(built, c.base(), c.n2(), c.m2(), c.p2(), built.r));
		}
		return Optional.empty();
	}

	/** A viable factorization with its predicted rank + chosen decomposition. */
	private record Candidate(NonCubicBilinearAlgorithm base, int n2, int m2, int p2, long predicted,
			SerendipitousBudProduct.BudDecomposition dec) {}

	private static long predictRank(SerendipitousBudProduct.BudDecomposition dec,
			FieldAwareLookup lookup, int n2, int m2, int p2) {
		long total = 0;
		long triv = rank(lookup, n2, m2, p2);
		if (triv < 0) return -1;
		total += (long) dec.trivial().length * triv;
		for (SerendipitousBudProduct.Bud bud : dec.buds()) {
			int k = bud.terms().length;
			long r = switch (bud.type()) {
				case U -> rank(lookup, n2, m2, k * p2);
				case V -> rank(lookup, k * n2, m2, p2);
				case W -> rank(lookup, n2, k * m2, p2);
			};
			if (r < 0) return -1;
			total += r;
		}
		return total;
	}

	private static long rank(FieldAwareLookup lookup, int n, int m, int p) {
		return lookup.findWithSource(n, m, p).map(ws -> (long) ws.alg().r).orElse(-1L);
	}
}
