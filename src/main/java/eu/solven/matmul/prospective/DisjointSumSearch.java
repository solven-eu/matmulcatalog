package eu.solven.matmul.prospective;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;

import eu.solven.matmul.isotropy.PairedSubProducts;
import eu.solven.matmul.recombination.Recombination;

/**
 * τ-theorem-style disjoint-MM-sum decomposition search (Pan 1980 /
 * Schönhage 1981; DIS09 §3; Schwartz-Zwecher 2025). For a target
 * matmul {@code ⟨n,m,p⟩}, enumerate ways to express it as a sum of
 * smaller matmul sub-tensors that together cover the target tensor,
 * optionally augmented with same-cubic-shape "trilinear-aggregation"
 * legs (Pan's pair-cost formula {@code a³ + 3a²} for two same-shape
 * cubic products).
 *
 * <h2>Why this exists</h2>
 * <p>Our existing {@link BlockSplitSearch} explores Strassen-like outer
 * recombinations — one specific algebraic identity (Strassen ⟨2,2,2⟩=7,
 * Laderman ⟨3,3,3⟩=23, etc.) parametrised by block allocations. That's
 * one corner of the larger trilinear-aggregation framework. The general
 * mechanism is: pick any collection of sub-MMs whose A/B/C supports cover
 * the target tensor, optionally share computation across same-shape
 * pairs. FMM-Lille's {@code ⟨17,17,17⟩=2934} is exactly this — 5 disjoint
 * sub-MMs + 1 TA pair, not a Strassen recombination at all.</p>
 *
 * <h2>Search strategy (beam search by ω_eff)</h2>
 * <p>The asymptotic exponent of matmul comes from cascading low-ω
 * sub-products (Strassen 1986, the "laser method"). We prioritise leaves
 * by their effective exponent
 * {@code ω_eff(⟨n,m,p⟩) = 3·log(R) / log(n·m·p)}, smaller-first. At
 * each search step we branch into the top-K candidates by ω_eff that
 * fit the remaining cover; the search explores a beam of partial
 * decompositions and returns the lowest-total-rank completion.</p>
 *
 * <h2>Cover constraints (and why they're not enough)</h2>
 * <p>A decomposition {@code Σᵢ ⟨nᵢ,mᵢ,pᵢ⟩} of {@code ⟨n,m,p⟩} must satisfy
 * the area-based cover bounds {@code Σ nᵢ·mᵢ ≥ n·m},
 * {@code Σ mᵢ·pᵢ ≥ m·p}, {@code Σ nᵢ·pᵢ ≥ n·p}. These are NECESSARY
 * (entry-counting arguments) but FAR FROM SUFFICIENT — they admit shape
 * collections that cannot actually be assembled into an exact-rank
 * algorithm. <strong>This implementation enforces only these loose
 * conditions, so the predicted rank is an OVER-OPTIMISTIC HEURISTIC,
 * not a constructive bound.</strong></p>
 *
 * <p>Concrete demonstrations of the false-positive flaw:</p>
 * <ul>
 *   <li>{@code ⟨4,4,4⟩} the search returns {@code 4·⟨2,2,2⟩=28} —
 *       but the actual exact rank is 47/F₂, 48/Q. 28 is unrealisable.</li>
 *   <li>{@code ⟨17,17,17⟩} returns ~2395, well below FMM-Lille's 2934.
 *       The shape collection passes the area cover but no Pan-style
 *       construction realises it at that rank.</li>
 * </ul>
 *
 * <p>The TRUE constraint is algebraic — there must exist linear-form
 * embeddings of each sub-tensor into the target trilinear form such
 * that the embeddings sum (with possible cancellation) to the target.
 * Verifying this requires constructive feasibility — out of scope here.</p>
 *
 * <p><strong>This search is currently useful only as a shape-exploration
 * tool</strong>: candidate shape multisets ranked by their area-cover-
 * passing total rank, to be HAND-VERIFIED against known constructive
 * algorithms. It is NOT a reliable rank predictor in its current form.
 * See follow-up tasks for tighter feasibility checks.</p>
 *
 * <h2>Same-shape TA post-pass</h2>
 * <p>After greedy/beam fill, scan for clusters of identical-cubic-shape
 * sub-products. For each cluster of {@code k ≥ 2} same-cubic ⟨a,a,a⟩
 * sub-products, replace {@code k·R(⟨a,a,a⟩)} with the Pan pair-cost
 * {@code (k/2)·(a³ + 3·a²) + (k mod 2)·R(⟨a,a,a⟩)} when profitable.</p>
 */
public final class DisjointSumSearch {

	private DisjointSumSearch() {}

	/** A single sub-MM term in the decomposition. */
	public record SubMM(int n, int m, int p, long rank) {
		public long aSupport() { return (long) n * m; }
		public long bSupport() { return (long) m * p; }
		public long cSupport() { return (long) n * p; }
		public boolean isCubic() { return n == m && m == p; }
		@Override
		public String toString() { return "⟨" + n + "," + m + "," + p + "⟩=" + rank; }
	}

	/** A TA leg replacing two cubic same-shape sub-products with a Pan pair-fuse. */
	public record TALeg(int firstIdx, int secondIdx, long fusedCost) {}

	public record DisjointSumPrediction(
			int targetN, int targetM, int targetP,
			List<SubMM> children,
			List<TALeg> taLegs,
			long totalRank) {

		public String label() {
			StringBuilder sb = new StringBuilder("disjoint-sum[");
			for (int i = 0; i < children.size(); i++) {
				if (i > 0) sb.append(" + ");
				sb.append(children.get(i).toString());
			}
			if (!taLegs.isEmpty()) {
				sb.append("; TA-legs=").append(taLegs.size());
			}
			sb.append("] = ").append(totalRank);
			return sb.toString();
		}
	}

	/**
	 * Public entry: beam-search the disjoint-sum decomposition space.
	 *
	 * @param n,m,p target shape
	 * @param sota  rank oracle for sub-product shapes
	 * @param beamK beam width (≥1; 1 = pure greedy)
	 * @param maxTerms cap on number of sub-MMs in the decomposition
	 * @param minAxis  minimum per-axis size for sub-MMs (≥ 1)
	 */
	public static Optional<DisjointSumPrediction> findBest(int n, int m, int p,
			Recombination.SotaResolver sota, int beamK, int maxTerms, int minAxis) {
		if (n < 1 || m < 1 || p < 1) return Optional.empty();
		if (beamK < 1) beamK = 1;
		if (maxTerms < 1) maxTerms = 8;
		if (minAxis < 1) minAxis = 1;

		// Build the candidate pool: every shape ⟨a,b,c⟩ with axes
		// in [minAxis..target] for which sota has a rank below the
		// trivial a·b·c (else we'd use the trivial bound directly).
		List<SubMM> pool = buildPool(n, m, p, sota, minAxis);
		if (pool.isEmpty()) return Optional.empty();
		// Sort by ω_eff ascending (lowest first), with tie-break on
		// support-area descending (use bigger pieces first when ω_eff ties).
		pool.sort(Comparator
				.<SubMM>comparingDouble(s -> omegaEff(s.n, s.m, s.p, s.rank))
				.thenComparingLong(s -> -(s.n * (long) s.m + s.m * (long) s.p + s.n * (long) s.p)));

		final long targetCoverA = (long) n * m;
		final long targetCoverB = (long) m * p;
		final long targetCoverC = (long) n * p;

		// Beam state: (partial decomposition, cover-remaining-A, B, C, rank-so-far).
		// Start with empty decomposition.
		List<BeamState> beam = new ArrayList<>();
		beam.add(new BeamState(List.of(), targetCoverA, targetCoverB, targetCoverC, 0L));

		// Track best-finished overall so we can prune.
		long bestFinished = Long.MAX_VALUE;
		List<SubMM> bestChildren = null;

		for (int depth = 0; depth < maxTerms && !beam.isEmpty(); depth++) {
			List<BeamState> nextBeam = new ArrayList<>();
			for (BeamState state : beam) {
				// For this state, evaluate top-K candidates by ω_eff that fit.
				// "fit" = (a, b, c) ≤ (n, m, p) AND adds support to at least
				// one axis that still has cover-remaining > 0.
				int picked = 0;
				for (SubMM cand : pool) {
					if (picked >= beamK) break;
					if (cand.n > n || cand.m > m || cand.p > p) continue;
					boolean helps = (state.remainA > 0 && cand.aSupport() > 0)
							|| (state.remainB > 0 && cand.bSupport() > 0)
							|| (state.remainC > 0 && cand.cSupport() > 0);
					if (!helps) continue;
					long newRank = state.rankSoFar + cand.rank;
					if (newRank >= bestFinished) continue;  // prune
					long newA = Math.max(0, state.remainA - cand.aSupport());
					long newB = Math.max(0, state.remainB - cand.bSupport());
					long newC = Math.max(0, state.remainC - cand.cSupport());
					List<SubMM> newChildren = appendOne(state.children, cand);
					BeamState ns = new BeamState(newChildren, newA, newB, newC, newRank);
					if (newA == 0 && newB == 0 && newC == 0) {
						// Cover complete — record candidate finish, no need to
						// extend further.
						long finalRank = applySameShapeTaFusion(newChildren).totalRank;
						if (finalRank < bestFinished) {
							bestFinished = finalRank;
							bestChildren = newChildren;
						}
					} else {
						nextBeam.add(ns);
					}
					picked++;
				}
			}
			// Trim nextBeam to the top-beamK·beamK states by rank-so-far (rough
			// upper bound on memory).
			if (nextBeam.size() > 16L * beamK) {
				nextBeam.sort(Comparator.comparingLong(s -> s.rankSoFar));
				beam = new ArrayList<>(nextBeam.subList(0, (int) (16L * beamK)));
			} else {
				beam = nextBeam;
			}
		}

		if (bestChildren == null) return Optional.empty();
		FusedResult fused = applySameShapeTaFusion(bestChildren);
		return Optional.of(new DisjointSumPrediction(
				n, m, p, bestChildren, fused.taLegs, fused.totalRank));
	}

	/** Convenience: defaults beamK=5, maxTerms=8, minAxis=2. */
	public static Optional<DisjointSumPrediction> findBest(int n, int m, int p,
			Recombination.SotaResolver sota) {
		return findBest(n, m, p, sota, 5, 8, 2);
	}

	private record BeamState(List<SubMM> children, long remainA, long remainB, long remainC,
			long rankSoFar) {}

	private static List<SubMM> appendOne(List<SubMM> prefix, SubMM tail) {
		List<SubMM> out = new ArrayList<>(prefix.size() + 1);
		out.addAll(prefix);
		out.add(tail);
		return out;
	}

	/**
	 * Build the candidate pool. Enumerate every {@code ⟨a,b,c⟩} with axes
	 * in {@code [minAxis, max(n,m,p)]} and ask sota for its rank. Only
	 * include shapes where {@code sota.rank} is strictly below the trivial
	 * {@code a·b·c} (else the trivial bound would dominate the decomposition).
	 *
	 * <p>The pool is bounded to keep enumeration tractable — at most a few
	 * thousand entries for typical targets.</p>
	 */
	private static List<SubMM> buildPool(int n, int m, int p,
			Recombination.SotaResolver sota, int minAxis) {
		List<SubMM> pool = new ArrayList<>();
		int maxAxis = Math.max(Math.max(n, m), p);
		// Cubic shapes first: ⟨a,a,a⟩ for a in [minAxis, maxAxis].
		// Then small rectangular shapes — full grid is too big at large
		// target, so cap.
		int rectCap = Math.min(maxAxis, 9);  // exhaustive non-cubic only up to 9
		for (int a = minAxis; a <= maxAxis; a++) {
			addIfBetter(pool, sota, a, a, a);
		}
		for (int a = minAxis; a <= rectCap; a++) {
			for (int b = minAxis; b <= rectCap; b++) {
				for (int c = minAxis; c <= rectCap; c++) {
					if (a == b && b == c) continue;  // cubic added above
					addIfBetter(pool, sota, a, b, c);
				}
			}
		}
		return pool;
	}

	private static void addIfBetter(List<SubMM> pool,
			Recombination.SotaResolver sota, int a, int b, int c) {
		long rank = sota.getRank(a, b, c);
		long trivial = (long) a * b * c;
		// Skip the trivial bound — it's the worst possible, never beneficial
		// to include in the cover.
		if (rank <= 0 || rank >= trivial) return;
		pool.add(new SubMM(a, b, c, rank));
	}

	/**
	 * ω_eff(⟨n,m,p⟩) = 3·log(R) / log(n·m·p). Smaller is better.
	 *
	 * <p>For cubic ⟨n,n,n⟩ this reduces to log_n(R), matching the
	 * classical matmul exponent of the recursion.</p>
	 */
	private static double omegaEff(int n, int m, int p, long rank) {
		double denom = Math.log((double) n * m * p);
		if (denom <= 0) return 3.0;
		return 3.0 * Math.log(rank) / denom;
	}

	/** Same-shape TA fusion (cubic-only first cut). */
	private record FusedResult(List<TALeg> taLegs, long totalRank) {}

	private static FusedResult applySameShapeTaFusion(List<SubMM> children) {
		long total = 0;
		List<TALeg> legs = new ArrayList<>();
		// Group cubic children by ⟨a,a,a⟩; pair off; non-cubic kept individual.
		// Pair cost: Pan's a³+3a², profitable iff < 2·rank(⟨a,a,a⟩).
		boolean[] paired = new boolean[children.size()];
		// Look for same-cubic-shape pairs.
		for (int i = 0; i < children.size(); i++) {
			if (paired[i] || !children.get(i).isCubic()) continue;
			SubMM ci = children.get(i);
			for (int j = i + 1; j < children.size(); j++) {
				if (paired[j] || !children.get(j).isCubic()) continue;
				SubMM cj = children.get(j);
				if (ci.n != cj.n) continue;  // same cubic shape required
				long pairCost = PairedSubProducts.pairCost(ci.n, ci.m, ci.p);
				long indep = ci.rank + cj.rank;
				if (pairCost < indep) {
					paired[i] = paired[j] = true;
					legs.add(new TALeg(i, j, pairCost));
					total += pairCost;
					break;
				}
			}
		}
		// Add unpaired children.
		for (int i = 0; i < children.size(); i++) {
			if (!paired[i]) total += children.get(i).rank;
		}
		return new FusedResult(legs, total);
	}
}
