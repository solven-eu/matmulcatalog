package eu.solven.matmul.search;

import eu.solven.matmul.catalog.Compose;

import eu.solven.matmul.catalog.Recombination;

import java.util.List;
import java.util.Optional;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * One-axis additive decomposition of a target {@code ⟨n,m,p⟩} matmul.
 * All three tensor modes can be split — each is a direct sum of the
 * matmul tensor along one mode, hence rank {@code r1 + r2}:
 * <ul>
 *   <li>{@code p = p1 + p2} — output-column tile ({@link Compose#concatRight}),</li>
 *   <li>{@code n = n1 + n2} — output-row tile ({@link Compose#concatBelow}),</li>
 *   <li>{@code m = m1 + m2} — inner/contraction <em>sum</em>
 *       ({@link Compose#concatInner}): {@code C = A1·B1 + A2·B2}.</li>
 * </ul>
 *
 * <p><b>Note (corrected 2026-06-04):</b> an earlier version of this
 * docstring claimed the {@code m} axis "cannot be split this way... needs
 * a real matmul algorithm". That conflated <em>tiling</em> (placing
 * disjoint output blocks side-by-side) with additive composition. The
 * m-split does not tile — it splits the shared inner dimension and
 * <em>accumulates</em> the two sub-products — but it is a perfectly valid
 * rank-{@code r1+r2} construction (a direct sum along the middle mode).
 * See {@link Compose#concatInner}.</p>
 *
 * <p>The rank is always {@code r1 + r2} — strictly additive — which sets
 * a baseline that Strassen-style recombination must beat to be preferred.
 * For "narrow" shapes (any axis = 2) the additive concat often *is* the
 * best constructive recipe.</p>
 */
public final class ConcatSplitSearch {

	private ConcatSplitSearch() {}

	public record ConcatSplit(
			int n, int m, int p,
			int axis,           // 0 = n-axis (concatBelow), 1 = m-axis (concatInner), 2 = p-axis (concatRight)
			int leftSize,
			int rightSize,
			long totalRank) {}

	/**
	 * Search all three axes; return the smaller-rank one. If none has a
	 * usable rank in {@code sota}, returns empty.
	 */
	public static Optional<ConcatSplit> findBest(int n, int m, int p,
			Recombination.SotaResolver sota) {
		Optional<ConcatSplit> best = Optional.empty();
		for (Optional<ConcatSplit> cand : List.of(
				findBestRight(n, m, p, sota),
				findBestBelow(n, m, p, sota),
				findBestInner(n, m, p, sota))) {
			if (cand.isEmpty()) continue;
			if (best.isEmpty() || cand.get().totalRank < best.get().totalRank) {
				best = cand;
			}
		}
		return best;
	}

	/** Split p axis only (concatRight). */
	public static Optional<ConcatSplit> findBestRight(int n, int m, int p,
			Recombination.SotaResolver sota) {
		long bestRank = Long.MAX_VALUE;
		int bestL = -1;
		for (int p1 = 1; p1 < p; p1++) {
			int p2 = p - p1;
			long r1 = sota.getRank(n, m, p1);
			long r2 = sota.getRank(n, m, p2);
			if (r1 <= 0 || r2 <= 0) continue;
			long tot = r1 + r2;
			if (tot < bestRank) { bestRank = tot; bestL = p1; }
		}
		if (bestL < 0) return Optional.empty();
		return Optional.of(new ConcatSplit(n, m, p, 2, bestL, p - bestL, bestRank));
	}

	/** Split n axis only (concatBelow). */
	public static Optional<ConcatSplit> findBestBelow(int n, int m, int p,
			Recombination.SotaResolver sota) {
		long bestRank = Long.MAX_VALUE;
		int bestL = -1;
		for (int n1 = 1; n1 < n; n1++) {
			int n2 = n - n1;
			long r1 = sota.getRank(n1, m, p);
			long r2 = sota.getRank(n2, m, p);
			if (r1 <= 0 || r2 <= 0) continue;
			long tot = r1 + r2;
			if (tot < bestRank) { bestRank = tot; bestL = n1; }
		}
		if (bestL < 0) return Optional.empty();
		return Optional.of(new ConcatSplit(n, m, p, 0, bestL, n - bestL, bestRank));
	}

	/**
	 * Split the inner/contraction m axis only (concatInner): the two halves
	 * share {@code n} and {@code p}, split {@code m = m1 + m2}, and their
	 * sub-products are <em>summed</em> ({@code C = A1·B1 + A2·B2}). Rank
	 * {@code r1 + r2}.
	 */
	public static Optional<ConcatSplit> findBestInner(int n, int m, int p,
			Recombination.SotaResolver sota) {
		long bestRank = Long.MAX_VALUE;
		int bestL = -1;
		for (int m1 = 1; m1 < m; m1++) {
			int m2 = m - m1;
			long r1 = sota.getRank(n, m1, p);
			long r2 = sota.getRank(n, m2, p);
			if (r1 <= 0 || r2 <= 0) continue;
			long tot = r1 + r2;
			if (tot < bestRank) { bestRank = tot; bestL = m1; }
		}
		if (bestL < 0) return Optional.empty();
		return Optional.of(new ConcatSplit(n, m, p, 1, bestL, m - bestL, bestRank));
	}

	/**
	 * Materialise the scheme implied by a concat split. Looks up the
	 * actual sub-algorithms via {@code lookup} and applies the relevant
	 * {@link Compose} operator.
	 */
	public static NonCubicBilinearAlgorithm materialise(ConcatSplit split,
			Recombination.AlgorithmLookup lookup) {
		if (split.axis == 2) {
			// p-axis tile (concatRight): share n,m; split p.
			NonCubicBilinearAlgorithm left = lookup.find(split.n, split.m, split.leftSize)
					.orElseThrow(() -> new IllegalStateException(
							"missing leaf ⟨" + split.n + "," + split.m + "," + split.leftSize + "⟩"));
			NonCubicBilinearAlgorithm right = lookup.find(split.n, split.m, split.rightSize)
					.orElseThrow(() -> new IllegalStateException(
							"missing leaf ⟨" + split.n + "," + split.m + "," + split.rightSize + "⟩"));
			return Compose.concatRight(left, right);
		} else if (split.axis == 1) {
			// m-axis contraction sum (concatInner): share n,p; split m; C = A1·B1 + A2·B2.
			NonCubicBilinearAlgorithm left = lookup.find(split.n, split.leftSize, split.p)
					.orElseThrow(() -> new IllegalStateException(
							"missing leaf ⟨" + split.n + "," + split.leftSize + "," + split.p + "⟩"));
			NonCubicBilinearAlgorithm right = lookup.find(split.n, split.rightSize, split.p)
					.orElseThrow(() -> new IllegalStateException(
							"missing leaf ⟨" + split.n + "," + split.rightSize + "," + split.p + "⟩"));
			return Compose.concatInner(left, right);
		} else {
			// n-axis tile (concatBelow): share m,p; split n.
			NonCubicBilinearAlgorithm top = lookup.find(split.leftSize, split.m, split.p)
					.orElseThrow(() -> new IllegalStateException(
							"missing leaf ⟨" + split.leftSize + "," + split.m + "," + split.p + "⟩"));
			NonCubicBilinearAlgorithm bot = lookup.find(split.rightSize, split.m, split.p)
					.orElseThrow(() -> new IllegalStateException(
							"missing leaf ⟨" + split.rightSize + "," + split.m + "," + split.p + "⟩"));
			return Compose.concatBelow(top, bot);
		}
	}
}
