package eu.solven.matmul.search;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.Recombination.SotaResolver;
import eu.solven.matmul.search.AnalyticalMaskSearch.SchemeSupports;

/**
 * The "support-pattern histogram" of a base scheme — the generic-formula
 * primitive behind fast recombination costing (2026-06-04 discussion).
 *
 * <p>Each of the {@code r} products of a base reads a <b>fixed</b> 6-tuple of
 * block index-sets: {@code (uRow, wRow)} on axis A, {@code (uCol, vRow)} on
 * axis B, {@code (vCol, wCol)} on axis C. Given an allocation {@code (a,b,c)},
 * {@link AnalyticalMaskSearch#shapesAt} computes that product's sub-shape as
 * {@code (min(max_{uRow}a, max_{wRow}a), min(max_{uCol}b, max_{vRow}b),
 * min(max_{vCol}c, max_{wCol}c))}. Two products with the <b>same</b> 6-tuple
 * therefore induce the same sub-shape for <em>every</em> allocation.</p>
 *
 * <p>So the cost depends on a product only through its 6-tuple. Grouping the
 * {@code r} products into the {@code D} distinct tuples (the
 * {@link Pattern}s, with multiplicities) turns per-allocation costing from
 * {@code O(r)} into {@code O(D)} — and {@code D ≪ r} for structured bases.
 * This single object underpins:</p>
 * <ul>
 *   <li>fast per-allocation / per-permutation multiset listing,</li>
 *   <li>a fast cost evaluator ({@link #cost}),</li>
 *   <li>the size→position assignment optimiser (the optimal ordering of a
 *       partition only sees the per-pattern maxes).</li>
 * </ul>
 *
 * <p>A base whose pattern-set is invariant under permuting its blocks is
 * "block-symmetric" (Strassen); one whose patterns break that symmetry
 * (Winograd, AlphaEvolve) makes the allocation <em>ordering</em> matter.</p>
 */
public final class BaseFingerprint {

	/** One distinct product support, with how many products share it. The six
	 *  arrays are sorted block-index sets; {@code count} is the multiplicity. */
	public record Pattern(int[] uRow, int[] wRow, int[] uCol, int[] vRow,
			int[] vCol, int[] wCol, int count) {}

	public final int n, m, p, r;
	public final List<Pattern> patterns;

	private BaseFingerprint(int n, int m, int p, int r, List<Pattern> patterns) {
		this.n = n; this.m = m; this.p = p; this.r = r; this.patterns = patterns;
	}

	public static BaseFingerprint of(NonCubicBilinearAlgorithm base) {
		return of(SchemeSupports.extract(base), base.n, base.m, base.p, base.r);
	}

	public static BaseFingerprint of(SchemeSupports sup, int n, int m, int p, int r) {
		// Group products by their full 6-tuple of index sets.
		LinkedHashMap<String, int[]> firstOf = new LinkedHashMap<>();
		Map<String, Integer> counts = new LinkedHashMap<>();
		for (int k = 0; k < r; k++) {
			String key = key(sup.uRowSupport[k]) + "|" + key(sup.wRowSupport[k]) + "|"
					+ key(sup.uColSupport[k]) + "|" + key(sup.vRowSupport[k]) + "|"
					+ key(sup.vColSupport[k]) + "|" + key(sup.wColSupport[k]);
			counts.merge(key, 1, Integer::sum);
			firstOf.putIfAbsent(key, new int[] { k });
		}
		List<Pattern> pats = new ArrayList<>();
		for (Map.Entry<String, int[]> e : firstOf.entrySet()) {
			int k = e.getValue()[0];
			pats.add(new Pattern(sup.uRowSupport[k], sup.wRowSupport[k], sup.uColSupport[k],
					sup.vRowSupport[k], sup.vColSupport[k], sup.wColSupport[k], counts.get(e.getKey())));
		}
		return new BaseFingerprint(n, m, p, r, pats);
	}

	/** Distinct support patterns {@code D} (≤ r). */
	public int distinctPatterns() { return patterns.size(); }

	/** Per-pattern sub-shape under allocation {@code (a,b,c)} — same formula as
	 *  {@link AnalyticalMaskSearch#shapesAt}, but evaluated once per pattern. */
	public int[] shapeOf(Pattern q, int[] a, int[] b, int[] c) {
		int subA = Math.min(maxIndexed(a, q.uRow()), maxIndexed(a, q.wRow()));
		int subB = Math.min(maxIndexed(b, q.uCol()), maxIndexed(b, q.vRow()));
		int subC = Math.min(maxIndexed(c, q.vCol()), maxIndexed(c, q.wCol()));
		return new int[] { subA, subB, subC };
	}

	/** Sub-shape multiset under {@code (a,b,c)} as {@code "n,m,p" → count},
	 *  computed in {@code O(D)} pattern evaluations rather than {@code O(r)}. */
	public Map<String, Integer> multiset(int[] a, int[] b, int[] c) {
		TreeMap<String, Integer> m = new TreeMap<>();
		for (Pattern q : patterns) {
			int[] s = shapeOf(q, a, b, c);
			m.merge(s[0] + "," + s[1] + "," + s[2], q.count(), Integer::sum);
		}
		return m;
	}

	/** Total recombination rank under {@code (a,b,c)} via the histogram. */
	public long cost(int[] a, int[] b, int[] c, SotaResolver sota) {
		long tot = 0;
		for (Pattern q : patterns) {
			int[] s = shapeOf(q, a, b, c);
			tot += (long) q.count() * sota.getRank(s[0], s[1], s[2]);
		}
		return tot;
	}

	private static int maxIndexed(int[] alloc, int[] indices) {
		int max = 0;
		for (int i : indices) if (alloc[i] > max) max = alloc[i];
		return max;
	}

	private static String key(int[] set) {
		StringBuilder sb = new StringBuilder();
		for (int v : set) sb.append(v).append('.');
		return sb.toString();
	}
}
