package eu.solven.matmul.search;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.recombination.Recombination.SotaResolver;
import eu.solven.matmul.recombination.AnalyticalMaskSearch.SchemeSupports;

/**
 * Per-axis support fingerprint — the granularity that actually compresses
 * (2026-06-04 follow-up to {@link BaseFingerprint}, where grouping by the full
 * 6-tuple gave {@code D = r} for dense schemes).
 *
 * <p>A product's sub-dim on a single axis depends on that axis ONLY: on axis A
 * it is {@code min(max_{uRow}a, max_{wRow}a)} — a function of the pair
 * {@code (uRow, wRow)} alone. Many products share that pair even when their full
 * tuples differ, so the per-axis distinct-group count {@code Dₐ} is typically
 * {@code ≪ r}. Grouping per axis lets us compute every product's axis sub-dim in
 * {@code O(Dₐ · |support|)} instead of {@code O(r · |support|)} — and, crucially,
 * it is the substrate for the size→position <b>assignment</b> optimiser: when we
 * permute one axis (other axes fixed) only that axis's groups re-evaluate.</p>
 */
public final class PerAxisFingerprint {

	/** Distinct support pair on one axis: {@code (s1, s2)} are the two index
	 *  sets whose maxes are min-combined, {@code members} the products sharing
	 *  this pair. The axis sub-dim is {@code min(max_{s1}alloc, max_{s2}alloc)}. */
	public record AxisGroup(int[] s1, int[] s2, int[] members) {}

	public final int n, m, p, r;
	public final List<AxisGroup> axisA, axisB, axisC;

	private PerAxisFingerprint(int n, int m, int p, int r,
			List<AxisGroup> axisA, List<AxisGroup> axisB, List<AxisGroup> axisC) {
		this.n = n; this.m = m; this.p = p; this.r = r;
		this.axisA = axisA; this.axisB = axisB; this.axisC = axisC;
	}

	public static PerAxisFingerprint of(NonCubicBilinearAlgorithm base) {
		return of(SchemeSupports.extract(base), base.n, base.m, base.p, base.r);
	}

	public static PerAxisFingerprint of(SchemeSupports sup, int n, int m, int p, int r) {
		return new PerAxisFingerprint(n, m, p, r,
				group(sup.uRowSupport, sup.wRowSupport, r),   // axis A: rows of A from U and W
				group(sup.uColSupport, sup.vRowSupport, r),   // axis B: cols of A (U) & rows of B (V)
				group(sup.vColSupport, sup.wColSupport, r));  // axis C: cols of B (V) & cols of C (W)
	}

	/** Group the r products by their (s1[k], s2[k]) index-set pair on one axis. */
	private static List<AxisGroup> group(int[][] s1, int[][] s2, int r) {
		LinkedHashMap<String, List<Integer>> byKey = new LinkedHashMap<>();
		for (int k = 0; k < r; k++) {
			String key = key(s1[k]) + "|" + key(s2[k]);
			byKey.computeIfAbsent(key, x -> new ArrayList<>()).add(k);
		}
		List<AxisGroup> out = new ArrayList<>();
		for (List<Integer> members : byKey.values()) {
			int k0 = members.get(0);
			int[] mem = members.stream().mapToInt(Integer::intValue).toArray();
			out.add(new AxisGroup(s1[k0], s2[k0], mem));
		}
		return out;
	}

	/** Fill {@code out[k]} = axis sub-dim of product k under {@code alloc},
	 *  evaluating each distinct group once. */
	public void subDims(List<AxisGroup> groups, int[] alloc, int[] out) {
		for (AxisGroup g : groups) {
			int d = Math.min(maxIndexed(alloc, g.s1()), maxIndexed(alloc, g.s2()));
			for (int k : g.members()) out[k] = d;
		}
	}

	/** Total recombination rank under {@code (a,b,c)}: per-axis sub-dims via the
	 *  grouped fingerprint, then one rank lookup per product. Identical result to
	 *  {@link AnalyticalMaskSearch#shapesAt}-based costing. */
	public long cost(int[] a, int[] b, int[] c, SotaResolver sota) {
		int[] sa = new int[r], sb = new int[r], sc = new int[r];
		subDims(axisA, a, sa);
		subDims(axisB, b, sb);
		subDims(axisC, c, sc);
		long tot = 0;
		for (int k = 0; k < r; k++) tot += sota.getRank(sa[k], sb[k], sc[k]);
		return tot;
	}

	public int distinctA() { return axisA.size(); }
	public int distinctB() { return axisB.size(); }
	public int distinctC() { return axisC.size(); }

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
