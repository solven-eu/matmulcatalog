package eu.solven.matmul.docs.explore;

import eu.solven.matmul.catalog.Recombination;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeMap;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Recombination.SotaResolver;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.search.AnalyticalMaskSearch;
import eu.solven.matmul.search.AnalyticalMaskSearch.SchemeSupports;

/**
 * Anti-symmetry illustration (2026-06-04): does the recombination cost depend on
 * the <em>ordering</em> of an allocation, for a fixed base scheme?
 *
 * <p>If the base were invariant under permuting its blocks, every ordering of a
 * given partition would induce the SAME multiset of per-product sub-shapes, so
 * the same cost — and the {@code compositions → partitions} reduction would be
 * valid. We measure the spread across orderings directly.</p>
 *
 * <p>The free row/col relabeling argument says permuting an allocation is the
 * same as permuting the base's blocks the other way:
 * {@code cost(B, σ·a) = cost(σ⁻¹·B, a)}. So the spread we print below is ALSO
 * the spread of cost over the base's block-relabelings — i.e. exactly the
 * choice a free I/O permutation gets to make. A non-zero spread means the
 * ordering matters for a single fixed base, and the "free" cost of a partition
 * is the MIN over its orderings.</p>
 */
public final class AllocationSymmetryProbe {

	private AllocationSymmetryProbe() {}

	private static SotaResolver sota(FieldAwareLookup lk) {
		return (p, q, r) -> {
			if (p == 0 || q == 0 || r == 0) return 0;
			if (p == 1) return q * r;
			if (q == 1) return p * r;
			if (r == 1) return p * q;
			int v = lk.findRank(p, q, r);
			return v >= Recombination.SotaResolver.UNKNOWN_RANK ? p * q * r : v;
		};
	}

	private static long cost(SchemeSupports sup, SotaResolver sota, int[] a, int[] b, int[] c) {
		int[][] shapes = AnalyticalMaskSearch.shapesAt(sup, a, b, c);
		long tot = 0;
		for (int[] s : shapes) tot += sota.getRank(s[0], s[1], s[2]);
		return tot;
	}

	/** Sorted "shape×count" summary of the per-product sub-shape multiset. */
	private static String multiset(SchemeSupports sup, int[] a, int[] b, int[] c) {
		int[][] shapes = AnalyticalMaskSearch.shapesAt(sup, a, b, c);
		TreeMap<String, Integer> m = new TreeMap<>();
		for (int[] s : shapes) m.merge(s[0] + "x" + s[1] + "x" + s[2], 1, Integer::sum);
		StringBuilder sb = new StringBuilder();
		m.forEach((k, v) -> sb.append(v).append("·⟨").append(k.replace('x', ',')).append("⟩  "));
		return sb.toString().trim();
	}

	/** Distinct orderings (permutations) of a multiset of part sizes. */
	private static List<int[]> orderings(int[] parts) {
		List<int[]> out = new ArrayList<>();
		permute(parts.clone(), 0, out);
		// dedup
		List<int[]> uniq = new ArrayList<>();
		java.util.Set<String> seen = new java.util.HashSet<>();
		for (int[] p : out) {
			String k = java.util.Arrays.toString(p);
			if (seen.add(k)) uniq.add(p);
		}
		return uniq;
	}

	private static void permute(int[] a, int k, List<int[]> out) {
		if (k == a.length) { out.add(a.clone()); return; }
		for (int i = k; i < a.length; i++) {
			int t = a[k]; a[k] = a[i]; a[i] = t;
			permute(a, k + 1, out);
			t = a[k]; a[k] = a[i]; a[i] = t;
		}
	}

	private static String fmt(int[] a) {
		StringBuilder sb = new StringBuilder("(");
		for (int i = 0; i < a.length; i++) { if (i > 0) sb.append(','); sb.append(a[i]); }
		return sb.append(')').toString();
	}

	public static void main(String[] args) throws Exception {
		FieldAwareLookup lk = new FieldAwareLookup("R");
		SotaResolver sota = sota(lk);

		// ============================================================
		// CASE 1 — small, eyeball-able: Laderman ⟨3,3,3⟩=23 on ⟨4,4,4⟩.
		// 4 = 2+1+1, so the ONLY partition into 3 parts is {2,1,1}, with
		// exactly 3 orderings. We vary axis A's ordering, hold B=C=(2,1,1).
		// ============================================================
		NonCubicBilinearAlgorithm lad =
				SchemeIO.read(eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section3/laderman_1976-3x3x3_m23_a98.json"));
		SchemeSupports supL = SchemeSupports.extract(lad);
		System.out.println("================ CASE 1: Laderman ⟨3,3,3⟩=23  on target ⟨4,4,4⟩ ================");
		System.out.println("partition of axis A = {2,1,1} (the only one); B=C fixed at (2,1,1).");
		System.out.println("If the base were block-symmetric, all 3 orderings would give the SAME multiset.\n");
		int[] bBal = { 2, 1, 1 }, cBal = { 2, 1, 1 };
		long minL = Long.MAX_VALUE, maxL = Long.MIN_VALUE;
		for (int[] a : orderings(new int[] { 2, 1, 1 })) {
			long cst = cost(supL, sota, a, bBal, cBal);
			minL = Math.min(minL, cst); maxL = Math.max(maxL, cst);
			System.out.printf("  A=%-9s  cost=%-4d  shapes: %s%n", fmt(a), cst, multiset(supL, a, bBal, cBal));
		}
		System.out.printf("%n  → spread over orderings: min=%d  max=%d  Δ=%d%n", minL, maxL, maxL - minL);
		System.out.println("  (Δ>0 ⇒ ordering matters ⇒ base is block-ASYMMETRIC ⇒ partitions DON'T collapse cost.)");
		System.out.println("  (the 'free I/O permutation' lets us pick the MIN ordering = " + minL
				+ " — but to FIND it we still scan orderings.)\n");

		// ============================================================
		// CASE 2 — the real target: AlphaEvolve ⟨5,5,5⟩=93 on ⟨16,16,16⟩.
		// Measure spread over orderings of two partitions, B=C balanced.
		// ============================================================
		NonCubicBilinearAlgorithm ae =
				SchemeIO.read(eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section5/alphaevolve-5x5x5_m93_a846.json"));
		SchemeSupports supA = SchemeSupports.extract(ae);
		int[] bal16 = { 4, 3, 3, 3, 3 };
		System.out.println("============ CASE 2: AlphaEvolve ⟨5,5,5⟩=93  on target ⟨16,16,16⟩ ============");
		System.out.println("vary axis A over orderings of a partition; B=C fixed balanced (4,3,3,3,3).\n");
		for (int[] part : new int[][] { { 4, 3, 3, 3, 3 }, { 8, 4, 2, 1, 1 }, { 12, 1, 1, 1, 1 } }) {
			List<int[]> ords = orderings(part);
			long mn = Long.MAX_VALUE, mx = Long.MIN_VALUE; int[] best = null;
			for (int[] a : ords) {
				long cst = cost(supA, sota, a, bal16, bal16);
				if (cst < mn) { mn = cst; best = a; }
				mx = Math.max(mx, cst);
			}
			System.out.printf("  partition %-13s  %3d orderings  cost min=%-5d max=%-5d Δ=%-4d  best A=%s%n",
					fmt(part), ords.size(), mn, mx, mx - mn, fmt(best));
		}
		System.out.println();
		System.out.println("Reading: Δ>0 confirms AlphaEvolve is block-asymmetric (settles the earlier caveat),");
		System.out.println("and that a free I/O permutation only buys the BEST ordering of a partition — the");
		System.out.println("'best ordering' is an assignment problem (size→position), the real lever for speed.");
	}
}
