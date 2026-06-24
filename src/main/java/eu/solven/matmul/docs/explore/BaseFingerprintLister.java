package eu.solven.matmul.docs.explore;

import eu.solven.matmul.recombination.Recombination;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.recombination.Recombination.SotaResolver;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.search.BaseFingerprint;
import eu.solven.matmul.search.PerAxisFingerprint;

/**
 * Lists the {@link BaseFingerprint} (support-pattern histogram) of several
 * bases, so the (a)symmetry that makes allocation-ordering matter is visible:
 * <ul>
 *   <li>{@code r} products → {@code D} distinct support patterns (compression);</li>
 *   <li>the pattern multiplicity histogram;</li>
 *   <li>a permutation-sensitivity probe: over all orderings of one balanced
 *       partition (axis A; B,C balanced), how many DISTINCT sub-shape multisets
 *       (and cost values) arise. 1 ⇒ block-symmetric on that partition; >1 ⇒
 *       ordering matters.</li>
 * </ul>
 */
public final class BaseFingerprintLister {

	private BaseFingerprintLister() {}

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

	private static List<int[]> orderings(int[] parts) {
		List<int[]> out = new ArrayList<>();
		permute(parts.clone(), 0, out);
		List<int[]> uniq = new ArrayList<>();
		Set<String> seen = new HashSet<>();
		for (int[] x : out) if (seen.add(java.util.Arrays.toString(x))) uniq.add(x);
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

	private static int[] balanced(int total, int parts) {
		int[] a = new int[parts];
		int base = total / parts, rem = total % parts;
		for (int i = 0; i < parts; i++) a[i] = base + (i < rem ? 1 : 0);
		return a;
	}

	private static void report(String label, String path, int target, SotaResolver sota) throws Exception {
		NonCubicBilinearAlgorithm base = SchemeIO.read(eu.solven.matmul.catalog.SchemeResolver.byHint(path));
		BaseFingerprint fp = BaseFingerprint.of(base);
		System.out.printf("%n=== %s : ⟨%d,%d,%d⟩=%d ===%n", label, fp.n, fp.m, fp.p, fp.r);
		System.out.printf("products r = %d   distinct support patterns D = %d   (compression r/D = %.2f×)%n",
				fp.r, fp.distinctPatterns(), fp.r / (double) fp.distinctPatterns());

		// Per-axis fingerprint: the granularity that actually compresses.
		PerAxisFingerprint pa = PerAxisFingerprint.of(base);
		System.out.printf("per-axis distinct groups: Dₐ=%d  D_b=%d  D_c=%d  (of r=%d → axis compression %.2f/%.2f/%.2f×)%n",
				pa.distinctA(), pa.distinctB(), pa.distinctC(), pa.r,
				pa.r / (double) pa.distinctA(), pa.r / (double) pa.distinctB(), pa.r / (double) pa.distinctC());

		// Multiplicity histogram: how many patterns occur with each count.
		TreeMap<Integer, Integer> mult = new TreeMap<>();
		for (BaseFingerprint.Pattern q : fp.patterns) mult.merge(q.count(), 1, Integer::sum);
		StringBuilder mh = new StringBuilder();
		mult.forEach((cnt, howMany) -> mh.append(howMany).append("×[count=").append(cnt).append("] "));
		System.out.println("pattern multiplicities: " + mh.toString().trim());

		// Permutation sensitivity at the balanced partition of `target`.
		int[] bal = balanced(target, fp.n);
		int[] balB = balanced(target, fp.m), balC = balanced(target, fp.p);
		Set<String> distinctMultisets = new HashSet<>();
		Set<Long> distinctCosts = new HashSet<>();
		long mn = Long.MAX_VALUE, mx = Long.MIN_VALUE;
		List<int[]> ords = orderings(bal);
		for (int[] a : ords) {
			Map<String, Integer> ms = fp.multiset(a, balB, balC);
			distinctMultisets.add(ms.toString());
			long c = fp.cost(a, balB, balC, sota);
			distinctCosts.add(c);
			mn = Math.min(mn, c); mx = Math.max(mx, c);
		}
		System.out.printf("perm-sensitivity @ balanced partition %s of ⟨%d³⟩ (B,C balanced):%n",
				java.util.Arrays.toString(bal), target);
		System.out.printf("  %d axis-A orderings → %d distinct multisets, %d distinct costs (min=%d max=%d Δ=%d)%n",
				ords.size(), distinctMultisets.size(), distinctCosts.size(), mn, mx, mx - mn);
		System.out.printf("  verdict: %s%n", distinctMultisets.size() == 1
				? "BLOCK-SYMMETRIC on this partition (ordering is free)"
				: "block-ASYMMETRIC (ordering matters → assignment problem)");
	}

	public static void main(String[] args) throws Exception {
		FieldAwareLookup lk = new FieldAwareLookup("R");
		SotaResolver sota = sota(lk);
		// label, path, target-cube for the permutation probe
		report("Strassen",   "src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json",       9, sota);
		report("Winograd",   "src/main/resources/schemes/known/section2/winograd_1971-2x2x2_m7_a24.json",  9, sota);
		report("Laderman",   "src/main/resources/schemes/known/section3/laderman_1976-3x3x3_m23_a98.json", 7, sota);
		report("AlphaEvolve","src/main/resources/schemes/known/section5/alphaevolve-5x5x5_m93_a846.json", 16, sota);
	}
}
