package eu.solven.matmul.docs.explore;

import org.apache.commons.math3.fraction.BigFraction;

/**
 * Attribute the SZ aggregation families' {@code 2.0·N²} tail to index-coincidence
 * classes. Families (a) Ŝ and (b) Ṡ\Ṡ1 emit one product per index-triple
 * {@code (i,j,k) ∈ [d]³}, {@code d=N/2+1}. Each product's count is a cubic in d;
 * bucketing by how many of {i,j,k} are DISTINCT (3 = generic, 2 = an edge, 1 = the
 * all-equal diagonal) and fitting each bucket exactly (finite differences) tells us
 * which triples carry the {@code N²} coefficient — the only lever not blocked by a
 * rank lower bound (correction, {@code 1.75}) or additivity.
 *
 * <p>Prints, per (family, #distinct) bucket, the exact d-polynomial and its
 * contribution to the {@code N²} coefficient of {@code tNew} (via d=N/2+1). The
 * generic 3-distinct bucket carries the {@code N³/3} leading term; its own {@code N²}
 * sub-leading and the 2-distinct bucket's leading are what sum to {@code 2.0}.</p>
 */
public final class ProbeAggTail {
	private ProbeAggTail() {}

	public static void main(String[] args) {
		// bucket[family][#distinct] as a function of d; family 0=(a), 1=(b).
		// Fit each with 4 sample points d=0..3 (cubic), verify with d=4,5.
		String[] fam = { "(a) Ŝ-symm ", "(b) Ṡ-barred" };
		BigFraction n2total = BigFraction.ZERO;
		System.out.println("bucket                     d-polynomial (a0+a1 d+a2 d²+a3 d³)         N²-coeff");
		for (int f = 0; f < 2; f++) {
			for (int nd = 3; nd >= 1; nd--) {
				long[] vals = new long[6];
				for (int d = 0; d < 6; d++) vals[d] = count(f, nd, d);
				BigFraction[] poly = cubicFit(vals);          // [a0,a1,a2,a3]
				// verify exactness at d=4,5
				for (int d = 4; d < 6; d++) {
					if (!eval(poly, d).equals(new BigFraction(vals[d]))) {
						throw new IllegalStateException("non-cubic bucket f=" + f + " nd=" + nd);
					}
				}
				// N²-coeff after d = N/2+1 : a2·(1/4) + a3·(3/4)
				BigFraction n2 = poly[2].multiply(new BigFraction(1, 4))
						.add(poly[3].multiply(new BigFraction(3, 4)));
				n2total = n2total.add(n2);
				System.out.printf("%s  %d-distinct : %-40s  %s%n",
						fam[f], nd, polyStr(poly), n2);
			}
		}
		System.out.println("--------");
		System.out.printf("families (a)+(b) total N² coefficient = %s  (expected 2.0 = 8/4)%n", n2total);
		System.out.println("[correction family (c) adds 7/4=1.75, family (d) is O(N); total tNew N² = 15/4]");
	}

	/** #products family f (0=a,1=b) emits for triples with exactly {@code nd} distinct indices, at dim d. */
	private static long count(int f, int nd, int d) {
		long c = 0;
		for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) for (int k = 0; k < d; k++) {
			int dist = distinct(i, j, k);
			if (dist != nd) continue;
			if (f == 0) {                                  // family (a): Ŝ condition, emits triple + bar
				if ((i <= j && j < k) || (k < j && j <= i)) c += 2;
			} else {                                       // family (b): barred always; unbarred unless all-equal
				c += 1;                                     // barred partner (same #distinct)
				if (dist > 1) c += 1;                       // unbarred, excluded only for i=j=k
			}
		}
		return c;
	}

	private static int distinct(int i, int j, int k) {
		int c = 1;
		if (j != i) c++;
		if (k != i && k != j) c++;
		return c;
	}

	/** Exact cubic through (0,v0)..(3,v3) via Newton forward differences → [a0,a1,a2,a3]. */
	private static BigFraction[] cubicFit(long[] v) {
		BigFraction f0 = new BigFraction(v[0]), f1 = new BigFraction(v[1]),
				f2 = new BigFraction(v[2]), f3 = new BigFraction(v[3]);
		BigFraction d1 = f1.subtract(f0);
		BigFraction d2 = f2.subtract(f1.multiply(2)).add(f0);
		BigFraction d3 = f3.subtract(f2.multiply(3)).add(f1.multiply(3)).subtract(f0);
		// falling-factorial → standard: a3=Δ³/6, a2=Δ²/2−Δ³/2, a1=Δ¹−Δ²/2+Δ³/3, a0=f0
		BigFraction a3 = d3.divide(6);
		BigFraction a2 = d2.divide(2).subtract(d3.divide(2));
		BigFraction a1 = d1.subtract(d2.divide(2)).add(d3.divide(3));
		return new BigFraction[] { f0, a1, a2, a3 };
	}

	private static BigFraction eval(BigFraction[] p, int d) {
		BigFraction x = new BigFraction(d), r = BigFraction.ZERO, pw = BigFraction.ONE;
		for (BigFraction a : p) { r = r.add(a.multiply(pw)); pw = pw.multiply(x); }
		return r;
	}

	private static String polyStr(BigFraction[] p) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < p.length; i++) {
			if (p[i].equals(BigFraction.ZERO)) continue;
			if (sb.length() > 0) sb.append(" + ");
			sb.append(p[i]).append(i == 0 ? "" : i == 1 ? " d" : " d^" + i);
		}
		return sb.length() == 0 ? "0" : sb.toString();
	}
}
