package eu.solven.matmul.recombination;

import java.util.Arrays;

/**
 * Pan's "pair simultaneous sub-products" trick (DIS09 §3, Layer 4).
 *
 * <p>Two sub-products of shapes {@code ⟨a,b,c⟩} and {@code ⟨b,c,a⟩}
 * (cyclic rotation) can be computed together using
 * {@code abc + ab + bc + ca} multiplications, instead of the
 * {@code 2·R(⟨a,b,c⟩)} needed when computed independently. Saving
 * applies when {@code 2·R(⟨a,b,c⟩) > abc + ab + bc + ca}, which
 * happens for "thick" shapes where the catalog can't drive the rank
 * far below {@code abc}.</p>
 *
 * <p>Reference: DIS09 page 15:</p>
 * <blockquote>
 * For any pair {@code ℓ &lt; ℓ' ≤ L}, Pan's technique is applicable
 * if either {@code (μ_ℓ, ν_ℓ, π_ℓ) = (ν_ℓ', π_ℓ', μ_ℓ')} or
 * {@code (μ_ℓ, ν_ℓ, π_ℓ) = (π_ℓ', μ_ℓ', ν_ℓ')}. In this case, we can
 * compute both products using {@code μ_ℓ ν_ℓ π_ℓ + μ_ℓ ν_ℓ + μ_ℓ π_ℓ
 * + ν_ℓ π_ℓ} base ring multiplications.
 * </blockquote>
 */
public final class PairedSubProducts {

	private PairedSubProducts() {}

	/**
	 * Pan's pair-product cost for two sub-products with shapes cyclically
	 * related ({@code ⟨a,b,c⟩} + {@code ⟨b,c,a⟩}).
	 */
	public static long pairCost(int a, int b, int c) {
		return (long) a * b * c + (long) a * b + (long) b * c + (long) c * a;
	}

	/**
	 * Greedy max-savings matching over the multiset of sub-product
	 * shapes. Returns the total rank after applying as many pair
	 * substitutions as profitable.
	 *
	 * @param subShapes each row is a {@code [n, m, p]} sub-product shape
	 * @param sota      per-shape rank lookup (used for both individual
	 *                  and savings calculations)
	 * @return total rank after pairing optimisation
	 */
	public static long applyPairing(int[][] subShapes, Recombination.SotaResolver sota) {
		int n = subShapes.length;
		long[] individualCost = new long[n];
		for (int i = 0; i < n; i++) {
			individualCost[i] = sota.getRank(subShapes[i][0], subShapes[i][1], subShapes[i][2]);
		}

		// Enumerate candidate pairs (i, j) where shapes are cyclically related and
		// pairing actually saves mults.
		int[][] cand = new int[n * (n - 1) / 2][3]; // (i, j, savings)
		int candCount = 0;
		for (int i = 0; i < n; i++) {
			int[] si = subShapes[i];
			if (si[0] == 0 || si[1] == 0 || si[2] == 0) continue; // degenerate
			for (int j = i + 1; j < n; j++) {
				int[] sj = subShapes[j];
				if (sj[0] == 0 || sj[1] == 0 || sj[2] == 0) continue;
				if (!panPairable(si, sj)) continue;
				long pair = pairCost(si[0], si[1], si[2]);
				long indiv = individualCost[i] + individualCost[j];
				long savings = indiv - pair;
				if (savings <= 0) continue;
				cand[candCount++] = new int[] { i, j, (int) Math.min(savings, Integer.MAX_VALUE) };
			}
		}

		// Greedy: pick pairs by largest savings first, skipping already-matched indices.
		int[][] active = Arrays.copyOf(cand, candCount);
		Arrays.sort(active, (a, b) -> Integer.compare(b[2], a[2]));
		boolean[] paired = new boolean[n];
		long total = 0;
		for (int[] p : active) {
			int i = p[0], j = p[1];
			if (paired[i] || paired[j]) continue;
			paired[i] = paired[j] = true;
			total += pairCost(subShapes[i][0], subShapes[i][1], subShapes[i][2]);
		}
		for (int i = 0; i < n; i++) {
			if (!paired[i]) total += individualCost[i];
		}
		return total;
	}

	/**
	 * {@code true} iff {@code b} is one of the three cyclic rotations of
	 * {@code a}: {@code (x,y,z)}, {@code (y,z,x)}, or {@code (z,x,y)}. Pure
	 * rotation relation (includes the identity); for the Pan-<em>fusion</em>
	 * predicate use {@link #panPairable}.
	 */
	public static boolean cyclicallyEquivalent(int[] a, int[] b) {
		return (a[0] == b[0] && a[1] == b[1] && a[2] == b[2])
				|| (a[0] == b[1] && a[1] == b[2] && a[2] == b[0])
				|| (a[0] == b[2] && a[1] == b[0] && a[2] == b[1]);
	}

	/**
	 * {@code true} iff two products of shapes {@code s} and {@code t} can be
	 * fused by Pan's pair trick. Pan computes one {@code ⟨a,b,c⟩} <em>and</em>
	 * one {@code ⟨b,c,a⟩} (a shape + a NON-trivial cyclic rotation). Crucially,
	 * <strong>transpose is a free isotropy</strong> — {@code ⟨a,b,c⟩} can be
	 * computed as {@code ⟨c,b,a⟩} via {@code Cᵀ = BᵀAᵀ} at zero arithmetic cost —
	 * so each product may be transposed first. Consequences:
	 * <ul>
	 *   <li>cubic {@code ⟨n,n,n⟩}: always pairable;</li>
	 *   <li>exactly two equal dims, e.g. {@code ⟨8,9,9⟩}: pairable — transpose
	 *       turns it into its rotation {@code ⟨9,9,8⟩};</li>
	 *   <li>all dims distinct, e.g. {@code ⟨8,9,10⟩}: NOT pairable — transpose
	 *       gives {@code ⟨10,9,8⟩}, still not a rotation (and cyclic rotation is
	 *       not a free isotropy of a bilinear product).</li>
	 * </ul>
	 * {@link #pairCost} is fully symmetric in {@code (a,b,c)}, so the fused cost
	 * is orientation-independent.
	 */
	public static boolean panPairable(int[] s, int[] t) {
		int[][] sv = { s, { s[2], s[1], s[0] } };
		int[][] tv = { t, { t[2], t[1], t[0] } };
		for (int[] a : sv) {
			for (int[] b : tv) {
				// b is a non-trivial rotation of a (identity matches only when cubic)
				if ((b[0] == a[1] && b[1] == a[2] && b[2] == a[0])
						|| (b[0] == a[2] && b[1] == a[0] && b[2] == a[1])) {
					return true;
				}
			}
		}
		return false;
	}
}
