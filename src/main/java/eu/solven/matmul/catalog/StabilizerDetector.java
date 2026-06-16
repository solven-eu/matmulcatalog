package eu.solven.matmul.catalog;

import java.util.HashMap;
import java.util.Map;
import java.util.Optional;

import eu.solven.matmul.BilinearAlgorithm;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;

/**
 * Detect non-trivial stabilizer subgroups of a matmul algorithm under the
 * tensor's symmetry group. Currently focuses on the <b>cyclic Z/3</b>
 * stabilizer — the most informative non-trivial discrete symmetry, present in
 * Strassen's `⟨2,2,2⟩=7` algorithm.
 *
 * <p><b>Theory recap.</b> The matmul tensor for cubic {@code ⟨n,n,n⟩} has a
 * Z/3 cyclic action — rotating the three tensor slots (A, B, C) — once we
 * work in the <b>trilinear</b> convention (i.e. with W's row index transposed
 * relative to the matmul-tensor convention; see
 * {@link Verifier#transposeW(BilinearAlgorithm)}).</p>
 *
 * <p>A rank-{@code r} decomposition {@code Σ_k u_k ⊗ v_k ⊗ w_k} of the
 * trilin tensor is Z/3-stabilized iff there exists a permutation
 * {@code σ : [0, r) → [0, r)} with {@code σ³ = id} and signs
 * {@code ε_k ∈ {±1}} such that</p>
 *
 * <pre>
 *   ε_k · u_{σ(k)} = v_k,    ε_k · v_{σ(k)} = w_k,    ε_k · w_{σ(k)} = u_k.
 * </pre>
 *
 * <p>Equivalently: the multiset of columns
 * {@code (u_k, v_k, w_k)} of the decomposition is invariant under cyclic
 * rotation {@code (u, v, w) → (v, w, u)} (up to overall sign on each
 * matched column).</p>
 *
 * <p><b>Detection algorithm.</b> For each base column k compute a "shifted
 * signature" — the canonical form of {@code (v_k, w_k, u_k)}. For each base
 * column k' compute its "original signature" {@code (u_{k'}, v_{k'}, w_{k'})}.
 * If every shifted signature matches exactly one original signature and the
 * induced permutation σ satisfies σ³ = id, the algorithm has a Z/3 stabilizer.</p>
 *
 * <p>Canonical form for sign-agnostic matching: multiply the triple
 * {@code (u, v, w)} by {@code -1} for ALL THREE simultaneously if needed,
 * so the first non-zero entry of the concatenation is positive. This
 * captures the unique sign flip that preserves the outer product
 * {@code u ⊗ v ⊗ w = (-u) ⊗ (-v) ⊗ w}'s symmetric variants.</p>
 *
 * <p><b>Limitations.</b> Only Z/3 is detected today. The full stabilizer
 * group can include Z/2 (transpose-like) and higher; computing those
 * exhaustively is open work — see {@code SYMMETRIES.md}.</p>
 */
public final class StabilizerDetector {

	private StabilizerDetector() {}

	/**
	 * Result of a stabilizer check: either {@code "Z/3"} (cyclic-3-symmetric),
	 * {@code "trivial"} (no Z/3 found), or {@code "n/a"} (non-cubic or empty
	 * scheme — Z/3 check doesn't apply).
	 */
	public static String detectZ3(NonCubicBilinearAlgorithm alg) {
		if (alg.n != alg.m || alg.m != alg.p) return "n/a";
		if (alg.r == 0) return "n/a";

		// Put the algorithm in the trilinear convention so cyclic-Z/3 is the
		// natural slot rotation.
		BilinearAlgorithm cubic = alg.asCubic();
		BilinearAlgorithm trilin = Verifier.transposeW(cubic);

		int r = trilin.r;
		int n2 = trilin.n * trilin.n;

		// Build a map from canonical-signature → column index for the original
		// (u, v, w) triples.
		Map<String, Integer> originalIndex = new HashMap<>(r * 2);
		for (int k = 0; k < r; k++) {
			String sig = canonicalSignature(trilin.U, trilin.V, trilin.W, k, n2);
			Integer prev = originalIndex.putIfAbsent(sig, k);
			if (prev != null) {
				// Two columns share a canonical signature — multiset duplicates
				// imply ambiguity. We keep the lowest index; the assignment
				// step below may still succeed, but won't necessarily map
				// shifted ↔ original 1-1 with a unique σ.
			}
		}

		// For each column k, compute the canonical signature of (v_k, w_k, u_k);
		// the column that matches it (in the original index) is σ(k).
		int[] sigma = new int[r];
		for (int k = 0; k < r; k++) {
			String shiftedSig = canonicalSignature(trilin.V, trilin.W, trilin.U, k, n2);
			Integer target = originalIndex.get(shiftedSig);
			if (target == null) return "trivial";
			sigma[k] = target;
		}

		// Check σ is a permutation (every index in [0,r) is hit exactly once)
		// and that σ³ = id.
		boolean[] seen = new boolean[r];
		for (int k = 0; k < r; k++) {
			if (sigma[k] < 0 || sigma[k] >= r || seen[sigma[k]]) return "trivial";
			seen[sigma[k]] = true;
		}
		for (int k = 0; k < r; k++) {
			if (sigma[sigma[sigma[k]]] != k) return "trivial";
		}

		return "Z/3";
	}

	/**
	 * Canonicalize the column triple {@code (u, v, w)} so two triples that
	 * differ by an overall sign flip on all three vectors share a signature.
	 * Returns the string of concatenated coefficients, optionally negated.
	 */
	private static String canonicalSignature(double[][] U, double[][] V, double[][] W,
			int k, int n2) {
		StringBuilder raw = new StringBuilder(3 * n2 * 4);
		appendCol(raw, U, k, n2);
		raw.append('|');
		appendCol(raw, V, k, n2);
		raw.append('|');
		appendCol(raw, W, k, n2);
		String s = raw.toString();

		// Resolve global sign: pick the variant whose first non-zero
		// coefficient is positive.
		double firstNonZero = firstNonZeroAcrossUVW(U, V, W, k, n2);
		if (firstNonZero < 0) {
			return negateNumericString(s);
		}
		return s;
	}

	private static void appendCol(StringBuilder sb, double[][] M, int k, int n2) {
		for (int i = 0; i < n2; i++) {
			sb.append(formatCoef(M[i][k]));
			sb.append(',');
		}
	}

	private static String formatCoef(double v) {
		if (v == 0.0) return "0";
		double rounded = Math.rint(v);
		if (Math.abs(v - rounded) < 1e-9) {
			return Long.toString((long) rounded);
		}
		return Double.toString(v);
	}

	private static double firstNonZeroAcrossUVW(double[][] U, double[][] V, double[][] W,
			int k, int n2) {
		for (int i = 0; i < n2; i++) if (U[i][k] != 0.0) return U[i][k];
		for (int i = 0; i < n2; i++) if (V[i][k] != 0.0) return V[i][k];
		for (int i = 0; i < n2; i++) if (W[i][k] != 0.0) return W[i][k];
		return 0.0;
	}

	private static String negateNumericString(String s) {
		StringBuilder out = new StringBuilder(s.length() + 16);
		int i = 0;
		while (i < s.length()) {
			char c = s.charAt(i);
			if (c == '|' || c == ',') {
				out.append(c);
				i++;
				continue;
			}
			// Parse a numeric token from i
			int j = i;
			while (j < s.length() && s.charAt(j) != ',' && s.charAt(j) != '|') j++;
			String token = s.substring(i, j);
			if (token.equals("0")) {
				out.append("0");
			} else if (token.startsWith("-")) {
				out.append(token.substring(1));
			} else {
				out.append('-').append(token);
			}
			i = j;
		}
		return out.toString();
	}

	/**
	 * Convenience: returns {@link Optional#empty()} for "n/a" (non-cubic or empty),
	 * or a present {@code "Z/3"} / {@code "trivial"} for genuine checks.
	 */
	public static Optional<String> stabilizerOrEmpty(NonCubicBilinearAlgorithm alg) {
		String r = detectZ3(alg);
		return "n/a".equals(r) ? Optional.empty() : Optional.of(r);
	}
}
