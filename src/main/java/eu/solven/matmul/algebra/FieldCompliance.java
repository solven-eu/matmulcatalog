package eu.solven.matmul.algebra;

import java.util.ArrayList;
import java.util.List;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Checks whether the coefficients of a scheme's factor matrices
 * actually live in the claimed {@link Field}. Complements
 * {@code Verifier.passesRandomMatmulSpotCheck} (which checks the
 * scheme COMPUTES matmul) by validating the FIELD CLAIM independently.
 *
 * <p>This is the missing check that would have prevented the
 * AT-Z mis-tagging at import time: a scheme using {@code 1/2}
 * coefficients claiming to be {@code Z} would fail
 * {@link #checkAllInField} with a {@link Discrepancy} pointing at the
 * specific (matrix, row, column) where the violation lives.</p>
 *
 * <p>Caveat: a scheme that USES integer coefficients but is
 * algorithmically wrong (the brokenness our denylist catches) will
 * <em>pass</em> the field-compliance check — that's only caught by
 * the matmul verifier. The two checks are complementary.</p>
 */
public final class FieldCompliance {

	private FieldCompliance() {}

	/** Coefficient floating-point tolerance for "is integer". */
	private static final double EPS = 1e-12;

	/**
	 * Where in the scheme a field-compliance violation occurs.
	 *
	 * @param factor   "U", "V", or "W"
	 * @param row      row index in the flat factor matrix
	 * @param col      column index (= product index)
	 * @param value    the offending coefficient
	 * @param reason   short human-readable reason
	 */
	public record Discrepancy(String factor, int row, int col, double value, String reason) {
		@Override public String toString() {
			return String.format("%s[%d][%d] = %s — %s", factor, row, col, value, reason);
		}
	}

	/**
	 * Returns the list of coefficient-field violations in
	 * {@code scheme} for the claimed {@code field}. Empty list means
	 * every coefficient is consistent with the claim.
	 *
	 * @param scheme       the algorithm to inspect
	 * @param claimedField what the scheme claims to be valid over
	 * @param maxReport    cap the result list at this many discrepancies
	 *                     (set Integer.MAX_VALUE for full enumeration)
	 */
	public static List<Discrepancy> checkAllInField(
			NonCubicBilinearAlgorithm scheme, Field claimedField, int maxReport) {
		List<Discrepancy> out = new ArrayList<>();
		scan(scheme.denseU(), "U", claimedField, out, maxReport);
		if (out.size() < maxReport) scan(scheme.denseV(), "V", claimedField, out, maxReport);
		if (out.size() < maxReport) scan(scheme.denseW(), "W", claimedField, out, maxReport);
		return out;
	}

	/** Convenience: first 20 violations. */
	public static List<Discrepancy> checkAllInField(
			NonCubicBilinearAlgorithm scheme, Field claimedField) {
		return checkAllInField(scheme, claimedField, 20);
	}

	/** True iff the scheme has zero coefficient-field violations. */
	public static boolean isCompliant(NonCubicBilinearAlgorithm scheme, Field claimedField) {
		return checkAllInField(scheme, claimedField, 1).isEmpty();
	}

	private static void scan(double[][] M, String label, Field f,
			List<Discrepancy> out, int maxReport) {
		for (int row = 0; row < M.length && out.size() < maxReport; row++) {
			for (int col = 0; col < M[row].length && out.size() < maxReport; col++) {
				double v = M[row][col];
				String reason = violation(v, f);
				if (reason != null) out.add(new Discrepancy(label, row, col, v, reason));
			}
		}
	}

	/** Returns a violation reason if {@code v} is not in {@code f}, else null. */
	private static String violation(double v, Field f) {
		if (Double.isNaN(v) || Double.isInfinite(v)) return "NaN or infinite";
		switch (f) {
			case Z:
				if (Math.abs(v - Math.rint(v)) > EPS) {
					return "not integer (rounding error " + Math.abs(v - Math.rint(v)) + ")";
				}
				return null;
			case Q:
			case R:
				// Any finite value is acceptable — we trust double precision to
				// represent reasonable rationals/reals.
				return null;
			case C:
				// Real-valued component check — complex schemes use a separate
				// {@code ComplexNonCubicBilinearAlgorithm} type; if you reach
				// here with a non-finite value it's a corrupt file.
				return null;
			case F2:
				// In F₂, only 0 and 1 are allowed.
				// -1 and 1 both reduce to 1 mod 2 BUT canonical F₂ schemes use {0, 1}.
				// We accept {0, ±1} since some imports leave the sign in.
				double a = Math.abs(v);
				if (a > EPS && Math.abs(a - 1.0) > EPS) {
					return "value " + v + " is not in {-1, 0, 1} (would not reduce cleanly mod 2)";
				}
				return null;
			case F3:
				// In F₃, allowed values reduce to {0, 1, 2} mod 3. Accept {-2, -1, 0, 1, 2}.
				double b = Math.abs(v);
				if (b > EPS && Math.abs(b - 1.0) > EPS && Math.abs(b - 2.0) > EPS) {
					return "value " + v + " is not in {-2, -1, 0, 1, 2}";
				}
				return null;
			default:
				return "unsupported field " + f;
		}
	}
}
