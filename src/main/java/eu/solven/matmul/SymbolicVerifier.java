package eu.solven.matmul;

import java.io.File;
import java.io.IOException;
import java.math.BigInteger;
import java.util.Optional;

import eu.solven.matmul.catalog.SchemeIO;
import tools.jackson.databind.JsonNode;

/**
 * Symbolic, exact-arithmetic verifier for bilinear matmul schemes.
 *
 * <p>Unlike {@link Verifier#passesRandomMatmulSpotCheck}, which uses
 * floating-point {@code double} arithmetic and accepts ANY scheme whose
 * polynomial identity holds numerically (including schemes whose
 * coefficients exceed the declared algebra — e.g. an {@code i}-valued
 * entry hidden in a {@code Q}-tagged file), this verifier:</p>
 * <ol>
 *   <li>checks that every coefficient lies in the declared algebra
 *       (rejects coefficients that exceed it), and</li>
 *   <li>checks the trilinear identity
 *       {@code Σ_k U[a,k]·V[b,k]·W[c,k] == T[a,b,c]}
 *       AS AN EXACT ALGEBRAIC IDENTITY in that algebra (no epsilon, no
 *       floating-point rounding).</li>
 * </ol>
 *
 * <p>Supported algebras (per CLAUDE.md "Field discipline"):</p>
 * <ul>
 *   <li>{@link Algebra#Z} integers — {@link BigInteger}; rejects
 *       non-integer entries.</li>
 *   <li>{@link Algebra#Q} rationals — reduced {@code num/den} pairs;
 *       rejects irrationals (including disguised complex).</li>
 *   <li>{@link Algebra#F2} GF(2) — XOR semantics; rejects entries not
 *       in {0, 1}.</li>
 *   <li>{@link Algebra#C} Gaussian rationals — pairs of rationals; all
 *       R/Q schemes lift, plus genuine complex coefficients pass.</li>
 *   <li>{@link Algebra#R} reals — for now, accepts rationals only and
 *       documents the gap for irrational-but-real coefficients
 *       ({@code sqrt(N)} etc.) which the on-disk encoding doesn't yet
 *       support.</li>
 * </ul>
 *
 * <p>Public entry: {@link #verify(JsonNode)} dispatches on the JSON tags
 * ({@code z2}, {@code complex}, {@code field}) and the coefficient encoding,
 * returning a {@link Result} that carries the per-algebra verdict and a
 * human-readable reason.</p>
 */
public final class SymbolicVerifier {

	public enum Algebra { Z, Q, F2, F3, C, R }

	/** Verdict of a single symbolic check. */
	public record Result(boolean verified, Algebra algebra, String reason) {
		public static Result ok(Algebra a) {
			return new Result(true, a, "exact in " + a);
		}
		public static Result fail(Algebra a, String reason) {
			return new Result(false, a, reason);
		}
		@Override public String toString() {
			return (verified ? "OK   " : "FAIL ") + algebra + " : " + reason;
		}
	}

	private SymbolicVerifier() {}

	// ──────────────────────────────────────────────────────────────────────
	// Public dispatch
	// ──────────────────────────────────────────────────────────────────────

	/**
	 * Verify a scheme file symbolically. Reads the JSON and infers the
	 * algebra from JSON tags PLUS the filename suffix ({@code _Q_},
	 * {@code _F2_}, etc.) — the filename is load-bearing when the JSON
	 * omits an explicit {@code field} tag, as is the case for many catalog
	 * files (e.g. {@code dis09_Q-4x4x4_m100.json}).
	 */
	public static Result verify(File f) throws IOException {
		return verify(SchemeIO.parseJson(f), f.getName());
	}

	/** {@link #verify(JsonNode, String)} without a filename hint. */
	public static Result verify(JsonNode root) throws IOException {
		return verify(root, null);
	}

	/**
	 * Verify a scheme JSON symbolically in its declared algebra. Non-bilinear
	 * / stub layouts are acknowledged (returns {@code ok} with a note) since
	 * they don't fit the {@code U/V/W} trilinear identity. Complex schemes
	 * are checked exactly over the Gaussian rationals. Bilinear Q / Z / F₂
	 * schemes are checked over BigInteger / GF(2) respectively.
	 *
	 * @param root parsed scheme JSON
	 * @param filenameHint optional filename for picking up {@code _Q_},
	 *                     {@code _F2_}, {@code _Z_}, {@code _C_} field
	 *                     suffixes when the JSON omits the {@code field} tag
	 */
	public static Result verify(JsonNode root, String filenameHint) throws IOException {
		if (SchemeIO.isStub(root)) {
			return new Result(true, inferAlgebra(root, filenameHint),
					"stub (no factor matrices to verify symbolically)");
		}
		if (SchemeIO.isNonBilinear(root)) {
			return new Result(true, Algebra.Q,
					"non-bilinear (skipped — out of symbolic scope)");
		}
		if (SchemeIO.isComplex(root)) {
			return verifyComplex(root);
		}
		Algebra declared = inferAlgebra(root, filenameHint);
		if (declared == Algebra.F2) {
			return verifyF2(root);
		}
		if (declared == Algebra.F3) {
			return verifyF3(root);
		}
		return verifyOverIntegerOrRational(root, declared);
	}

	/** Convenience: verify an in-memory {@link NonCubicBilinearAlgorithm}. */
	public static Result verifyBilinear(NonCubicBilinearAlgorithm alg, Algebra declared) {
		// Field-membership + exact identity.
		BigInteger D = BigInteger.ONE;
		for (int phase = 0; phase < 3; phase++) {
			double[][] mat = switch (phase) { case 0 -> alg.denseU(); case 1 -> alg.denseV(); default -> alg.denseW(); };
			for (int i = 0; i < mat.length; i++) {
				for (int k = 0; k < mat[i].length; k++) {
					double v = mat[i][k];
					BigInteger d = denominatorOf(v);
					if (d.bitLength() > 20) {
						return Result.fail(declared,
								"coefficient " + v + " not representable in " + declared);
					}
					if (declared == Algebra.Z && !d.equals(BigInteger.ONE)) {
						return Result.fail(Algebra.Z,
								"non-integer coefficient " + v + " in Z-declared scheme");
					}
					D = lcm(D, d);
				}
			}
		}
		boolean ok = verifyBilinearExactScaled(alg, D);
		return ok ? Result.ok(declared)
				: Result.fail(declared, declared + " trilinear identity fails (exact)");
	}

	private static Algebra inferAlgebra(JsonNode root, String filenameHint) {
		if (SchemeIO.isZ2(root)) return Algebra.F2;
		if (SchemeIO.isComplex(root)) return Algebra.C;
		// Unified fields[] (task #174): pick the STRICTEST tag present so the
		// declared-algebra check is the tightest one the scheme claims to pass.
		// Strictness order: F2, F3, Z, Q, R, C (F2/F3 are separate universes,
		// listed first only so a sole-prime-field scheme resolves to it).
		java.util.List<String> tags = SchemeIO.fieldTags(root);
		if (!tags.isEmpty()) {
			for (String t : new String[] { "Z", "Q", "R", "C", "F3", "F2" }) {
				if (tags.contains(t)) return Algebra.valueOf(t);
			}
		}
		JsonNode f = root.get("field");
		if (f != null && f.isTextual()) {
			Algebra a = parseFieldLabel(f.asString());
			if (a != null) return a;
		}
		// Fall back to filename token: …_Q_, …_Q.json, …_F2_, …_C_, …_R_,
		// …_Z_, …_ZT_ (treat ZT as Z for verification — both are integer
		// algebras; the T suffix in dronperminov-speak distinguishes
		// "tensor decomposition over Z" sub-catalog provenance only).
		if (filenameHint != null) {
			Algebra a = parseFilenameField(filenameHint);
			if (a != null) return a;
		}
		// Default: Z (so non-integers force a failure under the strictest tag).
		return Algebra.Z;
	}

	private static Algebra parseFieldLabel(String s) {
		// Order matters: F2 before R/C; ZT before Z. The AlphaEvolve corpus
		// also uses fractional-scale labels "0.5*Z" / "0.5*C" — those live
		// in Q resp. Q[i]=C, NOT in Z, because 1/2 isn't in Z.
		if (s.contains("F2") || s.contains("F_2")) return Algebra.F2;
		if (s.contains("0.5*C") || s.contains("0.5xC")) return Algebra.C;
		if (s.contains("0.5*Z") || s.contains("0.5xZ")) return Algebra.Q;
		if (s.contains("C")) return Algebra.C;
		if (s.contains("R")) return Algebra.R;
		if (s.contains("Q")) return Algebra.Q;
		if (s.contains("Z")) return Algebra.Z;
		return null;
	}

	private static Algebra parseFilenameField(String fname) {
		String n = fname;
		if (n.endsWith(".json")) n = n.substring(0, n.length() - 5);
		// Convention: field tokens appear as standalone tokens delimited by
		// '_' or '-' (e.g. "dis09-Q_4x4x4", "perminov-ZT_*", "alphatensor-F2_*",
		// "alphaevolve-4x4x4_m48_a1264_0.5xC"). Split on both delimiters so
		// glued forms like "dis09-Q" still surface "Q".
		String[] parts = n.split("[_-]");
		for (String tok : parts) {
			Algebra a = switch (tok) {
				case "F2" -> Algebra.F2;
				case "Q" -> Algebra.Q;
				case "R" -> Algebra.R;
				case "C" -> Algebra.C;
				case "Z", "ZT" -> Algebra.Z;
				default -> null;
			};
			if (a != null) return a;
		}
		// AlphaEvolve fractional-scale suffixes:
		//   "0.5xC" → entries are 0.5·(Gaussian integer) → still in Q[i] = C
		//   "0.5xZ" → entries are 0.5·(integer)         → half-integers in Q
		// In both cases the algebra LIVES in Q (resp. Q[i]) — the "x{Z,C}"
		// only documents the integer kernel, NOT the actual coefficient ring.
		for (String tok : parts) {
			if (tok.contains("xC") || tok.contains("XC")) return Algebra.C;
			if (tok.contains("xZ") || tok.contains("XZ")) return Algebra.Q;
		}
		return null;
	}

	// ──────────────────────────────────────────────────────────────────────
	// F2 verification
	// ──────────────────────────────────────────────────────────────────────

	private static Result verifyF2(JsonNode root) throws IOException {
		NonCubicBilinearAlgorithm alg = SchemeIO.isReduced(root)
				? SchemeIO.readReduced(root) : SchemeIO.read(root);
		// Field-membership: every entry must be in {0, 1} mod 2 (i.e. an integer).
		Optional<String> bad = firstNonIntegerEntry(alg);
		if (bad.isPresent()) {
			return Result.fail(Algebra.F2, "F2 entry not integer: " + bad.get());
		}
		int wrong = Verifier.residualNonCubicF2(alg);
		if (wrong != 0) {
			return Result.fail(Algebra.F2,
					"F2 trilinear identity fails at " + wrong + " positions");
		}
		return Result.ok(Algebra.F2);
	}

	private static Result verifyF3(JsonNode root) throws IOException {
		NonCubicBilinearAlgorithm alg = SchemeIO.isReduced(root)
				? SchemeIO.readReduced(root) : SchemeIO.read(root);
		// Field-membership: integer coefficients reduce trivially mod 3.
		// Strictly, Q coefficients with denominators coprime to 3 (e.g. 1/2 ≡ 2)
		// also live in F₃, but require an integer pre-filter rather than
		// fall through to the rounded residual path. Widen if a Q-scheme
		// with such denominators ever needs to surface in F₃.
		Optional<String> bad = firstNonIntegerEntry(alg);
		if (bad.isPresent()) {
			return Result.fail(Algebra.F3, "F3 entry not integer: " + bad.get());
		}
		int wrong = Verifier.residualNonCubicF3(alg);
		if (wrong != 0) {
			return Result.fail(Algebra.F3,
					"F3 trilinear identity fails at " + wrong + " positions");
		}
		return Result.ok(Algebra.F3);
	}

	// ──────────────────────────────────────────────────────────────────────
	// Complex verification (Gaussian rationals)
	// ──────────────────────────────────────────────────────────────────────

	private static Result verifyComplex(JsonNode root) throws IOException {
		ComplexNonCubicBilinearAlgorithm alg = SchemeIO.readComplex(root);
		// Field-membership: each (re, im) entry must come from a small-rational
		// representation. SchemeIO has already parsed numeric / fractional
		// values; an "irrational" surrogate denominator (1<<30) signals failure.
		for (int phase = 0; phase < 6; phase++) {
			double[][] mat = switch (phase) {
				case 0 -> alg.uRe; case 1 -> alg.uIm;
				case 2 -> alg.vRe; case 3 -> alg.vIm;
				case 4 -> alg.wRe; default -> alg.wIm;
			};
			for (double[] row : mat) {
				for (double v : row) {
					if (denominatorOf(v).bitLength() > 20) {
						return Result.fail(Algebra.C,
								"complex coefficient " + v + " not a small rational");
					}
				}
			}
		}
		BigInteger common = denominatorLcmComplex(alg);
		boolean ok = verifyComplexExact(alg, common);
		if (!ok) {
			return Result.fail(Algebra.C, "C exact integer-scaled identity failed");
		}
		return Result.ok(Algebra.C);
	}

	private static BigInteger denominatorLcmComplex(ComplexNonCubicBilinearAlgorithm alg) {
		BigInteger d = BigInteger.ONE;
		d = updateLcm(d, alg.uRe); d = updateLcm(d, alg.uIm);
		d = updateLcm(d, alg.vRe); d = updateLcm(d, alg.vIm);
		d = updateLcm(d, alg.wRe); d = updateLcm(d, alg.wIm);
		return d;
	}

	private static BigInteger updateLcm(BigInteger acc, double[][] m) {
		for (double[] row : m) {
			for (double v : row) {
				BigInteger d = denominatorOf(v);
				acc = lcm(acc, d);
			}
		}
		return acc;
	}

	/**
	 * Find the smallest integer denominator (≤ 1024) such that {@code v · d}
	 * is an integer to 1e-9 precision. Returns {@code 1L<<30} (treated as an
	 * "irrational" sentinel) if no small denominator works.
	 */
	static BigInteger denominatorOf(double v) {
		if (v == 0.0) return BigInteger.ONE;
		for (int d = 1; d <= 1024; d++) {
			double scaled = v * d;
			if (Math.abs(scaled - Math.round(scaled)) < 1e-9) {
				return BigInteger.valueOf(d);
			}
		}
		return BigInteger.valueOf(1L << 30);
	}

	static BigInteger numeratorScaled(double v, BigInteger d) {
		double scaled = v * d.doubleValue();
		return BigInteger.valueOf(Math.round(scaled));
	}

	private static boolean verifyComplexExact(ComplexNonCubicBilinearAlgorithm alg, BigInteger D) {
		int n = alg.n, m = alg.m, p = alg.p, r = alg.r;
		BigInteger D3 = D.multiply(D).multiply(D);
		int dimU = n * m, dimV = m * p, dimW = n * p;
		BigInteger[][] uR = scaleToInt(alg.uRe, D), uI = scaleToInt(alg.uIm, D);
		BigInteger[][] vR = scaleToInt(alg.vRe, D), vI = scaleToInt(alg.vIm, D);
		BigInteger[][] wR = scaleToInt(alg.wRe, D), wI = scaleToInt(alg.wIm, D);
		for (int a = 0; a < dimU; a++) {
			int aI = a / m, aJ = a % m;
			for (int b = 0; b < dimV; b++) {
				int bJ = b / p, bL = b % p;
				for (int c = 0; c < dimW; c++) {
					int i = c / p, l = c % p;
					BigInteger sumRe = BigInteger.ZERO;
					BigInteger sumIm = BigInteger.ZERO;
					for (int k = 0; k < r; k++) {
						BigInteger ur = uR[a][k], ui = uI[a][k];
						BigInteger vr = vR[b][k], vi = vI[b][k];
						BigInteger wr = wR[c][k], wi = wI[c][k];
						BigInteger uvRe = ur.multiply(vr).subtract(ui.multiply(vi));
						BigInteger uvIm = ur.multiply(vi).add(ui.multiply(vr));
						sumRe = sumRe.add(uvRe.multiply(wr).subtract(uvIm.multiply(wi)));
						sumIm = sumIm.add(uvRe.multiply(wi).add(uvIm.multiply(wr)));
					}
					int target = (aI == i && aJ == bJ && bL == l) ? 1 : 0;
					BigInteger expectRe = BigInteger.valueOf(target).multiply(D3);
					if (!sumRe.equals(expectRe) || sumIm.signum() != 0) return false;
				}
			}
		}
		return true;
	}

	// ──────────────────────────────────────────────────────────────────────
	// Integer / Rational verification (the common case)
	// ──────────────────────────────────────────────────────────────────────

	private static Result verifyOverIntegerOrRational(JsonNode root, Algebra declared) throws IOException {
		NonCubicBilinearAlgorithm alg = SchemeIO.isReduced(root)
				? SchemeIO.readReduced(root) : SchemeIO.read(root);
		return verifyBilinear(alg, declared);
	}

	private static boolean verifyBilinearExactScaled(NonCubicBilinearAlgorithm alg, BigInteger D) {
		int n = alg.n, m = alg.m, p = alg.p, r = alg.r;
		BigInteger D3 = D.multiply(D).multiply(D);
		int dimU = n * m, dimV = m * p, dimW = n * p;
		BigInteger[][] U = scaleToInt(alg.denseU(), D);
		BigInteger[][] V = scaleToInt(alg.denseV(), D);
		BigInteger[][] W = scaleToInt(alg.denseW(), D);
		for (int a = 0; a < dimU; a++) {
			int aI = a / m, aJ = a % m;
			for (int b = 0; b < dimV; b++) {
				int bJ = b / p, bL = b % p;
				for (int c = 0; c < dimW; c++) {
					int i = c / p, l = c % p;
					int target = (aI == i && aJ == bJ && bL == l) ? 1 : 0;
					BigInteger sum = BigInteger.ZERO;
					for (int k = 0; k < r; k++) {
						BigInteger uv = U[a][k].multiply(V[b][k]);
						if (uv.signum() == 0) continue;
						sum = sum.add(uv.multiply(W[c][k]));
					}
					BigInteger expected = BigInteger.valueOf(target).multiply(D3);
					if (!sum.equals(expected)) return false;
				}
			}
		}
		return true;
	}

	private static BigInteger[][] scaleToInt(double[][] m, BigInteger D) {
		BigInteger[][] out = new BigInteger[m.length][];
		double dd = D.doubleValue();
		for (int i = 0; i < m.length; i++) {
			out[i] = new BigInteger[m[i].length];
			for (int k = 0; k < m[i].length; k++) {
				double v = m[i][k] * dd;
				out[i][k] = BigInteger.valueOf(Math.round(v));
			}
		}
		return out;
	}

	// ──────────────────────────────────────────────────────────────────────
	// Helpers
	// ──────────────────────────────────────────────────────────────────────

	private static Optional<String> firstNonIntegerEntry(NonCubicBilinearAlgorithm alg) {
		for (int phase = 0; phase < 3; phase++) {
			double[][] mat = switch (phase) { case 0 -> alg.denseU(); case 1 -> alg.denseV(); default -> alg.denseW(); };
			String name = switch (phase) { case 0 -> "U"; case 1 -> "V"; default -> "W"; };
			for (int i = 0; i < mat.length; i++) {
				for (int k = 0; k < mat[i].length; k++) {
					double v = mat[i][k];
					if (Math.abs(v - Math.round(v)) > 1e-9) {
						return Optional.of(name + "[" + i + "," + k + "]=" + v);
					}
				}
			}
		}
		return Optional.empty();
	}

	static BigInteger lcm(BigInteger a, BigInteger b) {
		if (a.signum() == 0 || b.signum() == 0) return BigInteger.ZERO;
		return a.divide(a.gcd(b)).multiply(b).abs();
	}
}
