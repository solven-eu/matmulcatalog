package eu.solven.matmul.docs.verify;

import java.io.File;
import java.util.Optional;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.Compose;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Settle whether ⟨3,3,14⟩=95 and ⟨3,3,15⟩=103 (flagged "better than catalog"
 * by the band 11-16 sweep) are REAL, by explicitly building the concat of
 * on-disk NC schemes and running the Verifier:
 *   ⟨3,3,14⟩ = ⟨3,3,2⟩(=⟨2,3,3⟩=15) ⊞ ⟨3,3,12⟩(=80)  → 95   (catalog 98)
 *   ⟨3,3,15⟩ = ⟨3,3,3⟩(=23)         ⊞ ⟨3,3,12⟩(=80)  → 103  (catalog 105)
 * The materialiser missed these because it couldn't orient ⟨2,3,3⟩→⟨3,3,2⟩
 * and fell back to naive ⟨3,3,2⟩=18 (giving 18+80=98).
 */
public final class VerifyConcat3314 {

	private VerifyConcat3314() {}

	private static NonCubicBilinearAlgorithm orient(NonCubicBilinearAlgorithm a, int n, int m, int p) {
		if (a.n == n && a.m == m && a.p == p) return a;
		Optional<NonCubicBilinearAlgorithm> o = a.orientAs(n, m, p);
		if (o.isEmpty()) throw new IllegalStateException(
				"cannot orient ⟨" + a.n + "," + a.m + "," + a.p + "⟩ → ⟨" + n + "," + m + "," + p + "⟩");
		return o.get();
	}

	private static void tryConcat(String label, NonCubicBilinearAlgorithm left,
			NonCubicBilinearAlgorithm right, int catalog) {
		NonCubicBilinearAlgorithm c = Compose.concatRight(left, right);
		boolean ok = Verifier.isExactNonCubic(c);
		System.out.printf("%s : built ⟨%d,%d,%d⟩ rank=%d  verified=%s  catalog=%d  → %s%n",
				label, c.n, c.m, c.p, c.r, ok, catalog,
				ok && c.r < catalog ? "✔ REAL WIN (−" + (catalog - c.r) + ")"
						: ok ? "verified but not < catalog" : "✗ VERIFIER FAILED");
	}

	public static void main(String[] args) throws Exception {
		// ⟨2,3,3⟩=15 NC (AlphaTensor over Z), oriented to ⟨3,3,2⟩.
		NonCubicBilinearAlgorithm s233 = SchemeIO.read(eu.solven.matmul.catalog.SchemeResolver.byHint(
				"src/main/resources/schemes/known/section3/alphatensor_Z-2x3x3_m15_a58.json"));
		NonCubicBilinearAlgorithm s332 = orient(s233, 3, 3, 2);
		// ⟨3,3,12⟩=80 NC (FMM-Lille).
		NonCubicBilinearAlgorithm s3312 = SchemeIO.read(eu.solven.matmul.catalog.SchemeResolver.byHint(
				"src/main/resources/schemes/known/section12/fmm_lille-3x3x12_m80_a1724.json"));
		// ⟨3,3,3⟩=23 NC (Laderman).
		NonCubicBilinearAlgorithm s333 = SchemeIO.read(eu.solven.matmul.catalog.SchemeResolver.byHint(
				"src/main/resources/schemes/known/section3/laderman_1976-3x3x3_m23_a98.json"));

		System.out.printf("atoms: ⟨3,3,2⟩=%d  ⟨3,3,12⟩=%d  ⟨3,3,3⟩=%d%n", s332.r, s3312.r, s333.r);
		tryConcat("⟨3,3,14⟩ = ⟨3,3,2⟩ ⊞ ⟨3,3,12⟩", s332, s3312, 98);
		tryConcat("⟨3,3,15⟩ = ⟨3,3,3⟩ ⊞ ⟨3,3,12⟩", s333, s3312, 105);
	}
}
