package eu.solven.matmul.papers.hopcroftkerr1971;

import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;

/**
 * Diagnostic survey (task #7, HK71 constructive closure): for each ⟨2,n,p⟩ in
 * range, drive the existing emitters and report rank-vs-formula, exactness, and
 * which same-method cases throw. Establishes the TRUE current frontier of the
 * constructive HK implementation (the README's narrative predates the
 * sympy-derived (2,2,b3)/(1,1,b3) emitters and the FMM sync).
 */
class SurveyHk2npEmitterGaps {

	private static int hk(int n, int p) {
		return (int) Math.ceil((3.0 * n * p + Math.max(n, p)) / 2.0);
	}

	@Test
	@Tag("slow")
	void survey_square() {
		System.out.println("\n=== square ⟨2,n,n⟩: emitter vs formula ===");
		for (int n = 3; n <= 16; n++) {
			String res;
			try {
				NonCubicBilinearAlgorithm alg = HopcroftKerr2bc.buildSquare(n);
				boolean exact = Verifier.isExactNonCubic(alg.orientAs(2, n, n).orElse(alg))
						|| Verifier.isExactNonCubic(alg);
				res = String.format("r=%d (formula %d, delta %+d) exact=%b shape=%dx%dx%d",
						alg.r, hk(n, n), alg.r - hk(n, n), exact, alg.n, alg.m, alg.p);
			} catch (Throwable e) {
				res = "THROWS: " + e.getClass().getSimpleName() + ": "
						+ String.valueOf(e.getMessage()).replaceAll("\n.*", "");
			}
			System.out.printf("  n=%-3d %s%n", n, res);
		}
	}

	@Test
	@Tag("slow")
	void survey_rectangular() {
		System.out.println("\n=== rectangular ⟨2,n,p⟩ via Asymmetric.build: emitter vs formula ===");
		for (int n = 3; n <= 12; n++) {
			for (int p = n + 1; p <= Math.min(n + 6, 16); p++) {
				String res;
				try {
					NonCubicBilinearAlgorithm alg = HopcroftKerr2bcAsymmetric.build(n, p); // (smaller, larger)
					boolean exact = Verifier.isExactNonCubic(alg);
					res = String.format("r=%d (formula %d, delta %+d) exact=%b shape=%dx%dx%d",
							alg.r, hk(n, p), alg.r - hk(n, p), exact, alg.n, alg.m, alg.p);
				} catch (Throwable e) {
					res = "THROWS: " + e.getClass().getSimpleName() + ": "
							+ String.valueOf(e.getMessage()).replaceAll("\n.*", "");
				}
				System.out.printf("  n=%-3d p=%-3d %s%n", n, p, res);
			}
		}
	}
}
