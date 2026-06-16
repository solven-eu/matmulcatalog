package eu.solven.matmul.docs;

import java.io.File;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.SchemeIO;

public final class DebugF2Verify {
	public static void main(String[] a) throws Exception {
		NonCubicBilinearAlgorithm s = SchemeIO.read(eu.solven.matmul.catalog.SchemeResolver.byHint(a[0]));
		double[][] srcW = s.denseW();
		System.out.printf("⟨%d,%d,%d⟩ rank=%d%n", s.n, s.m, s.p, s.r);
		System.out.printf("F2 verify: %s%n", Verifier.isExactNonCubicF2(s) ? "PASS" : "FAIL");
		System.out.printf("Real spot-check: %s%n", Verifier.passesRandomMatmulSpotCheck(s) ? "PASS" : "FAIL");
		// Look at W shape and entry types
		System.out.printf("W shape: [%d][%d]%n", srcW.length, srcW[0].length);
	}
}
