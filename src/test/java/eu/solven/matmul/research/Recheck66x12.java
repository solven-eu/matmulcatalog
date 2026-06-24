package eu.solven.matmul.research;

import java.io.File;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.recombination.Recombination;
import eu.solven.matmul.catalog.SchemeIO;

public final class Recheck66x12 {
	public static void main(String[] args) throws Exception {
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		NonCubicBilinearAlgorithm scheme = Recombination.constructWithAllocation(
				strassen, lookup, new int[]{3,3}, new int[]{3,3}, new int[]{6,6});
		System.out.printf("⟨6,6,12⟩ via Strassen[3,3]/[3,3]/[6,6] + Q-lookup: rank=%d%n", scheme.r);
		System.out.printf("  spot-check: %s%n", Verifier.passesRandomMatmulSpotCheck(scheme) ? "PASS" : "FAIL");
	}
}
