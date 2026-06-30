package eu.solven.matmul.research;

import java.io.File;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.recombination.Recombination;
import eu.solven.matmul.catalog.SchemeIO;

/** Probe: does our Strassen[3,4]³ recombination for ⟨7,7,7⟩ now use
 *  the freshly-imported ⟨4,4,4⟩=48 leaf, hitting rank 249 (FMM's value)? */
public final class Recheck777ViaStrassen34 {
	public static void main(String[] args) throws Exception {
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		NonCubicBilinearAlgorithm scheme = Recombination.constructWithAllocation(
				strassen, lookup, new int[]{3,4}, new int[]{3,4}, new int[]{3,4});
		System.out.printf("⟨7,7,7⟩ via Strassen[3,4]³ + Q-lookup: rank=%d%n", scheme.r);
		boolean ok = Verifier.passesRandomMatmulSpotCheck(scheme);
		System.out.printf("  spot-check: %s%n", ok ? "PASS" : "FAIL");
		System.out.printf("  additions: %d%n", Verifier.additionCount(scheme));
	}
}
