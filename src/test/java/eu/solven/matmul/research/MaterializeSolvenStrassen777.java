package eu.solven.matmul.research;

import java.io.File;
import java.nio.file.Path;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.LineageTrackingLookup;
import eu.solven.matmul.recombination.Recombination;
import eu.solven.matmul.catalog.SchemeIO;

/** Cascade-materialise small cubics that benefit from FMM ⟨4,4,4⟩=48 leaf. */
public final class MaterializeSolvenStrassen777 {
	public static void main(String[] args) throws Exception {
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		FieldAwareLookup lookup = new FieldAwareLookup("Q");

		record Target(int n, int[] alloc) {}
		Target[] targets = {
				new Target(7, new int[]{3,4}),  // Strassen[3,4]³: 48+29·3+38·3 = 249
				new Target(8, new int[]{4,4}),  // Strassen[4,4]³: 7·48 = 336
				new Target(9, new int[]{4,5}),  // mix; need ⟨5,5,5⟩ + ⟨4,*,*⟩
				new Target(10, new int[]{4,6}), // mix; depends on ⟨6,*,*⟩
				new Target(11, new int[]{4,7}), // uses fresh ⟨7,7,7⟩=249
				new Target(12, new int[]{4,8}), // uses fresh ⟨8,8,8⟩
		};

		for (Target t : targets) {
			LineageTrackingLookup tracker = new LineageTrackingLookup(lookup);
			NonCubicBilinearAlgorithm scheme = Recombination.constructWithAllocation(
					strassen, tracker, t.alloc, t.alloc, t.alloc);
			boolean ok = Verifier.passesRandomMatmulSpotCheck(scheme);
			int adds = Verifier.additionCount(scheme);
			System.out.printf("⟨%d,%d,%d⟩ via Strassen%s³: rank=%d add=%d  %s%n",
					t.n, t.n, t.n, java.util.Arrays.toString(t.alloc),
					scheme.r, adds, ok ? "PASS" : "FAIL");
			if (!ok) { System.err.println("  ABORTING write"); continue; }
			Path out = Path.of(String.format(
					"src/main/resources/schemes/known/section%d/solven_strassen_2026-%dx%dx%d_m%d_a%d.json",
					t.n, t.n, t.n, t.n, scheme.r, adds));
			Lineage.Node lineage = new Lineage.RecombinationN(
					new Lineage.Atom("Strassen<2,2,2>=7"),
					t.alloc.clone(), t.alloc.clone(), t.alloc.clone(),
					tracker.leaves());
			SchemeIO.write(scheme, out.toFile(), lineage);
			System.out.println("  wrote " + out);
		}
	}
}
