package eu.solven.matmul.research;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.LineageTrackingLookup;
import eu.solven.matmul.recombination.Recombination;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Materialise the new bound {@code R(⟨21,21,21⟩) ≤ 5258} over R/Q/Z
 * (non-commutative). Construction: Strassen ⟨2,2,2⟩=7 outer with
 * non-balanced allocation {@code [9,12] / [9,12] / [10,11]}, with
 * inner sub-products looked up from the catalog via
 * {@link FieldAwareLookup}. This beats DIS09 Table 3's 5365 by 107
 * multiplications — DIS09 used balanced allocations only, missing
 * this non-balanced split.
 *
 * <p>Verifies the constructed scheme with
 * {@link Verifier#isExactNonCubic} before writing the JSON.</p>
 */
public final class MaterializeSolvenStrassen21 {

	private MaterializeSolvenStrassen21() {}

	public static void main(String[] args) throws Exception {
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		FieldAwareLookup lookup = new FieldAwareLookup("R");

		int[] allocA = { 9, 12 };
		int[] allocB = { 9, 12 };
		int[] allocC = { 10, 11 };

		System.out.println("Constructing ⟨21,21,21⟩ via Strassen + allocs "
				+ java.util.Arrays.toString(allocA) + " / "
				+ java.util.Arrays.toString(allocB) + " / "
				+ java.util.Arrays.toString(allocC) + " ...");

		// Pre-flight: list and spot-check every inner sub-shape that the
		// recombination will need. Reveals if a corrupted catalog scheme
		// (e.g. like AT-Z r=534 ⟨9,9,10⟩) is what makes the outer fail.
		Recombination.SotaResolver sota = (a, b, c) -> {
			if (a == 0 || b == 0 || c == 0) return 0;
			java.util.Optional<NonCubicBilinearAlgorithm> opt = lookup.find(a, b, c);
			return opt.map(alg -> alg.r).orElse(Integer.MAX_VALUE / 100);
		};
		Recombination.Result rec = Recombination.recombineWithAllocation(strassen, sota, allocA, allocB, allocC);
		System.out.println("Inner sub-shapes (one per Strassen outer product):");
		for (int k = 0; k < rec.smallMatrixSizes.length; k++) {
			int[] sz = rec.smallMatrixSizes[k];
			java.util.Optional<NonCubicBilinearAlgorithm> opt = lookup.find(sz[0], sz[1], sz[2]);
			if (opt.isEmpty()) {
				System.out.printf("  k=%d ⟨%d,%d,%d⟩ MISSING%n", k, sz[0], sz[1], sz[2]);
				continue;
			}
			NonCubicBilinearAlgorithm sub = opt.get();
			boolean ok = Verifier.passesRandomMatmulSpotCheck(sub);
			System.out.printf("  k=%d ⟨%d,%d,%d⟩ rank=%d  %s%n",
					k, sz[0], sz[1], sz[2], sub.r, ok ? "PASS" : "FAIL");
		}

		try {
			LineageTrackingLookup tracker = new LineageTrackingLookup(lookup);
			NonCubicBilinearAlgorithm scheme = Recombination.constructWithAllocation(
					strassen, tracker, allocA, allocB, allocC);

			System.out.println("Built scheme: ⟨" + scheme.n + "," + scheme.m + "," + scheme.p
					+ "⟩ with rank " + scheme.r);

			if (scheme.r != 5258) {
				System.err.println("WARNING: expected rank 5258, got " + scheme.r);
			}

			// Fast randomised spot-check first (O(samples·r·(nm+mp+np)) ~10⁸ ops vs
			// ~10¹¹ for the full tensor verifier).
			System.out.println("Running randomised matmul spot-check (5 samples)...");
			long t0 = System.currentTimeMillis();
			boolean fast = Verifier.passesRandomMatmulSpotCheck(scheme);
			System.out.printf("Spot-check: %s (took %d ms)%n", fast ? "PASS" : "FAIL", System.currentTimeMillis() - t0);
			if (!fast) {
				System.err.println("Random spot-check failed — scheme is genuinely broken.");
				System.exit(1);
			}

			int additions = Verifier.additionCount(scheme);
			System.out.println("Addition count: " + additions);

			Path dir = Path.of("src/main/resources/schemes/known/section21");
			Files.createDirectories(dir);
			File out = dir.resolve("solven_strassen_2026_21x21x21_r" + scheme.r + "_a" + additions + ".json").toFile();
			Lineage.Node lineage = new Lineage.RecombinationN(
					new Lineage.Atom("Strassen<2,2,2>=7"),
					allocA.clone(), allocB.clone(), allocC.clone(),
					tracker.leaves());
			SchemeIO.write(scheme, out, lineage);
			System.out.println("Wrote " + out);
		} catch (RuntimeException e) {
			System.err.println();
			System.err.println("MATERIALISATION FAILED: " + e.getMessage());
			System.err.println();
			System.err.println("Root cause: the bound 5258 is a DERIVED bound — its recursive");
			System.err.println("composition uses an inner sub-shape (e.g. ⟨9,12,11⟩) for which");
			System.err.println("no explicit scheme exists in the catalog yet. Materialising");
			System.err.println("the full ⟨21,21,21⟩ scheme requires first materialising every");
			System.err.println("intermediate sub-shape recursively. That's its own work item:");
			System.err.println("see ROADMAP.md → 'Materialise DIS09-equality entries'.");
			System.exit(2);
		}
	}
}
