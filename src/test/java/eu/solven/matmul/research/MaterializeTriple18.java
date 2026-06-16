package eu.solven.matmul.research;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Optional;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.Compose;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Materialise the {@code ⟨18,18,18⟩=3200} prediction surfaced by the
 * Kronecker-augmented strategy survey:
 *
 * <pre>
 *   ⟨18,18,18⟩  =  concat-p[9, 9] ( ⟨18,18,9⟩, ⟨18,18,9⟩ )
 *   ⟨18,18,9⟩   =  Kronecker ( ⟨3,6,3⟩=40 ⊗ ⟨6,3,3⟩=40 )
 * </pre>
 *
 * <p>Both leaf schemes are in the on-disk catalog. The result is a real
 * (non-phantom) constructive scheme at rank 3200, matching Sedoglavic's
 * FMM-digest value and beating our previous catalog rank 3402.</p>
 */
public final class MaterializeTriple18 {

	public static void main(String[] args) throws Exception {
		FieldAwareLookup lookup = new FieldAwareLookup("Q");

		Optional<NonCubicBilinearAlgorithm> outer = lookup.find(3, 6, 3);
		Optional<NonCubicBilinearAlgorithm> inner = lookup.find(6, 3, 3);
		if (outer.isEmpty() || inner.isEmpty()) {
			System.err.println("missing leaf scheme: ⟨3,6,3⟩=" + outer.isPresent()
					+ "  ⟨6,3,3⟩=" + inner.isPresent());
			return;
		}
		System.out.printf("Leaves: ⟨3,6,3⟩=%d  ⟨6,3,3⟩=%d%n", outer.get().r, inner.get().r);

		NonCubicBilinearAlgorithm leftHalf = Compose.kroneckerGeneral(outer.get(), inner.get());
		System.out.printf("Kronecker product: ⟨%d,%d,%d⟩  rank=%d  (expected ⟨18,18,9⟩=1600)%n",
				leftHalf.n, leftHalf.m, leftHalf.p, leftHalf.r);
		if (leftHalf.n != 18 || leftHalf.m != 18 || leftHalf.p != 9 || leftHalf.r != 1600) {
			System.err.println("shape/rank mismatch — aborting");
			return;
		}
		if (!Verifier.passesRandomMatmulSpotCheck(leftHalf)) {
			System.err.println("⟨18,18,9⟩ spot-check FAILED");
			return;
		}
		System.out.println("  ✓ ⟨18,18,9⟩=1600 spot-checks");

		NonCubicBilinearAlgorithm cube = Compose.concatRight(leftHalf, leftHalf);
		System.out.printf("Concat-right: ⟨%d,%d,%d⟩  rank=%d  (expected ⟨18,18,18⟩=3200)%n",
				cube.n, cube.m, cube.p, cube.r);
		if (cube.n != 18 || cube.m != 18 || cube.p != 18 || cube.r != 3200) {
			System.err.println("shape/rank mismatch — aborting");
			return;
		}
		if (!Verifier.passesRandomMatmulSpotCheck(cube)) {
			System.err.println("⟨18,18,18⟩ spot-check FAILED");
			return;
		}
		int adds = Verifier.additionCount(cube);
		System.out.printf("  ✓ ⟨18,18,18⟩=3200 spot-checks  adds=%d%n", adds);

		// Build lineage tree. ⟨3,6,3⟩ and ⟨6,3,3⟩ both resolve via
		// FieldAwareLookup to the canonical ⟨3,3,6⟩=40 leaf, oriented via
		// tensor symmetry (the file on disk is sorted-axes canonical).
		// Shared sub-tree for the ⟨18,18,9⟩ Kronecker — JSON dedups; pretty
		// form repeats inline.
		Lineage.Node leaf336 = new Lineage.Atom("fmm_lille-3x3x6_m40_a862");
		Lineage.Node leafOuter = new Lineage.Transpose(leaf336, "NMP->NPM"); // ⟨3,3,6⟩ → ⟨3,6,3⟩
		Lineage.Node leafInner = new Lineage.Transpose(leaf336, "NMP->PMN"); // ⟨3,3,6⟩ → ⟨6,3,3⟩
		Lineage.Node kronHalf = new Lineage.KronProduct(leafOuter, leafInner);
		Lineage.Node lineageCube = new Lineage.ConcatCols(kronHalf, kronHalf);  // same object reused

		Path out = Path.of(String.format(
				"src/main/resources/schemes/known/section18/derived_kron_concat-%dx%dx%d_m%d_a%d.json",
				18, 18, 18, cube.r, adds));
		Files.createDirectories(out.getParent());
		SchemeIO.write(cube, out.toFile(), lineageCube);
		System.out.println("wrote " + out);
		System.out.println("  lineage: " + Lineage.prettyString(lineageCube));

		// Also write the intermediate ⟨18,18,9⟩=1600 — it doesn't exist on disk
		// and is reusable for other compositions (e.g. concat-n with ⟨18,18,9⟩
		// for ⟨36,18,9⟩, or further Kronecker products).
		Path mid = Path.of(String.format(
				"src/main/resources/schemes/known/section18/derived_kron-%dx%dx%d_m%d_a%d.json",
				18, 18, 9, leftHalf.r, Verifier.additionCount(leftHalf)));
		Files.createDirectories(mid.getParent());
		SchemeIO.write(leftHalf, mid.toFile(), kronHalf);
		System.out.println("wrote " + mid);
		System.out.println("  lineage: " + Lineage.prettyString(kronHalf));
	}
}
