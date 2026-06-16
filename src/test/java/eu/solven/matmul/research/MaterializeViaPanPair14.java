package eu.solven.matmul.research;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;

import java.util.List;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.RecombinationWithPair;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * First end-to-end materialisation using {@link
 * eu.solven.matmul.catalog.PanPairProduct}: build ⟨14,14,14⟩ via
 * Strassen[7,7]³ with 3 of the 7 Strassen sub-products fused into
 * Pan pair-products, leaving 1 solo.
 *
 * <p>Predicted rank: {@code 3·490 + R(⟨7,7,7⟩)}. With our catalog's
 * solven-strassen-2026 ⟨7,7,7⟩=249 (built via Strassen[3,4]³ over
 * DPS 2025 ⟨4,4,4⟩=48), that's {@code 1470 + 249 = 1719}, vs the
 * current solven-strassen-2026 ⟨14,14,14⟩=1743 (no pairing) and the
 * unverified target ⟨14,14,14⟩=1720 (with Perminov-original ⟨7,7,7⟩=250).</p>
 */
public final class MaterializeViaPanPair14 {

	public static void main(String[] args) throws Exception {
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		FieldAwareLookup lookup = new FieldAwareLookup("Q");

		// Pair Strassen sub-products 0+1, 2+3, 4+5; leave 6 solo.
		RecombinationWithPair.Pairing pairing = new RecombinationWithPair.Pairing(
				new int[][] { {0, 1}, {2, 3}, {4, 5} },
				new int[]  { 6 });

		NonCubicBilinearAlgorithm scheme = RecombinationWithPair.constructWithPairing(
				strassen, lookup, new int[]{7, 7}, new int[]{7, 7}, new int[]{7, 7},
				pairing);

		System.out.printf("⟨14,14,14⟩ via Strassen[7,7]³ + Pan pair(3 pairs + 1 solo): rank=%d%n", scheme.r);
		boolean ok = Verifier.passesRandomMatmulSpotCheck(scheme);
		System.out.printf("  spot-check: %s%n", ok ? "PASS" : "FAIL");
		if (!ok) { System.err.println("  ABORT — scheme is broken"); System.exit(1); }
		int adds = Verifier.additionCount(scheme);
		System.out.printf("  additions: %d%n", adds);

		Path out = Path.of(String.format(
				"src/main/resources/schemes/known/section14/solven_strassen_2026-14x14x14_m%d_a%d.json",
				scheme.r, adds));
		Files.createDirectories(out.getParent());

		// Lineage: Pan-pair recombination of Strassen<2,2,2>=7 with allocs
		// [7,7]³ and a single ⟨7,7,7⟩ inner leaf. Each pair {0,1},{2,3},{4,5}
		// fuses two Strassen sub-products into one Pan pair-product, leaving
		// product 6 as a solo sub-product → 4 + 3 = 7 outer products at
		// rank 1 + 2 inner + 1 inner = 1719 total when ⟨7,7,7⟩ = 249.
		Lineage.Node lineage = new Lineage.RecombinationWithPairN(
				new Lineage.Atom("Strassen<2,2,2>=7"),
				pairing.pairs(),
				pairing.solo(),
				List.of(new Lineage.Atom("7x7x7-direct")));
		SchemeIO.write(scheme, out.toFile(), lineage);
		System.out.println("Wrote " + out);
	}
}
