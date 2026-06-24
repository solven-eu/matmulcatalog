package eu.solven.matmul.research;

import eu.solven.matmul.search.CitedBound;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.LineageTrackingLookup;
import eu.solven.matmul.recombination.Recombination;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Materialises every cubic ⟨n,n,n⟩ bound from
 * {@code NEW_BOUNDS.md} that this repo claims as strict improvement
 * over DIS09 and that is constructable without Layer-4 paired
 * sub-products. For each (n, allocation), runs
 * {@link Recombination#constructWithAllocation}, spot-checks via
 * {@link Verifier#passesRandomMatmulSpotCheck}, and writes the scheme.
 *
 * <p>Skipped: n=14 — the new bound 1720 requires Layer 4 paired
 * sub-products which isn't a constructive recipe yet; un-paired
 * Strassen[7,7]³ gives 1750 which doesn't beat DIS09's 1728.</p>
 */
public final class MaterializeSolvenStrassenNew {

	private MaterializeSolvenStrassenNew() {}

	private record Bound(int n, int[] allocA, int[] allocB, int[] allocC, int dis09) {}

	public static void main(String[] args) throws Exception {
		Bound[] bounds = {
				// n=14: Strassen[7,7]³ — depends on best ⟨7,7,7⟩ available
				new Bound(14, new int[]{7, 7}, new int[]{7, 7}, new int[]{7, 7}, 1728),
				// n=17: Strassen[8,9]³ — base[8,9]/[8,9]/[8,9]
				new Bound(17, new int[]{8, 9}, new int[]{8, 9}, new int[]{8, 9}, 2898),
				// n=19: Strassen[9,10]³ — without paired gives ~4044 (DIS09 4073)
				new Bound(19, new int[]{9, 10}, new int[]{9, 10}, new int[]{9, 10}, 4073),
				// n=21: Strassen[9,12]/[9,12]/[10,11] non-balanced (done previously but re-runs OK)
				new Bound(21, new int[]{9, 12}, new int[]{9, 12}, new int[]{10, 11}, 5365),
				// n=23: Strassen[11,12]³ — without paired ~6738 (DIS09 6806)
				new Bound(23, new int[]{11, 12}, new int[]{11, 12}, new int[]{11, 12}, 6806),
		};

		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		FieldAwareLookup lookup = new FieldAwareLookup("R");
		// Search uses formula-aware bounds (Pan TA + Hopcroft-Kerr) on top of the
		// catalog so unmatched leaves get a tighter estimate. Materialisation
		// still uses `lookup` to find real schemes.
		eu.solven.matmul.recombination.Recombination.SotaResolver sota =
				new eu.solven.matmul.search.CitedBound(lookup, false);

		for (Bound b : bounds) {
			System.out.printf("%n=== ⟨%d,%d,%d⟩ via Strassen + allocs %s / %s / %s ===%n",
					b.n, b.n, b.n,
					java.util.Arrays.toString(b.allocA),
					java.util.Arrays.toString(b.allocB),
					java.util.Arrays.toString(b.allocC));

			NonCubicBilinearAlgorithm scheme;
			LineageTrackingLookup tracker = new LineageTrackingLookup(lookup);
			try {
				scheme = Recombination.constructWithAllocation(strassen, tracker, sota,
						b.allocA, b.allocB, b.allocC);
			} catch (RuntimeException e) {
				System.err.println("  CONSTRUCT FAILED: " + e.getMessage());
				continue;
			}
			System.out.printf("  Built scheme: ⟨%d,%d,%d⟩ rank=%d%n",
					scheme.n, scheme.m, scheme.p, scheme.r);

			long t0 = System.currentTimeMillis();
			boolean ok = Verifier.passesRandomMatmulSpotCheck(scheme);
			System.out.printf("  Spot-check: %s (took %d ms)%n",
					ok ? "PASS" : "FAIL", System.currentTimeMillis() - t0);
			if (!ok) {
				System.err.println("  ABORTING write (scheme is broken)");
				continue;
			}

			int additions = Verifier.additionCount(scheme);
			System.out.printf("  Additions: %d  |  vs DIS09: %s (Δ=%d)%n",
					additions,
					scheme.r < b.dis09 ? "WIN"
							: (scheme.r == b.dis09 ? "TIE" : "LOSE"),
					scheme.r - b.dis09);

			Path dir = Path.of("src/main/resources/schemes/known/section" + b.n);
			Files.createDirectories(dir);
			File out = dir.resolve(String.format(
					"solven_strassen_2026_%dx%dx%d_r%d_a%d.json",
					b.n, b.n, b.n, scheme.r, additions)).toFile();
			Lineage.Node lineage = new Lineage.RecombinationN(
					new Lineage.Atom("Strassen<2,2,2>=7"),
					b.allocA.clone(), b.allocB.clone(), b.allocC.clone(),
					tracker.leaves());
			SchemeIO.write(scheme, out, lineage);
			System.out.println("  Wrote " + out);
		}
	}
}
