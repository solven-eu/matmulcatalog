package eu.solven.matmul.research;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;

import java.util.List;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.papers.pan1978.PanPairProduct;
import eu.solven.matmul.recombination.RecombinationWithPair;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Generalises {@link MaterializeViaPanPair14}: for cubic
 * {@code ⟨2k, 2k, 2k⟩} targets via Strassen[k,k]³, pair 3 of the 7
 * Strassen sub-products into Pan pair-products (rank {@code k³+3k²}
 * each) and leave 1 solo (rank {@code R(⟨k,k,k⟩)}). Net effect: rank
 * drops from {@code 7·R(⟨k,k,k⟩)} to {@code 3·(k³+3k²) + R(⟨k,k,k⟩)},
 * a save of {@code 4·R(⟨k,k,k⟩) − 3·(k³+3k²)} whenever positive.
 *
 * <p>Pan pair-product saves multiplications iff
 * {@code 2·R(⟨k,k,k⟩) > k³ + 3k²}, i.e. when the leaf rank exceeds
 * {@code (k³ + 3k²)/2}. For our catalog at the time of writing:</p>
 *
 * <table>
 *   <caption>Per-cubic-size profitability of Pan pair fusion in Strassen[k,k]³</caption>
 *   <tr><th>k</th><th>k³+3k²</th><th>R(⟨k,k,k⟩)</th><th>2R</th><th>pair saves?</th><th>net save/pair</th></tr>
 *   <tr><td>7</td><td>490</td><td>249</td><td>498</td><td>YES</td><td>8</td></tr>
 *   <tr><td>8</td><td>704</td><td>336</td><td>672</td><td>NO</td><td>-32</td></tr>
 *   <tr><td>9</td><td>972</td><td>486</td><td>972</td><td>BREAK-EVEN</td><td>0</td></tr>
 *   <tr><td>10</td><td>1300</td><td>?</td><td>?</td><td>?</td><td>?</td></tr>
 *   <tr><td>11</td><td>1694</td><td>896</td><td>1792</td><td>YES</td><td>98</td></tr>
 *   <tr><td>12</td><td>2160</td><td>1040</td><td>2080</td><td>NO</td><td>-80</td></tr>
 * </table>
 *
 * <p>So this materialises ⟨14,14,14⟩, ⟨22,22,22⟩ and similar — wherever
 * the inequality holds.</p>
 */
public final class MaterializeViaPanPair {

	public static void main(String[] args) throws Exception {
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		FieldAwareLookup lookup = new FieldAwareLookup("Q");

		// Probe each candidate k.
		int[] kCandidates = { 7, 8, 9, 10, 11, 12 };
		for (int k : kCandidates) {
			int leafRank = lookup.find(k, k, k).map(alg -> alg.r).orElse(-1);
			int pairCost = PanPairProduct.rank(k, k, k);
			int viaPair = 3 * pairCost + leafRank;
			int unpaired = 7 * leafRank;
			int net = unpaired - viaPair;
			System.out.printf("k=%-2d  R(⟨%d,%d,%d⟩)=%-5d  pair-cost=%-5d  unpaired=%-5d  via-pair=%-5d  net=%+d%n",
					k, k, k, k, leafRank, pairCost, unpaired, viaPair, net);
			if (net <= 0 || leafRank < 0) {
				System.out.println("    skip (pair not profitable)");
				continue;
			}

			int n = 2 * k;
			RecombinationWithPair.Pairing pairing = new RecombinationWithPair.Pairing(
					new int[][] { {0, 1}, {2, 3}, {4, 5} }, new int[] { 6 });
			NonCubicBilinearAlgorithm scheme = RecombinationWithPair.constructWithPairing(
					strassen, lookup, new int[]{k, k}, new int[]{k, k}, new int[]{k, k},
					pairing);

			if (scheme.r != viaPair) {
				System.err.printf("    UNEXPECTED rank mismatch: predicted %d, got %d%n", viaPair, scheme.r);
				continue;
			}
			if (!Verifier.passesRandomMatmulSpotCheck(scheme)) {
				System.err.println("    SPOT-CHECK FAILED");
				continue;
			}
			int adds = Verifier.additionCount(scheme);
			Path out = Path.of(String.format(
					"src/main/resources/schemes/known/section%d/solven_strassen_2026-%dx%dx%d_m%d_a%d.json",
					n, n, n, n, scheme.r, adds));
			Files.createDirectories(out.getParent());
			Lineage.Node lineage = new Lineage.RecombinationWithPairN(
					new Lineage.Atom("Strassen<2,2,2>=7"),
					pairing.pairs(),
					pairing.solo(),
					List.of(new Lineage.Atom(k + "x" + k + "x" + k + "-direct")));
			SchemeIO.write(scheme, out.toFile(), lineage);
			System.out.printf("    wrote %s  (adds=%d)%n", out, adds);
		}
	}
}
