package eu.solven.matmul.search;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.AxisSplitBases;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.PairedSubProducts;
import eu.solven.matmul.catalog.Recombination;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.papers.laderman1976.Laderman23;

/**
 * Diagnostic for the {@code Q⟨17,17,17⟩} gap: FMM-Lille reports rank 2934
 * via a Strassen ⟨2,2,2⟩ recombination on the [9,8]³ allocation, with the
 * diagonal ⟨9,9,9⟩ products pair-fused via Pan's pair-product trick. Our
 * {@link BlockSplitSearch#findBestMultiBaseSplit} currently returns 2940 —
 * 6 too many — even though
 * {@link PairedSubProducts#applyPairing} is wired in.
 *
 * <p>This test:</p>
 * <ol>
 *   <li>Builds the standard "root pool" (Strassen, Laderman, AT 2x2x3,
 *       AT 2x3x3, mul211, mul121, mul112; 7 entries — close enough to the
 *       "8 entries" target the task description quotes).</li>
 *   <li>Runs the search for ⟨17,17,17⟩ with a fresh {@code CitedBound}
 *       over field {@code Q}.</li>
 *   <li>Prints the winning (base, allocation, rank) plus the raw sub-shape
 *       list, then manually re-runs {@code applyPairing} on the [9,8]³
 *       Strassen allocation to compare.</li>
 * </ol>
 */
public class TestPairFusingDiagonal17 {

	private static List<BlockSplitSearch.NamedBase> rootPool() throws Exception {
		List<BlockSplitSearch.NamedBase> pool = new ArrayList<>();

		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		pool.add(new BlockSplitSearch.NamedBase("Strassen ⟨2,2,2⟩=7", strassen));

		NonCubicBilinearAlgorithm laderman = NonCubicBilinearAlgorithm.fromCubic(Laderman23.get());
		pool.add(new BlockSplitSearch.NamedBase("Laderman ⟨3,3,3⟩=23", laderman));

		tryAdd(pool, "AT-Z ⟨2,2,3⟩=11",
				"src/main/resources/schemes/known/section3/alphatensor_Z-2x2x3_m11_a25.json");
		tryAdd(pool, "AT-Z ⟨2,3,3⟩=15",
				"src/main/resources/schemes/known/section3/alphatensor_Z-2x3x3_m15_a58.json");

		pool.add(new BlockSplitSearch.NamedBase("mul211 ⟨2,1,1⟩=2", AxisSplitBases.mul211()));
		pool.add(new BlockSplitSearch.NamedBase("mul121 ⟨1,2,1⟩=2", AxisSplitBases.mul121()));
		pool.add(new BlockSplitSearch.NamedBase("mul112 ⟨1,1,2⟩=2", AxisSplitBases.mul112()));
		return pool;
	}

	private static void tryAdd(List<BlockSplitSearch.NamedBase> pool, String label, String path) {
		try {
			File f = new File(path);
			if (f.exists()) {
				pool.add(new BlockSplitSearch.NamedBase(label, SchemeIO.read(f)));
			}
		} catch (Exception ignored) {
			// best-effort
		}
	}

	@Test
	public void diagnose_17x17x17_pair_fuse_gap() throws Exception {
		List<BlockSplitSearch.NamedBase> pool = rootPool();
		System.out.printf("Pool size: %d%n", pool.size());

		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		Recombination.SotaResolver sota = new CitedBound(lookup, false);

		// --- (1) What does the balanced-only search return? (Cheap.)
		Optional<BlockSplitSearch.MultiBaseSplitCandidate> best =
				BlockSplitSearch.findBestMultiBaseSplit(17, 17, 17, pool, sota,
						true /* balancedOnly */, Long.MAX_VALUE);
		assertThat(best).isPresent();
		BlockSplitSearch.MultiBaseSplitCandidate c = best.get();
		System.out.println("=== Best from balanced-only search ===");
		System.out.printf("  base   = %s%n", c.baseLabel());
		System.out.printf("  alloc  = A=%s B=%s C=%s%n",
				Arrays.toString(c.allocA()), Arrays.toString(c.allocB()), Arrays.toString(c.allocC()));
		System.out.printf("  rank   = %d (FMM target: 2934, gap=%d)%n",
				c.rank(), c.rank() - 2934);

		// --- (2) Manually inspect the [9,8]³ Strassen allocation: what sub-shapes does it
		//     produce, what's the unpaired total, what does applyPairing return?
		NonCubicBilinearAlgorithm strassen = pool.get(0).base();
		assertThat(strassen.n).isEqualTo(2);

		int[] a98 = { 9, 8 };
		Recombination.Result r = Recombination.recombineWithAllocation(strassen, sota, a98, a98, a98);
		System.out.println();
		System.out.println("=== Strassen [9,8]³ recombination ===");
		System.out.printf("  unpaired total rank = %d%n", r.totalRank);
		System.out.println("  sub-product shapes:");
		long sum = 0;
		for (int k = 0; k < r.smallMatrixSizes.length; k++) {
			int[] s = r.smallMatrixSizes[k];
			int rk = sota.getRank(s[0], s[1], s[2]);
			sum += rk;
			System.out.printf("    M%-2d: ⟨%d,%d,%d⟩ → rank=%d%n",
					k + 1, s[0], s[1], s[2], rk);
		}
		System.out.printf("  (sum check = %d)%n", sum);

		// --- (3) Pairing analysis on those sub-shapes
		long paired = PairedSubProducts.applyPairing(r.smallMatrixSizes, sota);
		System.out.println();
		System.out.println("=== applyPairing on [9,8]³ Strassen sub-shapes ===");
		System.out.printf("  paired total rank = %d (vs unpaired %d, savings %d)%n",
				paired, r.totalRank, r.totalRank - paired);

		// Look for cyclically-equivalent pairs
		int[][] subs = r.smallMatrixSizes;
		System.out.println("  pairwise cyclic-equivalence:");
		for (int i = 0; i < subs.length; i++) {
			for (int j = i + 1; j < subs.length; j++) {
				boolean cyc = PairedSubProducts.cyclicallyEquivalent(subs[i], subs[j]);
				if (cyc) {
					long pcost = PairedSubProducts.pairCost(subs[i][0], subs[i][1], subs[i][2]);
					long indiv = sota.getRank(subs[i][0], subs[i][1], subs[i][2])
							+ sota.getRank(subs[j][0], subs[j][1], subs[j][2]);
					System.out.printf("    M%d ⟨%d,%d,%d⟩ ~ M%d ⟨%d,%d,%d⟩  pairCost=%d  indiv=%d  saving=%d%n",
							i + 1, subs[i][0], subs[i][1], subs[i][2],
							j + 1, subs[j][0], subs[j][1], subs[j][2],
							pcost, indiv, indiv - pcost);
				}
			}
		}

		// --- (4) FMM's claimed recipe: explicitly assemble the 5 solos + 2 diag ⟨9,9,9⟩
		System.out.println();
		System.out.println("=== FMM's claimed recipe (reference) ===");
		int[][] fmmShapes = {
				{ 9, 9, 8 }, { 8, 8, 8 }, { 8, 8, 8 }, { 8, 9, 9 }, { 9, 8, 9 },
				{ 9, 9, 9 }, { 9, 9, 9 }
		};
		long fmmUnpaired = 0;
		for (int[] s : fmmShapes) {
			int rk = sota.getRank(s[0], s[1], s[2]);
			System.out.printf("  ⟨%d,%d,%d⟩ → rank=%d%n", s[0], s[1], s[2], rk);
			fmmUnpaired += rk;
		}
		long fmmPaired = PairedSubProducts.applyPairing(fmmShapes, sota);
		System.out.printf("  fmm-unpaired = %d, fmm-paired = %d%n", fmmUnpaired, fmmPaired);

		// --- (5) Try other ⟨2,2,2⟩=7 variants on [9,8]³ — Winograd, fmm-lille
		System.out.println();
		System.out.println("=== Other ⟨2,2,2⟩=7 schemes on [9,8]³ ===");
		String[] variants = {
				"src/main/resources/schemes/known/section2/winograd_1971-2x2x2_m7_a24.json",
		};
		for (String pth : variants) {
			File f = new File(pth);
			if (!f.exists()) continue;
			NonCubicBilinearAlgorithm alg = SchemeIO.read(f);
			Recombination.Result rr = Recombination.recombineWithAllocation(alg, sota, a98, a98, a98);
			long pp = PairedSubProducts.applyPairing(rr.smallMatrixSizes, sota);
			System.out.printf("  %s: unpaired=%d, paired=%d, shapes=%s%n",
					f.getName(), rr.totalRank, pp,
					Arrays.deepToString(rr.smallMatrixSizes));
		}

		// --- (6) Try transposed allocation [8,9]³
		System.out.println();
		System.out.println("=== Strassen [8,9]³ recombination ===");
		int[] a89 = { 8, 9 };
		Recombination.Result rT = Recombination.recombineWithAllocation(strassen, sota, a89, a89, a89);
		long pT = PairedSubProducts.applyPairing(rT.smallMatrixSizes, sota);
		System.out.printf("  unpaired=%d, paired=%d, shapes=%s%n",
				rT.totalRank, pT, Arrays.deepToString(rT.smallMatrixSizes));

		// --- (7) Mixed-orientation [9,8] / [8,9] / [9,8] (a 2-mode CW-like alloc)
		System.out.println();
		System.out.println("=== Strassen mixed [9,8]/[8,9]/[9,8] ===");
		Recombination.Result rM = Recombination.recombineWithAllocation(strassen, sota, a98, a89, a98);
		long pM = PairedSubProducts.applyPairing(rM.smallMatrixSizes, sota);
		System.out.printf("  unpaired=%d, paired=%d, shapes=%s%n",
				rM.totalRank, pM, Arrays.deepToString(rM.smallMatrixSizes));

		// --- (8) DIS09 S_{X,Y,Z} internal-permutation orbit on Strassen.
		// For 2×2×2, X,Y,Z ∈ S_2 = {I, J} (J = row-swap). 8 candidate variants
		// (some collide under signature deduplication). Currently NOT enumerated
		// by SymmetryTransforms.s3Orbit (which only does slot-permutations).
		System.out.println();
		System.out.println("=== DIS09 S_{X,Y,Z} (internal perms) on Strassen, [9,8]³ ===");
		long bestVariantTotal = Long.MAX_VALUE;
		String bestVariantLabel = null;
		java.util.Set<String> seenSigs = new java.util.HashSet<>();
		for (int x = 0; x < 2; x++) {
			for (int y = 0; y < 2; y++) {
				for (int z = 0; z < 2; z++) {
					NonCubicBilinearAlgorithm var = applyInternalPerm(strassen, x == 1, y == 1, z == 1);
					String sig = signature(var);
					if (!seenSigs.add(sig)) continue;
					Recombination.Result rv = Recombination.recombineWithAllocation(var, sota, a98, a98, a98);
					long pv = PairedSubProducts.applyPairing(rv.smallMatrixSizes, sota);
					long bestV = Math.min(rv.totalRank, pv);
					int[] mult = countMultiset(rv.smallMatrixSizes);
					System.out.printf("  X=%s Y=%s Z=%s : unpaired=%d, paired=%d, multiset[⟨888⟩,⟨889⟩-c,⟨899⟩-c,⟨999⟩]=%s%n",
							x == 1 ? "J" : "I", y == 1 ? "J" : "I", z == 1 ? "J" : "I",
							rv.totalRank, pv, Arrays.toString(mult));
					if (bestV < bestVariantTotal) {
						bestVariantTotal = bestV;
						bestVariantLabel = "X=" + (x == 1 ? "J" : "I")
								+ " Y=" + (y == 1 ? "J" : "I")
								+ " Z=" + (z == 1 ? "J" : "I");
					}
				}
			}
		}
		System.out.printf("  BEST internal-perm variant: %s → %d (FMM target 2934, our base 2940)%n",
				bestVariantLabel, bestVariantTotal);
	}

	/**
	 * Returns the multiset count [c⟨8,8,8⟩, c⟨8,8,9⟩-cyclic, c⟨8,9,9⟩-cyclic, c⟨9,9,9⟩]
	 * for an array of sub-shapes (only shapes drawn from {8,9}³ are counted).
	 */
	private static int[] countMultiset(int[][] shapes) {
		int c888 = 0, c889 = 0, c899 = 0, c999 = 0;
		for (int[] s : shapes) {
			int n8 = 0, n9 = 0;
			for (int v : s) { if (v == 8) n8++; else if (v == 9) n9++; }
			if (n8 == 3) c888++;
			else if (n8 == 2 && n9 == 1) c889++;
			else if (n8 == 1 && n9 == 2) c899++;
			else if (n9 == 3) c999++;
		}
		return new int[] { c888, c889, c899, c999 };
	}

	/**
	 * DIS09 §3 transform S_{X,Y,Z}(U,V,W) with X,Y,Z restricted to 2×2 permutation
	 * matrices. {@code swapA,swapB,swapC} toggle whether to apply the row-swap J on the
	 * A/B/C axis. For 2×2×2 Strassen this gives 8 candidate variants (some collide
	 * under signature deduplication).
	 */
	private static NonCubicBilinearAlgorithm applyInternalPerm(NonCubicBilinearAlgorithm alg,
			boolean swapA, boolean swapB, boolean swapC) {
		int a = alg.n, b = alg.m, c = alg.p, r = alg.r;
		double[][] srcU = alg.denseU();
		double[][] srcV = alg.denseV();
		double[][] srcW = alg.denseW();
		double[][] U2 = new double[a * b][r];
		double[][] V2 = new double[b * c][r];
		double[][] W2 = new double[a * c][r];
		for (int i = 0; i < a; i++) for (int j = 0; j < b; j++) {
			int iP = swapA ? (a - 1 - i) : i;
			int jP = swapB ? (b - 1 - j) : j;
			for (int k = 0; k < r; k++) U2[iP * b + jP][k] = srcU[i * b + j][k];
		}
		for (int j = 0; j < b; j++) for (int l = 0; l < c; l++) {
			int jP = swapB ? (b - 1 - j) : j;
			int lP = swapC ? (c - 1 - l) : l;
			for (int k = 0; k < r; k++) V2[jP * c + lP][k] = srcV[j * c + l][k];
		}
		for (int i = 0; i < a; i++) for (int l = 0; l < c; l++) {
			int iP = swapA ? (a - 1 - i) : i;
			int lP = swapC ? (c - 1 - l) : l;
			for (int k = 0; k < r; k++) W2[iP * c + lP][k] = srcW[i * c + l][k];
		}
		return new NonCubicBilinearAlgorithm(a, b, c, U2, V2, W2);
	}

	private static String signature(NonCubicBilinearAlgorithm a) {
		double[][] srcU = a.denseU();
		double[][] srcV = a.denseV();
		double[][] srcW = a.denseW();
		StringBuilder sb = new StringBuilder();
		for (double[] row : srcU) for (double v : row) sb.append(v).append(',');
		sb.append('|');
		for (double[] row : srcV) for (double v : row) sb.append(v).append(',');
		sb.append('|');
		for (double[] row : srcW) for (double v : row) sb.append(v).append(',');
		return sb.toString();
	}
}
