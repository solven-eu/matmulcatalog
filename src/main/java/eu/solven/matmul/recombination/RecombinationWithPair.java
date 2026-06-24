package eu.solven.matmul.recombination;

import eu.solven.matmul.papers.pan1978.PanPairProduct;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Variant of {@link Recombination#constructWithAllocation} that fuses
 * cyclically-paired sub-products using {@link PanPairProduct}'s joint
 * scheme. For each pair {@code (kA, kB)} the inner sub-algorithm is
 * replaced by a single {@link PanPairProduct.PairScheme} that computes
 * both sub-products simultaneously in {@code abc + ab + bc + ca}
 * multiplications, vs {@code 2·R(⟨a,b,c⟩)} for two independent
 * lookups.
 *
 * <p>For {@code ⟨14,14,14⟩} via Strassen[7,7]³, pairing 3 of the 7
 * sub-products into 3 pairs + 1 solo cuts the rank from
 * {@code 7·R(⟨7,7,7⟩)} (e.g. {@code 7·249 = 1743}) down to
 * {@code 3·R_pair(7,7,7) + R(⟨7,7,7⟩) = 3·490 + 249 = 1719}.</p>
 *
 * <p>Currently restricted to cubic same-shape pairs (so the cyclic
 * permutation is trivial). General cyclic pairs (⟨a,b,c⟩+⟨b,c,a⟩
 * with {@code a≠b≠c}) are handled by {@link PanPairProduct} itself,
 * but in practice the cubic case is what matters for Strassen-style
 * outers with balanced allocations.</p>
 */
public final class RecombinationWithPair {

	private RecombinationWithPair() {}

	/**
	 * Pairing decision: which outer base slots are fused, and which
	 * are solo. Slots = {@code [0, base.r)}; every slot must appear
	 * exactly once (either in a pair or as solo).
	 */
	public record Pairing(int[][] pairs, int[] solo) {
		public Pairing {
			java.util.Set<Integer> seen = new java.util.HashSet<>();
			for (int[] p : pairs) {
				if (p.length != 2) throw new IllegalArgumentException("pair must have length 2");
				for (int k : p) if (!seen.add(k))
					throw new IllegalArgumentException("slot " + k + " appears more than once");
			}
			for (int k : solo) if (!seen.add(k))
				throw new IllegalArgumentException("slot " + k + " appears more than once");
		}
		public int totalSlots() {
			return pairs.length * 2 + solo.length;
		}
	}

	/**
	 * Build the composed ⟨targetA, targetB, targetC⟩ algorithm via
	 * outer {@code base} with allocations {@code allocA, allocB, allocC},
	 * fusing the specified pairs through {@link PanPairProduct}.
	 *
	 * <p>For paired slots, both sub-shapes must equal the same cubic
	 * ⟨a,b,c⟩ (so the cyclic-permuted pair shape ⟨b,c,a⟩ matches).
	 * For solo slots, {@code lookup} provides the sub-algorithm.</p>
	 */
	public static NonCubicBilinearAlgorithm constructWithPairing(
			NonCubicBilinearAlgorithm base, Recombination.AlgorithmLookup lookup,
			int[] allocA, int[] allocB, int[] allocC, Pairing pairing) {
		if (pairing.totalSlots() != base.r) {
			throw new IllegalArgumentException(
					"pairing covers " + pairing.totalSlots() + " slots but base.r = " + base.r);
		}
		Recombination.SotaResolver sota = (a, b, c) ->
				lookup.find(a, b, c).map(alg -> alg.r).orElse(Recombination.SotaResolver.UNKNOWN_RANK);
		Recombination.Result rec = Recombination.recombineWithAllocation(
				base, sota, allocA, allocB, allocC);

		int targetA = sum(allocA), targetB = sum(allocB), targetC = sum(allocC);
		int[] cumA = cumulative(allocA);
		int[] cumB = cumulative(allocB);
		int[] cumC = cumulative(allocC);

		// Resolve sub-algorithms / pair schemes per slot.
		record SlotPlan(int[] slots, int rank, NonCubicBilinearAlgorithm solo, PanPairProduct.PairScheme pair) {}
		List<SlotPlan> plan = new ArrayList<>();
		for (int[] p : pairing.pairs()) {
			int[] szA = rec.smallMatrixSizes[p[0]];
			int[] szB = rec.smallMatrixSizes[p[1]];
			if (!Arrays.equals(szA, szB)) {
				throw new IllegalArgumentException(
						"pair " + Arrays.toString(p) + " has mismatched sub-shapes: "
								+ Arrays.toString(szA) + " vs " + Arrays.toString(szB));
			}
			if (szA[0] == 0 || szA[1] == 0 || szA[2] == 0) {
				throw new IllegalArgumentException("degenerate pair slot " + Arrays.toString(p));
			}
			PanPairProduct.PairScheme ps = PanPairProduct.build(szA[0], szA[1], szA[2]);
			plan.add(new SlotPlan(p, ps.rank(), null, ps));
		}
		for (int k : pairing.solo()) {
			int[] sz = rec.smallMatrixSizes[k];
			if (sz[0] == 0 || sz[1] == 0 || sz[2] == 0) continue;
			Optional<NonCubicBilinearAlgorithm> hit = lookup.find(sz[0], sz[1], sz[2]);
			if (hit.isEmpty()) {
				throw new IllegalStateException(
						"missing solo leaf ⟨" + sz[0] + "," + sz[1] + "," + sz[2] + "⟩ for slot " + k);
			}
			plan.add(new SlotPlan(new int[]{k}, hit.get().r, hit.get(), null));
		}

		int totalRank = plan.stream().mapToInt(SlotPlan::rank).sum();
		double[][] U = new double[targetA * targetB][totalRank];
		double[][] V = new double[targetB * targetC][totalRank];
		double[][] W = new double[targetA * targetC][totalRank];

		int kStart = 0;
		for (SlotPlan sp : plan) {
			if (sp.solo() != null) {
				embedSolo(U, V, W, base, sp.solo(), sp.slots()[0], kStart,
						cumA, cumB, cumC, allocA, allocB, allocC, targetB, targetC);
			} else {
				embedPair(U, V, W, base, sp.pair(), sp.slots()[0], sp.slots()[1], kStart,
						cumA, cumB, cumC, allocA, allocB, allocC, targetB, targetC);
			}
			kStart += sp.rank();
		}
		return new NonCubicBilinearAlgorithm(targetA, targetB, targetC, U, V, W);
	}

	private static void embedSolo(double[][] U, double[][] V, double[][] W,
			NonCubicBilinearAlgorithm base, NonCubicBilinearAlgorithm sub,
			int kBase, int kStart,
			int[] cumA, int[] cumB, int[] cumC,
			int[] allocA, int[] allocB, int[] allocC,
			int targetB, int targetC) {
		double[][] baseU = base.denseU();
		double[][] baseV = base.denseV();
		double[][] baseW = base.denseW();
		double[][] subU = sub.denseU();
		double[][] subV = sub.denseV();
		double[][] subW = sub.denseW();
		embedFactor(U, baseU, subU, kBase, kStart,
				base.n, base.m, sub.n, sub.m, cumA, cumB, targetB);
		embedFactor(V, baseV, subV, kBase, kStart,
				base.m, base.p, sub.m, sub.p, cumB, cumC, targetC);
		embedFactor(W, baseW, subW, kBase, kStart,
				base.n, base.p, sub.n, sub.p, cumA, cumC, targetC);
	}

	private static void embedPair(double[][] U, double[][] V, double[][] W,
			NonCubicBilinearAlgorithm base, PanPairProduct.PairScheme pair,
			int kA, int kB, int kStart,
			int[] cumA, int[] cumB, int[] cumC,
			int[] allocA, int[] allocB, int[] allocC,
			int targetB, int targetC) {
		double[][] baseU = base.denseU();
		double[][] baseV = base.denseV();
		double[][] baseW = base.denseW();
		int a = pair.a(), b = pair.b(), c = pair.c();
		// U[a_global][kStart + k] = base.U[abOuter][kA]·αA[innerAB][k]
		//                         + base.U[abOuter][kB]·αU[innerBC][k]
		embedPairFactor(U, baseU, kA, kB, pair.alphaA(), pair.alphaU(),
				base.n, base.m, a, b, b, c, cumA, cumB, targetB, kStart);
		// V[ab_global][kStart + k] = base.V[abOuter][kA]·βB[innerBC][k]
		//                          + base.V[abOuter][kB]·βV[innerCA][k]
		embedPairFactor(V, baseV, kA, kB, pair.betaB(), pair.betaV(),
				base.m, base.p, b, c, c, a, cumB, cumC, targetC, kStart);
		// W[ab_global][kStart + k] = base.W[abOuter][kA]·W_C[innerAC][k]
		//                          + base.W[abOuter][kB]·W_Cp[innerBA][k]
		embedPairFactor(W, baseW, kA, kB, pair.W_C(), pair.W_Cp(),
				base.n, base.p, a, c, b, a, cumA, cumC, targetC, kStart);
	}

	/**
	 * For a pair (kA, kB), each inner product k carries TWO outer-coefficient
	 * contributions, one for each Strassen slot. {@code innerA} maps the "C"
	 * (first MM) input/output dims via slot kA; {@code innerB} maps the "C'"
	 * (second MM) via slot kB.
	 */
	private static void embedPairFactor(double[][] dst, double[][] base,
			int kA, int kB,
			double[][] innerA, double[][] innerB,
			int baseRows, int baseCols,
			int subRowsA, int subColsA, int subRowsB, int subColsB,
			int[] cumRows, int[] cumCols, int targetCols, int kStart) {
		int subRank = innerA[0].length;
		int baseDim = baseRows * baseCols;
		// Process slot kA with innerA dimensions.
		for (int aBase = 0; aBase < baseDim; aBase++) {
			double cA = base[aBase][kA];
			if (cA == 0.0) continue;
			int iBase = aBase / baseCols;
			int jBase = aBase % baseCols;
			int rowOff = cumRows[iBase], colOff = cumCols[jBase];
			int blockRowSize = cumRows[iBase + 1] - cumRows[iBase];
			int blockColSize = cumCols[jBase + 1] - cumCols[jBase];
			for (int iSub = 0; iSub < subRowsA; iSub++) {
				if (iSub >= blockRowSize) continue;
				for (int jSub = 0; jSub < subColsA; jSub++) {
					if (jSub >= blockColSize) continue;
					int aSub = iSub * subColsA + jSub;
					int aGlob = (rowOff + iSub) * targetCols + (colOff + jSub);
					double[] dstRow = dst[aGlob];
					double[] srcRow = innerA[aSub];
					for (int kSub = 0; kSub < subRank; kSub++) {
						double s = srcRow[kSub];
						if (s == 0.0) continue;
						dstRow[kStart + kSub] += cA * s;
					}
				}
			}
		}
		// Process slot kB with innerB dimensions.
		for (int aBase = 0; aBase < baseDim; aBase++) {
			double cB = base[aBase][kB];
			if (cB == 0.0) continue;
			int iBase = aBase / baseCols;
			int jBase = aBase % baseCols;
			int rowOff = cumRows[iBase], colOff = cumCols[jBase];
			int blockRowSize = cumRows[iBase + 1] - cumRows[iBase];
			int blockColSize = cumCols[jBase + 1] - cumCols[jBase];
			for (int iSub = 0; iSub < subRowsB; iSub++) {
				if (iSub >= blockRowSize) continue;
				for (int jSub = 0; jSub < subColsB; jSub++) {
					if (jSub >= blockColSize) continue;
					int aSub = iSub * subColsB + jSub;
					int aGlob = (rowOff + iSub) * targetCols + (colOff + jSub);
					double[] dstRow = dst[aGlob];
					double[] srcRow = innerB[aSub];
					for (int kSub = 0; kSub < subRank; kSub++) {
						double s = srcRow[kSub];
						if (s == 0.0) continue;
						dstRow[kStart + kSub] += cB * s;
					}
				}
			}
		}
	}

	// ── helpers ──
	private static int sum(int[] a) {
		int s = 0;
		for (int x : a) s += x;
		return s;
	}
	private static int[] cumulative(int[] alloc) {
		int[] cum = new int[alloc.length + 1];
		for (int i = 0; i < alloc.length; i++) cum[i + 1] = cum[i] + alloc[i];
		return cum;
	}

	/** Standard embedFactor (mirrored from Recombination, package-private accessor). */
	private static void embedFactor(double[][] dst, double[][] base, double[][] sub,
			int kBase, int kStart,
			int baseRows, int baseCols, int subRows, int subCols,
			int[] cumRows, int[] cumCols, int targetCols) {
		int subRank = sub[0].length;
		int baseDim = baseRows * baseCols;
		int subDim = subRows * subCols;
		for (int aBase = 0; aBase < baseDim; aBase++) {
			double c = base[aBase][kBase];
			if (c == 0.0) continue;
			int iBase = aBase / baseCols;
			int jBase = aBase % baseCols;
			int rowOff = cumRows[iBase];
			int colOff = cumCols[jBase];
			int blockRowSize = cumRows[iBase + 1] - cumRows[iBase];
			int blockColSize = cumCols[jBase + 1] - cumCols[jBase];
			for (int aSub = 0; aSub < subDim; aSub++) {
				int iSub = aSub / subCols;
				int jSub = aSub % subCols;
				if (iSub >= blockRowSize || jSub >= blockColSize) continue;
				int aGlob = (rowOff + iSub) * targetCols + (colOff + jSub);
				double[] srcRow = sub[aSub];
				double[] dstRow = dst[aGlob];
				for (int kSub = 0; kSub < subRank; kSub++) {
					double s = srcRow[kSub];
					if (s == 0.0) continue;
					dstRow[kStart + kSub] += c * s;
				}
			}
		}
	}
}
