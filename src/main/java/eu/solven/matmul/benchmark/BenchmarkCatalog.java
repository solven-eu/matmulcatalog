package eu.solven.matmul.benchmark;

import java.util.ArrayList;
import java.util.List;

import eu.solven.matmul.f2.sat.SatMatmulPipeline;

/**
 * Curated staircase of decomposition problems, organized from trivial to
 * research-grade. The intent is to enable diff-tracking the impact of each
 * strategy change across the whole ladder, so a tweak that helps the hard
 * cases doesn't silently regress the easy ones.
 *
 * <p>Problems are organized in tiers:</p>
 * <ul>
 *   <li>{@link #trivial()} — milliseconds. Always run.</li>
 *   <li>{@link #fast()} — seconds. Run by default.</li>
 *   <li>{@link #slow()} — minutes to an hour. Run on demand.</li>
 *   <li>{@link #research()} — likely unsolved. For overnight runs only.</li>
 * </ul>
 */
public final class BenchmarkCatalog {

	private BenchmarkCatalog() {}

	public static List<BenchmarkProblem> trivial() {
		List<BenchmarkProblem> ps = new ArrayList<>();
		// ⟨2,2,2⟩ upper-triangular: tight at r=4 (Strassen-style triangular).
		ps.add(triangular222Z2(4));
		ps.add(triangular222Z2(3));
		return ps;
	}

	public static List<BenchmarkProblem> fast() {
		List<BenchmarkProblem> ps = new ArrayList<>();
		// ⟨2,2,2⟩ dense Z/2: Strassen tight at r=7.
		ps.add(dense222Z2(7));
		ps.add(dense222Z2(6));
		// ⟨3,3,3⟩ diagonal-plus-(0,1): tight at r=5.
		ps.add(diagPlusOne333Z2(5));
		ps.add(diagPlusOne333Z2(4));
		return ps;
	}

	public static List<BenchmarkProblem> slow() {
		List<BenchmarkProblem> ps = new ArrayList<>();
		// ⟨2,2,3⟩ dense Z/2: known upper bound r=11. SAT4J doesn't scale well
		// beyond ⟨2,2,2⟩, so this tier really only works for kissat-based strategies.
		ps.add(dense223Z2(11));
		// ⟨2,3,3⟩ dense Z/2: believed-tight 15; r=14 is the research target.
		ps.add(dense233Z2(15));
		return ps;
	}

	public static List<BenchmarkProblem> research() {
		List<BenchmarkProblem> ps = new ArrayList<>();
		// Phase 1.6 research target: would prove R_{Z/2}(⟨2,3,3⟩) ≥ 15 if UNSAT.
		ps.add(dense233Z2(14));
		// Phase 2: Laderman-equivalent Z/2 algorithm.
		ps.add(dense333Z2(23));
		return ps;
	}

	// ───────────────────────────────────────────────────────────────────────────
	// Z/2 problem builders
	// ───────────────────────────────────────────────────────────────────────────

	public static BenchmarkProblem dense222Z2(int r) {
		int[][][] t = SatMatmulPipeline.z2DenseMatmulTensor(2, 2, 2);
		return new BenchmarkProblem("222_dense_z2_r" + r, Alphabet.Z2, 2, 2, 2,
				"dense", t, null, r);
	}

	public static BenchmarkProblem dense223Z2(int r) {
		int[][][] t = SatMatmulPipeline.z2DenseMatmulTensor(2, 2, 3);
		return new BenchmarkProblem("223_dense_z2_r" + r, Alphabet.Z2, 2, 2, 3,
				"dense", t, null, r);
	}

	public static BenchmarkProblem dense233Z2(int r) {
		int[][][] t = SatMatmulPipeline.z2DenseMatmulTensor(2, 3, 3);
		return new BenchmarkProblem("233_dense_z2_r" + r, Alphabet.Z2, 2, 3, 3,
				"dense", t, null, r);
	}

	public static BenchmarkProblem dense333Z2(int r) {
		int[][][] t = SatMatmulPipeline.z2DenseMatmulTensor(3, 3, 3);
		return new BenchmarkProblem("333_dense_z2_r" + r, Alphabet.Z2, 3, 3, 3,
				"dense", t, null, r);
	}

	public static BenchmarkProblem triangular222Z2(int r) {
		int[][][] t = SatMatmulPipeline.z2UpperTriangularMatmulTensor(2);
		return new BenchmarkProblem("222_uppertri_z2_r" + r, Alphabet.Z2, 2, 2, 2,
				"upper-tri", t, null, r);
	}

	/** ⟨3,3,3⟩ restricted: only the diagonal (0,0), (1,1), (2,2) plus (0,1) is non-zero. */
	public static BenchmarkProblem diagPlusOne333Z2(int r) {
		boolean[] mask = SatMatmulPipeline.diagonalPlusOne(3, 0, 1);
		int[][][] t = SatMatmulPipeline.z2RestrictedMatmulTensor(3, mask);
		boolean[] forceZero = new boolean[9];
		for (int i = 0; i < 9; i++) forceZero[i] = !mask[i];
		boolean[][] perSlot = { forceZero, forceZero, forceZero };
		return new BenchmarkProblem("333_diag+1_z2_r" + r, Alphabet.Z2, 3, 3, 3,
				"diag+1", t, perSlot, r);
	}
}
