package eu.solven.matmul.additions;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * CSE-minimised additive complexity of a bilinear matmul scheme — the number of
 * additions in a (heuristically) minimal linear straight-line program that
 * evaluates it, together with the full replayable {@link Slp}s so the
 * construction can be reconstructed. Splits the scheme into its three
 * independent linear-form systems and minimises each with
 * {@link LinearCircuitMinimizer}:
 *
 * <ul>
 *   <li><b>A-side</b>: the {@code r} columns of {@code U} (over {@code n·m}
 *       A-variables) → the left factors {@code L_k}.</li>
 *   <li><b>B-side</b>: the {@code r} columns of {@code V} (over {@code m·p}
 *       B-variables) → the right factors {@code R_k}.</li>
 *   <li><b>Output</b>: the {@code n·p} rows of {@code W} (over the {@code r}
 *       product variables {@code m_k}) → the outputs {@code C_i}.</li>
 * </ul>
 *
 * <p>The three systems share no variables, so the minimum total is their sum.
 * {@link Result#reconstructs} confirms all three SLPs reproduce the scheme's
 * factor matrices exactly.</p>
 */
public final class SchemeAdditiveComplexity {

	private SchemeAdditiveComplexity() {}

	public record Result(int minimal, int naive, int aSide, int bSide, int outSide, int scalarMults,
			Slp aSlp, Slp bSlp, Slp outSlp) {
		public int saved() {
			return naive - minimal;
		}

		/** True iff every side's SLP reproduces the scheme's factor matrices. */
		public boolean reconstructs(NonCubicBilinearAlgorithm alg) {
			return aSlp.reconstructs(columnsAsForms(alg.denseU(), alg.r))
					&& bSlp.reconstructs(columnsAsForms(alg.denseV(), alg.r))
					&& outSlp.reconstructs(alg.denseW());
		}

		@Override
		public String toString() {
			return String.format(
					"additions: minimal=%d (naive=%d, saved=%d) [A=%d, B=%d, out=%d]%s",
					minimal, naive, saved(), aSide, bSide, outSide,
					scalarMults > 0 ? " +" + scalarMults + " scalar-mults" : "");
		}
	}

	public static Result analyse(NonCubicBilinearAlgorithm alg) {
		LinearCircuitMinimizer.Result a = LinearCircuitMinimizer.minimize(columnsAsForms(alg.denseU(), alg.r));
		LinearCircuitMinimizer.Result b = LinearCircuitMinimizer.minimize(columnsAsForms(alg.denseV(), alg.r));
		LinearCircuitMinimizer.Result w = LinearCircuitMinimizer.minimize(alg.denseW()); // rows already forms
		int minimal = a.additions() + b.additions() + w.additions();
		int naive = a.naiveAdditions() + b.naiveAdditions() + w.naiveAdditions();
		int scalar = a.scalarMults() + b.scalarMults() + w.scalarMults();
		return new Result(minimal, naive, a.additions(), b.additions(), w.additions(), scalar,
				a.slp(), b.slp(), w.slp());
	}

	/** Transpose a {@code [vars][r]} factor matrix into {@code r} forms, each over
	 *  the {@code vars} input variables (one form per product/column). */
	static double[][] columnsAsForms(double[][] mat, int r) {
		int vars = mat.length;
		double[][] forms = new double[r][vars];
		for (int v = 0; v < vars; v++) {
			for (int k = 0; k < r; k++) {
				forms[k][v] = mat[v][k];
			}
		}
		return forms;
	}
}
