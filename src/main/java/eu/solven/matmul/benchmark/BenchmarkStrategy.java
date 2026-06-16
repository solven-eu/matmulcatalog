package eu.solven.matmul.benchmark;

/**
 * A combination of encoder options + symmetry-breaking preprocessor + solver.
 * Keep {@code id} filesystem-safe — it's used as a directory name.
 */
public final class BenchmarkStrategy {

	/** Default wall-time budget per run (10 minutes). */
	public static final long DEFAULT_TIMEOUT_SEC = 600;

	public final String id;
	/** Hand-coded column lex-ordering in {@code Z2CnfEncoder}. Essential for SAT4J on hard ranks. */
	public final boolean columnLex;
	/** Run BreakID as a CNF→CNF preprocessor before invoking the solver. */
	public final boolean useBreakId;
	public final SolverChoice solver;
	/** Wall-time budget in seconds. Solvers honour this where they can ({@code --time}, {@code --maxtime}, SAT4J {@code setTimeout}). */
	public final long timeoutSec;

	public BenchmarkStrategy(String id, boolean columnLex, boolean useBreakId, SolverChoice solver) {
		this(id, columnLex, useBreakId, solver, DEFAULT_TIMEOUT_SEC);
	}

	public BenchmarkStrategy(String id, boolean columnLex, boolean useBreakId, SolverChoice solver,
			long timeoutSec) {
		this.id = id;
		this.columnLex = columnLex;
		this.useBreakId = useBreakId;
		this.solver = solver;
		this.timeoutSec = timeoutSec;
	}

	@Override
	public String toString() {
		return String.format("%s{lex=%b,breakid=%b,solver=%s,timeoutSec=%d}",
				id, columnLex, useBreakId, solver, timeoutSec);
	}
}
