package eu.solven.matmul.benchmark;

import java.util.List;

import eu.solven.matmul.f2.sat.CryptoMiniSatSolver;
import eu.solven.matmul.f2.sat.KissatSolver;
import eu.solven.matmul.f2.sat.Sat4jSolver;
import eu.solven.matmul.f2.sat.SolverOutcome;

/**
 * Which SAT solver to drive. {@link #KISSAT} and {@link #CRYPTOMINISAT} require
 * the respective binary on {@code $PATH}; {@link #SAT4J} is pure-Java and
 * always available.
 */
public enum SolverChoice {
	SAT4J {
		@Override
		public boolean isAvailable() { return true; }

		@Override
		public SolverOutcome solveBounded(int varCount, List<int[]> clauses, long timeoutSec) {
			return Sat4jSolver.solveBounded(varCount, clauses, timeoutSec);
		}
	},
	KISSAT {
		@Override
		public boolean isAvailable() { return KissatSolver.isAvailable(); }

		@Override
		public SolverOutcome solveBounded(int varCount, List<int[]> clauses, long timeoutSec) {
			return KissatSolver.solveBounded(varCount, clauses, timeoutSec);
		}
	},
	CRYPTOMINISAT {
		@Override
		public boolean isAvailable() { return CryptoMiniSatSolver.isAvailable(); }

		@Override
		public SolverOutcome solveBounded(int varCount, List<int[]> clauses, long timeoutSec) {
			return CryptoMiniSatSolver.solveBounded(varCount, clauses, timeoutSec);
		}
	};

	public abstract boolean isAvailable();

	/**
	 * @param timeoutSec wall-time budget in seconds; {@code 0} = no limit.
	 *                   On timeout the outcome will be {@link SolverOutcome#timeout()}.
	 */
	public abstract SolverOutcome solveBounded(int varCount, List<int[]> clauses, long timeoutSec);
}
