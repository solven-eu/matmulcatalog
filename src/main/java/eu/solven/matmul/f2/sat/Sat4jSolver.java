package eu.solven.matmul.f2.sat;

import java.util.List;
import java.util.Optional;

import org.sat4j.core.VecInt;
import org.sat4j.minisat.SolverFactory;
import org.sat4j.specs.ContradictionException;
import org.sat4j.specs.ISolver;
import org.sat4j.specs.TimeoutException;

/**
 * Thin wrapper around SAT4J. Takes a CNF formula in DIMACS-style integer
 * representation and returns an assignment if satisfiable, {@link Optional#empty()}
 * otherwise.
 */
public final class Sat4jSolver {

	private Sat4jSolver() {}

	/**
	 * @param varCount number of distinct boolean variables (1-indexed in clauses)
	 * @param clauses  each clause is an int[] of signed literals
	 * @return         the satisfying assignment (length {@code varCount},
	 *                 {@code assignment[i]} = value of variable {@code i+1}) if SAT,
	 *                 {@code Optional.empty()} if UNSAT
	 */
	public static Optional<boolean[]> solve(int varCount, List<int[]> clauses) {
		SolverOutcome out = solveBounded(varCount, clauses, 0);
		switch (out.state) {
		case SAT:
			return out.assignment;
		case UNSAT:
			return Optional.empty();
		default:
			throw new RuntimeException("SAT4J: " + out.state
					+ (out.message == null ? "" : ": " + out.message));
		}
	}

	/**
	 * Solve with an explicit wall-time budget. Pass {@code timeoutSec = 0} for
	 * unlimited. Note: SAT4J's {@code setTimeout(int)} is in CPU-seconds, which
	 * approximates wall-time on single-threaded workloads. {@link SolverOutcome#timeout()}
	 * is returned if the limit is hit.
	 */
	public static SolverOutcome solveBounded(int varCount, List<int[]> clauses, long timeoutSec) {
		ISolver solver = SolverFactory.newDefault();
		solver.newVar(varCount);
		solver.setExpectedNumberOfClauses(clauses.size());
		if (timeoutSec > 0) {
			solver.setTimeout((int) Math.min(timeoutSec, Integer.MAX_VALUE));
		}

		try {
			for (int[] clause : clauses) {
				solver.addClause(new VecInt(clause));
			}
		} catch (ContradictionException e) {
			return SolverOutcome.unsat();
		}

		boolean satisfiable;
		try {
			satisfiable = solver.isSatisfiable();
		} catch (TimeoutException e) {
			return SolverOutcome.timeout();
		}

		if (!satisfiable) return SolverOutcome.unsat();

		int[] model = solver.model();
		boolean[] assignment = new boolean[varCount];
		for (int lit : model) {
			int absLit = Math.abs(lit);
			if (absLit >= 1 && absLit <= varCount) {
				assignment[absLit - 1] = lit > 0;
			}
		}
		return SolverOutcome.sat(assignment);
	}
}
