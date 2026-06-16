package eu.solven.matmul.f2.sat;

import java.util.Optional;

/**
 * Tri-state solver result: SAT (with assignment), UNSAT, or TIMEOUT — the
 * latter being what you get when {@code solveBounded} was given a wall-time
 * budget that the solver couldn't beat.
 */
public final class SolverOutcome {

	public enum State { SAT, UNSAT, TIMEOUT, ERROR }

	public final State state;
	public final Optional<boolean[]> assignment;
	public final String message;

	private SolverOutcome(State state, Optional<boolean[]> assignment, String message) {
		this.state = state;
		this.assignment = assignment;
		this.message = message;
	}

	public static SolverOutcome sat(boolean[] assignment) {
		return new SolverOutcome(State.SAT, Optional.of(assignment), null);
	}

	public static SolverOutcome unsat() {
		return new SolverOutcome(State.UNSAT, Optional.empty(), null);
	}

	public static SolverOutcome timeout() {
		return new SolverOutcome(State.TIMEOUT, Optional.empty(), null);
	}

	public static SolverOutcome error(String message) {
		return new SolverOutcome(State.ERROR, Optional.empty(), message);
	}
}
