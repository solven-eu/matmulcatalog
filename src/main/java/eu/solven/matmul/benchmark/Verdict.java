package eu.solven.matmul.benchmark;

/** Outcome of a single benchmark run. */
public enum Verdict {
	SAT,
	UNSAT,
	/** Solver started but did not return within the wall-time budget. */
	TIMEOUT,
	/** Solver returned an error / SAT result that failed verification. */
	ERROR,
	/** No encoder is implemented for the alphabet ({@code Z3} or {@code Z}); row recorded as placeholder. */
	NO_ENCODER,
	/** A required external binary (e.g. kissat, BreakID) was not on PATH. */
	UNAVAILABLE
}
