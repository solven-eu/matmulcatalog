package eu.solven.matmul.benchmark;

/**
 * Result of one ({@link BenchmarkProblem}, {@link BenchmarkStrategy}) run.
 * All timing fields are wall-clock ms.
 */
public final class BenchmarkResult {

	public final Verdict verdict;
	/** Total wall time from encode-start to verdict (or timeout). */
	public final long wallMs;
	public final long encodeMs;
	public final long breakIdMs;
	public final long solveMs;

	/** CNF size handed to the solver (post-BreakID if used). */
	public final int finalVars;
	public final int finalClauses;
	/** CNF size before BreakID; equals the final size when BreakID is off. */
	public final int preBreakIdVars;
	public final int preBreakIdClauses;

	/** Optional error/exception message for {@link Verdict#ERROR} / {@link Verdict#UNAVAILABLE}. */
	public final String message;

	public BenchmarkResult(Verdict verdict, long wallMs, long encodeMs, long breakIdMs,
			long solveMs, int finalVars, int finalClauses, int preBreakIdVars,
			int preBreakIdClauses, String message) {
		this.verdict = verdict;
		this.wallMs = wallMs;
		this.encodeMs = encodeMs;
		this.breakIdMs = breakIdMs;
		this.solveMs = solveMs;
		this.finalVars = finalVars;
		this.finalClauses = finalClauses;
		this.preBreakIdVars = preBreakIdVars;
		this.preBreakIdClauses = preBreakIdClauses;
		this.message = message;
	}

	static BenchmarkResult unavailable(String message) {
		return new BenchmarkResult(Verdict.UNAVAILABLE, 0, 0, 0, 0, 0, 0, 0, 0, message);
	}

	static BenchmarkResult noEncoder(String message) {
		return new BenchmarkResult(Verdict.NO_ENCODER, 0, 0, 0, 0, 0, 0, 0, 0, message);
	}
}
