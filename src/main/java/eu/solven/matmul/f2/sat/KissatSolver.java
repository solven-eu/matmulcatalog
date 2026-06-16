package eu.solven.matmul.f2.sat;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.concurrent.TimeUnit;

/**
 * Bridge to the {@code kissat} SAT solver (Biere et al. — winner of multiple
 * SAT Competitions). Production-grade CDCL solver; routinely orders of
 * magnitude faster than SAT4J on large structured instances.
 *
 * <h3>Install on macOS</h3>
 * <pre>brew install kissat</pre>
 *
 * <h3>Exit-code contract</h3>
 * <ul>
 *   <li>10 = SATISFIABLE (model on stdout in {@code v} lines)</li>
 *   <li>20 = UNSATISFIABLE</li>
 *   <li>other = error / timeout</li>
 * </ul>
 */
public final class KissatSolver {

	private static final String BINARY = System.getProperty("kissat.binary", "kissat");

	private KissatSolver() {}

	/** Returns true iff a {@code kissat} binary is reachable on {@code $PATH}. */
	public static boolean isAvailable() {
		try {
			Process p = new ProcessBuilder(BINARY, "--version").redirectErrorStream(true).start();
			boolean done = p.waitFor(2, TimeUnit.SECONDS);
			if (!done) {
				p.destroyForcibly();
				return false;
			}
			return p.exitValue() == 0;
		} catch (IOException | InterruptedException e) {
			return false;
		}
	}

	/**
	 * Solve via kissat with an explicit wall-time budget. Passes {@code --time=N}
	 * (real-time seconds) to kissat; if the solver exhausts its budget, kissat
	 * exits with a status indicating UNKNOWN and this method returns
	 * {@link SolverOutcome#timeout()}.
	 *
	 * @param timeoutSec wall-clock seconds; pass {@code 0} for no limit.
	 */
	public static SolverOutcome solveBounded(int varCount, List<int[]> clauses, long timeoutSec) {
		File inputFile;
		try {
			inputFile = File.createTempFile("kissat-in-", ".cnf");
			Cnf.writeDimacs(inputFile, varCount, clauses);
		} catch (IOException e) {
			return SolverOutcome.error("Failed to write CNF temp file: " + e.getMessage());
		}

		try {
			List<String> cmd = new ArrayList<>();
			cmd.add(BINARY);
			cmd.add("-q");
			if (timeoutSec > 0) cmd.add("--time=" + timeoutSec);
			cmd.add(inputFile.getAbsolutePath());

			ProcessBuilder pb = new ProcessBuilder(cmd);
			pb.redirectErrorStream(false);
			Process p = pb.start();

			List<String> lines = new ArrayList<>();
			try (BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()))) {
				String line;
				while ((line = br.readLine()) != null) lines.add(line);
			}
			int exit = p.waitFor();

			if (exit == 20) return SolverOutcome.unsat();
			if (exit == 0 || exit == 1) {
				// kissat returns 0 or 1 when it gives up (UNKNOWN). Either time
				// limit or other resource exhaustion. Treat as TIMEOUT.
				return SolverOutcome.timeout();
			}
			if (exit != 10) {
				return SolverOutcome.error("kissat exit code " + exit
						+ "; stdout: " + String.join("\\n", lines));
			}

			boolean[] assignment = new boolean[varCount];
			for (String line : lines) {
				if (line.startsWith("v ")) {
					String[] tokens = line.substring(2).trim().split("\\s+");
					for (String tok : tokens) {
						if (tok.isEmpty()) continue;
						int lit = Integer.parseInt(tok);
						if (lit == 0) break;
						int abs = Math.abs(lit);
						if (abs >= 1 && abs <= varCount) {
							assignment[abs - 1] = lit > 0;
						}
					}
				}
			}
			return SolverOutcome.sat(assignment);
		} catch (IOException | InterruptedException e) {
			Thread.currentThread().interrupt();
			return SolverOutcome.error("kissat invocation failed: " + e.getMessage());
		} finally {
			inputFile.delete();
		}
	}

	/**
	 * Solve via kissat. Writes CNF to a temp file, invokes kissat, parses the
	 * model from {@code v} lines on stdout. No wall-time budget — for benchmarks
	 * with explicit timeouts use {@link #solveBounded}.
	 *
	 * @return empty for UNSAT, the satisfying assignment for SAT.
	 * @throws RuntimeException on solver error (exit code not 10 or 20).
	 */
	public static Optional<boolean[]> solve(int varCount, List<int[]> clauses) {
		File inputFile;
		try {
			inputFile = File.createTempFile("kissat-in-", ".cnf");
			Cnf.writeDimacs(inputFile, varCount, clauses);
		} catch (IOException e) {
			throw new RuntimeException("Failed to write CNF temp file", e);
		}

		try {
			ProcessBuilder pb = new ProcessBuilder(BINARY, "-q", inputFile.getAbsolutePath());
			pb.redirectErrorStream(false); // ignore stderr verbosity; stdout has the model
			Process p = pb.start();

			List<String> lines = new ArrayList<>();
			try (BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()))) {
				String line;
				while ((line = br.readLine()) != null) lines.add(line);
			}
			int exit = p.waitFor();

			if (exit == 20) return Optional.empty(); // UNSAT
			if (exit != 10) {
				throw new RuntimeException("kissat exit code " + exit + "; stdout: "
						+ String.join("\\n", lines));
			}

			boolean[] assignment = new boolean[varCount];
			for (String line : lines) {
				if (line.startsWith("v ")) {
					String[] tokens = line.substring(2).trim().split("\\s+");
					for (String tok : tokens) {
						if (tok.isEmpty()) continue;
						int lit = Integer.parseInt(tok);
						if (lit == 0) break;
						int abs = Math.abs(lit);
						if (abs >= 1 && abs <= varCount) {
							assignment[abs - 1] = lit > 0;
						}
					}
				}
			}
			return Optional.of(assignment);
		} catch (IOException | InterruptedException e) {
			Thread.currentThread().interrupt();
			throw new RuntimeException("kissat invocation failed", e);
		} finally {
			inputFile.delete();
		}
	}
}
