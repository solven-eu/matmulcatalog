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
 * Bridge to the {@code cryptominisat5} SAT solver (Soos et al.). Production-
 * grade CDCL solver with **native XOR clause support** — particularly
 * relevant for our Z/2 matmul encoding, which is dominated by XOR-sum
 * constraints.
 *
 * <h3>Install on macOS</h3>
 * <pre>brew install cryptominisat</pre>
 *
 * <h3>Exit-code contract</h3>
 * <ul>
 *   <li>10 = SATISFIABLE (model on stdout as {@code v} lines)</li>
 *   <li>20 = UNSATISFIABLE</li>
 *   <li>other = error</li>
 * </ul>
 */
public final class CryptoMiniSatSolver {

	private static final String BINARY = System.getProperty("cryptominisat.binary", "cryptominisat5");

	private CryptoMiniSatSolver() {}

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

	public static Optional<boolean[]> solve(int varCount, List<int[]> clauses) {
		SolverOutcome out = solveBounded(varCount, clauses, 0);
		switch (out.state) {
		case SAT:
			return out.assignment;
		case UNSAT:
			return Optional.empty();
		default:
			throw new RuntimeException("cryptominisat5: " + out.state
					+ (out.message == null ? "" : ": " + out.message));
		}
	}

	/**
	 * Solve with an explicit wall-time budget (cryptominisat5 {@code --maxtime}).
	 * On timeout, returns {@link SolverOutcome#timeout()}.
	 */
	public static SolverOutcome solveBounded(int varCount, List<int[]> clauses, long timeoutSec) {
		File inputFile;
		try {
			inputFile = File.createTempFile("cms-in-", ".cnf");
			Cnf.writeDimacs(inputFile, varCount, clauses);
		} catch (IOException e) {
			return SolverOutcome.error("Failed to write CNF temp file: " + e.getMessage());
		}

		try {
			List<String> cmd = new ArrayList<>();
			cmd.add(BINARY);
			cmd.add("--verb=0");
			if (timeoutSec > 0) cmd.add("--maxtime=" + timeoutSec);
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
			if (exit == 15) return SolverOutcome.timeout(); // cryptominisat UNKNOWN
			if (exit != 10) {
				return SolverOutcome.error("cryptominisat5 exit code " + exit
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
			return SolverOutcome.error("cryptominisat5 invocation failed: " + e.getMessage());
		} finally {
			inputFile.delete();
		}
	}
}
