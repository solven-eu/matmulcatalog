package eu.solven.matmul.f2.sat;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Bridge to the BreakID SAT-preprocessor (Devriendt et al. 2016 —
 * "Improved Static Symmetry Breaking for SAT"). Takes a CNF formula and
 * returns the same formula augmented with automatically-detected symmetry-
 * breaking clauses.
 *
 * <h3>How it works</h3>
 * <ol>
 *   <li>Build the variable-incidence graph of the input CNF.</li>
 *   <li>Compute the graph automorphism group (via saucy or bliss).</li>
 *   <li>Derive lex-leader symmetry-breaking clauses from the generators.</li>
 *   <li>Emit the augmented CNF.</li>
 * </ol>
 *
 * <h3>Install</h3>
 * <pre>
 *   git clone https://bitbucket.org/krr/breakid.git
 *   cd breakid && mkdir build && cd build
 *   cmake .. && make && sudo make install
 *   # the binary should land at /usr/local/bin/breakid (or similar)
 * </pre>
 *
 * <p>If the binary isn't on {@code $PATH} when {@link #preprocess} is called,
 * the method falls back to returning the input unchanged. {@link #isAvailable}
 * lets callers check up-front whether BreakID can actually be invoked.</p>
 */
public final class BreakIdBridge {

	/**
	 * Override with {@code -Dbreakid.binary=/path/to/BreakID}. The default
	 * tries {@code BreakID} (the name of the binary produced by the upstream
	 * {@code make} build) and falls back to lowercase {@code breakid} if not
	 * found.
	 */
	private static final String BINARY = System.getProperty("breakid.binary",
			pickBinaryName());

	private static String pickBinaryName() {
		for (String candidate : new String[] { "BreakID", "breakid" }) {
			try {
				Process p = new ProcessBuilder(candidate, "-h")
						.redirectErrorStream(true).start();
				if (p.waitFor(2, java.util.concurrent.TimeUnit.SECONDS)) {
					return candidate;
				}
				p.destroyForcibly();
			} catch (IOException | InterruptedException ignored) {
				// try the next one
			}
		}
		return "BreakID"; // best guess; isAvailable() will return false anyway
	}

	private BreakIdBridge() {}

	/** Returns true iff a {@code BreakID} (or equivalent) binary is reachable on {@code $PATH}. */
	public static boolean isAvailable() {
		try {
			Process p = new ProcessBuilder(BINARY, "-h").redirectErrorStream(true).start();
			boolean done = p.waitFor(2, java.util.concurrent.TimeUnit.SECONDS);
			if (!done) {
				p.destroyForcibly();
				return false;
			}
			return true;
		} catch (IOException | InterruptedException e) {
			return false;
		}
	}

	/**
	 * Run BreakID on the input CNF. Returns the augmented CNF
	 * ({@code varCount} typically grows due to lex-leader auxiliary variables;
	 * the original variables 1..N keep their numbering).
	 *
	 * @throws IOException if BreakID isn't on PATH or fails. Callers should
	 *                     either check {@link #isAvailable} first or fall back
	 *                     to non-preprocessed solving on exception.
	 */
	/**
	 * BreakID's contract (per its README): reads DIMACS from a file given via
	 * {@code -f}, writes augmented DIMACS to <b>stdout</b>. Verbosity goes to
	 * stderr. We collect stdout into a temp file and pass to {@link Cnf#readDimacs}.
	 */
	public static Cnf.ReadResult preprocess(int varCount, List<int[]> clauses) throws IOException {
		File inputFile = File.createTempFile("breakid-in-", ".cnf");
		File outputFile = File.createTempFile("breakid-out-", ".cnf");
		try {
			Cnf.writeDimacs(inputFile, varCount, clauses);

			ProcessBuilder pb = new ProcessBuilder(
					BINARY,
					"-f", inputFile.getAbsolutePath(),
					"-v", "0",
					"-t", "60");
			pb.redirectOutput(outputFile);
			pb.redirectErrorStream(false); // keep stderr separate (verbosity)
			Process p = pb.start();
			try {
				int exit = p.waitFor();
				if (exit != 0) {
					throw new IOException("BreakID exited with code " + exit);
				}
			} catch (InterruptedException e) {
				Thread.currentThread().interrupt();
				throw new IOException("Interrupted while waiting for BreakID", e);
			}
			return Cnf.readDimacs(outputFile);
		} finally {
			inputFile.delete();
			outputFile.delete();
		}
	}
}
