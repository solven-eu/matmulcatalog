package eu.solven.matmul.f2.sat;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

/**
 * DIMACS-CNF read/write helpers for interoperating with external SAT tooling
 * (BreakID, Kissat, CryptoMiniSat, etc.).
 *
 * Format:
 * <pre>
 *   c optional comment lines
 *   p cnf &lt;numVars&gt; &lt;numClauses&gt;
 *   &lt;lit_1&gt; &lt;lit_2&gt; ... &lt;lit_k&gt; 0
 *   ...
 * </pre>
 */
public final class Cnf {

	private Cnf() {}

	public static void writeDimacs(File file, int varCount, List<int[]> clauses) throws IOException {
		try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(file)))) {
			pw.printf("p cnf %d %d%n", varCount, clauses.size());
			for (int[] clause : clauses) {
				StringBuilder sb = new StringBuilder();
				for (int lit : clause) {
					sb.append(lit).append(' ');
				}
				sb.append('0');
				pw.println(sb);
			}
		}
	}

	public static ReadResult readDimacs(File file) throws IOException {
		int varCount = 0;
		List<int[]> clauses = new ArrayList<>();
		try (BufferedReader br = new BufferedReader(new FileReader(file))) {
			String line;
			while ((line = br.readLine()) != null) {
				line = line.trim();
				if (line.isEmpty() || line.startsWith("c")) continue;
				if (line.startsWith("p ")) {
					String[] parts = line.split("\\s+");
					if (parts.length >= 4) {
						varCount = Integer.parseInt(parts[2]);
					}
					continue;
				}
				String[] parts = line.split("\\s+");
				int[] tmp = new int[parts.length];
				int n = 0;
				for (String p : parts) {
					if (p.isEmpty()) continue;
					int lit = Integer.parseInt(p);
					if (lit == 0) break;
					tmp[n++] = lit;
				}
				if (n > 0) {
					int[] clause = new int[n];
					System.arraycopy(tmp, 0, clause, 0, n);
					clauses.add(clause);
				}
			}
		}
		return new ReadResult(varCount, clauses);
	}

	public static final class ReadResult {
		public final int varCount;
		public final List<int[]> clauses;

		public ReadResult(int varCount, List<int[]> clauses) {
			this.varCount = varCount;
			this.clauses = clauses;
		}
	}
}
