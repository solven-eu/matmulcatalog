package eu.solven.matmul.catalog;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import eu.solven.matmul.NonBilinearAlgorithm;
import eu.solven.matmul.Verifier;

/**
 * Shared writer for non-bilinear scheme JSON files (Waksman 1970,
 * Rosowski 2019 Algorithm 1, Rosowski 2019 Thm 2/3, …). Centralises:
 *
 * <ul>
 *   <li>The sparse {@code [[row, value], …]} per-product encoding read
 *       back by {@link SchemeIO#readNonBilinear}.</li>
 *   <li>A {@link #formatCoef coefficient formatter} that preserves
 *       half-integer values (Waksman's W carries ±0.5) — naive
 *       {@code (int) v} casts would silently truncate those.</li>
 *   <li>Counting additions via {@link Verifier#additionCount(NonBilinearAlgorithm)}
 *       so the {@code "additions"} field is consistent across emitters.</li>
 * </ul>
 *
 * <p>Materialisers supply the algorithm + the {@link Metadata}; the
 * filename and target directory are caller-controlled (since
 * filename conventions vary slightly across sources).</p>
 */
public final class SparseNonBilinearWriter {

	private SparseNonBilinearWriter() {}

	/**
	 * Attribution / provenance bundle written into each scheme JSON.
	 * Captured as a record so callers can construct it once per
	 * materialiser run and reuse across many schemes.
	 */
	public record Metadata(
			String field,
			boolean commutative,
			String source,
			String attributionForRank,
			boolean discovery,
			int year,
			String reference) {}

	/** Write {@code alg} to {@code file} with attribution from {@code meta}. */
	public static void write(NonBilinearAlgorithm alg, File file, Metadata meta) throws IOException {
		write(alg, file, meta, null);
	}

	/**
	 * Write {@code alg} to {@code file} with attribution from {@code meta} and
	 * optional {@code lineage} tree. Pass {@code null} for {@code lineage} on
	 * generator-only schemes whose lineage is just "call the constructor"; the
	 * resulting file is then stub-incompatible (no enough info to regenerate)
	 * but small enough that it doesn't matter.
	 */
	public static void write(NonBilinearAlgorithm alg, File file, Metadata meta,
			Lineage.Node lineage) throws IOException {
		int adds = Verifier.additionCount(alg);
		try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(file)))) {
			pw.println("{");
			pw.println("  \"n\": [" + alg.n + ", " + alg.m + ", " + alg.p + "],");
			pw.println("  \"scheme_type\": \"non_bilinear\",");
			pw.println("  \"field\": \"" + esc(meta.field()) + "\",");
			if (meta.commutative()) pw.println("  \"commutative\": true,");
			pw.println("  \"rank\": " + alg.r + ",");
			pw.println("  \"additions\": " + adds + ",");
			emitSparse(pw, "u_a", alg.Ua, alg.r);
			emitSparse(pw, "u_b", alg.Ub, alg.r);
			emitSparse(pw, "v_a", alg.Va, alg.r);
			emitSparse(pw, "v_b", alg.Vb, alg.r);
			emitSparse(pw, "w", alg.W, alg.r);
			pw.println("  \"source\": \"" + esc(meta.source()) + "\",");
			pw.println("  \"year\": " + meta.year() + ",");
			pw.println("  \"reference\": \"" + esc(meta.reference()) + "\"");
			if (lineage != null) {
				pw.println(",");
				pw.println("  \"lineage\": " + Lineage.toJson(lineage));
			} else {
				pw.println();
			}
			pw.println("}");
		}
	}

	/** Compute the addition count without writing — same metric used by {@link #write}. */
	public static int additionCount(NonBilinearAlgorithm alg) {
		return Verifier.additionCount(alg);
	}

	private static void emitSparse(PrintWriter pw, String key, double[][] M, int r) {
		pw.print("  \"" + key + "\": [");
		for (int k = 0; k < r; k++) {
			if (k > 0) pw.print(", ");
			pw.print("[");
			boolean first = true;
			for (int row = 0; row < M.length; row++) {
				double v = M[row][k];
				if (v == 0.0) continue;
				if (!first) pw.print(", ");
				first = false;
				pw.print("[" + row + ", " + formatCoef(v) + "]");
			}
			pw.print("]");
		}
		pw.println("],");
	}

	/**
	 * Format a coefficient as JSON: integer-valued doubles emit as
	 * integers (compact), fractional values fall back to
	 * {@link Double#toString}. Preserves half-integers like the {@code ±0.5}
	 * coefficients in Waksman's output combination.
	 */
	private static String formatCoef(double v) {
		if (v == Math.rint(v) && !Double.isInfinite(v) && Math.abs(v) < 1e15) {
			return Long.toString((long) v);
		}
		return Double.toString(v);
	}

	private static String esc(String s) {
		return s.replace("\\", "\\\\").replace("\"", "\\\"");
	}
}
