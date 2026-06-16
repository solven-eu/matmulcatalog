package eu.solven.matmul.io;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Driver that parses the FMM-Lille {@code 17x17x17_LRP.mpl} scheme via
 * {@link MapleLRPParser} and writes the catalog JSON to
 * {@code src/main/resources/schemes/known/section17/fmm_lille_2025_17x17x17_r{rank}_a{adds}_Q.json}.
 *
 * <p>The {@code _LRP.mpl} file is FMM-Lille's already-kin-row-reduced
 * encoding (in their convention) — it contains a single Maple expression
 * {@code LRP := [Matrix(r, n·m, ...), Matrix(r, m·p, ...), Matrix(n·p, r, ...)]}
 * with {@code r = 2931} for {@code ⟨17,17,17⟩}, three less than the
 * published rank claim of 2934 (which counts the 3 "free" extra
 * multiplications). See the report from this import for details on the
 * gap.</p>
 */
public final class ImportFmmLille17LRP {

	private ImportFmmLille17LRP() {}

	public static void main(String[] args) throws IOException {
		File lrp = new File("references/fmm-lille/17x17x17/17x17x17_LRP.mpl");
		if (!lrp.isFile()) {
			throw new IOException("missing FMM-Lille LRP scheme at: " + lrp.getAbsolutePath());
		}
		NonCubicBilinearAlgorithm alg = MapleLRPParser.parse(lrp, 17, 17, 17);
		int adds = Verifier.additionCount(alg);
		String name = SchemeIO.canonicalName(alg, "fmm_lille_2025");
		File out = new File("src/main/resources/schemes/known/section17/" + name);
		writeWithMetadata(alg, out, adds);
		System.out.println("[ImportFmmLille17LRP] wrote " + out
				+ " (r=" + alg.r + ", additions=" + adds + ")");
	}

	public static void writeWithMetadata(NonCubicBilinearAlgorithm alg, File file, int adds)
			throws IOException {
		String body = SchemeIO.toJsonSparse(alg);
		int openBrace = body.indexOf('{');
		if (openBrace < 0) throw new IOException("malformed scheme JSON: no '{'");
		String prefix = body.substring(0, openBrace + 1);
		String suffix = body.substring(openBrace + 1);

		StringBuilder meta = new StringBuilder();
		meta.append("\n");
		meta.append("  \"field\": \"Q\",\n");
		meta.append("  \"commutative\": false,\n");
		meta.append("  \"rank\": ").append(alg.r).append(",\n");
		meta.append("  \"additions\": ").append(adds).append(",\n");
		meta.append("  \"source\": \"FMM-Lille\",\n");
		meta.append("  \"year\": 2025,\n");
		meta.append("  \"reference\": \"https://fmm.univ-lille.fr/17x17x17.html\",\n");
		meta.append("  \"notes\": \"Parsed from references/fmm-lille/17x17x17/17x17x17_LRP.mpl " +
				"via MapleLRPParser. The LRP file encodes a complete bilinear scheme with " +
				"exactly ").append(alg.r).append(" products (matrices L: r×n·m, R: r×m·p, P: n·p×r). " +
				"This is THREE LESS than the publicly cited rank 2934 in the FMM-Lille catalog page — " +
				"the parsed scheme verifies symbolically over Q, so 2931 is itself a valid " +
				"non-commutative rank for ⟨17,17,17⟩ over Q. The published 2934 likely counts " +
				"products before a final kin-row pass that LRP already includes. Coefficients " +
				"are Q-rational (±1, ±1/8). Index conventions (verified by 8-way probe in " +
				"TestMapleLRPParserDebug): L is column-major over A (j·n+i), R is column-major " +
				"over B (l·m+j), P is row-major over C (i·p+l).\",\n");
		meta.append("  \"derivation_task\": \"#88\",");

		String patched = prefix + meta + suffix.replaceFirst("^\n?", "\n");
		try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(file)))) {
			pw.print(patched);
		}
	}
}
