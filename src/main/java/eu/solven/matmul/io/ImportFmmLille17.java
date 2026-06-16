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
 * One-shot driver that parses the FMM-Lille {@code 17x17x17_raw.mpl}
 * scheme via {@link MapleSchemeParser} and writes the catalog JSON to
 * {@code src/main/resources/schemes/known/section17/fmm_lille_2025_17x17x17_r{rank}_a{adds}_Q.json}.
 *
 * <p>The output is the standard sparse-format JSON emitted by
 * {@link SchemeIO#toJsonSparse} with extra header metadata injected
 * (field, discovery, attribution, source, year, notes,
 * derivation_task). Importing this scheme is justified per the
 * "FMM-Lille construction not yet derivable" exception in
 * {@code CLAUDE.md} — the derivation side stays as a follow-up
 * task referenced by {@code "derivation_task"}.</p>
 */
public final class ImportFmmLille17 {

	private ImportFmmLille17() {}

	public static void main(String[] args) throws IOException {
		File raw = new File("references/fmm-lille/17x17x17/17x17x17_raw.mpl");
		if (!raw.isFile()) {
			throw new IOException("missing FMM-Lille raw scheme at: " + raw.getAbsolutePath());
		}
		NonCubicBilinearAlgorithm alg = MapleSchemeParser.parseRawFmmLille(raw, 17, 17, 17);
		int adds = Verifier.additionCount(alg);
		String name = SchemeIO.canonicalName(alg, "fmm_lille_2025");
		File out = new File("src/main/resources/schemes/known/section17/" + name);
		writeWithMetadata(alg, out, adds);
		System.out.println("[ImportFmmLille17] wrote " + out
				+ " (r=" + alg.r + ", additions=" + adds + ")");
	}

	/**
	 * Writes the scheme using {@link SchemeIO#toJsonSparse} and splices the
	 * extra metadata fields after the opening brace. Keeps SchemeIO's reader
	 * path unchanged (it ignores unknown top-level fields).
	 */
	public static void writeWithMetadata(NonCubicBilinearAlgorithm alg, File file, int adds)
			throws IOException {
		String body = SchemeIO.toJsonSparse(alg);
		int openBrace = body.indexOf('{');
		if (openBrace < 0) throw new IOException("malformed scheme JSON: no '{'");
		String prefix = body.substring(0, openBrace + 1);
		String suffix = body.substring(openBrace + 1);

		// Build the metadata block (each field followed by a comma, joined
		// to the existing JSON which itself starts with "\n  \"n\":...").
		StringBuilder meta = new StringBuilder();
		meta.append("\n");
		meta.append("  \"field\": \"Q\",\n");
		meta.append("  \"commutative\": false,\n");
		meta.append("  \"rank\": ").append(alg.r).append(",\n");
		meta.append("  \"additions\": ").append(adds).append(",\n");
		meta.append("  \"source\": \"FMM-Lille\",\n");
		meta.append("  \"year\": 2025,\n");
		meta.append("  \"reference\": \"https://fmm.univ-lille.fr/17x17x17.html\",\n");
		meta.append("  \"notes\": \"Parsed from references/fmm-lille/17x17x17/17x17x17_raw.mpl. " +
				"Coefficients are Q-rational (±1, ±1/8). The published rank is 2934 " +
				"(= 2931 post-kin-row-reduction + 3 free corrections); this import keeps " +
				"the raw 2945 products since kin-row unification is not implemented yet. " +
				"Scheme still computes ⟨17,17,17⟩ exactly. See also " +
				"references/fmm-lille/17x17x17/17x17x17_LRP.mpl and " +
				"references/fmm-lille/17x17x17/17x17x17_tensor.mpl.\",\n");
		meta.append("  \"derivation_task\": \"#88\",");

		String patched = prefix + meta + suffix.replaceFirst("^\n?", "\n");
		try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(file)))) {
			pw.print(patched);
		}
	}
}
