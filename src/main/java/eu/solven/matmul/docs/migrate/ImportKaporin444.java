package eu.solven.matmul.docs.migrate;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import eu.solven.matmul.ComplexNonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.SchemeIO;

import tools.jackson.databind.JsonNode;

/**
 * Import Kaporin's complex {@code C⟨4,4,4⟩=48} scheme (Doklady Mathematics 518,
 * 2024 — REFERENCES.md [87]). The explicit coefficients come from the author's
 * companion verification program {@code test444r48.for}
 * (references/papers/kaporin_2024_test444r48.for; the cyclic-symmetric
 * {@code r = p·q = 4·12} construction). The U/V/W factors were extracted from
 * that file (parsed + the {@code z(i,j,s)=(d_i/d_j)^s}, {@code x(i,j,it)},
 * {@code ip} permutation expanded into 48 rank-one terms) and laid out in the
 * catalog's storage convention (u/v row-major, w col-major) in
 * {@code /tmp/kaporin_raw.json}.
 *
 * <p>This driver reads that raw factorisation, <b>re-verifies it exactly</b> via
 * {@link Verifier#isExactComplex} (the same check {@code VerifyAllSchemes} runs),
 * and only then writes the canonical, content-hashed scheme JSON with full
 * metadata. Numerical scheme (complex floats, residual ~1e-15) — like
 * AlphaEvolve's C original; coefficients are NOT exact rationals.</p>
 *
 * <pre>mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.ImportKaporin444 \
 *     -Dexec.args="/tmp/kaporin_raw.json [--execute]"</pre>
 */
public final class ImportKaporin444 {
	private ImportKaporin444() {}

	private static final Path OUT_DIR = Path.of("src/main/resources/schemes/known/section4");
	private static final String SCHEME_URL =
			"https://cloud.mail.ru/public/Yfij/ErDxopqBh";   // test444r48.for (ref [8] of the paper)
	private static final String PAPER_URL =
			"https://journals.rcsi.science/2686-9543/article/view/269374";

	public static void main(String[] args) throws Exception {
		if (args.length < 1) {
			throw new IllegalArgumentException("usage: ImportKaporin444 <raw.json> [--execute]");
		}
		boolean execute = List.of(args).contains("--execute");
		File raw = new File(args[0]);

		ComplexNonCubicBilinearAlgorithm cx = SchemeIO.readComplex(raw);
		double residual = Verifier.residualComplex(cx);
		boolean exact = Verifier.isExactComplex(cx);
		System.out.printf("read ⟨%d,%d,%d⟩ r=%d  residual=%.3e  exact=%b%n",
				cx.n, cx.m, cx.p, cx.r, residual, exact);
		if (!exact) {
			throw new IllegalStateException("scheme does NOT verify as exact C matmul (residual "
					+ residual + ") — refusing to import");
		}
		if (!execute) {
			System.out.println("(DRY-RUN — verified OK; pass --execute to write the scheme)");
			return;
		}

		// Canonical body first (to obtain the content hash for the filename).
		Files.createDirectories(OUT_DIR);
		File tmp = OUT_DIR.resolve("4x4x4-r48-kaporin_2024-PENDING.json").toFile();
		SchemeIO.write(cx, tmp, "C");
		JsonNode root = SchemeIO.parseJson(tmp);
		String hash7 = SchemeIO.contentHashComplexJson(root).substring(0, 7);
		File out = OUT_DIR.resolve("4x4x4-r48-kaporin_2024-" + hash7 + ".json").toFile();
		SchemeIO.write(cx, out, "C");
		Files.deleteIfExists(tmp.toPath());

		Map<String, Object> meta = new LinkedHashMap<>();
		meta.put("fields", List.of("C"));
		meta.put("fields_not", List.of("F2", "F3", "Z", "Q", "R"));
		meta.put("commutative", false);
		meta.put("source", "Kaporin 2024");
		meta.put("year", 2024);
		meta.put("source_paper_url", PAPER_URL);
		meta.put("source_scheme_url", SCHEME_URL);
		// Discovery: an explicit, independently-constructed C⟨4,4,4⟩=48 published in
		// 2024 — predates AlphaEvolve 2025 (which is the catalog's other C=48). Mark
		// for audit against the existence claim of Li-Zhang-Ke 2023 (ref [12]).
		meta.put("discovery", "TBD");
		meta.put("verified", true);
		meta.put("notes", "Cyclic-symmetric semi-analytical Brent solution; complex floats "
				+ "(numerical, residual ~1e-15), not exact rationals. Extracted from the author's "
				+ "test444r48.for verification program. Differs from AlphaEvolve's C=48 original.");
		SchemeIO.addFields(out, meta, /* apply */ true);
		System.out.println("wrote " + out);
	}
}
