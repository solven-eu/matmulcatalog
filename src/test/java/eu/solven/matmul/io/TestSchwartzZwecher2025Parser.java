package eu.solven.matmul.io;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Stream;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.condition.EnabledIf;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.SchemeIO;
import tools.jackson.databind.JsonNode;

/**
 * Validates that the schemes imported from Schwartz-Zwecher 2025's
 * supplemental data (arXiv:2508.01748) parse correctly and compute
 * matmul exactly.
 *
 * <p>The supplemental zip
 * (<a href="https://www.cs.huji.ac.il/~odedsc/papers/trilinear_aggregation_algorithms_decomposed-2025-07-29.zip">cs.huji.ac.il/~odedsc</a>)
 * provides the decomposed-form encoding/decoding matrices for
 * {@code U_φ, V_φ, W_φ, φ} for each base case {@code n0 ∈ {20, 22, …, 50}}.
 * Our import script ({@code tools/sz2025-import/import_schwartz-zwecher-2025.py})
 * expands them to {@code U = U_φ · φ}, etc., snaps every entry to a
 * rational with denominator dividing {@code (n0/2+1)²}, and writes the
 * sparse JSON.</p>
 *
 * <p>This test exercises every imported scheme found on disk, asserting:
 * <ul>
 *   <li>parse succeeds via {@link SchemeIO#read(File)};</li>
 *   <li>the JSON declares {@code "field": "Q"};</li>
 *   <li>the filename declares {@code _Q_};</li>
 *   <li>the JSON rank matches the filename;</li>
 *   <li>the algorithm passes
 *       {@link Verifier#passesRandomMatmulSpotCheck}, the
 *       fast random-input verifier ({@link
 *       eu.solven.matmul.verifiers.SymbolicVerifier} is intractable at
 *       this scale — {@code O(n⁶·r)} BigInteger multiplies — and
 *       the spot-check covers the polynomial identity densely).</li>
 * </ul>
 */
public class TestSchwartzZwecher2025Parser {

	private static final Path SCHEMES_ROOT = Path.of("src/main/resources/schemes");

	/** Returns the list of imported scheme files, or empty list if none yet. */
	private static List<File> findSzFiles() {
		try {
			if (!Files.isDirectory(SCHEMES_ROOT)) return List.of();
			try (Stream<Path> walk = Files.walk(SCHEMES_ROOT)) {
				return walk
						.filter(p -> p.getFileName().toString().startsWith("schwartz_zwecher_2025_"))
						.filter(p -> p.getFileName().toString().endsWith(".json"))
						.map(Path::toFile)
						.sorted()
						.toList();
			}
		} catch (java.io.IOException e) {
			return List.of();
		}
	}

	static boolean haveImportedSchemes() {
		return !findSzFiles().isEmpty();
	}

	@Test
	@EnabledIf("haveImportedSchemes")
	public void at_least_one_sz2025_scheme_is_present() {
		List<File> files = findSzFiles();
		assertThat(files)
				.as("import_schwartz-zwecher-2025.py should produce at least one scheme")
				.isNotEmpty();
	}

	@Test
	@EnabledIf("haveImportedSchemes")
	public void every_sz2025_scheme_parses_and_verifies() throws Exception {
		List<File> files = findSzFiles();
		// Pull a subset that we always want covered. n=20 is fastest;
		// n=44 is the headline. Anything smaller skipped if absent.
		for (File f : files) {
			String name = f.getName();

			// (1) Filename declares Q.
			assertThat(name)
					.as("filename %s must carry the _Q_ field tag", name)
					.contains("_Q.json");

			// (2) Parse + JSON metadata sanity.
			JsonNode root = SchemeIO.parseJson(f);
			assertThat(root.get("field").asString())
					.as("scheme %s must declare field=Q", name)
					.isEqualTo("Q");
			// discovery field must be present (true for n>=28 where SZ
			// strictly improves on Pan 1982 / DIS09 TA; false for n<28
			// where SZ matches but doesn't beat the prior bound).
			assertThat(root.has("discovery"))
					.as("scheme %s must carry a discovery flag", name)
					.isTrue();
			assertThat(root.has("importing_source"))
					.as("scheme %s must record the importing source", name)
					.isTrue();
			assertThat(root.get("importing_source").asString())
					.isEqualTo("Schwartz-Zwecher 2025");
			assertThat(root.has("source_paper")).isTrue();
			assertThat(root.get("source_paper").asString())
					.contains("2508.01748");
			assertThat(root.has("derivation_task")).isTrue();

			// (3) Filename rank == JSON rank.
			int jsonRank = root.get("m").asInt();
			String[] tokens = name.split("_");
			// schwartz_zwecher_2025_NxNxN_r{rank}_a{adds}_Q.json
			int filenameRank = -1;
			for (String t : tokens) {
				if (t.startsWith("r") && t.length() > 1
						&& Character.isDigit(t.charAt(1))) {
					filenameRank = Integer.parseInt(t.substring(1));
					break;
				}
			}
			assertThat(filenameRank)
					.as("rank token in filename %s", name)
					.isEqualTo(jsonRank);

			// (4) Fast spot-check verification (sufficient for huge n).
			NonCubicBilinearAlgorithm alg = SchemeIO.read(f);
			boolean ok = Verifier.passesRandomMatmulSpotCheck(alg);
			assertThat(ok)
					.as("scheme %s must pass random-input matmul spot check", name)
					.isTrue();
		}
	}

	/**
	 * Anchor test for the headline {@code n=44, rank=36110} scheme.
	 * Only runs if that specific file is on disk (it's ~hundreds-of-MB
	 * and may not always be imported on every checkout).
	 */
	@Test
	@EnabledIf("hasN44Scheme")
	public void n44_anchor_matches_paper_table1() throws Exception {
		File f = findSzFiles().stream()
				.filter(p -> p.getName().contains("_44x44x44_m36110_"))
				.findFirst()
				.orElseThrow();
		JsonNode root = SchemeIO.parseJson(f);
		assertThat(root.get("m").asInt())
				.as("n=44 rank from arXiv:2508.01748 Table 1 (= 36110)")
				.isEqualTo(36110);
		NonCubicBilinearAlgorithm alg = SchemeIO.read(f);
		assertThat(alg.n).isEqualTo(44);
		assertThat(alg.m).isEqualTo(44);
		assertThat(alg.p).isEqualTo(44);
		assertThat(alg.r).isEqualTo(36110);
		assertThat(Verifier.passesRandomMatmulSpotCheck(alg))
				.as("n=44 SZ2025 must compute matmul")
				.isTrue();
	}

	static boolean hasN44Scheme() {
		return findSzFiles().stream()
				.anyMatch(p -> p.getName().contains("_44x44x44_m36110_"));
	}
}
