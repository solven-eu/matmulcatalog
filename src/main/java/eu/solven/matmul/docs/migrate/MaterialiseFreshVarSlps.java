package eu.solven.matmul.docs.migrate;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.additions.SlpJson;
import eu.solven.matmul.additions.SlpJson.SchemeSlp;
import eu.solven.matmul.catalog.SchemeIO;

import lombok.extern.slf4j.Slf4j;
import tools.jackson.databind.JsonNode;

/**
 * Register the programs the catalog already stores but never surfaced: the 168 legacy
 * Perminov {@code _reduced} files carry {@code u_fresh}/{@code v_fresh}/{@code w_fresh}
 * (an addition-reduced schedule, {@code complexity.reduced}) that {@code SchemeIO.readReduced}
 * only ever expanded into matrices — the manifest showed their naive {@code additions}
 * (e.g. 122 for the 58-addition ⟨3,3,3⟩=23 re-encoding) and nothing else.
 *
 * <p>Each program is read with {@link SlpJson#fromFreshVars}, replayed against the file's
 * own matrices (a program that does not reconstruct them is reported, never stamped) and
 * attached as {@code min_additions} + object-layout {@code slp} (tier {@code bound},
 * source {@code Perminov 2023}: the aggregator who published the reduced re-encoding),
 * unless the file already holds a program at a count at most as low — the mizopl programs
 * attached by {@code ImportUpstreamSchemeRepos} stay where they are lower. Idempotent,
 * re-runnable; dry by default.</p>
 *
 * <pre>mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.MaterialiseFreshVarSlps
 *     [-Dexec.args="--execute"]</pre>
 */
@Slf4j
public final class MaterialiseFreshVarSlps {
	private MaterialiseFreshVarSlps() {}

	private static final Path SCHEMES_ROOT = Path.of("src/main/resources/schemes");
	static final String METHOD = "Perminov's addition-reduced re-encoding (u_fresh/v_fresh/w_fresh of the _reduced"
			+ " layout), replayed exactly against the file's matrices by MaterialiseFreshVarSlps";
	static final String SOURCE = "Perminov 2023";

	public static void main(String[] args) throws Exception {
		boolean execute = List.of(args).contains("--execute");
		List<Path> files;
		try (Stream<Path> w = Files.walk(SCHEMES_ROOT)) {
			files = w.filter(f -> f.getFileName().toString().endsWith(".json")).toList();
		}
		int withProgram = 0, attached = 0, kept = 0, broken = 0;
		List<String> brokenFiles = new ArrayList<>();
		long t0 = System.nanoTime();
		for (Path p : files) {
			JsonNode root;
			try {
				root = SchemeIO.parseJson(p.toFile());
			} catch (Exception e) {
				continue;
			}
			if (!SlpJson.isFreshVarLayout(root) || !SchemeIO.isReduced(root)) {
				continue;
			}
			withProgram++;
			File f = p.toFile();
			try {
				NonCubicBilinearAlgorithm alg = SchemeIO.readBilinear(f);
				SchemeSlp program = SlpJson.fromFreshVars(root);
				if (!program.reconstructs(alg)) {
					broken++;
					brokenFiles.add(f.getName());
					log.warn("[BROKEN] {}: the fresh-variable program does not replay to the expanded matrices", f.getName());
					continue;
				}
				int declared = root.path("complexity").path("reduced").asInt(-1);
				if (declared >= 0 && declared != program.additions()) {
					log.info("{}: complexity.reduced={} but the program costs {} additions — the replay is what we stamp",
							f.getName(), declared, program.additions());
				}
				if (!execute) {
					log.info("[dry-run] {} would attach {} additions", f.getName(), program.additions());
					attached++;
					continue;
				}
				if (ImportUpstreamSchemeRepos.attachProgram(f, alg, program, METHOD, SOURCE)) {
					attached++;
				} else {
					kept++;
				}
			} catch (Exception e) {
				broken++;
				brokenFiles.add(f.getName());
				log.warn("[ERR] {}: {}", f.getName(), e.toString());
			}
			if ((attached + kept + broken) % 50 == 0) {
				log.info("[progress] {} programs seen, {} attached, {} kept, {} broken, {}ms", withProgram, attached, kept,
						broken, (System.nanoTime() - t0) / 1_000_000L);
			}
		}
		log.info("MaterialiseFreshVarSlps{}: {} files with a fresh-variable program, {} attached, {} already held a"
				+ " program at most as cheap, {} broken {} ({}ms). Next: GenerateCatalogManifest.",
				execute ? "" : " [dry-run]", withProgram, attached, kept, broken, brokenFiles,
				(System.nanoTime() - t0) / 1_000_000L);
	}
}
