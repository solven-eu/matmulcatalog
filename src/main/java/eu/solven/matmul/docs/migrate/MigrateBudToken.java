package eu.solven.matmul.docs.migrate;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.BudParetoSelection;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.catalog.SerendipitousBudProduct;

/**
 * One-shot migration: stamp the {@code _b{bud_score}} token onto existing
 * scheme filenames (user 2026-06-06), making bud-richness a first-class,
 * at-a-glance filename axis alongside {@code _m{rank}} / {@code _a{additions}}.
 *
 * <p>Per-file bud source, matching {@link GenerateCatalogManifest}'s priority:</p>
 * <ul>
 *   <li><b>Has matrices</b> (materialised, maxDim ≤ {@code MATERIALISE_MAX_DIM}) →
 *       {@link SerendipitousBudProduct#summarise} → exact score. A materialised
 *       scheme is always checkable, so a 0 score yields {@code _b0}
 *       ("checked, no buds").</li>
 *   <li><b>Lineage-only stub</b> (no matrices) → take the score the catalog
 *       inferred for it (from {@code docs/catalog.json}); if none, leave the file
 *       untouched — absence of {@code _b} honestly means "buds not computed".</li>
 * </ul>
 *
 * <p>Idempotent: files already carrying a {@code _b\d} token are skipped.
 * The token is inserted right after {@code _m{rank}[_a{additions}]} and before
 * any field suffix / {@code .json}.</p>
 */
public final class MigrateBudToken {

	private MigrateBudToken() {}

	/** Insert point: after the rank (and optional additions) tokens. */
	private static final Pattern INSERT_AT = Pattern.compile("(_[mr]\\d+)(_a\\d+)?");
	private static final Pattern HAS_BUD_TOKEN = Pattern.compile("_b\\d");

	public static void main(String[] args) throws Exception {
		// 1) bud_score the catalog inferred per file (covers lineage-only stubs).
		Map<String, Integer> catalogBuds = loadCatalogBuds(new File("docs/catalog.json"));
		System.out.printf("loaded %d catalog bud_scores%n", catalogBuds.size());

		File root = new File("src/main/resources/schemes");
		// sectionN dirs live under known/derived/curated since the 2026-06-09 split —
		// recurse to find them at any depth (not just at the schemes root).
		File[] sections;
		try (java.util.stream.Stream<java.nio.file.Path> w = java.nio.file.Files.walk(root.toPath())) {
			sections = w.filter(java.nio.file.Files::isDirectory)
					.filter(p -> p.getFileName().toString().startsWith("section"))
					.map(java.nio.file.Path::toFile).toArray(File[]::new);
		} catch (java.io.IOException e) { sections = null; }
		if (sections == null || sections.length == 0) { System.err.println("no sections"); return; }

		long start = System.currentTimeMillis();
		int seen = 0, renamed = 0, skippedToken = 0, skippedUnknown = 0;
		for (File sec : sections) {
			File[] files = sec.listFiles((d, n) -> n.endsWith(".json"));
			if (files == null) continue;
			for (File f : files) {
				seen++;
				String name = f.getName();
				if (HAS_BUD_TOKEN.matcher(name).find()) { skippedToken++; continue; }

				Integer score = scoreFor(f, catalogBuds);
				if (score == null) { skippedUnknown++; continue; }

				Matcher m = INSERT_AT.matcher(name);
				if (!m.find()) { skippedUnknown++; continue; } // no _m/_r anchor — leave it
				String newName = name.substring(0, m.end())
						+ "_b" + score
						+ name.substring(m.end());
				Files.move(f.toPath(), f.toPath().resolveSibling(newName));
				renamed++;

				if (renamed % 1000 == 0) {
					long ms = System.currentTimeMillis() - start;
					System.out.printf("[progress] %d renamed (%d seen) %dms%n", renamed, seen, ms);
				}
			}
		}
		System.out.printf("%n[done] seen=%d renamed=%d skipped(already-tokened)=%d skipped(buds-unknown)=%d %dms%n",
				seen, renamed, skippedToken, skippedUnknown, System.currentTimeMillis() - start);
	}

	/** Bud score for a file: from its own matrices when present, else the
	 *  catalog-inferred value for lineage-only stubs, else null (unknown). */
	private static Integer scoreFor(File f, Map<String, Integer> catalogBuds) {
		NonCubicBilinearAlgorithm alg = null;
		try {
			alg = SchemeIO.read(f);
		} catch (Exception ignored) {
			// stub / non-bilinear / unreadable
		}
		if (alg != null) {
			return BudParetoSelection.budScore(SerendipitousBudProduct.summarise(alg).summary());
		}
		return catalogBuds.get(f.getName());
	}

	/** Map basename(file) → bud_score for catalog entries where buds are known. */
	private static Map<String, Integer> loadCatalogBuds(File catalog) {
		Map<String, Integer> out = new HashMap<>();
		if (!catalog.isFile()) return out;
		JsonNode root = new JsonMapper().readTree(catalog);
		JsonNode schemes = root.get("schemes");
		if (schemes == null || !schemes.isArray()) return out;
		for (JsonNode s : schemes) {
			JsonNode file = s.get("file");
			if (file == null || !s.path("has_buds").asBoolean(false)) continue;
			String base = Path.of(file.asText()).getFileName().toString();
			out.put(base, s.path("bud_score").asInt(0));
		}
		return out;
	}
}
