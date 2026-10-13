package eu.solven.matmul.docs.migrate;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

import eu.solven.matmul.catalog.SchemeIO;

import lombok.extern.slf4j.Slf4j;
import tools.jackson.databind.JsonNode;

/**
 * Re-apply the provenance of a contribution's {@code --meta} file to the schemes it
 * was imported with — the audit step that turns {@code discovery: "TBD"} into a
 * verdict once the claim has been checked against the catalog's state at the release
 * date (CLAUDE.md, "Distinguish discoveries from re-discoveries").
 *
 * <p>Files are selected by CONTENT ({@code source} equal to the meta's
 * {@code source}; never by filename). Only the provenance keys are rewritten —
 * {@code discovery}, {@code attribution_for_rank}, {@code discovery_note},
 * {@code notes}, {@code source_url}, {@code year} — through the canonical writer
 * ({@link SchemeIO#updateFields}); matrices, hash and lineage are untouched.</p>
 *
 * <p>First use: issue #12 — the 27 Cloudroom Math 2026 schemes, imported with
 * {@code discovery: "TBD"} on 2026-10-11, audited against the dataset's own
 * {@code verification/baselines.md} (every source, 2026-10-08) and the issue's table
 * (this catalog at its 2026-10-05 data commit): all 27 below every published rank →
 * {@code true}.</p>
 *
 * <pre>mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.RestampContributedProvenance \
 *     -Dexec.args="--meta=references/contributions/cloudroom-math-2026.meta.json [--execute]"</pre>
 */
@Slf4j
public final class RestampContributedProvenance {
	private RestampContributedProvenance() {}

	private static final Path SCHEMES_ROOT = Path.of("src/main/resources/schemes");
	/** The provenance keys a meta file may carry, in the stamped order. */
	static final List<String> KEYS = List.of("year", "source_url", "discovery", "attribution_for_rank",
			"discovery_note", "notes");

	public static void main(String[] args) throws Exception {
		String metaPath = null;
		boolean execute = false;
		for (String a : args) {
			if (a.startsWith("--meta=")) metaPath = a.substring("--meta=".length());
			else if ("--execute".equals(a)) execute = true;
		}
		if (metaPath == null) {
			throw new IllegalArgumentException("usage: RestampContributedProvenance --meta=<file.json> [--execute]");
		}
		JsonNode meta = SchemeIO.parseJson(new File(metaPath));
		String source = meta.path("source").asString();
		if (source.isBlank()) throw new IllegalArgumentException(metaPath + ": no \"source\"");
		Map<String, Object> stamp = provenance(meta);
		log.info("source \"{}\" → re-stamping {}", source, stamp.keySet());

		List<Path> files;
		try (Stream<Path> w = Files.walk(SCHEMES_ROOT)) {
			files = w.filter(p -> p.getFileName().toString().endsWith(".json")).toList();
		}
		int matched = 0, changed = 0;
		List<String> report = new ArrayList<>();
		for (Path p : files) {
			JsonNode root;
			try {
				root = SchemeIO.parseJson(p.toFile());
			} catch (Exception e) {
				continue;
			}
			JsonNode src = root.get("source");
			if (src == null || !src.isTextual() || !source.equals(src.asString())) continue;
			matched++;
			boolean differs = false;
			for (Map.Entry<String, Object> e : stamp.entrySet()) {
				JsonNode current = root.get(e.getKey());
				if (current == null || !String.valueOf(unwrap(current)).equals(String.valueOf(e.getValue()))) {
					differs = true;
				}
			}
			String shape = root.path("n").toString() + " r=" + root.path("m").asInt();
			if (!differs) {
				report.add("  = " + p.getFileName() + " " + shape + " already stamped");
				continue;
			}
			changed++;
			report.add((execute ? "  ✓ " : "  ~ ") + p.getFileName() + " " + shape + "  discovery "
					+ root.path("discovery") + " → " + stamp.get("discovery"));
			if (execute) {
				// updateFields never clobbers an existing key: remove, then re-add (same keys, new values).
				SchemeIO.updateFields(p.toFile(), stamp, stamp.keySet(), true);
			}
		}
		report.forEach(log::info);
		log.info("RestampContributedProvenance: {} scheme(s) with source \"{}\", {} {} ({})", matched, source, changed,
				execute ? "re-stamped" : "would change", execute ? "--execute" : "DRY-RUN");
	}

	/** The meta's provenance keys as the values {@code ImportContributedSchemes} stamps. */
	static Map<String, Object> provenance(JsonNode meta) {
		Map<String, Object> out = new LinkedHashMap<>();
		for (String key : KEYS) {
			JsonNode v = meta.get(key);
			if (v == null || v.isNull()) continue;
			if (key.equals("discovery")) {
				String d = v.asString();
				out.put(key, "true".equals(d) ? Boolean.TRUE : "false".equals(d) ? Boolean.FALSE : d);
			} else if (v.isNumber()) {
				out.put(key, v.asInt());
			} else {
				out.put(key, v.asString());
			}
		}
		return out;
	}

	private static Object unwrap(JsonNode node) {
		if (node.isBoolean()) return node.booleanValue();
		if (node.isNumber()) return node.asInt();
		return node.isTextual() ? node.asString() : node.toString();
	}
}
