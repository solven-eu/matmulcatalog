package eu.solven.matmul.docs.migrate;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.MatrixJsonFormatter;
import eu.solven.matmul.catalog.SchemeIO;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.node.ObjectNode;

/**
 * Stamp the {@code "hash"} field onto every MATERIALISED scheme file (explicit
 * factor matrices, not a stub) that carries neither a JSON {@code hash} nor a
 * filename {@code -hash7} token — i.e. old-convention imports like
 * {@code bud-bases/section7/fmm-lille_5x7x7_r176_a3315.json}.
 *
 * <p>Why: {@code RecursiveMaterialiser.durableLeafRef} refuses to pin a lineage
 * leaf on a hash-less building block (a bare shape ref would be a cited bound,
 * not a precise scheme). So a hash-less bud-base import can WIN a serendipitous
 * search and then fail at persist time — the ⟨20,28,28⟩=8434 fmm-lille base hit
 * exactly this. Sibling of {@link StampStubHashes} (which handles lineage-only
 * stubs via replay); here the content is on disk, so we hash the parsed scheme
 * directly — no replay, no rank gate beyond declared-vs-parsed equality.</p>
 *
 * <p>Default DRY-RUN. {@code --execute} writes (via
 * {@link MatrixJsonFormatter#format}, the canonical scheme-JSON writer).</p>
 *
 * <pre>mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.StampImportHashes [-Dexec.args=--execute]</pre>
 */
public final class StampImportHashes {
	private StampImportHashes() {}

	private static final Pattern FILE_HASH = Pattern.compile("-([0-9a-f]{4,})\\.json$");

	public static void main(String[] args) throws Exception {
		boolean execute = List.of(args).contains("--execute");
		Path root = Path.of("src/main/resources/schemes");

		List<Path> targets;
		try (var s = Files.walk(root)) {
			targets = s.filter(p -> p.toString().endsWith(".json"))
					.filter(StampImportHashes::isUnstampedMaterialised)
					.sorted()
					.collect(Collectors.toList());
		}
		System.out.println("Found " + targets.size() + " hash-less materialised files (mode="
				+ (execute ? "EXECUTE" : "DRY-RUN") + ")…");

		int stamped = 0, mismatched = 0, errors = 0;
		for (Path f : targets) {
			try {
				JsonNode parsed = SchemeIO.parseJson(f.toFile());
				NonCubicBilinearAlgorithm alg = SchemeIO.read(f.toFile());
				int declaredRank = parsed.get("m").asInt();
				if (alg.r != declaredRank) {
					mismatched++;
					System.out.println("[RANK-MISMATCH] " + f.getFileName() + " parsed r=" + alg.r
							+ " vs declared m=" + declaredRank + " — NOT stamped, review by hand");
					continue;
				}
				String hash = SchemeIO.contentHash(alg);
				System.out.println((execute ? "[write] " : "[plan]  ") + f.getFileName()
						+ " -> hash=" + hash.substring(0, 7) + "…");
				if (execute) {
					Files.writeString(f, MatrixJsonFormatter.format(
							withHashAfterRank((ObjectNode) parsed, hash)));
				}
				stamped++;
			} catch (Exception e) {
				errors++;
				System.out.println("[skip]  " + f.getFileName() + " unreadable: "
						+ String.valueOf(e).replaceAll("\\s+", " "));
			}
		}

		System.out.println("\n=== " + (execute ? "STAMPED" : "PLAN") + " ===");
		System.out.println("stamped:             " + stamped);
		System.out.println("mismatched (review): " + mismatched);
		System.out.println("errors (skipped):    " + errors);
		if (!execute) System.out.println("\n(DRY-RUN — pass --execute to write)");
	}

	/** Materialised (has factors), no JSON hash, no filename hash token. */
	private static boolean isUnstampedMaterialised(Path p) {
		if (FILE_HASH.matcher(p.getFileName().toString()).find()) return false;
		try {
			JsonNode n = SchemeIO.parseJson(p.toFile());
			return (n.has("u") || n.has("u_sparse")) && !n.has("hash");
		} catch (Exception e) {
			return false;
		}
	}

	private static ObjectNode withHashAfterRank(ObjectNode obj, String hash) {
		ObjectNode out = obj.objectNode();
		boolean inserted = false;
		for (var it = obj.properties().iterator(); it.hasNext();) {
			var e = it.next();
			out.set(e.getKey(), e.getValue());
			if (!inserted && "m".equals(e.getKey())) {
				out.put("hash", hash);
				inserted = true;
			}
		}
		if (!inserted) {
			out.put("hash", hash);
		}
		return out;
	}
}
