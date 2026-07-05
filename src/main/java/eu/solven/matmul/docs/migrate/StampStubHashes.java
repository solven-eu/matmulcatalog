package eu.solven.matmul.docs.migrate;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.MatrixJsonFormatter;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.search.LineageReplayer;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.node.ObjectNode;

/**
 * Stamp the {@code "hash"} field onto every lineage-only stub that lacks one, by
 * replaying the stub and hashing the materialised content
 * ({@link SchemeIO#contentHash}).
 *
 * <p>Why this matters: pinned {@code shape@hash} lineage refs resolve against a
 * stub ONLY via its stamped {@code "hash"} (filenames are cosmetic — nothing reads
 * them). An unstamped stub therefore makes every inbound pin dangle, and since
 * 2026-07-05 {@link LineageReplayer} THROWS on a dangling pin instead of silently
 * substituting the shape-best sibling (see references/PURGE_REFCOUNT_POLICY.md).
 * Canonical case: {@code 31x31x31-r14878-dis09_Q-6ecc4a5.json} had no stamp, so
 * {@code 30x31x31-r14573-derived}'s pin {@code 31x31x31@6ecc4a5} dangled — and the
 * old fallback replayed it from the WRONG parent (the rank-best ⟨31,31,31⟩
 * sibling, r=14519) rather than the DIS09 cube its lineage recorded.</p>
 *
 * <p>Safety gate per stub: the replayed rank must equal the declared rank
 * ({@code m}) or nothing is stamped. A filename {@code -hash7} token that does
 * NOT prefix the computed hash is EXPECTED for rational-bearing stubs named
 * before the exact-token rehash ({@link RehashRationalsComplex} era — the whole
 * dis09_Q cube family): the stamp is content truth, the stale name is logged as
 * {@code [STALE-NAME]} and should be renamed to the new token afterwards (safe:
 * filenames are pure labels — but {@code GenerateCatalogManifest.hash7Of} reads
 * the filename, so leaving the divergence undercounts {@code used_as_base}).
 * Stubs whose lineage cannot replay (e.g. {@code RosowskiTheorem2} parametric
 * leaves not yet wired into the replayer) are logged and skipped. The lookup is
 * chosen per stub from its {@code fields[]} — F₂/F₃-only stubs replay against
 * the matching modular lookup, everything else against the broadest
 * ({@code C}).</p>
 *
 * <p>Default DRY-RUN. {@code --execute} writes. Re-runnable; re-run
 * {@code AuditLineageRefs} afterwards to confirm PINNED_DANGLING dropped.</p>
 *
 * <pre>mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.StampStubHashes [-Dexec.args=--execute]</pre>
 */
public final class StampStubHashes {
	private StampStubHashes() {}

	private static final Pattern FILE_HASH7 = Pattern.compile("-([0-9a-f]{7})\\.json$");

	public static void main(String[] args) throws Exception {
		boolean execute = List.of(args).contains("--execute");
		Path root = Path.of("src/main/resources/schemes");

		List<Path> stubs;
		try (var s = Files.walk(root)) {
			stubs = s.filter(p -> p.toString().endsWith(".json"))
					.filter(StampStubHashes::isUnstampedStub)
					.sorted()
					.collect(Collectors.toList());
		}
		System.out.println("Found " + stubs.size() + " unstamped stubs (mode="
				+ (execute ? "EXECUTE" : "DRY-RUN") + ")…");

		FieldAwareLookup lookupC = new FieldAwareLookup("C");
		FieldAwareLookup lookupF2 = new FieldAwareLookup("F2");
		FieldAwareLookup lookupF3 = new FieldAwareLookup("F3");

		AtomicInteger stamped = new AtomicInteger(), mismatched = new AtomicInteger();
		AtomicInteger unreplayable = new AtomicInteger();
		long start = System.currentTimeMillis();
		int done = 0;
		for (Path f : stubs) {
			try {
				JsonNode parsed = SchemeIO.parseJson(f.toFile());
				FieldAwareLookup lookup = lookupFor(parsed, lookupC, lookupF2, lookupF3);
				NonCubicBilinearAlgorithm alg =
						LineageReplayer.withDefaultPool(lookup).replayFromFile(f.toFile());

				int declaredRank = parsed.get("m").asInt();
				String hash = SchemeIO.contentHash(alg);
				Matcher fh = FILE_HASH7.matcher(f.getFileName().toString());
				String fileHash7 = fh.find() ? fh.group(1) : null;

				if (alg.r != declaredRank) {
					mismatched.incrementAndGet();
					System.out.println("[RANK-MISMATCH] " + f.getFileName() + " replay r=" + alg.r
							+ " vs declared m=" + declaredRank + " — NOT stamped, review by hand");
					continue;
				}
				boolean staleName = fileHash7 != null && !hash.startsWith(fileHash7);
				System.out.println((execute ? "[write] " : "[plan]  ") + f.getFileName()
						+ " -> hash=" + hash.substring(0, 7) + "…"
						+ (staleName ? "  [STALE-NAME: filename token " + fileHash7
								+ " predates the rational rehash — rename to the new token]" : ""));
				if (execute) {
					Files.writeString(f, MatrixJsonFormatter.format(withHashAfterRank(
							(ObjectNode) parsed, hash)));
				}
				stamped.incrementAndGet();
			} catch (Exception | StackOverflowError e) {
				unreplayable.incrementAndGet();
				System.out.println("[skip]  " + f.getFileName() + " un-replayable: "
						+ String.valueOf(e).replaceAll("\\s+", " "));
			}
			done++;
			if (done % 5 == 0 || done == stubs.size()) {
				System.out.println("[progress] " + done + "/" + stubs.size() + " stubs, "
						+ (System.currentTimeMillis() - start) + "ms elapsed");
			}
		}

		System.out.println("\n=== " + (execute ? "STAMPED" : "PLAN") + " ===");
		System.out.println("stamped:                 " + stamped.get());
		System.out.println("mismatched (review):     " + mismatched.get());
		System.out.println("un-replayable (skipped): " + unreplayable.get());
		if (!execute) System.out.println("\n(DRY-RUN — pass --execute to write)");
	}

	private static boolean isUnstampedStub(Path p) {
		try {
			JsonNode n = SchemeIO.parseJson(p.toFile());
			return n.has("lineage") && !n.has("u") && !n.has("u_sparse") && !n.has("hash");
		} catch (Exception e) {
			return false;
		}
	}

	/** F₂/F₃-only stubs replay against the matching modular lookup; everything
	 *  else (char-0, or no fields[] at all) against the broadest ({@code C}). */
	private static FieldAwareLookup lookupFor(JsonNode parsed,
			FieldAwareLookup c, FieldAwareLookup f2, FieldAwareLookup f3) {
		Set<String> tags = Set.copyOf(SchemeIO.fieldTags(parsed));
		boolean char0 = tags.contains("Z") || tags.contains("Q") || tags.contains("R")
				|| tags.contains("C");
		if (!char0 && tags.contains("F2")) return f2;
		if (!char0 && tags.contains("F3")) return f3;
		return c;
	}

	/** Rebuild the node with {@code hash} inserted right after the rank key
	 *  {@code m} — the slot stamped stubs already use — so canonical formatting
	 *  produces the same key order as sibling files. */
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
