package eu.solven.matmul.docs.migrate;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.search.LineageReplayer;
import lombok.extern.slf4j.Slf4j;
import tools.jackson.databind.JsonNode;

/**
 * Re-stamp lineage-only stubs whose recorded rank ({@code m}) no longer matches what
 * their lineage replays today.
 *
 * <p>A stub's {@code m} is a snapshot taken at materialisation time. Its lineage can
 * still move afterwards: a bare {@code shape} / legacy filename-stem ref resolves to
 * the catalog-best of the day, a {@code SerendipitousProduct} re-derives its
 * decomposition against today's inner ranks, and a pinned ref to another <em>stub</em>
 * pins that stub's lineage, not its matrices. A {@code Project} of such a parent then
 * dead-code-eliminates a different number of products — the 2026-10-14 audit found two
 * projections recorded 3–4 below their replay (⟨13,13,17⟩ 1876 → 1879, ⟨20,26,30⟩
 * 8590 → 8594), i.e. two false rank claims, and 510 stubs recorded <em>above</em> their
 * replay (stale, honest-but-loose). {@code VerifySchemes --check=rank} fails on the
 * first kind; this driver repairs either kind by rewriting the stub at its replayed
 * rank (same lineage, provenance keys carried over, file renamed to the new rank/hash).
 *
 * <pre>
 *   --direction=over|all   over (default): only recorded &lt; replayed (false claims);
 *                          all: also tighten recorded &gt; replayed
 *   --execute              write (default: dry-run, report only)
 *   --schemes-root=PATH    default src/main/resources/schemes
 * </pre>
 */
@Slf4j
public final class RestampStubRanks {
	private RestampStubRanks() {}

	/**
	 * Provenance keys carried over to the rewritten stub. Everything else is either
	 * owned by the stub writer (n, m, hash, fields, scheme_type, lineage*) or a metric
	 * computed from the OLD scheme (additions, buds, has_buds, projection_margin, zt, …)
	 * that would be stale at the new rank — the enrichment stampers recompute those.
	 */
	private static final Set<String> CARRY_KEYS = Set.of("source", "year", "source_url", "discovery",
			"attribution_for_rank", "discovery_note", "notes", "commutative", "verified", "category");

	public record Drift(Path file, int recorded, int replayed) {
		public boolean overClaim() {
			return replayed > recorded;
		}
	}

	public static void main(String[] args) throws Exception {
		List<String> argv = Arrays.asList(args);
		boolean execute = argv.contains("--execute");
		boolean all = argv.stream().anyMatch(a -> a.equals("--direction=all"));
		Path root = Path.of(argv.stream().filter(a -> a.startsWith("--schemes-root=")).map(a -> a.substring(15))
				.findFirst().orElse("src/main/resources/schemes"));
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		LineageReplayer replayer = LineageReplayer.withDefaultPool(lookup);
		List<Path> stubs = stubs(root);
		log.info("RestampStubRanks: {} stubs under {} (direction={}, {})", stubs.size(), root, all ? "all" : "over",
				execute ? "EXECUTE" : "dry-run");
		List<Drift> drifts = drifts(stubs, replayer);
		List<Drift> todo = drifts.stream().filter(d -> all || d.overClaim()).collect(Collectors.toList());
		log.info("{} stub(s) drifted ({} over-claims, {} loose); {} selected", drifts.size(),
				drifts.stream().filter(Drift::overClaim).count(), drifts.stream().filter(d -> !d.overClaim()).count(),
				todo.size());
		for (Drift d : todo) {
			log.info("  {} recorded m={} replays r={} ({})", root.relativize(d.file()), d.recorded(), d.replayed(),
					d.overClaim() ? "OVER-CLAIM" : "loose");
			if (execute) {
				Path out = restamp(d.file(), replayer);
				log.info("    -> {}", root.relativize(out));
			}
		}
		if (!execute && !todo.isEmpty()) log.info("dry-run: pass --execute to rewrite");
	}

	static List<Path> stubs(Path root) throws java.io.IOException {
		List<Path> out = new ArrayList<>();
		try (Stream<Path> s = Files.walk(root)) {
			for (Path p : s.filter(p -> p.toString().endsWith(".json")).collect(Collectors.toList())) {
				if (SchemeIO.isStub(SchemeIO.parseJson(p.toFile()))) out.add(p);
			}
		}
		return out;
	}

	/** Replays every stub (parallel) and returns those whose replayed rank differs from the recorded one. */
	static List<Drift> drifts(List<Path> stubs, LineageReplayer replayer) {
		List<Drift> out = java.util.Collections.synchronizedList(new ArrayList<>());
		AtomicInteger done = new AtomicInteger();
		long t0 = System.currentTimeMillis();
		stubs.parallelStream().forEach(p -> {
			try {
				int recorded = SchemeIO.parseJson(p.toFile()).get("m").asInt();
				int replayed = replayer.replayFromFile(p.toFile()).r;
				if (replayed != recorded) out.add(new Drift(p, recorded, replayed));
			} catch (RuntimeException | java.io.IOException e) {
				log.warn("  {} does not replay: {}", p.getFileName(), e.toString());
			}
			int d = done.incrementAndGet();
			if (d % 500 == 0) log.info("[progress] {}/{} replayed, {} drifted, {}ms", d, stubs.size(), out.size(),
					System.currentTimeMillis() - t0);
		});
		out.sort(java.util.Comparator.comparing(d -> d.file().toString()));
		return out;
	}

	/**
	 * Rewrite one stub at its replayed rank: same lineage and {@code fields[]}, provenance
	 * keys carried over, new name from the replayed scheme (old note token kept), old file
	 * deleted. Returns the new file.
	 */
	static Path restamp(Path stub, LineageReplayer replayer) throws java.io.IOException {
		JsonNode old = SchemeIO.parseJson(stub.toFile());
		Lineage.Node lineage = SchemeIO.readLineage(old)
				.orElseThrow(() -> new IllegalStateException(stub + " has no lineage"));
		NonCubicBilinearAlgorithm rep = replayer.replayFromFile(stub.toFile());
		List<String> fields = SchemeIO.fieldTags(old);
		File out = stub.resolveSibling(SchemeIO.canonicalName(rep, noteOf(stub.getFileName().toString()))).toFile();
		SchemeIO.writeStub(rep, out, lineage, fields.isEmpty() ? null : fields);
		Map<String, Object> carry = new LinkedHashMap<>();
		old.properties().forEach(e -> {
			if (!CARRY_KEYS.contains(e.getKey())) return;
			JsonNode v = e.getValue();
			if (v.isTextual()) carry.put(e.getKey(), v.asString());
			else if (v.isBoolean()) carry.put(e.getKey(), v.asBoolean());
			else if (v.isIntegralNumber()) carry.put(e.getKey(), v.asLong());
			else if (v.isNumber()) carry.put(e.getKey(), v.asDouble());
			else if (v.isArray()) {
				List<String> items = new ArrayList<>();
				v.forEach(x -> items.add(x.asString()));
				carry.put(e.getKey(), items);
			}
		});
		if (!carry.isEmpty()) SchemeIO.addFields(out, carry, true);
		if (!out.toPath().equals(stub)) Files.delete(stub);
		return out.toPath();
	}

	/** The note token of a {@code {shape}-r{rank}-{note}-{hash7}.json} name; {@code derived} when unparseable. */
	static String noteOf(String fileName) {
		java.util.regex.Matcher m = java.util.regex.Pattern.compile("^\\d+x\\d+x\\d+-r\\d+-(.+)-[0-9a-f]{7}\\.json$")
				.matcher(fileName);
		return m.matches() ? m.group(1) : "derived";
	}

	/** For tests: the drift of one stub, if any. */
	static Optional<Drift> driftOf(Path stub, LineageReplayer replayer) throws java.io.IOException {
		int recorded = SchemeIO.parseJson(stub.toFile()).get("m").asInt();
		int replayed = replayer.replayFromFile(stub.toFile()).r;
		return replayed == recorded ? Optional.empty() : Optional.of(new Drift(stub, recorded, replayed));
	}
}
