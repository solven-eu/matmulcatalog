package eu.solven.matmul.docs.migrate;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
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
import eu.solven.matmul.catalog.MatrixJsonFormatter;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.search.LineageReplayer;
import eu.solven.matmul.verifiers.Verifier;
import lombok.extern.slf4j.Slf4j;
import tools.jackson.databind.JsonNode;
import tools.jackson.databind.node.ObjectNode;

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
 * replay (stale, honest-but-loose — 114 of them below the catalog-best their shape
 * reported, 2026-10-18). {@code VerifySchemes --check=rank} fails on the first kind;
 * this driver repairs either kind by rewriting the stub at its replayed rank (same
 * lineage, provenance keys carried over, file renamed to the new rank/hash).
 *
 * <p>The bulk run is four phases: (1) <b>stage</b> — every drifted stub is replayed (in
 * parallel; the replay is the slow part) and its honest version written to a staging
 * directory outside the catalog, so nothing the other replays resolve through changes
 * under them; (2) <b>place</b> — sequentially, each staged file is moved to its canonical
 * name, unless the catalog already holds that content (the canonical name itself, or any
 * same-shape file carrying the content hash — a known import the derivation merely
 * re-wraps, {@code Project(ConcatCols(naive, X)) = X}): then the stale stub is a pure
 * duplicate and nothing is written; (3) <b>re-pin</b> — every lineage that pinned a
 * restamped stub's old stamped hash is re-pointed to the new hash (chains followed to
 * the end), exact by construction since the stub's lineage is unchanged; (4) the stale
 * files are deleted. A session kill between phases leaves the catalog consistent: stale
 * stubs stay loose and resolvable, honest files already placed are recognised as twins
 * on the next run.
 *
 * <pre>
 *   --direction=over|all   over (default): only recorded &lt; replayed (false claims);
 *                          all: also tighten recorded &gt; replayed
 *   --execute              write (default: dry-run, report only)
 *   --schemes-root=PATH    default src/main/resources/schemes
 *   --staging=PATH         default target/restamp-staging (wiped at start)
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

	/** One staged rewrite: the stale stub, its honest version in the staging dir, the canonical target. */
	public record Staged(Path stub, Path staged, Path out, int[] shape, String oldHash, String newHash) {}

	/** One placed rewrite: where it was, where it is now (an existing twin when {@code duplicate}),
	 *  and the stamped hash before/after. */
	public record Restamped(Path oldFile, Path newFile, String oldHash, String newHash, boolean duplicate) {
		public Restamped(Path oldFile, Path newFile, String oldHash, String newHash) {
			this(oldFile, newFile, oldHash, newHash, false);
		}
	}

	public static void main(String[] args) throws Exception {
		List<String> argv = Arrays.asList(args);
		boolean execute = argv.contains("--execute");
		boolean all = argv.stream().anyMatch(a -> a.equals("--direction=all"));
		Path root = Path.of(argv.stream().filter(a -> a.startsWith("--schemes-root=")).map(a -> a.substring(15))
				.findFirst().orElse("src/main/resources/schemes"));
		Path staging = Path.of(argv.stream().filter(a -> a.startsWith("--staging=")).map(a -> a.substring(10))
				.findFirst().orElse("target/restamp-staging"));
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
		}
		if (!execute) {
			if (!todo.isEmpty()) log.info("dry-run: pass --execute to rewrite");
			return;
		}

		// Phase 1: stage (parallel). Replays only read the catalog; the honest files land
		// outside it, so no replay sees a half-written catalog.
		wipe(staging);
		Files.createDirectories(staging);
		List<Staged> staged = Collections.synchronizedList(new ArrayList<>());
		List<String> failed = Collections.synchronizedList(new ArrayList<>());
		AtomicInteger n = new AtomicInteger();
		long t0 = System.currentTimeMillis();
		todo.parallelStream().forEach(d -> {
			try {
				staged.add(stage(d.file(), replayer, staging));
			} catch (Exception | StackOverflowError e) {
				failed.add(root.relativize(d.file()) + ": " + e);
				log.error("  {} not staged: {}", root.relativize(d.file()), e.toString());
			}
			int k = n.incrementAndGet();
			if (k % 50 == 0) log.info("[progress] staged {}/{}, {} failed, {}ms", k, todo.size(), failed.size(),
					System.currentTimeMillis() - t0);
		});
		log.info("staged {} honest stub(s) under {} ({} failed)", staged.size(), staging, failed.size());

		// Phase 2: place (sequential, path order).
		List<Staged> ordered = new ArrayList<>(staged);
		ordered.sort(Comparator.comparing(s -> s.stub().toString()));
		List<Restamped> done = new ArrayList<>();
		int duplicates = 0;
		for (Staged s : ordered) {
			Restamped r = place(s, lookup, false);
			done.add(r);
			if (r.duplicate()) duplicates++;
			log.info("  {} -> {}{}", root.relativize(s.stub()), root.relativize(r.newFile()),
					r.duplicate() ? "  (duplicate of an existing twin)" : "");
		}

		// Phase 3: re-pin; phase 4: delete the stale files.
		int repinned = rewritePins(root, done);
		int deleted = 0;
		for (Restamped r : done) {
			if (!r.newFile().equals(r.oldFile()) && Files.deleteIfExists(r.oldFile())) deleted++;
		}
		wipe(staging);
		log.info("restamped {} stub(s) ({} written, {} duplicates of an existing twin, {} failed); "
				+ "re-pinned {} file(s) that pinned a restamped stub's old hash; deleted {} stale file(s)",
				done.size(), done.size() - duplicates, duplicates, failed.size(), repinned, deleted);
		for (String f : failed) log.warn("  FAILED {}", f);
	}

	private static void wipe(Path dir) throws IOException {
		if (!Files.exists(dir)) return;
		try (Stream<Path> s = Files.walk(dir)) {
			for (Path p : s.sorted(Comparator.reverseOrder()).collect(Collectors.toList())) Files.deleteIfExists(p);
		}
	}

	static List<Path> stubs(Path root) throws IOException {
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
		List<Drift> out = Collections.synchronizedList(new ArrayList<>());
		AtomicInteger done = new AtomicInteger();
		long t0 = System.currentTimeMillis();
		stubs.parallelStream().forEach(p -> {
			try {
				int recorded = SchemeIO.parseJson(p.toFile()).get("m").asInt();
				int replayed = replayer.replayFromFile(p.toFile()).r;
				if (replayed != recorded) out.add(new Drift(p, recorded, replayed));
			} catch (RuntimeException | IOException e) {
				log.warn("  {} does not replay: {}", p.getFileName(), e.toString());
			}
			int d = done.incrementAndGet();
			if (d % 500 == 0) log.info("[progress] {}/{} replayed, {} drifted, {}ms", d, stubs.size(), out.size(),
					System.currentTimeMillis() - t0);
		});
		out.sort(Comparator.comparing(d -> d.file().toString()));
		return out;
	}

	/**
	 * Phase 1 for one stub: replay it, spot-check the result, and write its honest version
	 * (same lineage and {@code fields[]}, provenance keys carried over, old note token kept)
	 * into {@code stagingDir} under a name unique to the stub. Safe to run in parallel: the
	 * catalog is only read.
	 */
	static Staged stage(Path stub, LineageReplayer replayer, Path stagingDir) throws IOException {
		JsonNode old = SchemeIO.parseJson(stub.toFile());
		Lineage.Node lineage = SchemeIO.readLineage(old)
				.orElseThrow(() -> new IllegalStateException(stub + " has no lineage"));
		NonCubicBilinearAlgorithm rep = replayer.replayFromFile(stub.toFile());
		// The lineage operators are trusted (VerifyAllSchemes certifies them), but a stub
		// that is about to become a catalog record gets the cheap independent check too.
		if (!Verifier.passesRandomMatmulSpotCheck(rep, 2_000, 42L)) {
			throw new IllegalStateException(stub + ": the replayed scheme fails the matmul spot-check — not rewriting");
		}
		List<String> fields = SchemeIO.fieldTags(old);
		String canonical = SchemeIO.canonicalName(rep, noteOf(stub.getFileName().toString()));
		Path out = stub.resolveSibling(canonical);
		Path staged = stagingDir.resolve(stub.getFileName().toString() + ".new.json");
		SchemeIO.writeStub(rep, staged.toFile(), lineage, fields.isEmpty() ? null : fields);
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
		if (!carry.isEmpty()) SchemeIO.addFields(staged.toFile(), carry, true);
		return new Staged(stub, staged, out, new int[] { rep.n, rep.m, rep.p }, SchemeIO.readHash(old),
				SchemeIO.contentHash(rep));
	}

	/**
	 * Phase 2 for one staged rewrite: move it to its canonical name — or, when the catalog
	 * already holds that content (the canonical name itself, or any same-shape file the
	 * {@code lookup} knows that carries the content hash, e.g. a known import the derivation
	 * merely re-wraps), drop the staged copy: the stale stub is a pure duplicate, the twin
	 * keeps its lineage, provenance and computed metrics, and phase 3 re-points the
	 * duplicate's pins to the twin's hash (identical by construction). The stale stub is
	 * deleted now when {@code deleteOld}; the bulk run keeps it until the pins are rewritten.
	 */
	static Restamped place(Staged s, FieldAwareLookup lookup, boolean deleteOld) throws IOException {
		Optional<Path> twin = existingTwin(s.stub(), s.out(), s.shape(), s.newHash(), lookup);
		if (twin.isPresent()) {
			log.info("  {} is a duplicate of {} (same content) — stale copy removed, nothing overwritten",
					s.stub().getFileName(), twin.get().getFileName());
			Files.deleteIfExists(s.staged());
			if (deleteOld) Files.deleteIfExists(s.stub());
			return new Restamped(s.stub(), twin.get(), s.oldHash(), s.newHash(), true);
		}
		Files.move(s.staged(), s.out(), StandardCopyOption.REPLACE_EXISTING);
		FieldAwareLookup.onSchemeWritten(s.out().toFile());
		if (deleteOld && !s.out().equals(s.stub())) Files.deleteIfExists(s.stub());
		return new Restamped(s.stub(), s.out(), s.oldHash(), s.newHash());
	}

	/** Single-file convenience: rewrite the stub and delete the stale file. Returns the new file. */
	static Path restamp(Path stub, LineageReplayer replayer) throws IOException {
		return restamp(stub, replayer, null, true).newFile();
	}

	/** Single-file: stage into a temporary directory, then place. */
	static Restamped restamp(Path stub, LineageReplayer replayer, FieldAwareLookup lookup, boolean deleteOld)
			throws IOException {
		Path staging = Files.createTempDirectory("restamp-");
		try {
			return place(stage(stub, replayer, staging), lookup, deleteOld);
		} finally {
			wipe(staging);
		}
	}

	/**
	 * Phase 3: every lineage in the catalog that pins a restamped stub's OLD stamped hash
	 * ({@code shape@hash}, 7 or 64 hex digits) is re-pointed to the NEW hash. Exact by
	 * construction — the restamped stub carries the same lineage, so the pinning file
	 * already replayed through it to the same matrices; only the address changed. Returns
	 * the number of files rewritten (canonically, via {@link MatrixJsonFormatter}).
	 */
	static int rewritePins(Path root, List<Restamped> restamped) throws IOException {
		Map<String, String> next = new HashMap<>();
		for (Restamped r : restamped) {
			if (r.oldHash() != null && r.oldHash().length() >= 7 && !r.oldHash().equals(r.newHash())) {
				next.put(r.oldHash(), r.newHash());
			}
		}
		// Transitive closure: a duplicate's pins point at its twin's hash, and the twin may
		// itself have been restamped in the same run (old → twin-old → twin-new).
		Map<String, Restamped> byOldHash7 = new HashMap<>();
		for (Map.Entry<String, String> e : next.entrySet()) {
			String fin = e.getValue();
			for (int hop = 0; hop < 16 && next.containsKey(fin) && !next.get(fin).equals(fin); hop++) fin = next.get(fin);
			byOldHash7.put(e.getKey().substring(0, 7), new Restamped(null, null, e.getKey(), fin));
		}
		if (byOldHash7.isEmpty()) return 0;
		List<Path> files;
		try (Stream<Path> s = Files.walk(root)) {
			files = s.filter(p -> p.toString().endsWith(".json")).collect(Collectors.toList());
		}
		int rewritten = 0;
		for (Path p : files) {
			String raw = Files.readString(p);
			if (!raw.contains("\"lineage\"") || !raw.contains("@")) continue;
			JsonNode node = SchemeIO.parseJson(p.toFile());
			if (rewriteRefs(node.get("lineage"), byOldHash7) == 0) continue;
			Files.writeString(p, MatrixJsonFormatter.format(node));
			rewritten++;
		}
		return rewritten;
	}

	private static int rewriteRefs(JsonNode node, Map<String, Restamped> byOldHash7) {
		if (node == null) return 0;
		int n = 0;
		if (node.isObject()) {
			ObjectNode o = (ObjectNode) node;
			JsonNode op = o.get("op");
			JsonNode ref = o.get("ref");
			if (op != null && "Atom".equals(op.asString()) && ref != null && ref.isTextual()) {
				String r = ref.asString();
				int at = r.indexOf('@');
				if (at > 0 && r.length() - at - 1 >= 7) {
					String hash = r.substring(at + 1);
					Restamped hit = byOldHash7.get(hash.substring(0, 7));
					if (hit != null && (hit.oldHash().startsWith(hash) || hash.startsWith(hit.oldHash()))) {
						String newHash = hash.length() <= 7 ? hit.newHash().substring(0, 7) : hit.newHash();
						o.put("ref", r.substring(0, at + 1) + newHash);
						n++;
					}
				}
			}
			for (Map.Entry<String, JsonNode> e : o.properties()) n += rewriteRefs(e.getValue(), byOldHash7);
		} else if (node.isArray()) {
			for (JsonNode c : node) n += rewriteRefs(c, byOldHash7);
		}
		return n;
	}

	/**
	 * A file that already holds the replayed content: the stub's own canonical name at the
	 * new rank, or any file at the shape whose name carries the content hash7 and whose
	 * stamped hash (when present) agrees. Explicit files and stubs alike — the canonical
	 * name's hash7 is the content-hash prefix for both.
	 */
	static Optional<Path> existingTwin(Path stub, Path out, int[] shape, String newHash, FieldAwareLookup lookup) {
		if (Files.exists(out) && !out.equals(stub)) return Optional.of(out);
		if (lookup == null) return Optional.empty();
		String suffix = "-" + newHash.substring(0, 7) + ".json";
		for (Path c : lookup.findFiles(shape[0], shape[1], shape[2])) {
			if (c.toAbsolutePath().equals(stub.toAbsolutePath()) || !c.getFileName().toString().endsWith(suffix)) {
				continue;
			}
			try {
				String stamped = SchemeIO.readHash(SchemeIO.parseJson(c.toFile()));
				if (stamped == null || stamped.equals(newHash)) return Optional.of(c);
			} catch (Exception e) {
				// unreadable candidate — not a twin we can vouch for
			}
		}
		return Optional.empty();
	}

	/** The note token of a {@code {shape}-r{rank}-{note}-{hash7}.json} name; {@code derived} when unparseable. */
	static String noteOf(String fileName) {
		java.util.regex.Matcher m = java.util.regex.Pattern.compile("^\\d+x\\d+x\\d+-r\\d+-(.+)-[0-9a-f]{7}\\.json$")
				.matcher(fileName);
		return m.matches() ? m.group(1) : "derived";
	}

	/** For tests: the drift of one stub, if any. */
	static Optional<Drift> driftOf(Path stub, LineageReplayer replayer) throws IOException {
		int recorded = SchemeIO.parseJson(stub.toFile()).get("m").asInt();
		int replayed = replayer.replayFromFile(stub.toFile()).r;
		return replayed == recorded ? Optional.empty() : Optional.of(new Drift(stub, recorded, replayed));
	}
}
