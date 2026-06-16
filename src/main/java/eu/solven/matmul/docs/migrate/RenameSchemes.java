package eu.solven.matmul.docs.migrate;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.MatrixJsonFormatter;
import eu.solven.matmul.catalog.SchemeIO;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.node.ObjectNode;

/**
 * Catalog-wide filename normalisation to {@code <shape>-r<rank>-<note>-<hash7>.json}
 * (e.g. {@code 3x3x3-r23-ladermann71_Z-1a2b3c4.json}), AND — load-bearing — rewriting
 * every filename-stem lineage ref to a filename-independent {@code shape@<fullhash>}
 * ref so renames never dangle a lineage again.
 *
 * <p>Two phases, atomic-by-construction (build the full map, THEN rewrite+rename):</p>
 * <ol>
 *   <li><b>Map</b>: for each scheme file derive shape (JSON {@code n}), rank (JSON
 *       {@code m}), note (descriptive prefix of the OLD name, {@code derived_recursive
 *       → derived}), and full content hash (stamped {@code "hash"} field, else computed
 *       from matrices). Record {@code oldStem → "shape@fullhash"} for ref rewriting.</li>
 *   <li><b>Rewrite + rename</b>: replace each lineage {@code Atom.ref} that equals a
 *       known oldStem with its {@code shape@fullhash}; write the JSON to the new name.</li>
 * </ol>
 *
 * <p>Default is DRY-RUN (reports the plan, collisions, unresolvable refs — writes
 * nothing). Pass {@code --execute} to perform it. Re-run the manifest + DetectCyclicStubs
 * afterwards to confirm 0 corrupted / 0 cycles.</p>
 *
 * <pre>mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.RenameSchemes [-Dexec.args=--execute]</pre>
 */
public final class RenameSchemes {
	private RenameSchemes() {}

	private static final Pattern SHAPE = Pattern.compile("(\\d+)x(\\d+)x(\\d+)");
	/** Already-canonical stem: {@code {shape}-r{rank}-{note}-{hash7}} (note may contain
	 *  '_' but not '-'; hash is 4-7 hex). Group 1 is the note. */
	private static final Pattern CANONICAL =
			Pattern.compile("^\\d+x\\d+x\\d+-r\\d+-(.+)-[0-9a-f]{4,7}$");

	private record Entry(Path path, String oldStem, String shape, int rank, String note, String fullHash) {
		String newName() { return shape + "-r" + rank + "-" + note + "-" + fullHash.substring(0, 7) + ".json"; }
	}

	public static void main(String[] args) throws Exception {
		// --refs-only: rewrite lineage refs in place (filename-stem → shape@hash / bare
		// shape), keeping filenames UNTOUCHED. Decouples lineages from filenames and
		// makes every leaf resolvable (the enabler for field inference + the later
		// rename), without the risk of a full rename. Implies execute.
		boolean refsOnly = List.of(args).contains("--refs-only");
		boolean execute = refsOnly || List.of(args).contains("--execute");
		Path root = Path.of("src/main/resources/schemes");

		List<Path> files;
		try (var s = Files.walk(root)) {
			files = s.filter(p -> p.toString().endsWith(".json")).collect(Collectors.toList());
		}
		System.out.println("Scanning " + files.size() + " scheme files (mode=" + (execute ? "EXECUTE" : "DRY-RUN") + ")…");

		// Phase 1: build the entry list + oldStem → shape@hash map.
		List<Entry> entries = new ArrayList<>();
		Map<String, String> refMap = new HashMap<>();  // oldStem -> "shape@fullhash"
		Map<String, Path> oldStemToPath = new HashMap<>();
		int hashErrors = 0, processed = 0;
		for (Path f : files) {
			String stem = f.getFileName().toString().replaceFirst("\\.json$", "");
			Entry e = buildEntry(f, stem);
			if (e == null) { hashErrors++; continue; }
			entries.add(e);
			refMap.put(stem, e.shape() + "@" + e.fullHash().substring(0, 7));
			oldStemToPath.put(stem, f);
			if (++processed % 2000 == 0) System.out.println("[progress] " + processed + "/" + files.size() + " mapped");
		}

		// Collisions: two distinct files mapping to the same new filename.
		Map<String, List<Entry>> byNew = new LinkedHashMap<>();
		for (Entry e : entries) byNew.computeIfAbsent(e.newName(), k -> new ArrayList<>()).add(e);
		// Same new name from distinct files: benign if identical content (same full
		// hash → same file written, the extra old files are just dedup'd); a REAL
		// collision only if content differs (a 7-char-hash clash that would lose data).
		List<Map.Entry<String, List<Entry>>> collisions = byNew.entrySet().stream()
				.filter(en -> en.getValue().size() > 1)
				.filter(en -> en.getValue().stream().map(Entry::fullHash).distinct().count() > 1)
				.collect(Collectors.toList());
		int benignDupes = byNew.values().stream().filter(v -> v.size() > 1
				&& v.stream().map(Entry::fullHash).distinct().count() == 1)
				.mapToInt(v -> v.size() - 1).sum();

		// Phase 2 scan: count refs to rewrite + unresolvable filename-ish refs.
		int refsToRewrite = 0, toBareShape = 0, stuck = 0;
		List<String> stuckSample = new ArrayList<>();
		for (Entry e : entries) {
			for (String ref : lineageRefs(e.path())) {
				if (refMap.containsKey(ref)) refsToRewrite++;
				else if (looksLikeFilenameStem(ref) && !ref.contains("@")) {
					if (shapeOf(ref) != null) toBareShape++;
					else { stuck++; if (stuckSample.size() < 10) stuckSample.add(ref); }
				}
			}
		}

		// Only LEGACY-named files are renamed. An already-canonical name is kept
		// verbatim — never re-hashed — so a re-run can't churn the tree on hash-source
		// drift (sparse u_sparse files recompute a different hash than they were named
		// with). Refs are still normalised everywhere (in place).
		long actualRenames = entries.stream()
				.filter(e -> !CANONICAL.matcher(e.oldStem()).matches()).count();
		System.out.println("\n=== PLAN ===");
		System.out.println("files processed:        " + entries.size() + "  (" + hashErrors + " skipped: unreadable/no-hash)");
		System.out.println("ACTUAL name changes (legacy → canonical): " + actualRenames
				+ "  (already-canonical names kept verbatim)");
		System.out.println("lineage refs → shape@hash (pinned): " + refsToRewrite);
		System.out.println("lineage refs → bare shape (drifted stem, rename-proof): " + toBareShape);
		System.out.println("benign same-content dupes: " + benignDupes + "  (extra old files dedup'd to one)");
		System.out.println("REAL collisions (differing content, same name): " + collisions.size());
		for (var c : collisions.stream().limit(8).collect(Collectors.toList()))
			System.out.println("    " + c.getKey() + "  ← " + c.getValue().stream().map(Entry::oldStem).collect(Collectors.joining(", ")));
		System.out.println("truly stuck refs (no shape, left as-is): " + stuck + "  " + stuckSample);
		System.out.println("\nsample renames (legacy only):");
		for (Entry e : entries.stream()
				.filter(e -> !CANONICAL.matcher(e.oldStem()).matches())
				.limit(8).collect(Collectors.toList()))
			System.out.println("    " + e.oldStem() + ".json\n      → " + e.newName());

		if (!execute) {
			System.out.println("\n(DRY-RUN — nothing written. Pass --execute to perform; collisions MUST be 0 first.)");
			return;
		}
		if (!refsOnly && !collisions.isEmpty()) {
			System.out.println("\nABORT: " + collisions.size() + " new-name collisions — refusing to execute (would overwrite). Resolve first.");
			return;
		}

		// Phase 2 execute: rewrite refs; rename ONLY legacy-named files (canonical names
		// are kept verbatim — write in place only when their refs changed). UNLESS
		// --refs-only (then always write in place).
		int written = 0, refsChanged = 0, renamed = 0;
		for (Entry e : entries) {
			JsonNode json = SchemeIO.parseJson(e.path().toFile());
			if (!(json instanceof ObjectNode obj)) {
				continue;
			}
			boolean canonical = CANONICAL.matcher(e.oldStem()).matches();
			String before = obj.toString();
			rewriteRefs(obj.get("lineage"), refMap);
			boolean refsChangedThis = !obj.toString().equals(before);
			if (refsChangedThis) refsChanged++;
			// Idempotent: an already-canonical file with no ref change is left untouched.
			if (canonical && !refsChangedThis && !refsOnly) {
				continue;
			}
			boolean doRename = !refsOnly && !canonical;
			Path out = doRename ? e.path().getParent().resolve(e.newName()) : e.path();
			Files.writeString(out, MatrixJsonFormatter.format(obj));
			if (doRename && !out.equals(e.path())) {
				Files.deleteIfExists(e.path());
				renamed++;
			}
			if (++written % 2000 == 0) System.out.println("[progress] " + written + "/" + entries.size());
		}
		System.out.println("\nDONE: renamed " + renamed + " legacy files; wrote " + written
				+ " files total (" + refsChanged + " had ref changes). "
				+ "Re-run GenerateCatalogManifest + DetectCyclicStubs to validate.");
	}

	private static Entry buildEntry(Path f, String stem) {
		try {
			JsonNode root = SchemeIO.parseJson(f.toFile());
			String note = noteFromStem(stem);
			JsonNode nNode = root.get("n");
			if (nNode == null || !nNode.isArray() || nNode.size() != 3) return null;  // no shape
			String shape = nNode.get(0).asInt() + "x" + nNode.get(1).asInt() + "x" + nNode.get(2).asInt();
			int rank = root.has("m") ? root.get("m").asInt() : (root.has("rank") ? root.get("rank").asInt() : -1);
			if (rank < 0) return null;
			String hash;
			if (root.has("hash")) {
				// Stamped content hash (newer stubs + materialiser output) — the value
				// findByHash matches; use it directly.
				hash = root.get("hash").asString();
			} else if (root.has("u") && root.has("v") && root.has("w")) {
				// Loadable matrices (imports): compute the content hash findByHash will match.
				NonCubicBilinearAlgorithm alg = SchemeIO.isReduced(root) ? SchemeIO.readReduced(root) : SchemeIO.read(root);
				hash = SchemeIO.contentHash(alg);
			} else {
				// Hashless stub / cited entry (no matrices). findByHash SKIPS stubs, so a
				// ref to it resolves by SHAPE regardless — the hash here is only a stable,
				// unique FILENAME label. Derive it cheaply from the canonical lineage (no
				// replay), falling back to metadata for cited entries.
				String basis = root.has("lineage_compact") ? root.get("lineage_compact").asString()
						: root.has("lineage_str") ? root.get("lineage_str").asString()
						: root.has("lineage") ? root.get("lineage").toString()
						: shape + "r" + rank + note + (root.has("reference") ? root.get("reference").asString() : "");
				hash = sha256hex(shape + ";" + basis);
			}
			if (hash.length() < 7) hash = hash + "0000000";
			return new Entry(f, stem, shape, rank, note, hash);
		} catch (Exception ex) {
			return null;  // genuinely unreadable → skip (reported as count)
		}
	}

	private static String sha256hex(String s) {
		try {
			byte[] d = java.security.MessageDigest.getInstance("SHA-256").digest(s.getBytes(java.nio.charset.StandardCharsets.UTF_8));
			StringBuilder sb = new StringBuilder();
			for (byte b : d) sb.append(String.format("%02x", b));
			return sb.toString();
		} catch (Exception e) {
			return Integer.toHexString(s.hashCode());
		}
	}

	/** Descriptive prefix of the old filename (before the shape token); derived_recursive → derived. */
	private static String noteFromStem(String stem) {
		// IDEMPOTENCY: a stem that is ALREADY canonical ({shape}-r{rank}-{note}-{hash7})
		// keeps its note verbatim. Without this the shape-first canonical form fell to
		// the prefix branch below → note "" → "scheme", so EVERY canonical file looked
		// like a rename (non-idempotent; would churn the whole tree on a re-run).
		Matcher c = CANONICAL.matcher(stem);
		if (c.matches()) {
			return c.group(1);
		}
		Matcher m = SHAPE.matcher(stem);
		String note = m.find() ? stem.substring(0, m.start()) : stem;
		note = note.replaceAll("[-_]+$", "");
		if (note.isBlank()) note = "scheme";
		if (note.equals("derived_recursive")) note = "derived";
		return note;
	}

	/** A ref that is plausibly an old filename stem (has a shape and an author-ish prefix). */
	private static boolean looksLikeFilenameStem(String ref) {
		if (ref.startsWith("naive") || ref.contains("(")) return false;  // naive-/parametric
		Matcher m = SHAPE.matcher(ref);
		return m.find() && m.start() > 0;  // something before the shape → a named stem, not a bare shape
	}

	private static List<String> lineageRefs(Path f) {
		List<String> out = new ArrayList<>();
		try {
			collectRefs(SchemeIO.parseJson(f.toFile()).get("lineage"), out);
		} catch (Exception e) { /* none */ }
		return out;
	}

	private static void collectRefs(JsonNode node, List<String> out) {
		if (node == null) return;
		if (node.isObject()) {
			JsonNode ref = node.get("ref");
			if (ref != null && ref.isString()) out.add(ref.asString());
			for (var e : node.properties()) collectRefs(e.getValue(), out);
		} else if (node.isArray()) {
			for (JsonNode c : node) collectRefs(c, out);
		}
	}

	private static void rewriteRefs(JsonNode node, Map<String, String> refMap) {
		if (node == null) return;
		if (node instanceof ObjectNode obj) {
			JsonNode ref = obj.get("ref");
			if (ref != null && ref.isString()) {
				String r = ref.asString();
				if (refMap.containsKey(r)) {
					obj.put("ref", refMap.get(r));            // pinned → shape@hash7
				} else if (looksLikeFilenameStem(r) && !r.contains("@")) {
					int[] sh = shapeOf(r);                    // drifted filename-stem (its file
					if (sh != null) obj.put("ref", sh[0] + "x" + sh[1] + "x" + sh[2]); // doesn't exact-exist) → bare shape-ref (same resolution, rename-proof)
				}
			}
			for (var e : obj.properties()) rewriteRefs(e.getValue(), refMap);
		} else if (node.isArray()) {
			for (JsonNode c : node) rewriteRefs(c, refMap);
		}
	}

	private static int[] shapeOf(String ref) {
		Matcher m = SHAPE.matcher(ref);
		return m.find() ? new int[] { Integer.parseInt(m.group(1)),
				Integer.parseInt(m.group(2)), Integer.parseInt(m.group(3)) } : null;
	}
}
