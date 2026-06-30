package eu.solven.matmul.docs.migrate;

import eu.solven.matmul.docs.verify.AuditLineageRefs;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.MatrixJsonFormatter;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.search.LineageReplayer;

import lombok.extern.slf4j.Slf4j;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.node.ObjectNode;

/**
 * Repair pass for dangling pinned lineage refs (policy <b>b1: verified re-pin,
 * regenerate-in-place on content change</b>) — the {@link AuditLineageRefs}
 * follow-up.
 *
 * <p>A dangling ref is {@code shape@hash} whose hash resolves to nothing — in
 * practice because the 2026-06 purge replaced the pinned ingredient with an
 * equal-rank substitute, so bit-exact reproduction of the parent is impossible
 * <i>by design</i>. Per-file verification ladder (first tier that succeeds wins):</p>
 * <ol>
 *   <li><b>REPINNED_BITEXACT</b> — some candidate assignment reproduces the
 *       parent's stored content hash exactly (drift cases). Lineage updated,
 *       content untouched.</li>
 *   <li><b>REGENERATED</b> — a candidate assignment replays to a scheme of the
 *       <b>same shape and rank</b> that passes the randomized matmul spot-check.
 *       The parent's content legitimately changed (its ingredient was replaced):
 *       lineage re-pinned, stamped {@code hash} updated, stale per-content metrics
 *       ({@code additions}, {@code buds}, …) removed for later re-stamping by
 *       {@code EnrichSchemeMetrics}. Stub-only — materialised parents whose
 *       content would change are reported as {@code MATERIALISED_CHANGE} and left
 *       untouched (their on-disk matrices are still valid; rewriting them is a
 *       separate decision).</li>
 *   <li><b>QUARANTINE</b> — no candidate replays to the stored rank (incl. rank
 *       regressions — the safety net) → listed, untouched.</li>
 * </ol>
 *
 * <p><b>Fixpoint:</b> regenerating a stub changes its stamped hash, so files
 * pinning the OLD stamped hash dangle in turn. Each pass records
 * {@code oldHash → newHash} and re-scans; iteration ends when no file needs
 * repair (lineage depth bounds the pass count). Default DRY-RUN; {@code --apply}
 * writes canonically via {@link MatrixJsonFormatter}.</p>
 */
@Slf4j
public class RepinDanglingLineageRefs {

	private static final Path ROOT = Path.of("src/main/resources/schemes");
	private static final Pattern PINNED = Pattern.compile("^(\\d+)x(\\d+)x(\\d+)@([0-9a-f]{6,})$");
	private static final Pattern FILE_HASH7 = Pattern.compile("-([0-9a-f]{7})\\.json$");
	private static final int MAX_COMBOS = 64;
	private static final int MAX_PASSES = 6;
	private static final List<String> STALE_METRIC_FIELDS =
			List.of("additions", "min_additions", "has_buds", "buds", "bud_score", "projection_margin", "slp");

	private final FieldAwareLookup lookup = new FieldAwareLookup("C");
	private final LineageReplayer replayer = LineageReplayer.withDefaultPool(lookup);
	/** candidate scheme file → full content/stamped hash. */
	private final Map<Path, String> candidateHash = new LinkedHashMap<>();
	/** dangling ref → candidate hash that worked elsewhere (try-first hint). */
	private final Map<String, String> successHint = new LinkedHashMap<>();
	/** regenerated old stamped hash → new hash (fixpoint propagation). */
	private final Map<String, String> rewriteMap = new LinkedHashMap<>();
	/** quarantined file → rewriteMap size at the time; retry ONLY if the map has
	 *  grown since (a new substitute could rescue it). Without this, every fixpoint
	 *  pass re-chews the whole quarantine set with identical verdicts — the
	 *  2026-06-10 apply run burned 5 redundant passes (~30 min) doing exactly that. */
	private final Map<Path, Integer> quarantinedAt = new LinkedHashMap<>();

	public static void main(String[] args) throws Exception {
		boolean apply = false;
		int maxDim = 20;
		for (String a : args) {
			if (a.equals("--apply")) apply = true;
			else if (a.startsWith("--max-dim=")) maxDim = Integer.parseInt(a.substring("--max-dim=".length()));
		}
		new RepinDanglingLineageRefs().run(apply, maxDim);
	}

	void run(boolean apply, int maxDim) throws Exception {
		long start = System.nanoTime();
		List<Path> lineageFiles = new ArrayList<>();
		try (Stream<Path> walk = Files.walk(ROOT)) {
			for (Path p : walk.filter(q -> q.toString().endsWith(".json")).sorted().toList()) {
				String raw = Files.readString(p);
				if (raw.contains("\"lineage\"") && raw.contains("@")) lineageFiles.add(p);
			}
		}
		log.info("{} lineage-bearing files with pinned refs (policy b1, {}, max-dim {})",
				lineageFiles.size(), apply ? "APPLY" : "DRY-RUN", maxDim);

		Map<String, List<String>> outcome = new TreeMap<>();
		for (int pass = 1; pass <= MAX_PASSES; pass++) {
			List<Path> affected = new ArrayList<>();
			for (Path p : lineageFiles) {
				Integer qAt = quarantinedAt.get(p);
				if (qAt != null && qAt == rewriteMap.size()) {
					continue; // already quarantined and nothing new could rescue it
				}
				JsonNode root = SchemeIO.parseJson(p.toFile());
				JsonNode lineage = root.get("lineage");
				if (lineage == null || lineage.isNull()) continue;
				if (!needsRepair(lineage)) continue;
				affected.add(p);
			}
			if (affected.isEmpty()) {
				log.info("pass {}: fixpoint reached (no file needs repair)", pass);
				break;
			}
			affected.sort(java.util.Comparator.comparingInt(RepinDanglingLineageRefs::maxDimOf));
			log.info("pass {}: {} files need repair", pass, affected.size());

			int done = 0;
			for (Path p : affected) {
				done++;
				if (done % 50 == 0) {
					log.info("[progress] pass {}: {}/{} processed, {} ms — {}", pass, done,
							affected.size(), (System.nanoTime() - start) / 1_000_000, sizes(outcome));
				}
				String verdict;
				try {
					verdict = repairOne(p, apply, maxDim);
				} catch (Exception | StackOverflowError e) {
					verdict = "QUARANTINE";
					log.warn("  {} → exception: {}: {}", p.getFileName(),
							e.getClass().getSimpleName(), String.valueOf(e.getMessage()));
				}
				if (verdict.equals("QUARANTINE")) {
					quarantinedAt.put(p, rewriteMap.size());
				} else {
					quarantinedAt.remove(p);
				}
				outcome.computeIfAbsent(verdict, k -> new ArrayList<>()).add(p.toString());
			}
			if (!apply) {
				log.info("DRY-RUN: stopping after pass 1 (fixpoint propagation needs writes)");
				break;
			}
		}

		log.info("==================================================================");
		log.info(" Re-pin b1 {} — outcomes (a file may appear once per pass):", apply ? "APPLIED" : "DRY-RUN");
		log.info("==================================================================");
		for (var e : outcome.entrySet()) {
			log.info("  {}  ×{}", String.format("%-22s", e.getKey()), e.getValue().size());
		}
		Path quarantine = Path.of("target/lineage-quarantine.txt");
		List<String> qLines = new ArrayList<>();
		for (String k : List.of("QUARANTINE", "MATERIALISED_CHANGE", "DEFERRED_TOO_BIG")) {
			for (String f : outcome.getOrDefault(k, List.of())) qLines.add(k + "\t" + f);
		}
		Files.createDirectories(quarantine.getParent());
		Files.write(quarantine, qLines);
		log.info("Quarantine/deferred list ({} entries) → {}", qLines.size(), quarantine);
		log.info("Done in {} ms.", (System.nanoTime() - start) / 1_000_000);
	}

	private static String sizes(Map<String, List<String>> m) {
		StringBuilder sb = new StringBuilder();
		m.forEach((k, v) -> sb.append(k).append('=').append(v.size()).append(' '));
		return sb.toString();
	}

	private static int maxDimOf(Path p) {
		try {
			JsonNode n = SchemeIO.parseJson(p.toFile()).get("n");
			int max = 0;
			for (JsonNode d : n) max = Math.max(max, d.asInt());
			return max;
		} catch (Exception e) {
			return Integer.MAX_VALUE;
		}
	}

	/** A file needs repair iff it has a truly-dangling pin or pins a rewritten hash. */
	private boolean needsRepair(JsonNode lineage) {
		List<String> refs = new ArrayList<>();
		AuditLineageRefs.collectAtomRefs(lineage, refs);
		for (String ref : refs) {
			Matcher m = PINNED.matcher(ref);
			if (!m.matches()) continue;
			String hash = m.group(4);
			if (rewriteMap.keySet().stream().anyMatch(old -> old.startsWith(hash) || hash.startsWith(old))) {
				return true;
			}
			if (isDangling(m)) return true;
		}
		return false;
	}

	/** Dangling = neither a materialised candidate's content hash nor a stub's stamped hash matches. */
	private boolean isDangling(Matcher pinned) {
		int n = Integer.parseInt(pinned.group(1));
		int mm = Integer.parseInt(pinned.group(2));
		int p = Integer.parseInt(pinned.group(3));
		String hash = pinned.group(4);
		if (lookup.findByHash(n, mm, p, hash).isPresent()) return false;
		for (Path c : lookup.findFiles(n, mm, p)) {
			String h = hashOf(c);
			if (h != null && h.startsWith(hash)) return false;
		}
		return true;
	}

	/** @return outcome label for this file. */
	private String repairOne(Path p, boolean apply, int maxDim) throws Exception {
		if (maxDimOf(p) > maxDim) {
			return "DEFERRED_TOO_BIG";
		}
		ObjectNode root = (ObjectNode) SchemeIO.parseJson(p.toFile());
		JsonNode lineage = root.get("lineage");

		// Refs to fix: truly-dangling pins + pins to rewritten hashes.
		List<String> refs = new ArrayList<>();
		AuditLineageRefs.collectAtomRefs(lineage, refs);
		Set<String> toFix = new LinkedHashSet<>();
		Map<String, String> forcedAssignment = new LinkedHashMap<>(); // rewritten-hash pins: target known
		for (String ref : refs) {
			Matcher m = PINNED.matcher(ref);
			if (!m.matches()) continue;
			String hash = m.group(4);
			Optional<String> rewritten = rewriteMap.entrySet().stream()
					.filter(e -> e.getKey().startsWith(hash) || hash.startsWith(e.getKey()))
					.map(Map.Entry::getValue).findFirst();
			if (rewritten.isPresent()) {
				forcedAssignment.put(ref, rewritten.get());
			} else if (isDangling(m)) {
				toFix.add(ref);
			}
		}
		if (toFix.isEmpty() && forcedAssignment.isEmpty()) return "ALREADY_OK";

		boolean hasMatrices = root.has("u") || root.has("u_sparse");
		String storedHash = SchemeIO.readHash(root);
		String fileHash7 = null;
		Matcher fh = FILE_HASH7.matcher(p.getFileName().toString());
		if (fh.find()) fileHash7 = fh.group(1);
		String exactAnchor = hasMatrices ? SchemeIO.contentHash(SchemeIO.readBilinear(p.toFile()))
				: (storedHash != null ? storedHash : fileHash7);
		int storedRank = root.get("m") != null ? root.get("m").asInt() : -1;
		JsonNode shape = root.get("n");
		int sn = shape.get(0).asInt(), sm = shape.get(1).asInt(), sp = shape.get(2).asInt();

		// Candidate hashes per dangling ref.
		Map<String, List<String>> candidates = new LinkedHashMap<>();
		for (String ref : toFix) {
			Matcher m = PINNED.matcher(ref);
			m.matches();
			int n = Integer.parseInt(m.group(1)), mm = Integer.parseInt(m.group(2)), pp = Integer.parseInt(m.group(3));
			List<String> hashes = new ArrayList<>();
			String hint = successHint.get(ref);
			if (hint != null) hashes.add(hint);
			for (Path c : lookup.findFiles(n, mm, pp)) {
				String h = hashOf(c);
				if (h != null && !hashes.contains(h)) hashes.add(h);
			}
			if (hashes.isEmpty()) return "QUARANTINE";
			candidates.put(ref, hashes);
		}

		List<Map<String, String>> combos = cartesian(candidates, MAX_COMBOS);
		ObjectNode bestRegen = null;
		NonCubicBilinearAlgorithm bestRegenAlg = null;
		Map<String, String> bestRegenAssignment = null;
		for (Map<String, String> assignment : combos) {
			Map<String, String> full = new LinkedHashMap<>(forcedAssignment);
			full.putAll(assignment);
			ObjectNode trial = root.deepCopy();
			rewriteRefs(trial.get("lineage"), full);
			Optional<Lineage.Node> node = SchemeIO.readLineage(trial);
			if (node.isEmpty()) return "QUARANTINE";
			NonCubicBilinearAlgorithm replayed;
			try {
				replayed = replayer.replay(node.get());
			} catch (Exception | StackOverflowError e) {
				continue;
			}
			String replayedHash = SchemeIO.contentHash(replayed);
			// Tier 1: bit-exact (or filename-hash7 prefix when that's the only anchor).
			if (exactAnchor != null
					&& (replayedHash.equals(exactAnchor) || replayedHash.startsWith(exactAnchor))) {
				full.forEach(successHint::put);
				if (apply) {
					writeRepaired(p, trial, node.get(), null);
				}
				return "REPINNED_BITEXACT";
			}
			// Tier 2 candidate: same shape + rank, spot-check valid.
			if (bestRegen == null && replayed.n == sn && replayed.m == sm && replayed.p == sp
					&& replayed.r == storedRank
					&& Verifier.passesRandomMatmulSpotCheck(replayed, 2_000, 42L)) {
				bestRegen = trial;
				bestRegenAlg = replayed;
				bestRegenAssignment = full;
				// keep scanning — a later combo might still hit tier 1
			}
		}
		if (bestRegen != null) {
			if (hasMatrices) {
				return "MATERIALISED_CHANGE"; // valid matrices on disk; rewriting them is out of scope
			}
			bestRegenAssignment.forEach(successHint::put);
			String newHash = SchemeIO.contentHash(bestRegenAlg);
			if (storedHash != null) rewriteMap.put(storedHash, newHash);
			if (apply) {
				writeRepaired(p, bestRegen, SchemeIO.readLineage(bestRegen).orElseThrow(), newHash);
			}
			return "REGENERATED";
		}
		return "QUARANTINE";
	}

	private void writeRepaired(Path p, ObjectNode trial, Lineage.Node node, String newHash)
			throws Exception {
		trial.put("lineage_str", Lineage.prettyString(node));
		trial.put("lineage_compact", Lineage.prettyCompact(node));
		if (newHash != null) {
			trial.put("hash", newHash);
			// Content changed: per-content metrics are stale; EnrichSchemeMetrics re-stamps.
			for (String f : STALE_METRIC_FIELDS) trial.remove(f);
		}
		Files.writeString(p, MatrixJsonFormatter.format(trial));
	}

	private String hashOf(Path schemeFile) {
		return candidateHash.computeIfAbsent(schemeFile, f -> {
			try {
				JsonNode root = SchemeIO.parseJson(f.toFile());
				if (root.has("u") || root.has("u_sparse")) {
					return SchemeIO.contentHash(SchemeIO.readBilinear(f.toFile()));
				}
				String stamped = SchemeIO.readHash(root);
				if (stamped != null) return stamped;
				Optional<Lineage.Node> node = SchemeIO.readLineage(root);
				if (node.isEmpty()) return null;
				return SchemeIO.contentHash(replayer.replay(node.get()));
			} catch (Exception | StackOverflowError e) {
				return null;
			}
		});
	}

	/** Rewrite every Atom ref present in {@code assignment} to its new pinned value. */
	private static void rewriteRefs(JsonNode node, Map<String, String> assignment) {
		if (node == null) return;
		if (node.isObject()) {
			ObjectNode o = (ObjectNode) node;
			JsonNode op = o.get("op");
			JsonNode ref = o.get("ref");
			if (op != null && "Atom".equals(op.asString()) && ref != null && ref.isTextual()) {
				String newHash = assignment.get(ref.asString());
				if (newHash != null) {
					String shapePart = ref.asString().substring(0, ref.asString().indexOf('@'));
					o.put("ref", shapePart + "@" + newHash);
				}
			}
			o.properties().forEach(e -> rewriteRefs(e.getValue(), assignment));
		} else if (node.isArray()) {
			node.forEach(child -> rewriteRefs(child, assignment));
		}
	}

	/** Cartesian product of candidate assignments, capped at {@code cap} combos. */
	private static List<Map<String, String>> cartesian(Map<String, List<String>> candidates, int cap) {
		List<Map<String, String>> out = new ArrayList<>();
		out.add(new LinkedHashMap<>());
		for (var e : candidates.entrySet()) {
			List<Map<String, String>> next = new ArrayList<>();
			for (Map<String, String> base : out) {
				for (String h : e.getValue()) {
					if (next.size() >= cap) break;
					Map<String, String> ext = new LinkedHashMap<>(base);
					ext.put(e.getKey(), h);
					next.add(ext);
				}
			}
			out = next;
			if (out.size() >= cap) break;
		}
		return out;
	}
}
