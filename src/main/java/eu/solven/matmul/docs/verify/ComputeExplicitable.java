package eu.solven.matmul.docs.verify;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.SchemeIO;

import lombok.extern.slf4j.Slf4j;

import tools.jackson.databind.JsonNode;

/**
 * Computes the {@code explicitable} property for every catalog scheme: is the entry a
 * precise EXPLICIT scheme (matrices, or a lineage that bottoms out entirely in precise
 * leaves) — as opposed to a CITED BOUND (a best-at-shape {@code @sota}/{@code -direct}/
 * bare leaf, for which no exact construction is pinned)?
 *
 * <p>{@code explicitable = true} iff:</p>
 * <ul>
 *   <li>the file carries explicit matrices ({@code u}/{@code u_sparse}) — trivially explicit; OR</li>
 *   <li>it is a stub whose every lineage leaf is explicitable: {@code @naive} (elementary),
 *       or a {@code shape@hash} / named-canonical ref to an explicitable scheme (recursively).</li>
 * </ul>
 * <p>{@code explicitable = false} iff any leaf is a best-at-shape CITED BOUND
 * ({@code @sota} / {@code -direct} / bare {@code NxMxP}), or a pin dangles, or a cycle is hit.
 * Note this is NOT {@code buildable}: an HK71 cited bound has a buildable FORMULA yet is
 * not an explicit scheme — it would be {@code explicitable=false} until materialised.</p>
 *
 * <p>READ-ONLY: reports the distribution and writes {@code target/explicitable.tsv}
 * (relative-file &lt;tab&gt; true|false). Stamping JSON / the manifest is a separate step.</p>
 *
 * <pre>mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.verify.ComputeExplicitable</pre>
 */
@Slf4j
public final class ComputeExplicitable {

	private ComputeExplicitable() {}

	private static final Path SCHEMES_ROOT = Path.of("src/main/resources/schemes");
	private static final Pattern FN = Pattern.compile("^(\\d+x\\d+x\\d+)-r\\d+-.+-([0-9a-f]{4,})$");

	/**
	 * Reusable: compute {@code explicitable} for every file (DAG, memoised). Used by
	 * {@code GenerateCatalogManifest} to emit the attribute, and by {@link #main} to report.
	 */
	public static Map<Path, Boolean> computeAll(List<Path> files) {
		Map<String, Path> byPinnedKey = new HashMap<>();
		Map<String, Path> byCanonKey = new HashMap<>();
		for (Path p : files) {
			String stem = p.getFileName().toString().replaceFirst("\\.json$", "");
			byCanonKey.putIfAbsent(Lineage.canonicalKey(stem), p);
			Matcher m = FN.matcher(stem);
			if (m.matches()) {
				String hash7 = m.group(2).substring(0, Math.min(7, m.group(2).length()));
				byPinnedKey.put(m.group(1) + "@" + hash7, p);
			}
		}
		Explicitor ex = new Explicitor(byPinnedKey, byCanonKey, new HashMap<>(), new java.util.HashSet<>());
		Map<Path, Boolean> out = new HashMap<>();
		for (Path p : files) out.put(p, ex.explicitable(p));
		return out;
	}

	public static void main(String[] args) throws IOException {
		List<Path> files;
		try (Stream<Path> w = Files.walk(SCHEMES_ROOT)) {
			files = w.filter(p -> p.toString().endsWith(".json")).sorted().toList();
		}
		log.info("Computing explicitable over {} scheme files", files.size());

		Map<Path, Boolean> result = computeAll(files);
		int explicit = 0, citedBound = 0;
		StringBuilder tsv = new StringBuilder();
		for (Path p : files) {
			boolean v = result.getOrDefault(p, false);
			tsv.append(SCHEMES_ROOT.relativize(p)).append('\t').append(v).append('\n');
			if (v) explicit++; else citedBound++;
		}
		Path out = Path.of("target", "explicitable.tsv");
		Files.createDirectories(out.getParent());
		Files.writeString(out, tsv.toString());

		// Split by band: explicit schemes (≤16, matrices) vs stubs (>16, lineage-only).
		log.info("================ EXPLICITABLE ================");
		log.info("  total                 : {}", files.size());
		log.info("  explicitable = true   : {}  (precise: matrices, or fully-pinned lineage)", explicit);
		log.info("  explicitable = false  : {}  (cited bounds: @sota/-direct/bare or dangling leaf)", citedBound);
		log.info("  → target/explicitable.tsv");
	}

	private enum Reason { NONE, CITED, DANGLING, CYCLE }

	private static final class Explicitor {
		final Map<String, Path> byPinnedKey;
		final Map<String, Path> byCanonKey;
		final Map<Path, Boolean> memo;
		final java.util.Set<Path> onPath;
		Reason reason = Reason.NONE;

		Explicitor(Map<String, Path> byPinnedKey, Map<String, Path> byCanonKey,
				Map<Path, Boolean> memo, java.util.Set<Path> onPath) {
			this.byPinnedKey = byPinnedKey;
			this.byCanonKey = byCanonKey;
			this.memo = memo;
			this.onPath = onPath;
		}

		boolean explicitable(Path p) {
			Boolean cached = memo.get(p);
			if (cached != null) return cached;
			if (!onPath.add(p)) { reason = Reason.CYCLE; return false; } // cycle → not explicitable
			boolean result;
			try {
				JsonNode root = SchemeIO.parseJson(p.toFile());
				if (!SchemeIO.isStub(root)) {
					// Any concrete representation on disk is explicit — bilinear (u/v/w),
					// non_bilinear (Waksman/Rosowski u_a/u_b/…), complex, or reduced. Only a
					// scheme_type=="stub" is lineage-only.
					result = true;
				} else {
					JsonNode lineage = root.get("lineage");
					if (lineage == null || lineage.isNull()) {
						result = false;                              // stub with no lineage → a bound
					} else {
						List<String> refs = new ArrayList<>();
						collectRefs(lineage, refs);
						result = true;
						for (String ref : refs) {
							if (!leafExplicitable(ref)) { result = false; break; }
						}
					}
				}
			} catch (Exception e) {
				result = false;                                      // unreadable → treat as a bound
			}
			onPath.remove(p);
			memo.put(p, result);
			return result;
		}

		/** Is a single lineage leaf ref precise+explicitable? {@code @naive} yes; a pinned/named
		 *  ref recurses into its file; a bare / {@code -direct} / {@code @sota} ref is a cited bound. */
		private boolean leafExplicitable(String ref) {
			int at = ref.indexOf('@');
			if (at > 0) {
				String tail = ref.substring(at + 1);
				if (tail.equals("naive")) return true;
				if (tail.equals("sota")) { reason = Reason.CITED; return false; }
				String shape = ref.substring(0, at);
				String hash7 = tail.substring(0, Math.min(7, tail.length()));
				Path t = byPinnedKey.get(shape + "@" + hash7);
				if (t == null) { reason = Reason.DANGLING; return false; }
				return explicitable(t);
			}
			// marker-suffixed or prefixed naive (naive-NxMxP / NxMxP-naive)
			if (ref.contains("naive")) return true;
			// bare shape or "-direct"  → best-at-shape cited bound
			if (ref.endsWith("-direct") || ref.matches("\\d+x\\d+x\\d+")) {
				reason = Reason.CITED;
				return false;
			}
			// named canonical ref (perminov_…, alphatensor_…) → recurse into the exact file
			Path named = byCanonKey.get(Lineage.canonicalKey(ref));
			if (named != null) return explicitable(named);
			reason = Reason.CITED;                                   // unknown ref form → conservative
			return false;
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
	}
}
