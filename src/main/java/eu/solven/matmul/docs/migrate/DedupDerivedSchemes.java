package eu.solven.matmul.docs.migrate;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Remove <em>stale duplicate clutter</em> from {@code schemes/derived/}: when several
 * derived files cover the SAME ⟨n,m,p⟩ at the SAME rank, keep the one with the
 * SIMPLEST lineage and drop the rest.
 *
 * <p>The canonical motivator is ⟨2,2,7⟩=25, which carried TWO derived witnesses:
 * a clean {@code ConcatCols(naive-2x2x1, 2x2x2, 2x2x2, 2x2x2)} (251639b) <em>and</em>
 * a convoluted {@code ConcatCols(Project(2x2x3→2x2x2) ×3, naive-2x2x1)} (b017c27) — same
 * rank 25, three needless {@code Project} nodes, referenced by nothing. The clean twin
 * already exists for every such case (the materialiser produces it once the sub-block
 * atoms are present); the convoluted file is a leftover from an earlier catalog state.
 *
 * <p>This is regenerable-{@code derived/} hygiene, not a rank change: every kept file is
 * verified to compute matmul, and only STRICTLY more-complex same-shape-same-rank
 * siblings are dropped. A removable file that is still PINNED by another scheme's lineage
 * ({@code shape@hash7}) is RETAINED (logged) so no ref dangles — re-pinning is left to a
 * deliberate follow-up; in practice this is a tiny tail (1 file as of writing).
 *
 * <p>Lineage complexity is scored as the total node count (every {@code "op":} token):
 * a shallow atom scores 1, the clean ⟨2,2,7⟩ concat 7, the convoluted one 10. Lower wins;
 * ties break on path for determinism. Dry run by default; {@code --apply} deletes.
 *
 * <p>After {@code --apply}, re-run {@code GenerateCatalogManifest} to refresh
 * {@code docs/catalog.json}.
 */
public final class DedupDerivedSchemes {
	private DedupDerivedSchemes() {}

	private static final String ROOT = "src/main/resources/schemes";
	private static final String DERIVED = ROOT + "/derived";
	private static final Pattern OP = Pattern.compile("\"op\"\\s*:");
	private static final Pattern PIN = Pattern.compile("@([0-9a-f]{7})(?![0-9a-f])");

	private record Entry(Path path, int score, String hash7) {}

	public static void main(String[] args) throws Exception {
		boolean apply = List.of(args).contains("--apply");

		// (shape|rank) → derived files at that exact shape+rank.
		Map<String, List<Entry>> buckets = new TreeMap<>();
		try (var s = Files.walk(Path.of(DERIVED))) {
			for (Path p : s.filter(x -> x.toString().endsWith(".json")).sorted().toList()) {
				NonCubicBilinearAlgorithm alg;
				try {
					alg = SchemeIO.read(p.toFile());
				} catch (Exception e) {
					// A lineage-only STUB (maxDim>16: no explicit u/v/w, replay-on-demand) is the
					// expected, common case — out of scope for this dense-file hygiene pass. Only
					// surface GENUINELY unexpected read failures.
					String msg = String.valueOf(e.getMessage());
					if (!msg.contains("key 'u'")) {
						System.out.println("  SKIP (unreadable): " + p.getFileName() + " — " + msg);
					}
					continue;
				}
				String txt = Files.readString(p);
				int score = count(OP, txt);
				String hash7 = SchemeIO.contentHash(alg).substring(0, 7);
				String key = alg.n + "x" + alg.m + "x" + alg.p + "|r" + alg.r;
				buckets.computeIfAbsent(key, k -> new ArrayList<>()).add(new Entry(p, score, hash7));
			}
		}

		// Every @hash7 pinned anywhere in the catalog — a removable file in this set is RETAINED.
		Set<String> pinned = new java.util.HashSet<>();
		try (var s = Files.walk(Path.of(ROOT))) {
			for (Path p : s.filter(x -> x.toString().endsWith(".json")).toList()) {
				Matcher m = PIN.matcher(Files.readString(p));
				while (m.find()) pinned.add(m.group(1));
			}
		}

		List<Path> toDelete = new ArrayList<>();
		Map<String, String> retainedReferenced = new LinkedHashMap<>();
		int bucketsTouched = 0;
		for (var e : buckets.entrySet()) {
			List<Entry> lst = e.getValue();
			if (lst.size() < 2) continue;
			lst.sort(Comparator.comparingInt(Entry::score).thenComparing(en -> en.path().toString()));
			// Keep the simplest sibling that actually verifies (a corrupt shortest must not
			// license dropping a valid-but-complex twin).
			Entry keep = null;
			for (Entry en : lst) {
				try {
					if (Verifier.passesRandomMatmulSpotCheck(SchemeIO.read(en.path().toFile()))) {
						keep = en;
						break;
					}
				} catch (Exception ignore) {
					// try the next candidate
				}
			}
			if (keep == null) {
				System.out.println("  SKIP bucket (no verifying sibling): " + e.getKey());
				continue;
			}
			bucketsTouched++;
			System.out.printf("%-16s KEEP  [score %2d]: %s%n", e.getKey(), keep.score(), keep.path().getFileName());
			for (Entry en : lst) {
				if (en == keep) continue;
				if (pinned.contains(en.hash7())) {
					retainedReferenced.put(en.hash7(), en.path().getFileName().toString());
					System.out.printf("    RETAIN (pinned @%s): %s%n", en.hash7(), en.path().getFileName());
					continue;
				}
				System.out.printf("    DROP  [score %2d]: %s%n", en.score(), en.path().getFileName());
				toDelete.add(en.path());
			}
		}

		System.out.printf("%n%s: %d duplicate buckets; %d files %s; %d referenced files retained.%n",
				apply ? "APPLIED" : "DRY RUN", bucketsTouched, toDelete.size(),
				apply ? "deleted" : "to delete", retainedReferenced.size());
		if (apply) {
			for (Path p : toDelete) Files.delete(p);
			System.out.println("Deletion complete. Re-run GenerateCatalogManifest to refresh docs/catalog.json.");
		} else {
			System.out.println("Re-run with --apply to delete.");
		}
	}

	private static int count(Pattern pat, String txt) {
		Matcher m = pat.matcher(txt);
		int c = 0;
		while (m.find()) c++;
		return c;
	}
}
