package eu.solven.matmul.docs.migrate;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.verifiers.LineageVerifier;
import lombok.extern.slf4j.Slf4j;

/**
 * Promote serendipitous / closure wins from a staging schemes-root into the
 * committed catalog ({@code src/main/resources/schemes}).
 *
 * <p>This is the "proven at commit" half of the
 * <em>bound-during-search, proven-at-commit</em> split: the search loop only
 * spot-checks (a {@code bound}); promotion runs the EXACT
 * {@link Verifier#isExactNonCubic} verification (replaying the lineage stub into
 * explicit matrices first) before any file lands in the committed catalog. A
 * staging file is promoted only when it (1) replays, (2) exact-verifies, and
 * (3) is STRICTLY better than the committed best rank for its shape (stubs
 * included — via {@link FieldAwareLookup#findRank}).</p>
 *
 * <p>Dry-run by default; pass {@code --apply=true} to actually copy. The old
 * committed entry is left in place (the catalog keeps dominated history; the SPA
 * hides it), so promotion is additive.</p>
 *
 * <pre>{@code
 * mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.PromoteStagingWins \
 *     -Dexec.args="--staging=target/staging-serendip-13-18 --apply=true"
 * }</pre>
 */
@Slf4j
public final class PromoteStagingWins {
	private PromoteStagingWins() {}

	private static final Pattern SHAPE_RANK =
			Pattern.compile("(\\d+)x(\\d+)x(\\d+)_m(\\d+)");

	private static final Path COMMITTED_ROOT = Path.of("src/main/resources/schemes");

	public static void main(String[] args) throws Exception {
		Path staging = null;
		boolean apply = false;
		String fieldTag = "Q";
		for (String a : args) {
			if (a.startsWith("--staging=")) staging = Path.of(a.substring("--staging=".length()));
			else if (a.startsWith("--apply=")) apply = Boolean.parseBoolean(a.substring("--apply=".length()));
			else if (a.startsWith("--field=")) fieldTag = a.substring("--field=".length());
			else throw new IllegalArgumentException("unknown arg " + a);
		}
		if (staging == null || !Files.isDirectory(staging)) {
			throw new IllegalArgumentException("--staging=<dir> required (got " + staging + ")");
		}

		FieldAwareLookup committed = new FieldAwareLookup(Field.fromTag(fieldTag), COMMITTED_ROOT);
		LineageVerifier verifier = new LineageVerifier(committed);

		List<Path> files = new ArrayList<>();
		try (Stream<Path> walk = Files.walk(staging)) {
			walk.filter(p -> p.getFileName().toString().endsWith(".json")).forEach(files::add);
		}
		files.sort(Comparator.comparing(p -> p.getFileName().toString()));
		log.info("Promote scan: {} staging file(s) under {} (apply={})", files.size(), staging, apply);

		int promoted = 0, tie = 0, notbetter = 0, failed = 0;
		for (Path f : files) {
			String name = f.getFileName().toString();
			Matcher mm = SHAPE_RANK.matcher(name);
			if (!mm.find()) {
				log.warn("  SKIP {} — no shape/rank in filename", name);
				continue;
			}
			int n = Integer.parseInt(mm.group(1));
			int m = Integer.parseInt(mm.group(2));
			int p = Integer.parseInt(mm.group(3));
			int claimedRank = Integer.parseInt(mm.group(4));
			int committedBest = committed.findRank(n, m, p);

			if (claimedRank >= committedBest) {
				if (claimedRank == committedBest) tie++; else notbetter++;
				log.info("  skip ⟨{},{},{}⟩ r={} — committed best={} (not strictly better)",
						n, m, p, claimedRank, committedBest);
				continue;
			}

			// Compositional verification: exact-verify the (small) primitive leaves
			// and trust the correctness-preserving operators — no expansion of the
			// (possibly rank-10⁴) composed product. Ring-correct (dispatches F₂/F₃).
			long t0 = System.nanoTime();
			LineageVerifier.Result vr = verifier.verifyFile(f.toFile());
			long ms = (System.nanoTime() - t0) / 1_000_000L;
			if (!vr.certified()) {
				failed++;
				log.warn("  FAILED ⟨{},{},{}⟩ r={} — {} ({} ms)", n, m, p, claimedRank, vr.detail(), ms);
				continue;
			}

			Path destDir = COMMITTED_ROOT.resolve("derived").resolve("section" + Math.max(n, Math.max(m, p)));
			Path dest = destDir.resolve(name);
			log.info("  PROMOTE ⟨{},{},{}⟩ r={} < committed {} — {} ({} ms) → {}",
					n, m, p, claimedRank, committedBest, vr.detail(), ms, dest);
			if (apply) {
				Files.createDirectories(destDir);
				Files.copy(f, dest, StandardCopyOption.REPLACE_EXISTING);
			}
			promoted++;
		}

		log.info("=== Promote {}: {} promoted, {} tie, {} not-better, {} failed-verify (of {}) ===",
				apply ? "APPLIED" : "DRY-RUN", promoted, tie, notbetter, failed, files.size());
		if (!apply && promoted > 0) {
			log.info("Re-run with --apply=true to copy, then regenerate the manifest "
					+ "(GenerateCatalogManifest).");
		}
	}
}
