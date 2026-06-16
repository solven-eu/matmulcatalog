package eu.solven.matmul.docs.migrate;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Optional;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Stream;

import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.SchemeIO;

import lombok.extern.slf4j.Slf4j;

/**
 * One-shot backfill: read every scheme JSON under
 * {@code src/main/resources/schemes/}, compute the new
 * {@code lineage_compact} string from the existing {@code lineage} DAG,
 * and insert it into the file as a sibling of {@code lineage_str}.
 *
 * <p>Idempotent: files that already have {@code lineage_compact} are
 * skipped. Files without any lineage (older schemes pre-dating the
 * lineage layer) are skipped silently.
 *
 * <p>The on-disk edit is text-based to preserve the original formatting
 * — we look for the {@code "lineage_str": "…"} line and insert a new
 * {@code "lineage_compact": "…"} line immediately after it. Parsing +
 * re-serialising the whole JSON would rewrite the dense factor-matrix
 * block and produce a giant diff for no useful change.
 *
 * <p>Invoke:
 * <pre>
 *   mvn -q exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.BackfillLineageCompact
 *   # or with --dry-run to count without writing:
 *   mvn -q exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.BackfillLineageCompact \
 *       -Dexec.args="--dry-run"
 * </pre>
 */
@Slf4j
public final class BackfillLineageCompact {

	private static final Path SCHEMES_ROOT = Path.of("src/main/resources/schemes");

	private BackfillLineageCompact() {}

	public static void main(String[] args) throws IOException {
		boolean dryRunMutable = false;
		Path rootMutable = SCHEMES_ROOT;
		for (String a : args) {
			if ("--dry-run".equals(a)) dryRunMutable = true;
			else if (a.startsWith("--root=")) rootMutable = Path.of(a.substring("--root=".length()));
			else if ("--help".equals(a) || "-h".equals(a)) {
				System.out.println("Usage: BackfillLineageCompact [--dry-run] [--root=PATH]");
				return;
			} else throw new IllegalArgumentException("Unknown arg: " + a);
		}
		final boolean dryRun = dryRunMutable;
		final Path root = rootMutable;

		AtomicInteger seen = new AtomicInteger();
		AtomicInteger updated = new AtomicInteger();
		AtomicInteger alreadyHad = new AtomicInteger();
		AtomicInteger noLineage = new AtomicInteger();
		AtomicInteger failed = new AtomicInteger();

		try (Stream<Path> walk = Files.walk(root)) {
			walk.filter(p -> p.toString().endsWith(".json"))
					.filter(Files::isRegularFile)
					.forEach(p -> {
						seen.incrementAndGet();
						try {
							Outcome o = processOne(p, dryRun);
							switch (o) {
								case UPDATED -> updated.incrementAndGet();
								case ALREADY_HAD -> alreadyHad.incrementAndGet();
								case NO_LINEAGE -> noLineage.incrementAndGet();
							}
						} catch (RuntimeException | IOException e) {
							failed.incrementAndGet();
							log.warn("failed on {}: {}", p, e.getMessage());
						}
					});
		}

		log.info("Backfill done (dry-run={}): seen={}, updated={}, already-had={}, no-lineage={}, failed={}",
				dryRun, seen.get(), updated.get(), alreadyHad.get(),
				noLineage.get(), failed.get());
	}

	private enum Outcome { UPDATED, ALREADY_HAD, NO_LINEAGE }

	private static Outcome processOne(Path p, boolean dryRun) throws IOException {
		String content = Files.readString(p);
		if (content.contains("\"lineage_compact\"")) {
			return Outcome.ALREADY_HAD;
		}
		Optional<Lineage.Node> nodeOpt = SchemeIO.readLineage(p.toFile());
		if (nodeOpt.isEmpty()) {
			return Outcome.NO_LINEAGE;
		}
		String compact = Lineage.prettyCompact(nodeOpt.get());
		String escaped = compact.replace("\\", "\\\\").replace("\"", "\\\"");
		String updated = insertCompactLine(content, escaped);
		if (updated == null) {
			// Couldn't locate a place to insert — log and skip rather than corrupt.
			log.warn("couldn't locate insertion point in {} (no lineage_str line?)", p);
			return Outcome.NO_LINEAGE;
		}
		if (!dryRun) {
			// Route the string-surgery result through the single shared formatter
			// so the on-disk file stays canonical (and is validated as JSON).
			Files.writeString(p, eu.solven.matmul.catalog.MatrixJsonFormatter.format(updated));
		}
		return Outcome.UPDATED;
	}

	/**
	 * Insert {@code "lineage_compact": "…"} on its own line. Preferred
	 * placement is immediately AFTER the {@code "lineage_str"} line;
	 * fallback (for older files that have {@code "lineage"} but no
	 * {@code "lineage_str"}) is immediately BEFORE the {@code "lineage"}
	 * key line. Returns {@code null} if neither anchor is found. The
	 * insertion preserves the file's leading whitespace per line.
	 */
	private static String insertCompactLine(String content, String escapedCompact) {
		int strKeyIdx = content.indexOf("\"lineage_str\"");
		if (strKeyIdx >= 0) {
			int lineEnd = content.indexOf('\n', strKeyIdx);
			if (lineEnd < 0) return null;
			int lineStart = content.lastIndexOf('\n', strKeyIdx) + 1;
			String indent = content.substring(lineStart, strKeyIdx);
			String prefix = content.substring(0, lineEnd + 1);
			String suffix = content.substring(lineEnd + 1);
			String newLine = indent + "\"lineage_compact\": \"" + escapedCompact + "\",\n";
			return prefix + newLine + suffix;
		}
		// Fallback: insert before the "lineage" key. The key may be
		// printed as "lineage" :  with a leading space (Jackson default
		// pretty-print) or "lineage": without one — accept both.
		int linKeyIdx = content.indexOf("\"lineage\"");
		if (linKeyIdx < 0) return null;
		int lineStart = content.lastIndexOf('\n', linKeyIdx) + 1;
		String indent = content.substring(lineStart, linKeyIdx);
		String prefix = content.substring(0, lineStart);
		String suffix = content.substring(lineStart);
		String newLine = indent + "\"lineage_compact\": \"" + escapedCompact + "\",\n";
		return prefix + newLine + suffix;
	}
}
