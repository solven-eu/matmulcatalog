package eu.solven.matmul.docs.migrate;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.SchemeIO;
import lombok.extern.slf4j.Slf4j;
import tools.jackson.databind.JsonNode;

/**
 * Regenerates the {@code -<hash7>} filename token for every scheme whose content
 * hash drifted under the double→string coefficient migration — i.e. the
 * <strong>rational</strong> schemes (coefficients now stored exactly as
 * {@code "1/4"}) and the <strong>complex</strong> schemes ({@code [re,im]} pairs).
 *
 * <p>{@link SchemeIO#contentHash} now hashes the exact canonical token for
 * rational schemes and {@link SchemeIO#contentHashComplexJson} hashes the verbatim
 * {@code [re,im]} list for complex ones, so recomputing yields a hash stable
 * against the representation change. Pure-integer schemes keep the legacy
 * {@code round(v·1e6)} basis byte-for-byte, so their hashes (and filenames) do not
 * move — this driver leaves them untouched.</p>
 *
 * <p>Pure rename (no content edit), so git rename-detection keeps history. Default
 * is DRY-RUN; pass {@code --execute}. Aborts on any 7-hex collision. Re-run the
 * manifest afterwards.</p>
 *
 * <pre>mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.RehashRationalsComplex [-Dexec.args=--execute]</pre>
 */
@Slf4j
public final class RehashRationalsComplex {

	private RehashRationalsComplex() {}

	private static final Path ROOT = Path.of("src/main/resources/schemes");
	/** {@code <shape>-r<rank>-<note>-<hash7>.json}; note is greedy so it may contain '-'. */
	private static final Pattern NAME = Pattern.compile("^(\\d+x\\d+x\\d+)-r(\\d+)-(.+)-([0-9a-f]{7})\\.json$");

	public static void main(String[] args) throws Exception {
		boolean execute = List.of(args).contains("--execute");
		List<Path> files;
		try (Stream<Path> walk = Files.walk(ROOT)) {
			files = walk.filter(p -> p.toString().endsWith(".json")).sorted().toList();
		}

		AtomicInteger scanned = new AtomicInteger();
		AtomicInteger errors = new AtomicInteger();
		AtomicInteger complexRenamed = new AtomicInteger();
		AtomicInteger rationalRenamed = new AtomicInteger();
		List<Path[]> plan = new ArrayList<>();           // {old, new}
		Map<Path, Long> newNameCount = new HashMap<>();  // collision detection

		for (Path f : files) {
			scanned.incrementAndGet();
			String name = f.getFileName().toString();
			Matcher mt = NAME.matcher(name);
			if (!mt.matches()) {
				continue; // non-conforming name (handled elsewhere)
			}
			String shape = mt.group(1), note = mt.group(3), curH7 = mt.group(4);
			int rank = Integer.parseInt(mt.group(2));
			try {
				JsonNode root = SchemeIO.parseJson(f.toFile());
				if (SchemeIO.isStub(root) || SchemeIO.isNonBilinear(root)) {
					continue;
				}
				boolean complex = root.path("complex").asBoolean(false) || SchemeIO.isComplex(root);
				String newHash;
				if (complex) {
					newHash = SchemeIO.contentHashComplexJson(root);
				} else {
					NonCubicBilinearAlgorithm alg =
							SchemeIO.isReduced(root) ? SchemeIO.readReduced(root) : SchemeIO.read(root);
					newHash = SchemeIO.contentHash(alg);
				}
				String newH7 = newHash.substring(0, 7);
				if (!newH7.equals(curH7)) {
					Path dst = f.resolveSibling(shape + "-r" + rank + "-" + note + "-" + newH7 + ".json");
					plan.add(new Path[] { f, dst });
					newNameCount.merge(dst, 1L, Long::sum);
					if (complex) {
						complexRenamed.incrementAndGet();
					} else {
						rationalRenamed.incrementAndGet();
					}
				}
			} catch (Exception e) {
				errors.incrementAndGet();
				log.warn("failed on {}: {}", name, e.getMessage());
			}
		}

		// Collision guard: two distinct sources must never map to one new name.
		List<Path> collisions = newNameCount.entrySet().stream()
				.filter(e -> e.getValue() > 1).map(Map.Entry::getKey).toList();

		// Dump the full old→new plan (with shape@hash7 refs) for ref-overlap analysis.
		StringBuilder dump = new StringBuilder("oldName\tnewName\tshape\toldRef\tnewRef\n");
		for (Path[] pr : plan) {
			String on = pr[0].getFileName().toString(), nn = pr[1].getFileName().toString();
			Matcher om = NAME.matcher(on), nm = NAME.matcher(nn);
			om.matches();
			nm.matches();
			dump.append(on).append('\t').append(nn).append('\t').append(om.group(1)).append('\t')
					.append(om.group(1)).append('@').append(om.group(4)).append('\t')
					.append(nm.group(1)).append('@').append(nm.group(4)).append('\n');
		}
		Files.writeString(Path.of("target/rehash-plan.tsv"), dump.toString());

		log.info("=== PLAN (rehash rationals + complex) ===");
		log.info("scanned={} to-rename={} (rational={} complex={}) errors={} collisions={}",
				scanned.get(), plan.size(), rationalRenamed.get(), complexRenamed.get(),
				errors.get(), collisions.size());
		plan.stream().limit(12).forEach(pr ->
				log.info("  {}\n    -> {}", pr[0].getFileName(), pr[1].getFileName()));

		if (!collisions.isEmpty()) {
			collisions.forEach(c -> log.error("COLLISION target: {}", c.getFileName()));
			throw new IllegalStateException("aborting: " + collisions.size() + " hash-7 collisions");
		}

		// Lineage refs are `shape@<hash7>` content-addresses; a re-hashed scheme's
		// old ref would dangle. Build old→new ref map and rewrite every occurrence
		// across the catalog IN LOCKSTEP. (Rewriting a referrer edits only its
		// lineage text, never its U/V/W, so its own content hash — hence its
		// filename — is unaffected; no cascade.)
		Map<String, String> refMap = new HashMap<>();
		for (Path[] pr : plan) {
			Matcher om = NAME.matcher(pr[0].getFileName().toString());
			Matcher nm = NAME.matcher(pr[1].getFileName().toString());
			om.matches();
			nm.matches();
			refMap.put(om.group(1) + "@" + om.group(4), nm.group(1) + "@" + nm.group(4));
		}
		Pattern ref = Pattern.compile("\\d+x\\d+x\\d+@[0-9a-f]{7}");
		int refFilesChanged = 0, refsRewritten = 0;
		for (Path f : files) {
			String text = Files.readString(f);
			Matcher m = ref.matcher(text);
			StringBuilder out = null;
			int last = 0, local = 0;
			while (m.find()) {
				String repl = refMap.get(m.group());
				if (repl != null) {
					if (out == null) {
						out = new StringBuilder();
					}
					out.append(text, last, m.start()).append(repl);
					last = m.end();
					local++;
				}
			}
			if (out != null) {
				out.append(text, last, text.length());
				refsRewritten += local;
				refFilesChanged++;
				if (execute) {
					Files.writeString(f, out.toString());
				}
			}
		}
		log.info("lineage refs to rewrite: {} occurrences in {} files", refsRewritten, refFilesChanged);

		if (!execute) {
			log.info("(DRY-RUN — nothing written. Pass --execute to perform.)");
			return;
		}
		int done = 0;
		for (Path[] pr : plan) {
			Files.move(pr[0], pr[1]);
			if (++done % 500 == 0) {
				log.info("[progress] {}/{} renamed", done, plan.size());
			}
		}
		log.info("RehashRationalsComplex: renamed {} files, rewrote {} refs in {} files.",
				plan.size(), refsRewritten, refFilesChanged);
	}
}
