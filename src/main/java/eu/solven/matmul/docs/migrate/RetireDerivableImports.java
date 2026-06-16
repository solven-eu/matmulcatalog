package eu.solven.matmul.docs.migrate;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import eu.solven.matmul.catalog.SchemeIO;

/**
 * Retire the FMM-Lille imports that are now <em>derived</em> at the same rank
 * (the 4 serendipitous shapes ⟨8,9,9⟩=430, ⟨6,8,9⟩=296, ⟨4,8,12⟩=272,
 * ⟨8,8,12⟩=504, each reproduced by the bud-serendipitous engine). For each:
 *
 * <ol>
 *   <li><b>Re-pin</b> every {@code shape@<fmmHash7>} lineage ref to
 *       {@code shape@<derivedHash7>} — so no pinned ref dangles after the drop
 *       (the replayer <em>would</em> fall back to the bare shape, but an explicit
 *       re-pin is cleaner and avoids the "bud-poor sibling" footgun the resolver
 *       warns about). Pins verified to be 7-char content-hash prefixes; cascade
 *       depth is 0 (no dependent is itself pinned-to), so a single pass suffices —
 *       we still re-scan and fail loud if any {@code @<fmmHash7>} survives.</li>
 *   <li><b>Drop</b> the redundant fmm import file.</li>
 * </ol>
 *
 * <p>Filenames are pure labels (content-driven resolution via
 * {@link SchemeIO#contentHash}), so the edited dependents need no rename. Dry run
 * by default; {@code --apply} writes.</p>
 */
public final class RetireDerivableImports {
	private RetireDerivableImports() {}

	private static final String ROOT = "src/main/resources/schemes";

	/** {shape, fmm-filename-substring, derived-stub-filename-substring}. */
	private static final String[][] RETIRE = {
			{ "8x9x9", "8x9x9-r430-fmm_lille", "derived_recursive-8x9x9_m430" },
			{ "6x8x9", "6x8x9-r296-fmm_lille", "derived_recursive-6x8x9_m296" },
			{ "4x8x12", "4x8x12-r272-fmm_lille", "derived_recursive-4x8x12_m272" },
			{ "8x8x12", "8x8x12-r504-fmm_lille", "derived_recursive-8x8x12_m504" } };

	public static void main(String[] args) throws Exception {
		boolean apply = List.of(args).contains("--apply");

		Map<String, String> remap = new LinkedHashMap<>();   // "shape@old7" -> "shape@new7"
		List<Path> fmmToDelete = new ArrayList<>();
		for (String[] r : RETIRE) {
			Path fmm = findOne(r[1]);
			Path der = findOne(r[2]);
			String oldH = SchemeIO.contentHash(SchemeIO.read(fmm.toFile())).substring(0, 7);
			String newH = SchemeIO.contentHash(SchemeIO.read(der.toFile())).substring(0, 7);
			remap.put(r[0] + "@" + oldH, r[0] + "@" + newH);
			fmmToDelete.add(fmm);
			System.out.printf("%-7s re-pin %s@%s → %s@%s   (drop %s)%n",
					r[0], r[0], oldH, r[0], newH, fmm.getFileName());
		}

		List<Path> all;
		try (var s = Files.walk(Path.of(ROOT))) {
			all = s.filter(p -> p.toString().endsWith(".json")).sorted().toList();
		}
		int filesEdited = 0, refsRewritten = 0;
		for (Path p : all) {
			String txt = Files.readString(p), out = txt;
			for (var e : remap.entrySet()) {
				// match the 7-char pin only when NOT followed by another hex digit,
				// so a (hypothetical) longer hash is never partially clobbered.
				Pattern pat = Pattern.compile(Pattern.quote(e.getKey()) + "(?![0-9a-f])");
				Matcher m = pat.matcher(out);
				int c = 0;
				StringBuilder sb = new StringBuilder();
				while (m.find()) { m.appendReplacement(sb, Matcher.quoteReplacement(e.getValue())); c++; }
				m.appendTail(sb);
				if (c > 0) { out = sb.toString(); refsRewritten += c; }
			}
			if (!out.equals(txt)) {
				filesEdited++;
				if (apply) Files.writeString(p, out);
			}
		}
		System.out.printf("%n%s: %d refs re-pinned across %d files; %d fmm imports %s%n",
				apply ? "APPLIED" : "DRY RUN", refsRewritten, filesEdited, fmmToDelete.size(),
				apply ? "dropped" : "to drop");

		if (apply) {
			for (Path p : fmmToDelete) Files.delete(p);
			// fail loud if any old pin survived (would dangle)
			int dangling = 0;
			try (var s = Files.walk(Path.of(ROOT))) {
				for (Path p : s.filter(x -> x.toString().endsWith(".json")).toList()) {
					String t = Files.readString(p);
					for (String oldRef : remap.keySet()) {
						if (Pattern.compile(Pattern.quote(oldRef) + "(?![0-9a-f])").matcher(t).find()) {
							System.out.println("  STILL DANGLING: " + oldRef + " in " + p.getFileName());
							dangling++;
						}
					}
				}
			}
			if (dangling > 0) throw new IllegalStateException(dangling + " dangling refs remain — investigate");
			System.out.println("Post-check: 0 dangling refs. Retirement complete.");
		} else {
			System.out.println("Re-run with --apply to write + drop.");
		}
	}

	private static Path findOne(String filenameSubstring) throws Exception {
		try (var s = Files.walk(Path.of(ROOT))) {
			List<Path> hits = s.filter(p -> p.getFileName().toString().contains(filenameSubstring)
					&& p.toString().endsWith(".json")).toList();
			if (hits.size() != 1) {
				throw new IllegalStateException("expected exactly 1 file matching '" + filenameSubstring
						+ "', found " + hits.size() + ": " + hits);
			}
			return hits.get(0);
		}
	}
}
