package eu.solven.matmul.docs.explore;

import java.io.File;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.recombination.BlockSplitSearch;
import eu.solven.matmul.search.RecursiveClosureSota;
import eu.solven.matmul.search.RecursiveMaterialiser;
import eu.solven.matmul.util.ProgressMonitor;

import lombok.extern.slf4j.Slf4j;
import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;

/**
 * Targeted downward-projection sweep over exactly the shapes where FMM-Lille beats
 * our on-disk catalog. For each WORSE shape it runs
 * {@link RecursiveMaterialiser#projectInto} — restrict a slightly-larger held
 * parent (incl. the high-projection-margin PanTA cube) down to the shape + DCE —
 * persisting any strict improvement. Largest max-axis first, so a freshly-projected
 * parent is available when its smaller child projects in the same pass.
 *
 * <p>This is the actionable half of the FMM gap (see DISCOVERIES_PENDING 2026-06-08):
 * our coordinate projection reaches FMM's *real published* ranks (e.g. ⟨27,28,28⟩
 * 11718 → 10442); the closure simply had not been run over the 17–32 band.</p>
 *
 * <pre>MAVEN_OPTS="-Xmx4g -XX:+HeapDumpOnOutOfMemoryError -XX:HeapDumpPath=target/oom-dumps/" \
 *   mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.explore.ProjectFmmGaps -Dexec.args="2"</pre>
 * (arg = max passes; default 2.)
 */
@Slf4j
public final class ProjectFmmGaps {
	private ProjectFmmGaps() {}

	private static final int MAX_DIM = eu.solven.matmul.catalog.CatalogLimits.MAX_DIM;
	private static final Path SCHEMES_ROOT = Path.of("src/main/resources/schemes");

	public static void main(String[] args) throws Exception {
		// Args (positional [limit] [maxPasses] still accepted for back-compat):
		//   --limit=N          cap the shape count (WORSE mode only); ≤0 = all
		//   --passes=N         max projection passes (default 2)
		//   --shapes=NxMxP,…   project exactly these shapes (comma-separated; x or , inner sep)
		//   --shapes-file=PATH project shapes listed in a file (one "NxMxP" / "n,m,p" per line,
		//                      blank lines and #comments ignored)
		// With no --shapes/--shapes-file, the FMM-WORSE list (from fmm-cross-check.md) is used.
		int limit = 0, maxPasses = 2, pos = 0;
		List<int[]> explicit = null;
		for (String a : args) {
			if (a.startsWith("--limit=")) limit = Integer.parseInt(a.substring(8).trim());
			else if (a.startsWith("--passes=")) maxPasses = Integer.parseInt(a.substring(9).trim());
			else if (a.startsWith("--shapes=")) explicit = parseShapeList(a.substring(9));
			else if (a.startsWith("--shapes-file=")) explicit = readShapeFile(a.substring(14).trim());
			else if (a.matches("-?\\d+")) { if (pos++ == 0) limit = Integer.parseInt(a); else maxPasses = Integer.parseInt(a); }
		}

		List<int[]> worse; // each: {n,m,p,ours,fmm} — ours/fmm are 0 for explicit shapes
		if (explicit != null) {
			worse = explicit;
			log.info("explicit shape list: {} shape(s)", worse.size());
		} else {
			// WORSE shapes are already computed by FmmCrossCheck — parse its report rather
			// than recomputing "ours" (findRank over all FMM shapes is pathologically slow).
			worse = new ArrayList<>();
			java.util.regex.Pattern row = java.util.regex.Pattern.compile(
					"\\u27e8(\\d+),(\\d+),(\\d+)\\u27e9\\s*\\|\\s*(\\d+)\\s*\\|\\s*(\\d+)");
			boolean inWorse = false;
			for (String line : java.nio.file.Files.readAllLines(Path.of("references/fmm-cross-check.md"))) {
				if (line.startsWith("## WORSE")) { inWorse = true; continue; }
				if (line.startsWith("## ") && inWorse) break; // next section ends WORSE
				if (!inWorse) continue;
				java.util.regex.Matcher mt = row.matcher(line);
				if (mt.find()) worse.add(new int[] {
						Integer.parseInt(mt.group(1)), Integer.parseInt(mt.group(2)), Integer.parseInt(mt.group(3)),
						Integer.parseInt(mt.group(4)), Integer.parseInt(mt.group(5)) });
			}
		}
		// Best gap-per-second first: maximise (gap / projection-cost), where cost ≈
		// maxDim³ (PanTA-cube build + survivor count scale with the axis size). A big
		// gap on a small shape closes the most rank fastest, so the cross-check
		// improves quickest even if the sweep is interrupted (user 2026-06-08). Tie
		// secondary: larger shape first (parent-before-child cascade within a score).
		worse.sort(Comparator
				.<int[]>comparingDouble(r -> {
					long mx = Math.max(r[0], Math.max(r[1], r[2]));
					double gap = r[3] - r[4];
					return -gap / (double) (mx * mx * mx);
				})
				.thenComparingInt(r -> -Math.max(r[0], Math.max(r[1], r[2]))));
		final List<int[]> shapes = (limit > 0 && worse.size() > limit)
				? new ArrayList<>(worse.subList(0, limit)) : worse;
		log.info("shapes to project: {} (limit={}, maxPasses={})", shapes.size(), limit, maxPasses);

		log.info("building FieldAwareLookup(R)…");
		FieldAwareLookup lookup = new FieldAwareLookup("R");
		log.info("lookup ready; starting projection.");

		// Projection-only materialiser: empty pool/sota (projectInto uses neither).
		List<BlockSplitSearch.NamedBase> pool = List.of();
		RecursiveClosureSota sota = new RecursiveClosureSota(lookup, pool, true, true);
		RecursiveMaterialiser mat = new RecursiveMaterialiser(lookup, pool, sota, SCHEMES_ROOT, true, true, true);

		// Parent-centric (scatter) sweep: resolves each slightly-larger parent shape
		// EXACTLY ONCE per pass and fans it out to every target child it covers (vs the
		// old child-centric loop that re-loaded the same big stub once per child). The
		// gap/second ordering above is now moot — scatter derives its own largest-first
		// parent order internally — but is kept harmless for the explicit-shapes path.
		int totalWins = mat.projectScatter(shapes, maxPasses);
		log.info("ProjectFmmGaps done: {} projection win(s) over {} shape(s); {} parent replays skipped by margin prune.",
				totalWins, shapes.size(), RecursiveMaterialiser.prunedParents());
	}

	/** Parse {@code --shapes=27x28x28,31x32x32} (comma/semicolon outer, x inner). */
	private static List<int[]> parseShapeList(String s) {
		List<int[]> out = new ArrayList<>();
		for (String tok : s.split("[;,]+")) {
			tok = tok.trim();
			if (tok.isEmpty()) continue;
			out.add(parseShape(tok));
		}
		return out;
	}

	/** Read a shapes file: one shape per line ("NxMxP" or "n,m,p" or "n m p"),
	 *  blank lines and {@code #}-comments ignored. */
	private static List<int[]> readShapeFile(String path) throws java.io.IOException {
		List<int[]> out = new ArrayList<>();
		for (String line : java.nio.file.Files.readAllLines(Path.of(path))) {
			String t = line.trim();
			if (t.isEmpty() || t.startsWith("#")) continue;
			out.add(parseShape(t));
		}
		return out;
	}

	/** A single shape token split on x, comma or whitespace into {n,m,p,0,0}. */
	private static int[] parseShape(String tok) {
		String[] p = tok.trim().split("[xX*,\\s]+");
		if (p.length != 3) throw new IllegalArgumentException("bad shape token: '" + tok + "'");
		return new int[] { Integer.parseInt(p[0]), Integer.parseInt(p[1]), Integer.parseInt(p[2]), 0, 0 };
	}
}
