package eu.solven.matmul.docs.migrate;

import java.io.File;
import java.math.BigInteger;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.SchemeAnalysis;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.io.NpzSchemeReader;
import eu.solven.matmul.papers.khoruzhii2026.LitaNpzCubes;
import eu.solven.matmul.search.LineageReplayer;
import eu.solven.matmul.verifiers.Verifier;

import lombok.extern.slf4j.Slf4j;

/**
 * Import the Khoruzhii–Serafin–Gelß–Pokutta 2026 LITA cubes (REFERENCES.md [82];
 * <a href="https://github.com/khoruzhii/lita">khoruzhii/lita</a> {@code schemes/*.npz},
 * {@code ⟨N,N,N⟩} for {@code 13 ≤ N ≤ 32}) — the source of FMM-Lille's 2026-09
 * index refresh, where every large cube and (by projection) 509 more shapes moved
 * below our catalog. Our {@code LitaTaConstruction} ports their June generator
 * (N ≥ 19); the repo has since gone through "LITA3 … LITA7 odd" (Aug–Sep 2026)
 * with explicit rational schemes at {@code R_even = N³/3 + 3N² + 37N/6 + 5},
 * {@code R_odd = N³/3 + 7N²/2 + 14N/3 − 11/2} (32³ = 14197, ω ≈ 2.7587).
 *
 * <p>Per file, a SIZE-AWARE gate (these are the largest schemes in the catalog —
 * 32³ is 1024 × 14197 per factor):</p>
 * <ol>
 *   <li><b>Exact symbolic proof</b> ({@link Verifier#isExactNonCubic}, BigInteger over
 *       the lcm denominator) when the sparse term estimate fits {@code --max-exact-terms}
 *       (default 60M); otherwise a 20-sample random matmul spot-check — an operational
 *       check, NOT an algebraic proof, and recorded as such in {@code verification}.</li>
 *   <li><b>Fields from the coefficients</b>: Q, R, C always (Q-exact rationals). F₂/F₃
 *       only when (a) the lcm of denominators is coprime to p — otherwise the scheme is
 *       provably NOT reducible mod p as written → {@code fields_not} — and (b) the DENSE
 *       mod-p verification ({@code dimU·dimV·dimW·r} triple products) is affordable;
 *       when it is not, the field is left UNCLAIMED (neither list).</li>
 *   <li>Catalog comparison over Q; disk round-trip after the canonical write.</li>
 * </ol>
 *
 * <p>Metadata: {@code source "Khoruzhii 2026"}, {@code discovery: true} (these ranks
 * beat every catalog; FMM cites the repo), {@code source_paper_url} / {@code source_scheme_url}
 * pointing at the repo / the exact file. Dry-run by default; {@code --execute}
 * archives the {@code .npz} under {@link LitaNpzCubes#DIR} and writes a HASH-STAMPED
 * STUB {@code schemes/known/section{N}/{N}x{N}x{N}-r{R}-khoruzhii_2026-{hash7}.json}
 * whose lineage is the data-backed atom {@code TA_lita_npz(n=N)} (replayed by
 * {@code LineageReplayer} through {@link NpzSchemeReader}). Explicit catalog JSON is
 * not an option at this size: 32³ is 117 MB canonical (GitHub's file limit is
 * 100 MB) and the 20 cubes total 491 MB, against 12 MB of archives.</p>
 *
 * <pre>mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.ImportKhoruzhiiLita \
 *     -Dexec.args="path/to/npz-dir [--only=NxNxN] [--max-exact-terms=60000000] [--execute]"</pre>
 */
@Slf4j
public final class ImportKhoruzhiiLita {
	private ImportKhoruzhiiLita() {}

	private static final Path SCHEMES_ROOT = Path.of("src/main/resources/schemes");
	private static final String REPO_URL = "https://github.com/khoruzhii/lita";
	private static final String RAW_BASE = "https://github.com/khoruzhii/lita/blob/master/schemes/";
	/** Dense mod-p verification budget: dimU·dimV·dimW·r triple products. ~2·10¹¹ ≈ a few minutes. */
	private static final long MAX_DENSE_FP_TRIPLES = 200_000_000_000L;

	public static void main(String[] args) throws Exception {
		Path dir = null;
		String only = null;
		long maxExactTerms = 60_000_000L;
		boolean execute = false;
		for (String a : args) {
			if ("--execute".equals(a)) execute = true;
			else if (a.startsWith("--only=")) only = a.substring("--only=".length());
			else if (a.startsWith("--max-exact-terms=")) maxExactTerms = Long.parseLong(a.substring("--max-exact-terms=".length()));
			else if (a.startsWith("--")) throw new IllegalArgumentException("unknown option " + a);
			else dir = Path.of(a);
		}
		if (dir == null || !Files.isDirectory(dir)) {
			throw new IllegalArgumentException("usage: ImportKhoruzhiiLita <npz-dir> [--only=NxNxN] [--max-exact-terms=N] [--execute]");
		}
		List<Path> files;
		try (Stream<Path> ls = Files.list(dir)) {
			files = ls.filter(f -> f.toString().endsWith(".npz")).sorted(
					java.util.Comparator.comparingLong(f -> f.toFile().length())).toList();
		}
		FieldAwareLookup lookup = new FieldAwareLookup(Field.fromTag("Q"), SCHEMES_ROOT);
		List<String> report = new ArrayList<>();
		List<String> failures = new ArrayList<>();
		int written = 0;
		long t0 = System.nanoTime();
		for (Path f : files) {
			String stem = f.getFileName().toString().replace(".npz", "");
			if (only != null && !stem.startsWith(only + "_")) continue;
			long s = System.nanoTime();
			NpzSchemeReader.Loaded loaded;
			try {
				loaded = NpzSchemeReader.read(f.toFile());
			} catch (Exception e) {
				failures.add(stem + ": unreadable — " + e);
				log.error("{}: unreadable — {}", stem, e.toString());
				continue;
			}
			NonCubicBilinearAlgorithm alg = loaded.alg();
			long est = Verifier.estimateExactTerms(alg, Long.MAX_VALUE / 4);
			long dense = (long) alg.n * alg.m * (long) alg.m * alg.p * (long) alg.n * alg.p * alg.r;
			log.info("{} → ⟨{},{},{}⟩ r={} field={} fmt={} maxDen={} lcmDen={} exactTerms≈{} denseFpTriples≈{}",
					stem, alg.n, alg.m, alg.p, alg.r, loaded.coefficientField(), loaded.format(),
					loaded.maxDenominator(), loaded.lcmDenominator(), est, dense);

			// 1. verification tier
			boolean ok;
			String tier;
			if (est <= maxExactTerms) {
				ok = Verifier.isExactNonCubic(alg);
				tier = "exact-symbolic (BigInteger over lcm denominator, " + est + " terms)";
			} else {
				double eps = 1e-8 * Math.sqrt(alg.n * alg.p) * Math.sqrt(alg.r) * loaded.maxDenominator();
				ok = Verifier.passesRandomMatmulSpotCheck(alg, 20, eps);
				tier = "random matmul spot-check (20 samples, eps=" + eps + "; exact term estimate " + est
						+ " exceeds the " + maxExactTerms + " cap — NOT an algebraic proof)";
			}
			if (!ok) {
				failures.add(stem + ": VERIFICATION FAILED — " + tier);
				log.error("{}: VERIFICATION FAILED — {}", stem, tier);
				continue;
			}
			// 2. fields
			List<String> fields = new ArrayList<>();
			List<String> fieldsNot = new ArrayList<>();
			List<String> unclaimed = new ArrayList<>();
			for (int p : new int[] { 2, 3 }) {
				String tag = "F" + p;
				if (loaded.lcmDenominator().mod(BigInteger.valueOf(p)).signum() == 0) {
					fieldsNot.add(tag); // a denominator divisible by p: not reducible mod p as written
				} else if (dense <= MAX_DENSE_FP_TRIPLES) {
					if (Verifier.isExactNonCubicFp(alg, p)) fields.add(tag);
					else fieldsNot.add(tag);
				} else {
					unclaimed.add(tag);
				}
			}
			boolean integer = loaded.maxDenominator() == 1;
			if (integer) fields.add("Z");
			fields.add("Q");
			fields.add("R");
			fields.add("C");
			int catalog = lookup.findRank(alg.n, alg.m, alg.p);
			String vs = alg.r < catalog ? "BEATS " + catalog + " (−" + (catalog - alg.r) + ")"
					: alg.r == catalog ? "ties " + catalog : "ABOVE " + catalog;
			String line = String.format("%s → ⟨%d,%d,%d⟩ r=%d  %s  fields=%s not=%s unclaimed=%s  catalog(Q): %s  (%d ms)",
					stem, alg.n, alg.m, alg.p, alg.r, ok ? "OK" : "FAIL", fields, fieldsNot, unclaimed, vs,
					(System.nanoTime() - s) / 1_000_000L);
			log.info("{}  [{}]", line, tier);
			report.add(line);
			if (!execute) continue;

			// 3. register: archive under the resource dir + a hash-stamped STUB on the
			//    data-backed atom TA_lita_npz(n=N). NOT explicit JSON: 32³ is 117 MB canonical
			//    (over GitHub's 100 MB file limit; all 20 total 491 MB) vs 12 MB of archives —
			//    the catalog stores big formula cubes as stubs replayed on demand.
			if (alg.n != alg.m || alg.m != alg.p) {
				throw new IllegalStateException(stem + ": TA_lita_npz covers cubes only, got ⟨" + alg.n + "," + alg.m
						+ "," + alg.p + "⟩");
			}
			Files.createDirectories(LitaNpzCubes.DIR);
			Path archived = LitaNpzCubes.DIR.resolve(f.getFileName().toString());
			if (!archived.toAbsolutePath().equals(f.toAbsolutePath())) {
				Files.copy(f, archived, java.nio.file.StandardCopyOption.REPLACE_EXISTING);
			}
			int maxDim = Math.max(alg.n, Math.max(alg.m, alg.p));
			Path outDir = SCHEMES_ROOT.resolve("known").resolve("section" + maxDim);
			Files.createDirectories(outDir);
			File out = outDir.resolve(SchemeIO.canonicalName(alg, "khoruzhii_2026")).toFile();
			SchemeIO.writeStub(alg, out, new Lineage.Atom("TA_lita_npz(n=" + alg.n + ")"), fields);
			Map<String, Object> meta = new LinkedHashMap<>();
			meta.put("fields", fields);
			meta.put("fields_not", fieldsNot);
			meta.put("commutative", false);
			meta.put("source", "Khoruzhii 2026");
			meta.put("year", 2026);
			meta.put("source_paper_url", REPO_URL);
			meta.put("source_scheme_url", RAW_BASE + f.getFileName());
			meta.put("discovery", true);
			meta.put("discovery_note", "Khoruzhii, Serafin, Gelß & Pokutta 2026, Local Improvements to Trilinear "
					+ "Aggregation (LITA), github.com/khoruzhii/lita — explicit rational ⟨" + alg.n + "," + alg.m + ","
					+ alg.p + "⟩=" + alg.r + " (R_even = N³/3+3N²+37N/6+5, R_odd = N³/3+7N²/2+14N/3−11/2; the "
					+ "Aug–Sep 2026 LITA3…LITA7 constructions, superseding the June generator ported in "
					+ "LitaTaConstruction). Cited by fmm.univ-lille.fr for this cube and, by projection, for "
					+ "its whole ⟨≤N⟩ family. Imported from the repo's " + f.getFileName()
					+ " (fmm.q.csr.v1) via NpzSchemeReader; REFERENCES.md [82].");
			meta.put("verified", true);
			meta.put("verification", tier);
			if (!unclaimed.isEmpty()) {
				meta.put("notes", "F-p reducibility for " + unclaimed + " left UNCLAIMED: denominators are coprime to p "
						+ "but the dense mod-p verification (" + dense + " triple products) exceeds the import budget.");
			}
			Field field = Field.fromTag("Q");
			for (SchemeAnalysis an : SchemeAnalysis.defaults()) {
				long as = System.nanoTime();
				for (Map.Entry<String, Object> e : an.analyse(alg, field).entrySet()) {
					meta.putIfAbsent(e.getKey(), e.getValue());
				}
				log.info("  {}: {} ms", an.name(), (System.nanoTime() - as) / 1_000_000L);
			}
			SchemeIO.updateFields(out, meta, List.of(), true);
			// Round-trip the STUB: replaying its lineage atom must rebuild exactly the
			// verified matrices (the stamped hash is what projection stubs pin as parent).
			NonCubicBilinearAlgorithm back = LineageReplayer.withDefaultPool(lookup).replayFromFile(out);
			if (!SchemeIO.contentHash(back).equals(SchemeIO.contentHash(alg))) {
				throw new IllegalStateException("STUB ROUND-TRIP content hash drifted: " + out);
			}
			written++;
			log.info("wrote {}", SCHEMES_ROOT.relativize(out.toPath()));
		}
		log.info("==================================================================");
		for (String r : report) log.info("  OK  {}", r);
		for (String x : failures) log.error("  FAIL {}", x);
		log.info("ImportKhoruzhiiLita: {} files, {} passed, {} failed, {} written ({}) — {} ms", files.size(),
				report.size(), failures.size(), written, execute ? "--execute" : "DRY-RUN",
				(System.nanoTime() - t0) / 1_000_000L);
		if (!failures.isEmpty()) {
			throw new IllegalStateException(failures.size() + " file(s) failed the import gate");
		}
	}
}
