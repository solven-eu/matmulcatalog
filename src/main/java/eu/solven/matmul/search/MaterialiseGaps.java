package eu.solven.matmul.search;

import lombok.extern.slf4j.Slf4j;

import eu.solven.matmul.catalog.CatalogLimits;

import eu.solven.matmul.catalog.FieldAwareLookup;

import eu.solven.matmul.catalog.SchemeIO;

import eu.solven.matmul.catalog.Recombination;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.Function;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;

/**
 * Closes the loop: detect → construct → verify → save.
 *
 * <p>For each gap reported by {@link BlockSplitSearch#findBestSplitNonCubic}
 * (R-field, dim ≤ {@value #MAX_DIM}), materialise the algorithm via
 * {@link Recombination#constructWithAllocation}, verify exactness, and
 * save under {@code section{max_dim}/derived_strassen_recombine_NxMxP_rR.json}.</p>
 *
 * <p>Run:</p>
 * <pre>mvn exec:java -Dexec.mainClass=eu.solven.matmul.search.MaterialiseGaps</pre>
 *
 * <p>After running, re-generate the manifest + derived bounds:</p>
 * <pre>mvn exec:java -Dexec.mainClass=eu.solven.matmul.catalog.GenerateCatalogManifest
 *mvn exec:java -Dexec.mainClass=eu.solven.matmul.catalog.GenerateDerivedBounds</pre>
 */
@Slf4j
public final class MaterialiseGaps {

	private static final int MAX_DIM = CatalogLimits.MAX_DIM;
	private static final String SCHEMES_ROOT = "src/main/resources/schemes";

	private MaterialiseGaps() {}

	public static void main(String[] args) throws IOException {
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint(SCHEMES_ROOT + "/section2/strassen-2x2x2_m7_a18.json"));

		// R-field lookup (used both for SotaResolver in search AND for constructWithAllocation).
		FieldAwareLookup lookupR = new FieldAwareLookup("R");
		Map<String, Integer> ranks = BlockSplitSearch.loadCatalogBestRanksForField("R");
		Function<int[], Optional<Integer>> rankLookup = BlockSplitSearch.rankLookupFromMap(ranks);

		Recombination.SotaResolver sota = (a, b, c) -> {
			if (a == 0 || b == 0 || c == 0) return 0;
			if (a == 1) return b * c;
			if (b == 1) return a * c;
			if (c == 1) return a * b;
			return rankLookup.apply(new int[] { a, b, c }).orElse(Recombination.SotaResolver.UNKNOWN_RANK);
		};

		int materialised = 0, skipped = 0, failed = 0;
		List<String> logLines = new ArrayList<>();
		for (int n = 2; n <= MAX_DIM; n++) {
			for (int m = n; m <= MAX_DIM; m++) {
				for (int p = m; p <= MAX_DIM; p++) {
					int tn = n, tm = m, tp = p;
					Recombination.SotaResolver pure = (a, b, c) -> {
						if (a == tn && b == tm && c == tp) return Recombination.SotaResolver.UNKNOWN_RANK;
						return sota.getRank(a, b, c);
					};
					Optional<BlockSplitSearch.NonCubicSplitCandidate> best =
							BlockSplitSearch.findBestSplitNonCubic(n, m, p, strassen, pure);
					if (best.isEmpty()) continue;
					BlockSplitSearch.NonCubicSplitCandidate c = best.get();
					if (c.rank() >= Integer.MAX_VALUE / 200) continue;
					Optional<Integer> direct = rankLookup.apply(new int[] { n, m, p });
					if (direct.isPresent() && direct.get() <= c.rank()) continue;

					String label = String.format("⟨%d,%d,%d⟩", n, m, p);
					try {
						NonCubicBilinearAlgorithm alg = Recombination.constructWithAllocation(
								strassen, lookupR, c.allocA(), c.allocB(), c.allocC());
						if (alg.r != c.rank()) {
							logLines.add("[WARN] " + label + " rank mismatch: predicted=" + c.rank()
									+ " actual=" + alg.r);
						}
						if (!Verifier.isExactNonCubic(alg)) {
							logLines.add("[FAIL] " + label + " verification failed (predicted r=" + c.rank() + ")");
							failed++;
							continue;
						}
						int maxDim = Math.max(Math.max(n, m), p);
						String dir = SCHEMES_ROOT + "/derived/section" + maxDim;
						new File(dir).mkdirs();
						String fn = String.format("derived_strassen_recombine_%dx%dx%d_r%d.json",
								n, m, p, alg.r);
						File out = new File(dir, fn);
						if (out.exists()) {
							logLines.add("[SKIP] " + label + " already exists: " + out.getPath());
							skipped++;
							continue;
						}
						SchemeIO.write(alg, out);
						logLines.add(String.format("[OK]   %s r=%d (was %d, Δ=%d) → %s",
								label, alg.r, direct.orElse(null),
								direct.isPresent() ? direct.get() - alg.r : 0, out.getPath()));
						materialised++;
					} catch (Exception e) {
						logLines.add("[ERR]  " + label + ": " + e.getClass().getSimpleName() + " — " + e.getMessage());
						failed++;
					}
				}
			}
		}

		logLines.forEach(log::info);
		log.info("");
		log.info("Summary: " + materialised + " materialised, " + skipped + " skipped, " + failed + " failed.");
	}
}
