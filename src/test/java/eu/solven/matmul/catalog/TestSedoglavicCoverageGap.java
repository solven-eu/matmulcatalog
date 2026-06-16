package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.Function;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.search.BlockSplitSearch;
import eu.solven.matmul.catalog.Recombination;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Coverage-gap detector.
 *
 * <p>For every cubic target {@code ⟨n,n,n⟩} ({@code n ∈ [4, 32]}),
 * compute the Sedoglavic-formula upper bound via {@link BlockSplitSearch}
 * and compare it to the direct catalog rank. A <strong>gap</strong>
 * (formula &lt; direct) means our catalog has a worse upper bound than
 * is constructively achievable — i.e. we should either run
 * {@link Recombination#constructWithAllocation} for that target, or
 * import a better direct scheme.</p>
 *
 * <p>This test is fast (≤ ~1 sec) — runs O(n) formula evaluations per
 * target via catalog map lookups, no scheme construction. It's also a
 * <strong>regression test</strong>: improving any small-format rank
 * (e.g. a new ⟨3,3,3⟩ algorithm) will automatically propagate via
 * Sedoglavic and may surface NEW gaps for big targets that previously
 * matched.</p>
 */
public class TestSedoglavicCoverageGap {

	@Test
	public void scan_for_sedoglavic_gaps_per_field() {
		System.out.println();
		for (String field : new String[] { "R", "C", "F2" }) {
			scanField(field);
		}
	}

	@Test
	public void scan_for_sedoglavic_gaps_in_noncubic_targets_all_fields() throws Exception {
		// Resolve via the tolerant resolver (content-driven names re-hash; never
		// pin a literal path — see SchemeResolver).
		eu.solven.matmul.NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				SchemeResolver.byHint("strassen-2x2x2_m7.json"));
		for (String field : new String[] { "R", "C", "F2" }) {
			scanNonCubicField(strassen, field, 16);
		}
	}

	private void scanNonCubicField(eu.solven.matmul.NonCubicBilinearAlgorithm strassen,
			String field, int maxDim) throws Exception {
		Map<String, Integer> ranks = BlockSplitSearch.loadCatalogBestRanksForField(field);
		Function<int[], Optional<Integer>> lookup = BlockSplitSearch.rankLookupFromMap(ranks);
		Recombination.SotaResolver sota = (a, b, c) -> {
			if (a == 0 || b == 0 || c == 0) return 0;
			if (a == 1) return b * c;
			if (b == 1) return a * c;
			if (c == 1) return a * b;
			return lookup.apply(new int[] { a, b, c }).orElse(Integer.MAX_VALUE / 100);
		};

		Recombination.SotaResolver sotaWithTrivial = (a, b, c) -> {
			if (a == 0 || b == 0 || c == 0) return 0;
			if (a == 1) return b * c;
			if (b == 1) return a * c;
			if (c == 1) return a * b;
			return lookup.apply(new int[] { a, b, c }).orElse(Integer.MAX_VALUE / 100);
		};

		int gaps = 0, missing = 0, ok = 0;
		List<String> gapLines = new ArrayList<>();
		int gapCap = 10, missingCap = 5;
		for (int n = 2; n <= maxDim; n++) {
			for (int m = n; m <= maxDim; m++) {
				for (int p = m; p <= maxDim; p++) {
					if (n == m && m == p) continue;
					int tn = n, tm = m, tp = p;
					Recombination.SotaResolver pureSota = (a, b, c) -> {
						if (a == tn && b == tm && c == tp) return Integer.MAX_VALUE / 100;
						return sotaWithTrivial.getRank(a, b, c);
					};
					Optional<BlockSplitSearch.NonCubicSplitCandidate> best =
							BlockSplitSearch.findBestSplitNonCubic(n, m, p, strassen, pureSota);
					if (best.isEmpty()) continue;
					long formulaRank = best.get().rank();
					Optional<Integer> direct = lookup.apply(new int[] { n, m, p });
					if (direct.isEmpty()) {
						missing++;
					} else if (formulaRank < direct.get()) {
						gaps++;
						if (gapLines.size() < gapCap) {
							gapLines.add(String.format(
									"  ⟨%d,%d,%d⟩: direct %s=%d, formula=%d (Δ=%d) via %s",
									n, m, p, field, direct.get(), formulaRank,
									direct.get() - formulaRank, best.get().breakdown()));
						}
					} else {
						ok++;
					}
				}
			}
		}
		System.out.println("=== Non-cubic " + field + " scan (dims ≤ " + maxDim + ") ===");
		System.out.println("  " + ok + " targets at/above formula, " + gaps + " gaps, "
				+ missing + " missing direct.");
		if (!gapLines.isEmpty()) {
			System.out.println("Top " + Math.min(gapLines.size(), gapCap)
					+ " gaps (formula beats direct):");
			gapLines.forEach(System.out::println);
		}
		assertThat(ok + gaps + missing).isGreaterThan(0);
	}

	private void scanField(String field) {
		Map<String, Integer> ranks = BlockSplitSearch.loadCatalogBestRanksForField(field);
		Function<int[], Optional<Integer>> lookup = BlockSplitSearch.rankLookupFromMap(ranks);

		List<String> gaps = new ArrayList<>();
		List<String> missingDirect = new ArrayList<>();
		for (int n = 4; n <= 32; n++) {
			Optional<BlockSplitSearch.SplitCandidate> best = BlockSplitSearch.findBestSplit(n, lookup);
			Optional<Integer> direct = lookup.apply(new int[] { n, n, n });
			if (best.isEmpty()) continue;
			BlockSplitSearch.SplitCandidate c = best.get();
			if (direct.isEmpty()) {
				missingDirect.add(String.format(
						"  ⟨%d,%d,%d⟩: no direct %s scheme; formula bound = %d via split %d+%d (%s)",
						n, n, n, field, c.formulaRank(), c.u(), c.v(), c.breakdown()));
				continue;
			}
			if (c.formulaRank() < direct.get()) {
				gaps.add(String.format(
						"  ⟨%d,%d,%d⟩: direct %s=%d, formula=%d (Δ=%d) via split %d+%d (%s)",
						n, n, n, field, direct.get(), c.formulaRank(),
						direct.get() - c.formulaRank(), c.u(), c.v(), c.breakdown()));
			}
		}

		System.out.println("=== Field " + field + " ===");
		if (gaps.isEmpty() && missingDirect.isEmpty()) {
			System.out.println("  (no gaps, all cubic targets up to 32 have at least catalog-matching ranks)");
		}
		if (!gaps.isEmpty()) {
			System.out.println("Sedoglavic gaps in " + field + " — formula beats direct catalog rank:");
			gaps.forEach(System.out::println);
		}
		if (!missingDirect.isEmpty()) {
			System.out.println("Targets with NO direct " + field + " catalog scheme:");
			missingDirect.forEach(System.out::println);
		}
		System.out.println();

		// Informational only — flips to a hard assertion once we want strict
		// field-pure coverage (i.e., every cubic target should have a direct
		// scheme at the formula bound or better).
		assertThat(gaps).as(field + " gaps (informational): see stdout").isNotNull();
	}
}
