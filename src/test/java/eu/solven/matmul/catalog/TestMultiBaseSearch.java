package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.Function;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.AxisSplitBases;
import eu.solven.matmul.papers.laderman1976.Laderman23;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.search.BlockSplitSearch;
import eu.solven.matmul.catalog.Recombination;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Multi-base block-split: prove it recovers Sedoglavic's ⟨7,7,7⟩=250
 * via Strassen, and check if any other base in the pool gives lower
 * rank for selected cubic targets.
 */
public class TestMultiBaseSearch {

	private static List<BlockSplitSearch.NamedBase> buildPool() throws Exception {
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		NonCubicBilinearAlgorithm laderman = NonCubicBilinearAlgorithm.fromCubic(Laderman23.get());
		// Hopcroft-Kerr ⟨2,3,3⟩ = 15 (and axis perms via orientAs at lookup time)
		// Skip for now — use what's clearly in catalog
		return List.of(
				new BlockSplitSearch.NamedBase("Strassen ⟨2,2,2⟩=7", strassen),
				new BlockSplitSearch.NamedBase("Laderman ⟨3,3,3⟩=23", laderman),
				new BlockSplitSearch.NamedBase("mul211 ⟨2,1,1⟩=2", AxisSplitBases.mul211()),
				new BlockSplitSearch.NamedBase("mul121 ⟨1,2,1⟩=2", AxisSplitBases.mul121()),
				new BlockSplitSearch.NamedBase("mul112 ⟨1,1,2⟩=2", AxisSplitBases.mul112())
		);
	}

	private static Recombination.SotaResolver buildSotaR() {
		Map<String, Integer> ranks = BlockSplitSearch.loadCatalogBestRanksForField("R");
		Function<int[], Optional<Integer>> lookup = BlockSplitSearch.rankLookupFromMap(ranks);
		return (a, b, c) -> {
			if (a == 0 || b == 0 || c == 0) return 0;
			if (a == 1) return b * c;
			if (b == 1) return a * c;
			if (c == 1) return a * b;
			return lookup.apply(new int[] { a, b, c }).orElse(Integer.MAX_VALUE / 100);
		};
	}

	@Test
	public void multibase_777_recovers_sedoglavic_250_or_better() throws Exception {
		Recombination.SotaResolver sota = buildSotaR();
		List<BlockSplitSearch.NamedBase> pool = buildPool();
		Optional<BlockSplitSearch.MultiBaseSplitCandidate> best =
				BlockSplitSearch.findBestMultiBaseSplit(7, 7, 7, pool, sota, true /* balanced */);
		assertThat(best).isPresent();
		System.out.println("⟨7,7,7⟩ best: " + best.get().breakdown());
		assertThat(best.get().rank()).isLessThanOrEqualTo(250);
	}

	@Test
	public void multibase_scan_cubic_4_to_16() throws Exception {
		Recombination.SotaResolver sota = buildSotaR();
		List<BlockSplitSearch.NamedBase> pool = buildPool();
		System.out.printf("%7s | %5s | %s%n", "target", "best", "via");
		System.out.println("-".repeat(80));
		for (int n = 4; n <= 16; n++) {
			Optional<BlockSplitSearch.MultiBaseSplitCandidate> best =
					BlockSplitSearch.findBestMultiBaseSplit(n, n, n, pool, sota, true);
			if (best.isEmpty()) {
				System.out.printf("%7s | %5s | —%n", "⟨" + n + "⟩³", "?");
				continue;
			}
			System.out.printf("%7s | %5d | %s%n",
					"⟨" + n + "⟩³", best.get().rank(), best.get().breakdown());
		}
	}
}
