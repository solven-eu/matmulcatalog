package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.Function;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.AxisSplitBases;
import eu.solven.matmul.papers.laderman1976.Laderman23;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.SymmetryTransforms;
import eu.solven.matmul.search.BlockSplitSearch;
import eu.solven.matmul.catalog.Recombination;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Phase C integration: expand each outer base in the pool to its S₃ orbit
 * (cyclic shifts + transpose) and run the multi-base block-split search.
 * Equivalent algorithms have different U/V/W column-support multisets per
 * axis, so they can produce different sub-shape distributions under
 * min-reduction — the framework picks the orbit element giving min rank.
 */
public class TestMultiBaseSymmetrySearch {

	private static List<BlockSplitSearch.NamedBase> buildSymmetryExpandedPool() throws Exception {
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		NonCubicBilinearAlgorithm laderman = NonCubicBilinearAlgorithm.fromCubic(Laderman23.get());

		List<BlockSplitSearch.NamedBase> pool = new ArrayList<>();
		expandAndAdd(pool, "Strassen ⟨2,2,2⟩=7", strassen);
		expandAndAdd(pool, "Laderman ⟨3,3,3⟩=23", laderman);
		expandAndAdd(pool, "mul211", AxisSplitBases.mul211());
		expandAndAdd(pool, "mul121", AxisSplitBases.mul121());
		expandAndAdd(pool, "mul112", AxisSplitBases.mul112());
		return pool;
	}

	private static void expandAndAdd(List<BlockSplitSearch.NamedBase> pool,
			String label, NonCubicBilinearAlgorithm base) {
		List<NonCubicBilinearAlgorithm> orbit = SymmetryTransforms.s3Orbit(base);
		for (int i = 0; i < orbit.size(); i++) {
			NonCubicBilinearAlgorithm a = orbit.get(i);
			String tag = label + " [σ=" + i + " ⟨" + a.n + "," + a.m + "," + a.p + "⟩]";
			pool.add(new BlockSplitSearch.NamedBase(tag, a));
		}
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
	public void symmetry_expanded_pool_777_meets_or_beats_plain_strassen() throws Exception {
		Recombination.SotaResolver sota = buildSotaR();
		List<BlockSplitSearch.NamedBase> pool = buildSymmetryExpandedPool();
		Optional<BlockSplitSearch.MultiBaseSplitCandidate> best =
				BlockSplitSearch.findBestMultiBaseSplit(7, 7, 7, pool, sota, true);
		assertThat(best).isPresent();
		System.out.println("⟨7,7,7⟩ best (S₃-expanded): " + best.get().breakdown());
		// Sedoglavic via Strassen ⟨2,2,2⟩=7 gives 250 — any symmetry variant
		// must do at least as well.
		assertThat(best.get().rank()).isLessThanOrEqualTo(250);
	}

	@Test
	public void symmetry_scan_cubic_4_to_12_report() throws Exception {
		Recombination.SotaResolver sota = buildSotaR();
		List<BlockSplitSearch.NamedBase> pool = buildSymmetryExpandedPool();
		System.out.printf("pool size (with S₃ orbit expansion): %d%n", pool.size());
		System.out.printf("%7s | %5s | %s%n", "target", "best", "via");
		System.out.println("-".repeat(80));
		for (int n = 4; n <= 12; n++) {
			Optional<BlockSplitSearch.MultiBaseSplitCandidate> best =
					BlockSplitSearch.findBestMultiBaseSplit(n, n, n, pool, sota, true);
			if (best.isEmpty()) {
				System.out.printf("%7s | %5s | —%n", "⟨" + n + "⟩³", "?");
				continue;
			}
			System.out.printf("%7s | %5d | %s%n",
					"⟨" + n + "⟩³", best.get().rank(), best.get().breakdown());
			assertThat(best.get().rank()).isLessThan(Long.MAX_VALUE / 100);
		}
	}
}
