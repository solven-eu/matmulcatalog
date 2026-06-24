package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.Function;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Tag;

import eu.solven.matmul.AxisSplitBases;
import eu.solven.matmul.papers.laderman1976.Laderman23;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.SymmetryTransforms;
import eu.solven.matmul.recombination.BlockSplitSearch;
import eu.solven.matmul.recombination.Recombination;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Probe: does enumerating all per-axis non-degenerate allocations
 * (balancedOnly=false) close the small Strassen-vs-Strassen gaps at
 * n=14 (DIS09 1728 vs our balanced 1750) and n=21 (DIS09 5365 vs 5409)?
 *
 * If yes, the catalog should switch to full enumeration for these
 * mid-sized cubic targets.
 */
@Tag("catalog-iterating")
public class TestDIS09UnbalancedProbe {

	private static List<BlockSplitSearch.NamedBase> buildPool() throws Exception {
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		NonCubicBilinearAlgorithm laderman = NonCubicBilinearAlgorithm.fromCubic(Laderman23.get());
		List<BlockSplitSearch.NamedBase> pool = new ArrayList<>();
		expandS3(pool, "Strassen", strassen);
		expandS3(pool, "Laderman", laderman);
		expandS3(pool, "mul211", AxisSplitBases.mul211());
		expandS3(pool, "mul121", AxisSplitBases.mul121());
		expandS3(pool, "mul112", AxisSplitBases.mul112());
		return pool;
	}

	private static void expandS3(List<BlockSplitSearch.NamedBase> pool,
			String label, NonCubicBilinearAlgorithm base) {
		List<NonCubicBilinearAlgorithm> orbit = SymmetryTransforms.s3Orbit(base);
		for (int i = 0; i < orbit.size(); i++) {
			NonCubicBilinearAlgorithm a = orbit.get(i);
			pool.add(new BlockSplitSearch.NamedBase(
					label + "[σ" + i + " ⟨" + a.n + "," + a.m + "," + a.p + "⟩]", a));
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
	public void probe_n14_unbalanced() throws Exception {
		Recombination.SotaResolver sota = buildSotaR();
		List<BlockSplitSearch.NamedBase> pool = buildPool();
		Optional<BlockSplitSearch.MultiBaseSplitCandidate> balanced =
				BlockSplitSearch.findBestMultiBaseSplit(14, 14, 14, pool, sota, true);
		Optional<BlockSplitSearch.MultiBaseSplitCandidate> unbalanced =
				BlockSplitSearch.findBestMultiBaseSplit(14, 14, 14, pool, sota, false);
		System.out.println("⟨14,14,14⟩ DIS09 = 1728");
		System.out.println("  balanced:   " + balanced.map(c -> c.rank() + " via " + c.breakdown()).orElse("—"));
		System.out.println("  unbalanced: " + unbalanced.map(c -> c.rank() + " via " + c.breakdown()).orElse("—"));
		assertThat(balanced).isPresent();
		assertThat(unbalanced).isPresent();
		assertThat(unbalanced.get().rank()).isLessThanOrEqualTo(balanced.get().rank());
	}

	@Test
	public void probe_n21_unbalanced() throws Exception {
		Recombination.SotaResolver sota = buildSotaR();
		List<BlockSplitSearch.NamedBase> pool = buildPool();
		Optional<BlockSplitSearch.MultiBaseSplitCandidate> balanced =
				BlockSplitSearch.findBestMultiBaseSplit(21, 21, 21, pool, sota, true);
		Optional<BlockSplitSearch.MultiBaseSplitCandidate> unbalanced =
				BlockSplitSearch.findBestMultiBaseSplit(21, 21, 21, pool, sota, false);
		System.out.println("⟨21,21,21⟩ DIS09 = 5365");
		System.out.println("  balanced:   " + balanced.map(c -> c.rank() + " via " + c.breakdown()).orElse("—"));
		System.out.println("  unbalanced: " + unbalanced.map(c -> c.rank() + " via " + c.breakdown()).orElse("—"));
		assertThat(balanced).isPresent();
		assertThat(unbalanced).isPresent();
		assertThat(unbalanced.get().rank()).isLessThanOrEqualTo(balanced.get().rank());
	}
}
