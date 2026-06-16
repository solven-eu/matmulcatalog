package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.TreeMap;
import java.util.function.Function;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Tag;

import eu.solven.matmul.AxisSplitBases;
import eu.solven.matmul.papers.laderman1976.Laderman23;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.SymmetryTransforms;
import eu.solven.matmul.search.BlockSplitSearch;
import eu.solven.matmul.catalog.Recombination;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Full DIS09 ⟨n,n,n⟩ ranks (R/Q/Z) for n=4..30 vs our multi-base + S₃
 * symmetry expansion using the modern catalog (Sedoglavic, AlphaTensor,
 * AlphaEvolve, etc.) as the SOTA resolver for sub-shapes.
 */
@Tag("catalog-iterating")
public class TestDIS09FullScan {

	/** DIS09 Table 3 — non-commutative ⟨n,n,n⟩ over R/Q/Z, n=4..30. */
	private static final Map<Integer, Integer> DIS09_R = new TreeMap<>(Map.ofEntries(
			Map.entry(4, 49),
			Map.entry(5, 100),
			Map.entry(6, 161),
			Map.entry(7, 258),
			Map.entry(8, 343),
			Map.entry(9, 522),
			Map.entry(10, 700),
			Map.entry(11, 923),
			Map.entry(12, 1125),
			Map.entry(13, 1450),
			Map.entry(14, 1728),
			Map.entry(15, 2108),
			Map.entry(16, 2401),
			Map.entry(17, 2972),
			Map.entry(18, 3306),
			Map.entry(19, 4073),
			Map.entry(20, 4340),
			Map.entry(21, 5365),
			Map.entry(22, 5566),
			Map.entry(23, 6806),
			Map.entry(24, 7000),
			Map.entry(25, 8448),
			Map.entry(26, 8658),
			Map.entry(27, 10330),
			Map.entry(28, 10556),
			Map.entry(29, 12468),
			Map.entry(30, 12710)));

	private static List<BlockSplitSearch.NamedBase> buildPool() throws Exception {
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		NonCubicBilinearAlgorithm laderman = NonCubicBilinearAlgorithm.fromCubic(Laderman23.get());
		// Hopcroft-Kerr 1971 ⟨2,3,3⟩=15 (proved tight in HK71; same scheme found
		// computationally by AlphaTensor 2022 over Z). S₃ orbit covers ⟨3,3,2⟩
		// and ⟨3,2,3⟩ orientations DIS09 calls "Hopcroft332".
		NonCubicBilinearAlgorithm hopcroft = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section3/alphatensor_Z-2x3x3_m15_a58.json"));

		List<BlockSplitSearch.NamedBase> pool = new ArrayList<>();
		expandAndAdd(pool, "Strassen ⟨2,2,2⟩=7", strassen);
		expandAndAdd(pool, "Laderman ⟨3,3,3⟩=23", laderman);
		expandAndAdd(pool, "Hopcroft-Kerr ⟨2,3,3⟩=15", hopcroft);
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
			String tag = label + " [σ" + i + " ⟨" + a.n + "," + a.m + "," + a.p + "⟩]";
			pool.add(new BlockSplitSearch.NamedBase(tag, a));
		}
	}

	private static Recombination.SotaResolver buildSotaR() {
		return buildSotaForField("R");
	}

	private static Recombination.SotaResolver buildSotaForField(String field) {
		Map<String, Integer> ranks = BlockSplitSearch.loadCatalogBestRanksForField(field);
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
	public void dis09_full_comparison_R_4_to_30_balanced() throws Exception {
		runScan(true, "R", "balanced");
	}

	@Test
	public void dis09_full_comparison_R_4_to_30_unbalanced() throws Exception {
		runScan(false, "R", "FULL enumeration");
	}

	@Test
	public void dis09_full_comparison_C_4_to_30_unbalanced() throws Exception {
		runScan(false, "C", "FULL enumeration");
	}

	@Test
	public void dis09_full_comparison_F2_4_to_30_unbalanced() throws Exception {
		runScan(false, "F2", "FULL enumeration");
	}

	private void runScan(boolean balanced, String field, String label) throws Exception {
		Recombination.SotaResolver sota = buildSotaForField(field);
		List<BlockSplitSearch.NamedBase> pool = buildPool();

		System.out.println();
		System.out.println("DIS09 Table 3 (NON-COMMUTATIVE) vs Strassen-symmetric multi-base "
				+ "(field=" + field + ", " + label + ")");
		System.out.println("DIS09 baseline is over a generic non-commutative ring "
				+ "(applies to R/Q/Z; F₂ and C can use more catalog entries).");
		System.out.println("=".repeat(95));
		System.out.printf("%5s | %7s | %7s | %6s | %6s | %s%n",
				"n", "DIS09", "ours", "Δ", "Δ%", "via");
		System.out.println("-".repeat(95));

		int wins = 0, ties = 0, losses = 0;
		long totalDis09 = 0, totalOurs = 0;
		for (int n = 4; n <= 30; n++) {
			int dis09 = DIS09_R.get(n);
			Optional<BlockSplitSearch.MultiBaseSplitCandidate> best =
					BlockSplitSearch.findBestMultiBaseSplit(n, n, n, pool, sota, balanced, 1_000_000L);
			long recursive = best.map(BlockSplitSearch.MultiBaseSplitCandidate::rank).orElse(Long.MAX_VALUE);
			long direct = sota.getRank(n, n, n);
			long ours;
			String via;
			if (direct < recursive) {
				ours = direct;
				via = "DIRECT catalog/Pan-TA = " + direct;
			} else if (best.isPresent()) {
				ours = recursive;
				via = best.get().breakdown();
			} else {
				System.out.printf("%5d | %7d | %7s | %6s | %6s | —%n", n, dis09, "?", "?", "?");
				continue;
			}
			long delta = ours - dis09;
			double pct = 100.0 * delta / dis09;
			String marker;
			if (delta < 0) { wins++; marker = "✓"; }
			else if (delta == 0) { ties++; marker = "="; }
			else { losses++; marker = "✗"; }
			totalDis09 += dis09;
			totalOurs += ours;
			System.out.printf("%5d | %7d | %7d | %+6d | %+6.2f | %s %s%n",
					n, dis09, ours, delta, pct, marker, via);
		}
		System.out.println("=".repeat(95));
		System.out.printf("Summary: wins=%d ties=%d losses=%d  | totals: DIS09=%d ours=%d (%.2f%%)%n",
				wins, ties, losses, totalDis09, totalOurs, 100.0 * (totalOurs - totalDis09) / totalDis09);
		assertThat(wins + ties).isGreaterThanOrEqualTo(losses);
	}
}
