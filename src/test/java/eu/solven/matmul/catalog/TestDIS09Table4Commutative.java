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
import eu.solven.matmul.recombination.BlockSplitSearch;
import eu.solven.matmul.papers.dis2009.PanTrilinearAggregation;
import eu.solven.matmul.recombination.Recombination;
import eu.solven.matmul.papers.rosowski2019.RosowskiBound;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.papers.waksman1970.WaksmanBound;

/**
 * Reproduces DIS09 Table 4 (COMMUTATIVE ⟨n,n,n⟩ bounds, n=4..30).
 *
 * <p>Commutative SOTA = min of:</p>
 * <ul>
 *   <li>NC catalog (R-filtered): any NC scheme is also a valid
 *       commutative upper bound</li>
 *   <li>Waksman 1970 closed-form family ({@link WaksmanBound})</li>
 *   <li>Rosowski 2019 commutative bilinear formulas
 *       ({@link RosowskiBound#bestCommutativeBound}) — includes the
 *       {@code ⟨3,3,3⟩=21} non-bilinear bound via Theorems 2/3.
 *       The actual scheme (non-bilinear) is in
 *       {@code references/rosowski-algorithms.md}; the bound is
 *       represented here as a formula.</li>
 *   <li>Pan TA closed-form ({@link PanTrilinearAggregation#cubicBound})</li>
 * </ul>
 *
 * <p>DIS09 Table 4 baseline: commutative bilinear algorithms only —
 * the values are the best commutative bilinear rank known in 2009.</p>
 */
@Tag("catalog-iterating")
public class TestDIS09Table4Commutative {

	/** DIS09 Table 4 — COMMUTATIVE ⟨n,n,n⟩, n=2..30. */
	private static final Map<Integer, Integer> DIS09_CMT = new TreeMap<>(Map.ofEntries(
			Map.entry(2, 7),
			Map.entry(3, 22),
			Map.entry(4, 46),
			Map.entry(5, 93),
			Map.entry(6, 141),
			Map.entry(7, 235),
			Map.entry(8, 316),
			Map.entry(9, 472),
			Map.entry(10, 595),
			Map.entry(11, 825),
			Map.entry(12, 987),
			Map.entry(13, 1318),
			Map.entry(14, 1525),
			Map.entry(15, 1941),
			Map.entry(16, 2212),
			Map.entry(17, 2762),
			Map.entry(18, 3060),
			Map.entry(19, 3757),
			Map.entry(20, 4158),
			Map.entry(21, 4938),
			Map.entry(22, 5440),
			Map.entry(23, 6382),
			Map.entry(24, 6900),
			Map.entry(25, 8083),
			Map.entry(26, 8658),
			Map.entry(27, 9994),
			Map.entry(28, 10556),
			Map.entry(29, 12109),
			Map.entry(30, 12710)));

	private static List<BlockSplitSearch.NamedBase> buildPool() throws Exception {
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		NonCubicBilinearAlgorithm laderman = NonCubicBilinearAlgorithm.fromCubic(Laderman23.get());
		NonCubicBilinearAlgorithm hopcroft = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section3/alphatensor_Z-2x3x3_m15_a58.json"));
		List<BlockSplitSearch.NamedBase> pool = new ArrayList<>();
		expandS3(pool, "Strassen", strassen);
		expandS3(pool, "Laderman", laderman);
		expandS3(pool, "Hopcroft-Kerr ⟨2,3,3⟩=15", hopcroft);
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
					label + " [σ" + i + "]", a));
		}
	}

	/**
	 * Commutative SOTA resolver: returns the min over Waksman,
	 * Rosowski commutative, Pan TA, and the NC catalog.
	 */
	private static Recombination.SotaResolver buildCommutativeSota() {
		Map<String, Integer> ncRanks = BlockSplitSearch.loadCatalogBestRanksForField("R");
		Function<int[], Optional<Integer>> ncLookup = BlockSplitSearch.rankLookupFromMap(ncRanks);
		return (a, b, c) -> {
			if (a == 0 || b == 0 || c == 0) return 0;
			if (a == 1) return b * c;
			if (b == 1) return a * c;
			if (c == 1) return a * b;

			long best = Long.MAX_VALUE;

			// NC catalog (NC ⊆ commutative as upper bound)
			int[] key = { a, b, c };
			Optional<Integer> nc = ncLookup.apply(key);
			if (nc.isPresent()) best = Math.min(best, nc.get());

			// Waksman family — min over the three axis orientations
			best = Math.min(best, WaksmanBound.forShape(a, b, c));
			best = Math.min(best, WaksmanBound.forShape(b, a, c));
			best = Math.min(best, WaksmanBound.forShape(a, c, b));

			// Rosowski commutative — min over all 6 axis permutations
			Optional<Long> ros = RosowskiBound.bestCommutativeBound(a, b, c);
			if (ros.isPresent()) best = Math.min(best, ros.get());

			// Pan TA (NC, applies to commutative as upper bound too) — only cubic
			if (a == b && b == c) {
				best = Math.min(best, PanTrilinearAggregation.cubicBound(a));
			}

			return best > Integer.MAX_VALUE / 100 ? Integer.MAX_VALUE / 100 : (int) best;
		};
	}

	@Test
	public void dis09_table4_commutative_R_4_to_30() throws Exception {
		Recombination.SotaResolver sota = buildCommutativeSota();
		List<BlockSplitSearch.NamedBase> pool = buildPool();

		System.out.println();
		System.out.println("DIS09 Table 4 (COMMUTATIVE) vs our commutative SOTA "
				+ "(Waksman + Rosowski + Pan TA + NC catalog)");
		System.out.println("=".repeat(95));
		System.out.printf("%5s | %7s | %7s | %6s | %6s | %s%n",
				"n", "DIS09", "ours", "Δ", "Δ%", "via");
		System.out.println("-".repeat(95));

		int wins = 0, ties = 0, losses = 0;
		long totalDis09 = 0, totalOurs = 0;
		for (int n = 4; n <= 30; n++) {
			int dis09 = DIS09_CMT.get(n);
			Optional<BlockSplitSearch.MultiBaseSplitCandidate> best =
					BlockSplitSearch.findBestMultiBaseSplit(n, n, n, pool, sota, false, 1_000_000L);
			long recursive = best.map(BlockSplitSearch.MultiBaseSplitCandidate::rank).orElse(Long.MAX_VALUE);
			long direct = sota.getRank(n, n, n);
			long ours = Math.min(recursive, direct);
			String via;
			if (direct < recursive) via = "DIRECT formula = " + direct;
			else if (best.isPresent()) via = best.get().breakdown();
			else via = "—";

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
