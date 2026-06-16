package eu.solven.matmul.search;

import java.io.File;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.search.AnalyticalMaskSearch.MaskCandidate;

/**
 * Layer 1 from {@code docs/notes/enumerating-rank-7-2x2x2-schemes.md}:
 *
 * <p>For every rank-7 ⟨2,2,2⟩ scheme in the catalog × all 8 axis-flip masks ×
 * unbalanced cubic targets n ∈ {17, 19, 21, 23, 25, 27} × cubic-allocations
 * (n₀, n₁) per axis, score the analytical shape multiset cost. Report which
 * (scheme, mask, allocation) achieves the best per-target rank.
 *
 * <p>This is a diagnostic test: it does NOT assert specific values (those
 * would change as the catalog evolves). It dumps a table to stdout so the
 * user can audit what Layer 1 actually discovers.
 */
class TestLayer1OrbitSweep {

	private static final int[] TARGETS = {17, 19, 21, 23, 25, 27};
	private static final String SECTION2 = "src/main/resources/schemes/known/section2/";
	// Non-commutative ⟨2,2,2⟩=7 schemes, one per distinct discrete orbit
	// (see TestCatalog2x2x2DiscreteOrbits). solven-winograd-cousin-axflip1
	// is dropped — orbit-equivalent to winograd-1971; fmm-lille-2x2x2 (a=18)
	// likewise dropped — orbit-equivalent to strassen-1969 and deleted from
	// the catalog as redundant. AT-F2 (F2 only), Waksman (commutative) excluded.
	private static final String[] CATALOG_2X2X2_R7 = {
			"strassen-2x2x2_m7_a18.json",
			"winograd_1971-2x2x2_m7_a24.json",
			"alphatensor_Z-2x2x2_m7_a22.json",
			"perminov_cr15_cn24_ZT_reduced-2x2x2_m7_a24.json",
	};

	@Test
	void sweep_unbalanced_cubics_dumps_winners() throws Exception {
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		CitedBound sota = new CitedBound(lookup);

		System.out.println();
		System.out.println("==============================================================");
		System.out.println(" Layer 1 sweep — catalog rank-7 ⟨2,2,2⟩ × axis-flip × cubic n");
		System.out.println("==============================================================");
		System.out.println();

		TreeMap<Integer, TargetWinner> winners = new TreeMap<>();
		for (int n : TARGETS) {
			TargetWinner tw = new TargetWinner(n);
			for (String scheme : CATALOG_2X2X2_R7) {
				NonCubicBilinearAlgorithm alg = SchemeIO.readBilinear(new File(SECTION2 + scheme));
				for (int n0 = (n + 1) / 2; n0 <= n - 1; n0++) {
					int n1 = n - n0;
					if (n1 <= 0) continue;
					int[] alloc = {n0, n1};
					List<MaskCandidate> top = AnalyticalMaskSearch.topKMasks(
							alg, alloc, alloc, alloc, sota::getRank, 1);
					if (top.isEmpty()) continue;
					MaskCandidate best = top.get(0);
					tw.record(scheme, alloc, best.mask, best.cost, best.shapes);
				}
			}
			winners.put(n, tw);
		}

		System.out.printf("%-4s | %-8s | %-50s | %-10s | %-4s | %-8s | %s%n",
				"n", "sota", "scheme", "alloc", "mask", "cost", "shape-multiset");
		System.out.println("-----|----------|---------" + "-".repeat(45) + "|------------|------|----------|-------------");
		for (var e : winners.entrySet()) {
			int n = e.getKey();
			long sotaN = sota.getRank(n, n, n);
			TargetWinner tw = e.getValue();
			for (Winner w : tw.distinctTop(3)) {
				String marker = w.cost < sotaN ? "  *NEW*" : (w.cost == sotaN ? "  =" : "");
				System.out.printf("%-4d | %-8d | %-50s | (%d,%d)%-6s | %-4d | %-8d%s | %s%n",
						n, sotaN, trim(w.scheme),
						w.alloc[0], w.alloc[1], "",
						w.mask, w.cost, marker, multisetCounts(w.shapes));
			}
			System.out.println();
		}
		System.out.println();
		System.out.println("==============================================================");
	}

	private static String trim(String s) {
		String base = s.replace(".json", "").replace("_2x2x2_m7", "");
		if (base.length() > 50) base = base.substring(0, 47) + "...";
		return base;
	}

	private static String multisetCounts(int[][] shapes) {
		TreeMap<String, Integer> counts = new TreeMap<>();
		for (int[] s : shapes) {
			String key = "⟨" + s[0] + "," + s[1] + "," + s[2] + "⟩";
			counts.merge(key, 1, Integer::sum);
		}
		StringBuilder sb = new StringBuilder();
		boolean first = true;
		for (var e : counts.entrySet()) {
			if (!first) sb.append(" + ");
			sb.append(e.getValue()).append("·").append(e.getKey());
			first = false;
		}
		return sb.toString();
	}

	private static final class Winner {
		final String scheme;
		final int[] alloc;
		final int mask;
		final long cost;
		final int[][] shapes;

		Winner(String scheme, int[] alloc, int mask, long cost, int[][] shapes) {
			this.scheme = scheme;
			this.alloc = alloc;
			this.mask = mask;
			this.cost = cost;
			this.shapes = shapes;
		}
	}

	private static final class TargetWinner {
		final int n;
		final List<Winner> all = new ArrayList<>();

		TargetWinner(int n) {
			this.n = n;
		}

		void record(String scheme, int[] alloc, int mask, long cost, int[][] shapes) {
			all.add(new Winner(scheme, alloc, mask, cost, shapes));
		}

		/** Top-K distinct (cost, scheme) — we collapse mask/alloc redundancy. */
		List<Winner> distinctTop(int k) {
			all.sort(Comparator.comparingLong(w -> w.cost));
			List<Winner> result = new ArrayList<>();
			Set<String> seenSchemes = new HashSet<>();
			long bestCost = Long.MAX_VALUE;
			for (Winner w : all) {
				if (bestCost == Long.MAX_VALUE) bestCost = w.cost;
				// Always include all winners at the best cost
				if (w.cost == bestCost) {
					if (!seenSchemes.contains(w.scheme)) {
						result.add(w);
						seenSchemes.add(w.scheme);
					}
				} else if (result.size() < k && !seenSchemes.contains(w.scheme)) {
					result.add(w);
					seenSchemes.add(w.scheme);
				}
				if (result.size() >= k && w.cost > bestCost) break;
			}
			return result;
		}
	}
}
