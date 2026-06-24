package eu.solven.matmul.search;

import eu.solven.matmul.recombination.RecombinationMultisetOrbit;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.List;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Guards {@link RecombinationMultisetOrbit.Result#dominanceFrontier()}: it must
 * be a strict subset of the canonical multisets AND — the load-bearing property —
 * it must be <em>lossless</em>, i.e. the actual rank minimum (over a concrete
 * allocation, via {@link FieldAwareLookup#findRank}) is always attained on a
 * frontier member. A regression here would silently drop the optimal multiset.
 */
public class TestRecombinationDominance {

	private static NonCubicBilinearAlgorithm base(int n, int m, int p, String token) throws Exception {
		FieldAwareLookup lk = new FieldAwareLookup(Field.Q);
		for (var path : lk.findFiles(n, m, p)) {
			if (path.getFileName().toString().contains(token)) {
				return SchemeIO.read(path.toFile());
			}
		}
		throw new IllegalStateException("no " + token + " base at " + n + "x" + m + "x" + p);
	}

	@Test
	public void strassen_222_frontier_is_6_and_subset() throws Exception {
		var res = RecombinationMultisetOrbit.enumerate(base(2, 2, 2, "strassen"), 3);
		List<String> frontier = res.dominanceFrontier();

		assertThat(res.canonicalMultisets).hasSize(40);
		assertThat(frontier).hasSize(6);
		assertThat(res.canonicalMultisets).containsAll(frontier);
	}

	@Test
	public void frontier_is_seed_independent_222() throws Exception {
		// de Groote: rank-7 ⟨2,2,2⟩ is a single orbit, so the frontier must NOT depend
		// on which representative seeds the sweep (the earlier bug compared seed-dependent
		// representativeShapes and returned 8 vs 7).
		var s = RecombinationMultisetOrbit.enumerate(base(2, 2, 2, "strassen"), 3);
		var w = RecombinationMultisetOrbit.enumerate(base(2, 2, 2, "winograd_1971"), 3);
		assertThat(s.dominanceFrontier()).hasSize(6);
		assertThat(w.dominanceFrontier()).hasSize(6);
	}

	@Test
	public void frontier_contains_the_rank_minimum_222_at_9_8() throws Exception {
		// Concrete allocation: each axis split 17 = 9 + 8 (block index 0 = the 9).
		int[][] sizes = { { 9, 8 }, { 9, 8 }, { 9, 8 } };
		FieldAwareLookup lk = new FieldAwareLookup(Field.Q);

		var res = RecombinationMultisetOrbit.enumerate(base(2, 2, 2, "strassen"), 3);
		List<String> frontier = res.dominanceFrontier();

		long best = Long.MAX_VALUE;
		String argmin = null;
		for (String key : res.canonicalMultisets) {
			long total = 0;
			for (int[] t : res.representativeShapes.get(key)) {
				int a = sizes[0][t[0]], b = sizes[1][t[1]], c = sizes[2][t[2]];
				total += lk.findRank(a, b, c);
			}
			if (total < best) {
				best = total;
				argmin = key;
			}
		}
		// The rank-optimal multiset must survive the (rank-agnostic) dominance prune.
		assertThat(frontier).contains(argmin);
	}

	@Test
	public void frontier_strictly_prunes_2x2x3() throws Exception {
		// ⟨2,2,3⟩ has thousands of canonical multisets; the frontier must be a tiny
		// fraction — and never larger than the full set.
		var res = RecombinationMultisetOrbit.enumerate(base(2, 2, 3, "alphatensor_Z"), 2);
		List<String> frontier = res.dominanceFrontier();

		assertThat(res.canonicalMultisets.size()).isGreaterThan(1000);
		assertThat(frontier.size()).isLessThan(res.canonicalMultisets.size() / 10);
		assertThat(res.canonicalMultisets).containsAll(frontier);
	}
}
