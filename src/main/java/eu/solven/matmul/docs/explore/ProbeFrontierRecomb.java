package eu.solven.matmul.docs.explore;

import java.io.File;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.SymmetryTransforms;
import eu.solven.matmul.algebra.Algebra;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.recombination.Recombination;
import eu.solven.matmul.recombination.Recombination.SotaResolver;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.recombination.AllocationOptimizer;
import eu.solven.matmul.recombination.AnalyticalMaskSearch.SchemeSupports;
import eu.solven.matmul.isotropy.RecombinationMultisetOrbit;
import eu.solven.matmul.isotropy.RecombinationMultisetOrbit.Result;
import eu.solven.matmul.search.SearchBudget;

/**
 * Demonstrate the frontier-of-supports idea: for a base scheme recombined onto a target
 * ⟨N,M,P⟩, optimise the allocation over EVERY GL-orbit dominance-frontier support (not just the
 * base's native one) and keep the cheapest. Compares native-support best vs frontier best vs the
 * catalog's current ⟨N,M,P⟩ rank — so a win here is a concrete gap reduction.
 *
 * <p>Args: {@code <baseSchemeFile> <N> <M> <P> [glBound] [structural]}. Uses the exact GL frontier
 * by default (tractable for base max-dim ≤ 3); pass {@code structural} to use the partial d≥4
 * frontier instead.
 */
public final class ProbeFrontierRecomb {
	private ProbeFrontierRecomb() {}

	public static void main(String[] args) throws Exception {
		File f = new File(args[0]);
		int N = Integer.parseInt(args[1]), M = Integer.parseInt(args[2]), P = Integer.parseInt(args[3]);
		int glBound = args.length > 4 ? Integer.parseInt(args[4]) : 2;
		boolean structural = args.length > 5 && args[5].equals("structural");
		NonCubicBilinearAlgorithm base = SchemeIO.read(f);
		SotaResolver sota = Recombination.catalogResolver(Algebra.nonCommutative(Field.R));

		System.out.printf("base %s ⟨%d,%d,%d⟩ r=%d  →  target ⟨%d,%d,%d⟩%n",
				f.getName(), base.n, base.m, base.p, base.r, N, M, P);
		long catalog = sota.getRank(N, M, P);
		System.out.printf("  catalog R⟨%d,%d,%d⟩ = %d%n", N, M, P, catalog);

		// 1) Native support — the current single-support approach.
		long t0 = System.nanoTime();
		var nat = AllocationOptimizer.optimize(SchemeSupports.extract(base), sota, N, M, P, SearchBudget.EXACT, null);
		System.out.printf("  NATIVE support : rank=%d  alloc=%s/%s/%s  (%.1fs)%n",
				nat.rank(), str(nat.allocA()), str(nat.allocB()), str(nat.allocC()), (System.nanoTime() - t0) / 1e9);

		// 2) ORBIT × FRONTIER — scan every axis ORIENTATION (s3Orbit: which base axis maps to which
		// target axis — the "orbitFlip" lever the default search skips) AND, per orientation, every
		// GL-orbit frontier support. All cheap, in multiset space. Keep the global cheapest.
		t0 = System.nanoTime();
		long best = Long.MAX_VALUE;
		String bestKey = null, bestOrient = null;
		int orientCount = 0, totalSupports = 0;
		for (NonCubicBilinearAlgorithm ob : SymmetryTransforms.s3Orbit(base)) {
			orientCount++;
			Result orbit = structural ? RecombinationMultisetOrbit.enumerateStructuralFrontier(ob, 1)
					: RecombinationMultisetOrbit.enumerate(ob, glBound);
			var frontier = orbit.dominanceFrontier();
			// Descending partitions, sized to THIS orientation's per-axis slot counts.
			var partA = descendingPartitions(N, ob.n);
			var partB = descendingPartitions(M, ob.m);
			var partC = descendingPartitions(P, ob.p);
			for (String key : frontier) {
				int[][] ms = orbit.representativeShapes.get(key);
				if (ms == null) continue;
				totalSupports++;
				long r = bestSortedCost(ms, partA, partB, partC, sota);
				if (r < best) { best = r; bestKey = key; bestOrient = ob.n + "x" + ob.m + "x" + ob.p; }
			}
		}
		System.out.printf("  ORBIT×FRONTIER (%d orientations, %d supports): rank=%d  (%.1fs)%n",
				orientCount, totalSupports, best, (System.nanoTime() - t0) / 1e9);
		if (bestKey != null)
			System.out.printf("    via orientation ⟨%s⟩ multiset: %s%n", bestOrient,
					RecombinationMultisetOrbit.prettySymbolic(bestKey, "n", "m", "p"));
		System.out.printf("  >>> native=%d  orbit=%d  catalog=%d  %s%n",
				nat.rank(), best, catalog,
				best < nat.rank() ? "ORBIT WINS by " + (nat.rank() - best) : "no gain over native");
	}

	/** Min Σ R(allocA[idx_n], allocB[idx_m], allocC[idx_p]) over descending partitions per axis. */
	private static long bestSortedCost(int[][] ms, java.util.List<int[]> partA, java.util.List<int[]> partB,
			java.util.List<int[]> partC, SotaResolver sota) {
		long best = Long.MAX_VALUE;
		for (int[] a : partA) {
			long lbA = 0; // a-only lower bound: every product at its a-size, b=c=1 (R monotone)
			for (int[] s : ms) lbA += sota.getRank(a[s[0]], 1, 1);
			if (lbA >= best) continue;
			for (int[] b : partB) {
				for (int[] c : partC) {
					long tot = 0;
					for (int[] s : ms) {
						tot += sota.getRank(a[s[0]], b[s[1]], c[s[2]]);
						if (tot >= best) break;
					}
					if (tot < best) best = tot;
				}
			}
		}
		return best;
	}

	/** All non-increasing partitions of {@code total} into exactly {@code parts} positive parts. */
	private static java.util.List<int[]> descendingPartitions(int total, int parts) {
		java.util.List<int[]> out = new java.util.ArrayList<>();
		dp(total, parts, total, new int[parts], 0, out);
		return out;
	}

	private static void dp(int total, int parts, int maxPart, int[] cur, int idx, java.util.List<int[]> out) {
		if (parts == 0) { if (total == 0) out.add(cur.clone()); return; }
		int hi = Math.min(maxPart, total - (parts - 1)); // leave ≥1 for each remaining part
		int lo = (total + parts - 1) / parts; // each remaining part ≤ v ⇒ v ≥ ⌈total/parts⌉
		for (int v = hi; v >= lo; v--) { cur[idx] = v; dp(total - v, parts - 1, v, cur, idx + 1, out); }
	}

	private static String str(int[] a) {
		StringBuilder sb = new StringBuilder("[");
		for (int i = 0; i < a.length; i++) { if (i > 0) sb.append(','); sb.append(a[i]); }
		return sb.append(']').toString();
	}
}
