package eu.solven.matmul.recombination;

import java.util.ArrayList;
import java.util.List;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.SymmetryTransforms;
import eu.solven.matmul.recombination.Recombination.SotaResolver;
import eu.solven.matmul.recombination.RecombinationMultisetOrbit.Result;

/**
 * Score a target ⟨N,M,P⟩ against the whole GL-orbit + axis-orientation menu of recombination
 * supports for a set of bases, and return the cheapest. This is the productised
 * {@code ProbeFrontierRecomb} logic: for each base, each axis ORIENTATION ({@code s3Orbit} — which
 * base axis tiles which target axis, the lever the default {@code AXIS_FLIP} pool only partly
 * covers), and each GL-orbit dominance-frontier multiset, it minimises the allocation in multiset
 * space. The winner's concrete scheme is rebuilt on demand via
 * {@link RecombinationMultisetOrbit#materialise} (no per-orientation scheme rebuild during scoring).
 *
 * <p>The allocation must use <b>descending</b> partitions: the multiset's "index 0 = largest block"
 * convention requires sorted block sizes, else a big-index product could be placed (invalidly) on a
 * tiny block.
 */
public final class FrontierRecombination {
	private FrontierRecombination() {}

	/** The cheapest recombination found: rank, the oriented base, the GL transform realising the
	 *  winning multiset (null for a structural d≥4 winner — scoring only), and the per-axis allocation. */
	public record Best(long rank, NonCubicBilinearAlgorithm orientedBase, int[][][] transform,
			int[] allocA, int[] allocB, int[] allocC, String multiset) {
		/** Rebuild the concrete winning scheme (orient + GL transform); requires {@link #transform} != null. */
		public NonCubicBilinearAlgorithm materialiseWinner() {
			if (transform == null) {
				throw new IllegalStateException("structural winner has no captured transform; re-derive via GL");
			}
			return RecombinationMultisetOrbit.materialise(orientedBase, transform[0], transform[1], transform[2]);
		}
	}

	/** Best recombination of {@code target} over all (base, orientation, frontier-multiset, allocation). */
	public static Best bestFor(int N, int M, int P, List<NonCubicBilinearAlgorithm> bases, SotaResolver sota) {
		Best best = null;
		for (NonCubicBilinearAlgorithm base : bases) {
			int maxDim = Math.max(base.n, Math.max(base.m, base.p));
			for (NonCubicBilinearAlgorithm ob : SymmetryTransforms.s3Orbit(base)) {
				Result orbit = maxDim <= 3 ? RecombinationMultisetOrbit.enumerate(ob, 2)
						: RecombinationMultisetOrbit.enumerateStructuralFrontier(ob, 1);
				List<int[]> partA = descendingPartitions(N, ob.n);
				List<int[]> partB = descendingPartitions(M, ob.m);
				List<int[]> partC = descendingPartitions(P, ob.p);
				if (partA.isEmpty() || partB.isEmpty() || partC.isEmpty()) continue; // base too coarse for target
				for (String key : orbit.dominanceFrontier()) {
					int[][] ms = orbit.representativeShapes.get(key);
					if (ms == null) continue;
					AllocCost ac = bestAlloc(ms, partA, partB, partC, sota);
					if (best == null || ac.cost < best.rank) {
						best = new Best(ac.cost, ob, orbit.representativeTransforms.get(key), ac.a, ac.b, ac.c, key);
					}
				}
			}
		}
		return best;
	}

	/**
	 * Best rank of a PRECOMPUTED frontier over all 6 axis ORIENTATIONS, derived by permuting the
	 * multiset's columns (no re-enumeration) — the efficient, index-friendly path. {@code dims} =
	 * the base's {@code (bn,bm,bp)}; {@code frontier} = block-index {@code [r][3]} arrays, one per
	 * frontier multiset (from {@link RecombinationMultisetOrbit.Result#representativeShapes}). For
	 * orientation π, target axis j is tiled by base axis π[j], so the per-product index on axis j is
	 * {@code ms[k][π[j]]} and the allocation has {@code dims[π[j]]} parts.
	 */
	public static long bestRankOverOrientations(int N, int M, int P, int[] dims, List<int[][]> frontier, SotaResolver sota) {
		long best = Long.MAX_VALUE;
		int[][] allPerms = { { 0, 1, 2 }, { 0, 2, 1 }, { 1, 0, 2 }, { 1, 2, 0 }, { 2, 0, 1 }, { 2, 1, 0 } };
		// Quotient orientations by the base's shape stabilizer: two perms with the SAME per-target
		// block-count tuple (dims[π]) are equivalent for cost — ⟨2,3,3⟩ has 3 distinct, ⟨2,2,2⟩ has 1,
		// not the naive 6. (Distinct dims-tuple captures it because equal-dim axes are interchangeable.)
		java.util.Set<String> seenOrient = new java.util.HashSet<>();
		List<int[]> perms = new ArrayList<>();
		for (int[] pi : allPerms)
			if (seenOrient.add(dims[pi[0]] + "," + dims[pi[1]] + "," + dims[pi[2]])) perms.add(pi);
		for (int[] pi : perms) {
			List<int[]> pa = descendingPartitions(N, dims[pi[0]]);
			List<int[]> pb = descendingPartitions(M, dims[pi[1]]);
			List<int[]> pc = descendingPartitions(P, dims[pi[2]]);
			if (pa.isEmpty() || pb.isEmpty() || pc.isEmpty()) continue; // base too coarse for target on this orientation
			for (int[][] ms : frontier) {
				for (int[] a : pa) {
					long lbA = 0;
					for (int[] s : ms) lbA += sota.getRank(a[s[pi[0]]], 1, 1);
					if (lbA >= best) continue;
					for (int[] b : pb) for (int[] c : pc) {
						long tot = 0;
						for (int[] s : ms) {
							tot += sota.getRank(a[s[pi[0]]], b[s[pi[1]]], c[s[pi[2]]]);
							if (tot >= best) break;
						}
						if (tot < best) best = tot;
					}
				}
			}
		}
		return best;
	}

	private record AllocCost(long cost, int[] a, int[] b, int[] c) {}

	/** Min Σ R(a[idx_n], b[idx_m], c[idx_p]) over descending partitions, returning the winning allocation. */
	private static AllocCost bestAlloc(int[][] ms, List<int[]> partA, List<int[]> partB, List<int[]> partC, SotaResolver sota) {
		long best = Long.MAX_VALUE;
		int[] bestA = null, bestB = null, bestC = null;
		for (int[] a : partA) {
			long lbA = 0; // a-only lower bound (b=c=1), R monotone ⇒ admissible prune
			for (int[] s : ms) lbA += sota.getRank(a[s[0]], 1, 1);
			if (lbA >= best) continue;
			for (int[] b : partB) {
				for (int[] c : partC) {
					long tot = 0;
					for (int[] s : ms) {
						tot += sota.getRank(a[s[0]], b[s[1]], c[s[2]]);
						if (tot >= best) break;
					}
					if (tot < best) { best = tot; bestA = a; bestB = b; bestC = c; }
				}
			}
		}
		return new AllocCost(best, bestA, bestB, bestC);
	}

	/** All non-increasing partitions of {@code total} into exactly {@code parts} positive parts. */
	static List<int[]> descendingPartitions(int total, int parts) {
		List<int[]> out = new ArrayList<>();
		if (total < parts) return out; // cannot give each part ≥1
		dp(total, parts, total, new int[parts], 0, out);
		return out;
	}

	private static void dp(int total, int parts, int maxPart, int[] cur, int idx, List<int[]> out) {
		if (parts == 0) { if (total == 0) out.add(cur.clone()); return; }
		int hi = Math.min(maxPart, total - (parts - 1));
		int lo = (total + parts - 1) / parts;
		for (int v = hi; v >= lo; v--) { cur[idx] = v; dp(total - v, parts - 1, v, cur, idx + 1, out); }
	}
}
