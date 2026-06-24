package eu.solven.matmul.docs.explore;

import eu.solven.matmul.recombination.Recombination;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.Compose;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.ProjectionSearch;
import lombok.extern.slf4j.Slf4j;

/**
 * THE BRUTAL CASE (user 2026-06-12): build ⟨32,32,32⟩ parents with as much
 * projection margin as possible and cascade them by chained DCE-projection
 * onto EVERY shape ⟨n,m,p⟩, 2 ≤ n,m,p ≤ 32, comparing against the catalog.
 *
 * <p>Method: optimal multi-drop projection is exponential, but a GIVEN mask is
 * exact and O(rank) ({@link ProjectionSearch.MaskedProjector}). So we run a
 * DAG-greedy sweep over the drop-count lattice: node (n,m,p) extends the best
 * of its ≤3 parent nodes (n+1,m,p)/(n,m+1,p)/(n,m,p+1) by the single kept
 * index whose removal leaves the fewest survivors. Every node's count is the
 * EXACT rank of a concrete constructible projection (bound: greedy mask
 * choice, label accordingly); wins vs {@code findRank} are materializable via
 * {@code Compose.project} on the recorded masks.</p>
 *
 * <p>Parents: Strassen⁵ and Winograd⁵ towers (their margin profiles differ —
 * Winograd's single-block products are the padding-friendly ones), plus the
 * catalog's own ⟨32,32,32⟩ when readable.</p>
 *
 * <p>Args: {@code [minDim=2]}.</p>
 */
@Slf4j
public class ScanCubeProjectionCascade {

	public static void main(String[] args) {
		int minDim = args.length > 0 ? Integer.parseInt(args[0]) : 2;
		FieldAwareLookup q = new FieldAwareLookup(Field.Q);

		Map<String, NonCubicBilinearAlgorithm> parents = new LinkedHashMap<>();
		q.findByHash(2, 2, 2, "db11bcc")
				.ifPresent(ws -> parents.put("strassen^5", tower(ws.alg(), 5)));
		q.findByHash(2, 2, 2, "511df05")
				.ifPresent(ws -> parents.put("winograd^5", tower(ws.alg(), 5)));
		q.find(32, 32, 32).ifPresent(alg -> parents.put("catalog-32cube", alg));

		for (Map.Entry<String, NonCubicBilinearAlgorithm> e : parents.entrySet()) {
			cascade(e.getKey(), e.getValue(), q, minDim);
		}
	}

	private static NonCubicBilinearAlgorithm tower(NonCubicBilinearAlgorithm base, int k) {
		NonCubicBilinearAlgorithm t = base;
		for (int i = 1; i < k; i++) {
			t = Compose.kroneckerGeneral(t, base);
		}
		return t;
	}

	record Node(boolean[] mn, boolean[] mm, boolean[] mp, long survivors) {}

	private static void cascade(String label, NonCubicBilinearAlgorithm parent,
			FieldAwareLookup q, int minDim) {
		long t0 = System.currentTimeMillis();
		log.info("{}: rank={} axisMargins={}", label, parent.r,
				Arrays.toString(ProjectionSearch.axisMargins(parent)));
		ProjectionSearch.MaskedProjector proj = new ProjectionSearch.MaskedProjector(parent);
		int dim = parent.n;  // cubic parent
		Map<Long, Node> lattice = new LinkedHashMap<>();
		boolean[] full = new boolean[dim];
		Arrays.fill(full, true);
		// THREE DISTINCT arrays — sharing one `full` array across the axes let a
		// single-axis bit-clear corrupt all three masks (the 2026-06-12 phantom-win
		// bug; caught by the harvest's write-boundary shape check).
		lattice.put(key(dim, dim, dim), new Node(full.clone(), full.clone(), full.clone(),
				parent.r));

		int wins = 0;
		int ties = 0;
		int losses = 0;
		List<String> winLines = new ArrayList<>();
		// BFS layers by total drops: every node's ≤3 lattice parents are ready.
		for (int drops = 1; drops <= 3 * (dim - minDim); drops++) {
			for (int n = minDim; n <= dim; n++) {
				for (int m = minDim; m <= dim; m++) {
					int p = 3 * dim - drops - n - m;
					if (p < minDim || p > dim) {
						continue;
					}
					Node best = null;
					for (int ax = 0; ax < 3; ax++) {
						int pn = n + (ax == 0 ? 1 : 0);
						int pm = m + (ax == 1 ? 1 : 0);
						int pp = p + (ax == 2 ? 1 : 0);
						if (pn > dim || pm > dim || pp > dim) {
							continue;
						}
						Node from = lattice.get(key(pn, pm, pp));
						if (from == null) {
							continue;
						}
						boolean[] axMask = ax == 0 ? from.mn() : ax == 1 ? from.mm() : from.mp();
						for (int i = 0; i < dim; i++) {
							if (!axMask[i]) {
								continue;
							}
							axMask[i] = false;
							long ceiling = best == null ? Long.MAX_VALUE / 4 : best.survivors();
							long s = proj.survivors(from.mn(), from.mm(), from.mp(), ceiling);
							if (best == null || s < best.survivors()) {
								best = new Node(from.mn().clone(), from.mm().clone(),
										from.mp().clone(), s);
							}
							axMask[i] = true;
						}
					}
					if (best == null) {
						continue;
					}
					// Fail-loud lattice invariant (the phantom-win bug's guard).
					if (count(best.mn()) != n || count(best.mm()) != m
							|| count(best.mp()) != p) {
						throw new IllegalStateException(String.format(
								"lattice mask/key drift at ⟨%d,%d,%d⟩: masks %d/%d/%d",
								n, m, p, count(best.mn()), count(best.mm()), count(best.mp())));
					}
					lattice.put(key(n, m, p), best);
					int cat = q.findRank(n, m, p);
					if (cat < Recombination.SotaResolver.UNKNOWN_RANK) {
						if (best.survivors() < cat) {
							wins++;
							winLines.add(String.format("⟨%d,%d,%d⟩ cascade %d < catalog %d",
									n, m, p, best.survivors(), cat));
						} else if (best.survivors() == cat) {
							ties++;
						} else {
							losses++;
						}
					}
				}
			}
			if (drops % 10 == 0) {
				log.info("[progress] {} layer {}/{}: {} nodes, {} wins/{} ties/{} losses, {}s",
						label, drops, 3 * (dim - minDim), lattice.size(), wins, ties, losses,
						(System.currentTimeMillis() - t0) / 1000);
			}
		}
		for (String w : winLines.subList(0, Math.min(30, winLines.size()))) {
			log.info("  WIN {} — materializable via Compose.project on the recorded masks", w);
		}
		log.info("{} cascade done: {} nodes, {} WINS / {} ties / {} losses vs catalog, {}s "
				+ "(greedy-mask bound — wins are exact constructible schemes)",
				label, lattice.size(), wins, ties, losses,
				(System.currentTimeMillis() - t0) / 1000);
	}

	private static int count(boolean[] mask) {
		int c = 0;
		for (boolean b : mask) {
			if (b) {
				c++;
			}
		}
		return c;
	}

	private static long key(int n, int m, int p) {
		return ((long) n << 20) | ((long) m << 10) | p;
	}
}
