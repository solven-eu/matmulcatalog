package eu.solven.matmul.docs.migrate;

import eu.solven.matmul.recombination.Recombination;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.CatalogPolicy;
import eu.solven.matmul.catalog.Compose;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.ProjectionSearch;
import eu.solven.matmul.catalog.SchemeIO;
import lombok.extern.slf4j.Slf4j;

/**
 * HARVEST of the 32-cube projection cascade (user "Go", 2026-06-12; finder:
 * {@code docs.explore.ScanCubeProjectionCascade} — Strassen⁵ alone beat the
 * catalog on 8,155 of 29,791 shapes ≤32). Re-runs the DAG-greedy chained-DCE
 * lattice for ONE parent and persists every strict catalog win through the
 * standard pipeline: {@code Compose.project} materialisation, rank +
 * spot-check verification, {@code Lineage.Project(parentLineage, keeps)},
 * {@code SchemeIO.write}/{@code writeStub} per {@code MATERIALISE_MAX_DIM},
 * under {@code derived/sectionN} (deterministic replayable output).
 *
 * <p>The tower parent itself is persisted as a lineage stub under
 * {@code margin-bases/section32} — rank-worse than the catalog 32-cube but the
 * strongest projecting base we have (the folder's exact Pareto criterion),
 * stamped with its target-free {@code projection_margins} triple.</p>
 *
 * <p>One parent per JVM invocation (≈15 min each — also dodges the
 * environment's long-JVM reaper): {@code --parent=strassen|winograd|catalog32}.
 * Runs are order-independent and idempotent-ish: each strictly gates against
 * the CURRENT catalog rank, so later runs only write further improvements.</p>
 */
@Slf4j
public class MaterialiseCubeCascade {

	public static void main(String[] args) {
		String which = args.length > 0 ? args[0].replace("--parent=", "") : "strassen";
		int minDim = args.length > 1 ? Integer.parseInt(args[1]) : 2;
		FieldAwareLookup q = new FieldAwareLookup(Field.Q);

		NonCubicBilinearAlgorithm parentAlg;
		Lineage.Node parentNode;
		switch (which) {
		case "strassen": {
			NonCubicBilinearAlgorithm base = q.findByHash(2, 2, 2, "db11bcc").orElseThrow().alg();
			parentAlg = tower(base, 5);
			parentNode = new Lineage.KronChain(List.of(atom222("db11bcc"), atom222("db11bcc"),
					atom222("db11bcc"), atom222("db11bcc"), atom222("db11bcc")));
			break;
		}
		case "winograd": {
			NonCubicBilinearAlgorithm base = q.findByHash(2, 2, 2, "511df05").orElseThrow().alg();
			parentAlg = tower(base, 5);
			parentNode = new Lineage.KronChain(List.of(atom222("511df05"), atom222("511df05"),
					atom222("511df05"), atom222("511df05"), atom222("511df05")));
			break;
		}
		case "catalog32": {
			parentAlg = q.find(32, 32, 32).orElseThrow();
			parentNode = new Lineage.Atom(
					"32x32x32@" + SchemeIO.contentHash(parentAlg).substring(0, 7));
			break;
		}
		default:
			throw new IllegalArgumentException("--parent=strassen|winograd|catalog32");
		}
		log.info("parent {}: rank={} axisMargins={}", which, parentAlg.r,
				Arrays.toString(ProjectionSearch.axisMargins(parentAlg)));

		if (!"catalog32".equals(which)) {
			writeParentToMarginBases(which, parentAlg, parentNode);
		}

		harvest(which, parentAlg, parentNode, q, minDim);
	}

	private static Lineage.Atom atom222(String hash) {
		return new Lineage.Atom("2x2x2@" + hash);
	}

	private static NonCubicBilinearAlgorithm tower(NonCubicBilinearAlgorithm base, int k) {
		NonCubicBilinearAlgorithm t = base;
		for (int i = 1; i < k; i++) {
			t = Compose.kroneckerGeneral(t, base);
		}
		return t;
	}

	record Node(boolean[] mn, boolean[] mm, boolean[] mp, long survivors) {}

	private static void harvest(String label, NonCubicBilinearAlgorithm parent,
			Lineage.Node parentNode, FieldAwareLookup q, int minDim) {
		long t0 = System.currentTimeMillis();
		ProjectionSearch.MaskedProjector proj = new ProjectionSearch.MaskedProjector(parent);
		int dim = parent.n;
		Map<Long, Node> lattice = new LinkedHashMap<>();
		boolean[] full = new boolean[dim];
		Arrays.fill(full, true);
		// THREE DISTINCT arrays — sharing one `full` array across the axes let a
		// single-axis bit-clear corrupt all three masks (the 2026-06-12 phantom-win
		// bug: every lattice node drifted below its key's shape, so survivor counts
		// of SMALLER shapes were compared against BIGGER shapes' catalog ranks).
		lattice.put(key(dim, dim, dim), new Node(full.clone(), full.clone(), full.clone(),
				parent.r));

		int wins = 0;
		int written = 0;
		int failed = 0;
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
					// Fail-loud lattice invariant (the phantom-win bug's guard):
					// a node's mask popcounts MUST match its key's shape.
					if (count(best.mn()) != n || count(best.mm()) != m
							|| count(best.mp()) != p) {
						throw new IllegalStateException(String.format(
								"lattice mask/key drift at ⟨%d,%d,%d⟩: masks %d/%d/%d",
								n, m, p, count(best.mn()), count(best.mm()), count(best.mp())));
					}
					lattice.put(key(n, m, p), best);
					int cat = q.findRank(n, m, p);
					if (cat < Recombination.SotaResolver.UNKNOWN_RANK && best.survivors() < cat) {
						wins++;
						if (persistWin(parent, parentNode, best, n, m, p, cat)) {
							written++;
						} else {
							failed++;
						}
					}
				}
			}
			if (drops % 10 == 0) {
				log.info("[progress] {} layer {}/{}: {} wins, {} written, {} failed, {}s",
						label, drops, 3 * (dim - minDim), wins, written, failed,
						(System.currentTimeMillis() - t0) / 1000);
			}
		}
		log.info("{} harvest done: {} wins, {} written, {} failed, {}s",
				label, wins, written, failed, (System.currentTimeMillis() - t0) / 1000);
	}

	private static boolean persistWin(NonCubicBilinearAlgorithm parent, Lineage.Node parentNode,
			Node node, int n, int m, int p, int catalogRank) {
		try {
			int[] keepN = toIdx(node.mn());
			int[] keepM = toIdx(node.mm());
			int[] keepP = toIdx(node.mp());
			NonCubicBilinearAlgorithm alg = Compose.project(parent, keepN, keepM, keepP);
			if (alg.r != node.survivors() || alg.n != n || alg.m != m || alg.p != p) {
				log.warn("⟨{},{},{}⟩: projection mismatch (got ⟨{},{},{}⟩ r={} vs survivors={},"
						+ " keeps {}/{}/{}) — skipped",
						n, m, p, alg.n, alg.m, alg.p, alg.r, node.survivors(),
						keepN.length, keepM.length, keepP.length);
				return false;
			}
			if (!Verifier.passesRandomMatmulSpotCheck(alg)) {
				log.warn("⟨{},{},{}⟩: spot-check FAILED — skipped", n, m, p);
				return false;
			}
			Lineage.Node lineage = new Lineage.Project(parentNode, keepN, keepM, keepP);
			int maxDim = Math.max(n, Math.max(m, p));
			File dir = new File("src/main/resources/schemes/derived/section" + maxDim);
			dir.mkdirs();
			boolean stub = maxDim > CatalogPolicy.MATERIALISE_MAX_DIM;
			int adds = stub ? 0 : Verifier.additionCount(alg);
			String fn = String.format("derived_cascade-%dx%dx%d_m%d%s.json",
					n, m, p, alg.r, stub ? "" : "_a" + adds);
			File out = new File(dir, fn);
			if (out.exists()) {
				return true;  // an equal-rank cascade entry already landed
			}
			if (stub) {
				SchemeIO.writeStub(alg, out, lineage);
			} else {
				SchemeIO.write(alg, out, lineage);
			}
			log.info("WIN ⟨{},{},{}⟩ {} < catalog {} → {}{}", n, m, p, alg.r, catalogRank,
					out.getName(), stub ? " (stub)" : "");
			return true;
		} catch (Exception e) {
			log.warn("⟨{},{},{}⟩: persist failed: {}", n, m, p, e.toString());
			return false;
		}
	}

	private static void writeParentToMarginBases(String which, NonCubicBilinearAlgorithm alg,
			Lineage.Node node) {
		try {
			File dir = new File("src/main/resources/schemes/margin-bases/section32");
			dir.mkdirs();
			File out = new File(dir, String.format("%s_tower5-32x32x32_m%d.json", which, alg.r));
			if (out.exists()) {
				return;
			}
			SchemeIO.writeStub(alg, out, node);
			int[] mu = ProjectionSearch.axisMargins(alg);
			SchemeIO.addFields(out, Map.of(
					"fields", List.of("F2", "F3", "Z", "Q", "R", "C"),
					"source", "kron_tower",
					"projection_margins", List.of(mu[0], mu[1], mu[2])), true);
			log.info("wrote margin-base parent {} (r={}, μ={})", out.getName(), alg.r,
					Arrays.toString(mu));
		} catch (Exception e) {
			log.warn("margin-base parent write failed: {}", e.toString());
		}
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

	private static int[] toIdx(boolean[] mask) {
		List<Integer> idx = new ArrayList<>();
		for (int i = 0; i < mask.length; i++) {
			if (mask[i]) {
				idx.add(i);
			}
		}
		return idx.stream().mapToInt(Integer::intValue).toArray();
	}

	private static long key(int n, int m, int p) {
		return ((long) n << 20) | ((long) m << 10) | p;
	}
}
