package eu.solven.matmul.docs.explore;

import eu.solven.matmul.recombination.Recombination;

import java.util.ArrayList;
import java.util.List;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.Compose;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SerendipitousBudProduct;
import eu.solven.matmul.catalog.SerendipitousBudProduct.BudDecomposition;
import lombok.extern.slf4j.Slf4j;

/**
 * TWO-SIDED serendipity scan (user 2026-06-12: "can we kronecker 2 bases, each
 * with a bud-structure, in a way both bud-structures can be leveraged?"). In a
 * PLAIN Kronecker product S1⊗S2, SAME-axis bud classes MULTIPLY (k₁-class ×
 * k₂-class → k₁k₂-class); using the bud-rich plain Kron — deliberately
 * rank-suboptimal — as the OUTER base against a third factor S3 can fuse
 * enlarged blocks no hierarchical route sees. Validated: plain Kron ⟨2,4,3⟩²
 * (rank 400, 16 W-quads) ⊗ ⟨3,2,3⟩ = 5824 exact &lt; 5884 = same grouping
 * through the rank-398 catalog ⟨4,16,9⟩ scheme.
 *
 * <p><b>Two-phase pricing.</b> Phase 1 prunes with the ANALYTIC savings bound:
 * sum σ over all axes' multiplied classes as if every product fused on every
 * axis simultaneously. That is physically impossible (fusion consumes the
 * product — a product in an S1 U-pair AND an S2 V-pair can serve only one),
 * but it is a SOUND upper bound: σ is superadditive (R(enl a+b) ≤ R(enl a) +
 * R(enl b)), so any realizable bud decomposition saves no more. Triples whose
 * bound can't beat the catalog are pruned exactly. Phase 2 EXACT-confirms the
 * survivors ({@code findBuds} over all orderings + {@code costOf} on the
 * materialized Kron). Only same-axis structures compound; different-axis
 * structures collide per-product.</p>
 *
 * <p>Any WIN logged in phase 2 is a true constructive upper bound (plain Kron
 * + Smith fusion, both verified mechanisms) — materialize via
 * {@code productViaBuds}, verify, and register.</p>
 *
 * <p>Args: {@code [maxDim=32] [maxFactorVolume=120] [maxKronRank=2000]}.</p>
 */
@Slf4j
public class ScanTwoSidedKron {

	record Shape(int n, int m, int p) {}

	public static void main(String[] args) {
		int maxDim = args.length > 0 ? Integer.parseInt(args[0]) : 32;
		int maxFactorVolume = args.length > 1 ? Integer.parseInt(args[1]) : 120;
		int maxKronRank = args.length > 2 ? Integer.parseInt(args[2]) : 2_000;
		FieldAwareLookup q = new FieldAwareLookup(Field.Q);

		List<Shape> factors = new ArrayList<>();
		for (int n = 1; n <= 16; n++) {
			for (int m = 1; m <= 16; m++) {
				for (int p = 1; p <= 16; p++) {
					if (n * m * p < 2 || n * m * p > maxFactorVolume) {
						continue;
					}
					if (q.findRank(n, m, p) < Recombination.SotaResolver.UNKNOWN_RANK) {
						factors.add(new Shape(n, m, p));
					}
				}
			}
		}
		log.info("two-sided Kron scan (exact greedy pricing): {} candidate factors, maxDim={}",
				factors.size(), maxDim);

		long t0 = System.currentTimeMillis();
		int triples = 0;
		int candidates = 0;
		int wins = 0;
		java.util.Map<Shape, int[][]> profiles = new java.util.HashMap<>();
		for (int i1 = 0; i1 < factors.size(); i1++) {
			Shape s1 = factors.get(i1);
			int[][] prof1 = profiles.computeIfAbsent(s1, s -> loadProfile(q, s));
			if (prof1 == null) {
				continue;
			}
			// Kron(S1,S2) ≡ Kron(S2,S1) for this pricing (class products and the
			// target are symmetric) — ordered pairs halve the sweep losslessly.
			for (Shape s2 : factors.subList(i1, factors.size())) {
				if (s1.n() * s2.n() > maxDim || s1.m() * s2.m() > maxDim
						|| s1.p() * s2.p() > maxDim) {
					continue;
				}
				int[][] prof2 = profiles.computeIfAbsent(s2, s -> loadProfile(q, s));
				if (prof2 == null
						|| (long) q.findRank(s1.n(), s1.m(), s1.p())
								* q.findRank(s2.n(), s2.m(), s2.p()) > maxKronRank) {
					continue;
				}
				long r1 = q.findRank(s1.n(), s1.m(), s1.p());
				long r2 = q.findRank(s2.n(), s2.m(), s2.p());
				NonCubicBilinearAlgorithm kron = null;
				List<BudDecomposition> decs = null;
				for (Shape s3 : factors) {
					int bn = s1.n() * s2.n() * s3.n();
					int bm = s1.m() * s2.m() * s3.m();
					int bp = s1.p() * s2.p() * s3.p();
					if (bn > maxDim || bm > maxDim || bp > maxDim) {
						continue;
					}
					int best = q.findRank(bn, bm, bp);
					if (best >= Recombination.SotaResolver.UNKNOWN_RANK) {
						continue;
					}
					triples++;
					long r3 = q.findRank(s3.n(), s3.m(), s3.p());
					// Phase 1: analytic upper bound on savings (over-counts
					// cross-axis overlaps; sound by σ superadditivity).
					long bound = 0;
					for (int ax = 0; ax < 3; ax++) {
						int singles1 = (int) r1 - sum(prof1[ax]);
						int singles2 = (int) r2 - sum(prof2[ax]);
						for (int k1 : prof1[ax]) {
							for (int k2 : prof2[ax]) {
								bound += sigma(q, s3, ax, k1 * k2, r3);
							}
							bound += (long) singles2 * sigma(q, s3, ax, k1, r3);
						}
						for (int k2 : prof2[ax]) {
							bound += (long) singles1 * sigma(q, s3, ax, k2, r3);
						}
					}
					if (r1 * r2 * r3 - bound >= best) {
						continue;
					}
					candidates++;
					// Phase 2: exact greedy on the materialized Kron.
					if (decs == null) {
						NonCubicBilinearAlgorithm a1 = q.find(s1.n(), s1.m(), s1.p()).orElseThrow();
						NonCubicBilinearAlgorithm a2 = q.find(s2.n(), s2.m(), s2.p()).orElseThrow();
						kron = Compose.kroneckerGeneral(a1, a2);
						decs = new ArrayList<>();
						for (SerendipitousBudProduct.BudType[] order
								: SerendipitousBudProduct.ALL_ORDERINGS) {
							decs.add(SerendipitousBudProduct.findBuds(kron, order));
						}
					}
					long price = Long.MAX_VALUE;
					for (BudDecomposition dec : decs) {
						price = Math.min(price, SerendipitousBudProduct.costOf(dec, q,
								s3.n(), s3.m(), s3.p()));
					}
					if (price < best) {
						wins++;
						log.info("WIN ⟨{},{},{}⟩: Kron[{}⊗{}]⊗{} = {} < catalog {} — exact greedy; "
								+ "materialize via productViaBuds + verify before entry",
								bn, bm, bp, s1, s2, s3, price, best);
					} else {
						log.info("refuted ⟨{},{},{}⟩: Kron[{}⊗{}]⊗{} exact {} ≥ catalog {} "
								+ "(analytic bound predicted {})",
								bn, bm, bp, s1, s2, s3, price, best, r1 * r2 * r3 - bound);
					}
				}
			}
		}
		log.info("two-sided Kron scan done: {} triples bounded, {} exact-confirmed candidates, "
				+ "{} wins, {}s", triples, candidates, wins, (System.currentTimeMillis() - t0) / 1000);
	}

	/** Per-axis independent class sizes ≥2 of the catalog-best scheme. */
	private static int[][] loadProfile(FieldAwareLookup q, Shape s) {
		NonCubicBilinearAlgorithm alg = q.find(s.n(), s.m(), s.p()).orElse(null);
		if (alg == null) {
			return null;
		}
		int[][] raw = SerendipitousBudProduct.independentClassSizes(alg);
		int[][] out = new int[3][];
		for (int ax = 0; ax < 3; ax++) {
			out[ax] = java.util.Arrays.stream(raw[ax]).filter(k -> k >= 2).toArray();
		}
		return out;
	}

	/** σ_ax(k) for inner s3 — 0 when the enlarged block is missing or k==1. */
	private static long sigma(FieldAwareLookup q, Shape s3, int ax, int k, long r3) {
		if (k <= 1) {
			return 0;
		}
		int enlarged = switch (ax) {
			case 0 -> q.findRank(s3.n(), s3.m(), k * s3.p());
			case 1 -> q.findRank(k * s3.n(), s3.m(), s3.p());
			default -> q.findRank(s3.n(), k * s3.m(), s3.p());
		};
		return enlarged >= Recombination.SotaResolver.UNKNOWN_RANK ? 0 : Math.max(0, k * r3 - enlarged);
	}

	private static int sum(int[] xs) {
		int s = 0;
		for (int x : xs) {
			s += x;
		}
		return s;
	}
}
