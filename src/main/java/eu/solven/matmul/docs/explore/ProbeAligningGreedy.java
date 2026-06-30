package eu.solven.matmul.docs.explore;

import eu.solven.matmul.recombination.Recombination;

import java.util.Random;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SerendipitousBudProduct;
import eu.solven.matmul.catalog.SerendipityCeiling;
import eu.solven.matmul.search.flip.FlipGraphWalk;
import eu.solven.matmul.search.flip.FlipObjectives;
import eu.solven.matmul.search.flip.FlipScheme;
import lombok.extern.slf4j.Slf4j;

/**
 * DIRECTED bud-structure driver (user 2026-06-11: "bud-structure might be
 * easier to drive by applying meta-flips in a specific way"): instead of a
 * random walk, an exhaustive ONE-STEP LOOKAHEAD greedy — every legal flip
 * (slot, pair, sign, variant) is scored by the true priced objective
 * (serendipitousCost, tie-broken by W-axis σ-potential), and the best
 * improving move is taken; random kicks + retry when stuck. Each flip rewrites
 * exactly one factor vector, so this is precisely "which single vector rewrite
 * best aligns a W column with an existing class" — the constructive move the
 * blind walk kept missing on ⟨2,4,3⟩.
 *
 * <p>Also prints the {@link SerendipityCeiling} bracket so the result is
 * positioned between the catalog value and the within-scope optimum.</p>
 */
@Slf4j
public class ProbeAligningGreedy {

	public static void main(String[] args) {
		FieldAwareLookup zLookup = new FieldAwareLookup(Field.Z);
		FieldAwareLookup qLookup = new FieldAwareLookup(Field.Q);
		int n2 = 3;
		int m2 = 2;
		int p2 = 3;
		long current = qLookup.findRank(6, 8, 9);

		NonCubicBilinearAlgorithm seedAlg = zLookup.find(2, 4, 3).orElseThrow();
		FlipScheme seed = FlipScheme.of(seedAlg);
		for (SerendipitousBudProduct.BudType t : SerendipitousBudProduct.BudType.values()) {
			SerendipityCeiling.AxisCeiling c = SerendipityCeiling.forAxis(
					t, 2, 4, 3, seed.rank(), qLookup, n2, m2, p2);
			log.info("ceiling[{}-buds]: maxSavings={} → ⟨6,8,9⟩ floor {} (≥{} classes;"
					+ " optimal-within-scope, likely unreachable)",
					t, c.maxSavings(), c.productRankFloor(), c.minClasses());
		}
		long seedCost = cost(seed, qLookup, n2, m2, p2);
		log.info("seed: rank={} predicted ⟨6,8,9⟩ ≤ {} | catalog {} ", seed.rank(),
				Math.ceilDiv(seedCost, 1_000_000), current);

		FlipScheme best = seed.copy();
		long bestCost = seedCost;
		Random rng = new Random(42);
		FlipScheme cur = seed.copy();
		long curCost = seedCost;
		for (int round = 0; round < 400; round++) {
			FlipScheme nextScheme = null;
			long nextCost = curCost;
			for (FlipScheme.Slot slot : FlipScheme.Slot.values()) {
				for (int[] cls : cur.signClasses(slot)) {
					for (int a = 0; a < cls.length; a++) {
						for (int b = 0; b < cls.length; b++) {
							if (a == b) {
								continue;
							}
							int sign = FlipScheme.signRatio(
									cur.vec(slot, cls[a]), cur.vec(slot, cls[b]));
							for (boolean variantB : new boolean[] { false, true }) {
								FlipScheme trial = cur.copy();
								if (!trial.flipWithinCap(slot, cls[a], cls[b], sign,
										variantB, 2)) {
									continue;
								}
								trial.dropZeroProducts();
								long c = cost(trial, qLookup, n2, m2, p2);
								if (c < nextCost) {
									nextCost = c;
									nextScheme = trial;
								}
							}
						}
					}
				}
			}
			if (nextScheme != null) {
				cur = nextScheme;
				curCost = nextCost;
				if (curCost < bestCost) {
					bestCost = curCost;
					best = cur.copy();
					log.info("round {}: improved — predicted ⟨6,8,9⟩ ≤ {} (potential tie {})",
							round, Math.ceilDiv(bestCost, 1_000_000), bestCost % 1_000_000);
				}
				continue;
			}
			// Stuck at a local optimum of the 1-step landscape: random kick.
			for (int k = 0; k < 3 + rng.nextInt(4); k++) {
				FlipGraphWalk.applyRandomFlip(cur, rng, 2);
			}
			cur.dropZeroProducts();
			curCost = cost(cur, qLookup, n2, m2, p2);
		}
		long predicted = Math.ceilDiv(bestCost, 1_000_000);  // cost = direct*1e6 - potential
		log.info("DIRECTED best: rank={} predicted ⟨6,8,9⟩ ≤ {} (seed {}, catalog {},"
				+ " W-ceiling floor {})",
				best.rank(), predicted, Math.ceilDiv(seedCost, 1_000_000), current,
				SerendipityCeiling.forAxis(SerendipitousBudProduct.BudType.W, 2, 4, 3,
						best.rank(), qLookup, n2, m2, p2).productRankFloor());
		if (predicted < current) {
			NonCubicBilinearAlgorithm baseAlg = best.toAlgorithm();
			SerendipitousBudProduct.BudDecomposition bestDec = null;
			long bestDecCost = Long.MAX_VALUE / 4;
			for (SerendipitousBudProduct.BudType[] order : SerendipitousBudProduct.ALL_ORDERINGS) {
				SerendipitousBudProduct.BudDecomposition dec =
						SerendipitousBudProduct.findBuds(baseAlg, order);
				long c = SerendipitousBudProduct.costOf(dec, qLookup, n2, m2, p2);
				if (c < bestDecCost) {
					bestDecCost = c;
					bestDec = dec;
				}
			}
			NonCubicBilinearAlgorithm product = SerendipitousBudProduct.productFromDecomposition(
					baseAlg, bestDec, qLookup, n2, m2, p2,
					java.util.EnumSet.allOf(SerendipitousBudProduct.BudType.class));
			log.info("BUILT ⟨{},{},{}⟩ rank={} exact={} — vs catalog {}",
					product.n, product.m, product.p, product.r,
					eu.solven.matmul.verifiers.Verifier.isExactNonCubic(product), current);
		}
	}

	/** Priced cost: serendipitousCost (truth) ·1e6 − W-axis σ-potential (plateau tie). */
	private static long cost(FlipScheme s, FieldAwareLookup q, int n2, int m2, int p2) {
		long direct = SerendipitousBudProduct.serendipitousCost(s.toAlgorithm(), q, n2, m2, p2);
		long inner = q.findRank(n2, m2, p2);
		long wPotential = 0;
		int[] wSizes = SerendipitousBudProduct.independentClassSizes(s.toAlgorithm())[2];
		for (int size : wSizes) {
			if (size >= 2) {
				long enlarged = q.findRank(n2, size * m2, p2);
				if (enlarged < Recombination.SotaResolver.UNKNOWN_RANK) {
					wPotential += Math.max(0, size * inner - enlarged);
				}
			}
		}
		return direct * 1_000_000 - wPotential;
	}
}
