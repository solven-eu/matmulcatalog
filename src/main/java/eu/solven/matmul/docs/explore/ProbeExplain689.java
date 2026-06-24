package eu.solven.matmul.docs.explore;

import eu.solven.matmul.recombination.Recombination;

import java.util.Arrays;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SerendipitousBudProduct;
import lombok.extern.slf4j.Slf4j;

/**
 * Pedagogical dump for the ⟨6,8,9⟩=296 walkthrough (user 2026-06-11): the
 * base's identity and bud structure term-by-term, the product assembly
 * arithmetic, and what the NON-serendipitous mechanisms reach (plain Kronecker
 * over every factorization, single-axis concatenation).
 */
@Slf4j
public class ProbeExplain689 {

	public static void main(String[] args) {
		FieldAwareLookup z = new FieldAwareLookup(Field.Z);
		FieldAwareLookup q = new FieldAwareLookup(Field.Q);

		// 1. The base.
		var ws = z.findWithSource(2, 4, 3).orElseThrow();
		NonCubicBilinearAlgorithm base = ws.alg();
		log.info("base ⟨2,4,3⟩: rank={} file={}", base.r,
				ws.path() == null ? "?" : ws.path().getFileName());

		// 2. Bud structure: greedy decompositions (all orderings) + W classes.
		for (SerendipitousBudProduct.BudType[] order : SerendipitousBudProduct.ALL_ORDERINGS) {
			SerendipitousBudProduct.BudDecomposition dec =
					SerendipitousBudProduct.findBuds(base, order);
			long c = SerendipitousBudProduct.costOf(dec, q, 3, 2, 3);
			log.info("ordering {}: {} buds + {} trivial → cost {}",
					Arrays.toString(order), dec.buds().size(), dec.trivial().length, c);
		}
		SerendipitousBudProduct.BudDecomposition dec = SerendipitousBudProduct.findBuds(base);
		double[][] w = base.denseW();
		for (SerendipitousBudProduct.Bud bud : dec.buds()) {
			StringBuilder sb = new StringBuilder();
			for (int t : bud.terms()) {
				double[] col = new double[w.length];
				for (int i = 0; i < w.length; i++) {
					col[i] = w[i][t];
				}
				sb.append(" term").append(t).append("=w").append(Arrays.toString(col));
			}
			log.info("{}-bud k={}:{}", bud.type(), bud.terms().length, sb);
		}
		int[][] classSizes = SerendipitousBudProduct.independentClassSizes(base);
		log.info("independent class sizes U={} V={} W={}", Arrays.toString(classSizes[0]),
				Arrays.toString(classSizes[1]), Arrays.toString(classSizes[2]));

		// 3. Assembly arithmetic.
		log.info("assembly: {}·R⟨3,2,3⟩={} + {}·R⟨3,4,3⟩={} → total {}",
				dec.trivial().length, dec.trivial().length * q.findRank(3, 2, 3),
				dec.buds().size(), dec.buds().size() * q.findRank(3, 4, 3),
				dec.trivial().length * q.findRank(3, 2, 3)
						+ dec.buds().size() * q.findRank(3, 4, 3));

		// 4. Alternatives WITHOUT serendipity.
		long bestKron = Long.MAX_VALUE;
		String bestKronLabel = "";
		for (int a = 1; a <= 6; a++) {
			for (int b = 1; b <= 8; b++) {
				for (int c = 1; c <= 9; c++) {
					if (6 % a != 0 || 8 % b != 0 || 9 % c != 0) {
						continue;
					}
					if (a * b * c == 1 || (a == 6 && b == 8 && c == 9)) {
						continue;  // identity factorizations are circular
					}
					long r1 = q.findRank(a, b, c);
					long r2 = q.findRank(6 / a, 8 / b, 9 / c);
					if (r1 >= Recombination.SotaResolver.UNKNOWN_RANK || r2 >= Recombination.SotaResolver.UNKNOWN_RANK) {
						continue;
					}
					if (r1 * r2 < bestKron) {
						bestKron = r1 * r2;
						bestKronLabel = "⟨" + a + "," + b + "," + c + "⟩(" + r1 + ") ⊗ ⟨"
								+ (6 / a) + "," + (8 / b) + "," + (9 / c) + "⟩(" + r2 + ")";
					}
				}
			}
		}
		log.info("best PLAIN Kronecker: {} = {}", bestKronLabel, bestKron);
		long bestConcat = Long.MAX_VALUE;
		String bestConcatLabel = "";
		for (int k = 1; k < 9; k++) {
			long s = q.findRank(6, 8, k) + q.findRank(6, 8, 9 - k);
			if (s < bestConcat) {
				bestConcat = s;
				bestConcatLabel = "⟨6,8," + k + "⟩+⟨6,8," + (9 - k) + "⟩";
			}
		}
		for (int k = 1; k < 6; k++) {
			long s = q.findRank(k, 8, 9) + q.findRank(6 - k, 8, 9);
			if (s < bestConcat) {
				bestConcat = s;
				bestConcatLabel = "⟨" + k + ",8,9⟩+⟨" + (6 - k) + ",8,9⟩";
			}
		}
		for (int k = 1; k < 8; k++) {
			long s = q.findRank(6, k, 9) + q.findRank(6, 8 - k, 9);
			if (s < bestConcat) {
				bestConcat = s;
				bestConcatLabel = "⟨6," + k + ",9⟩+⟨6," + (8 - k) + ",9⟩";
			}
		}
		log.info("best CONCAT split: {} = {}", bestConcatLabel, bestConcat);
		log.info("catalog ⟨6,8,9⟩ = {} (serendipitous)", q.findRank(6, 8, 9));
	}
}
