package eu.solven.matmul.docs.explore;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.search.LineageReplayer;
import lombok.extern.slf4j.Slf4j;

/**
 * Throwaway probe (fmm-gap 2026-07-09, the ±1 family): does OUR replayed
 * recombination contain two multiplications that coincide as bilinear forms
 * over the TARGET's variables (identical U and V columns up to a common
 * scale)? Leaves are built independently, but their feeding forms overlap —
 * a coinciding pair merges into one product (W-columns add), rank −1 free.
 */
@Slf4j
public final class ProbeDuplicateProducts {
	private ProbeDuplicateProducts() {}

	public static void main(String[] args) {
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		LineageReplayer replayer = LineageReplayer.withDefaultPool(lookup);
		String[] files = args.length > 0 ? args : new String[] {
				"src/main/resources/schemes/derived/section19/17x17x19-r3267-derived-e21e78f.json" };
		for (String f : files) {
			NonCubicBilinearAlgorithm alg = replayer.replayFromFile(new File(f));
			log.info("{} → r={}", f.substring(f.lastIndexOf('/') + 1), alg.r);
			Map<String, Integer> seen = new HashMap<>();
			int dups = 0;
			for (int k = 0; k < alg.r; k++) {
				String key = directionKey(alg.u(), alg.dimU(), k) + "|" + directionKey(alg.v(), alg.dimV(), k);
				Integer prev = seen.putIfAbsent(key, k);
				if (prev != null) {
					dups++;
					log.info("  DUPLICATE u⊗v: products {} and {} — mergeable, rank −1", prev, k);
				}
			}
			log.info("  {} coinciding u⊗v pairs", dups);
		}
	}

	/** Direction-normalised sparse key of column k (first nonzero scaled to 1). */
	private static String directionKey(eu.solven.matmul.FactorMatrix m, int rows, int k) {
		StringBuilder sb = new StringBuilder();
		double[] col = new double[rows];
		m.forEachInColumn(k, (row, v) -> col[row] = v);
		double scale = 0;
		for (double v : col) if (v != 0) { scale = v; break; }
		if (scale == 0) return "zero";
		for (int i = 0; i < rows; i++) {
			if (col[i] != 0) sb.append(i).append(':').append((float) (col[i] / scale)).append(',');
		}
		return sb.toString();
	}
}
