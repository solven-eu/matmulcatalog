package eu.solven.matmul.docs.explore;

import java.io.File;
import java.util.List;
import java.util.TreeMap;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.catalog.SerendipitousBudProduct;
import eu.solven.matmul.catalog.SerendipitousBudProduct.Bud;
import eu.solven.matmul.catalog.SerendipitousBudProduct.BudDecomposition;
import eu.solven.matmul.catalog.SerendipitousBudProduct.GridBud;

import lombok.extern.slf4j.Slf4j;

/**
 * Print the bud structure the serendipitous engine picks for one base against one
 * inner — in the notation of Perminov's 17–32 table
 * ({@code references/perminov-serendipitous-17-32.json}: {@code "<1,2,2> + 2<2,1,1> +
 * 15<1,1,1>"}), so a recipe we fail to reach can be compared term class by term
 * class with the published one.
 *
 * <p>Args: {@code <base.json> n2 m2 p2 [--field=Q]}. Prints every member of
 * {@link SerendipitousBudProduct#candidateDecompositions} with its priced cost, best
 * first. The cost is a <b>bound</b> (cheapest of a structural family, not an optimum
 * over all bud partitions).</p>
 */
@Slf4j
public class ProbeSerendipStructure {
	public static void main(String[] args) throws Exception {
		File baseFile = new File(args[0]);
		int n2 = Integer.parseInt(args[1]), m2 = Integer.parseInt(args[2]), p2 = Integer.parseInt(args[3]);
		String field = "Q";
		for (String a : args) if (a.startsWith("--field=")) field = a.substring("--field=".length());
		FieldAwareLookup lookup = new FieldAwareLookup(Field.fromTag(field));
		NonCubicBilinearAlgorithm base = SchemeIO.readBilinear(baseFile);

		int[][] sizes = SerendipitousBudProduct.independentClassSizes(base);
		log.info("base ⟨{},{},{}⟩ r={} — class sizes ≥2: U {} V {} W {}", base.n, base.m, base.p, base.r,
				histogram(sizes[0]), histogram(sizes[1]), histogram(sizes[2]));
		long inner = lookup.findRank(n2, m2, p2);
		log.info("inner ⟨{},{},{}⟩ = {} → plain Kronecker {}", n2, m2, p2, inner, base.r * inner);

		List<BudDecomposition> family = SerendipitousBudProduct.candidateDecompositions(base);
		TreeMap<Long, String> byCost = new TreeMap<>();
		for (BudDecomposition dec : family) {
			long cost = SerendipitousBudProduct.costOf(dec, lookup, n2, m2, p2);
			byCost.putIfAbsent(cost, structure(dec));
		}
		log.info("{} candidate decompositions, {} distinct costs:", family.size(), byCost.size());
		byCost.forEach((cost, s) -> log.info("  {}  {}", cost, s));
		log.info("⟨{},{},{}⟩ ≤ {} (bound: cheapest of the structural family)", base.n * n2, base.m * m2, base.p * p2,
				byCost.firstKey());
	}

	/** {@code "<1,2,2> + 2<2,1,1> + 15<1,1,1>"}. */
	static String structure(BudDecomposition dec) {
		TreeMap<String, Integer> groups = new TreeMap<>();
		for (Bud b : dec.buds()) {
			int k = b.terms().length;
			groups.merge(switch (b.type()) {
				case U -> "<1,1," + k + ">";
				case V -> "<" + k + ",1,1>";
				case W -> "<1," + k + ",1>";
			}, 1, Integer::sum);
		}
		for (GridBud g : dec.grids()) groups.merge("<" + g.a() + "," + g.b() + "," + g.c() + ">", 1, Integer::sum);
		if (dec.trivial().length > 0) groups.merge("<1,1,1>", dec.trivial().length, Integer::sum);
		StringBuilder sb = new StringBuilder();
		groups.forEach((shape, count) -> sb.append(sb.length() > 0 ? " + " : "").append(count > 1 ? count : "").append(shape));
		return sb.toString();
	}

	private static String histogram(int[] classSizes) {
		TreeMap<Integer, Integer> h = new TreeMap<>();
		for (int s : classSizes) if (s >= 2) h.merge(s, 1, Integer::sum);
		return h.toString();
	}
}
