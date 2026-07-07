package eu.solven.matmul.docs.explore;

import java.io.File;
import java.util.Arrays;
import java.util.Map;
import java.util.TreeMap;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.SymmetryTransforms;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.catalog.SerendipitousBudProduct;
import eu.solven.matmul.catalog.SerendipitousBudProduct.Bud;
import eu.solven.matmul.catalog.SerendipitousBudProduct.BudDecomposition;
import lombok.extern.slf4j.Slf4j;

/**
 * Throwaway probe (fmm-react 2026-07-06, ⟨20,24,25⟩ gap): compare the two
 * rank-204 ⟨5,5,12⟩ bases — our catalog's {@code perminov_ZT-61a6cb7} vs the
 * FMM-Lille import ({@code bud-bases/section12/fmm-lille_5x5x12_r204_a2326}) —
 * as serendipitous bases for {@code ⟨5,12,5⟩ ⊗ ⟨4,2,5⟩ → ⟨20,24,25⟩}.
 *
 * <p>FMM's page states {@code ⟨20,24,25⟩:6466 = (⟨5,12,5⟩:204 − 42) ⊗ ⟨4,2,5⟩:32
 * + ⟨16,2,5⟩:126 + 2·⟨4,6,5⟩:90 + 16·⟨4,4,5⟩:61}, i.e. a bud partition
 * 16×(size-2) + 2×(size-3) + 1×(size-4) saving 62 vs plain Kron 6528. Our
 * catalog scheme only yields 18×(size-2) (saving 54 → 6474). This probe prints,
 * per base and per orientation reaching ⟨5,12,5⟩, the bud class-size histogram
 * under every ordering and the priced serendipitous cost.</p>
 */
@Slf4j
public final class ProbeFmmSerendip202425 {
	private ProbeFmmSerendip202425() {}

	public static void main(String[] args) throws Exception {
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		Map<String, String> bases = Map.of(
				"ours-61a6cb7", "src/main/resources/schemes/known/section12/5x5x12-r204-perminov_ZT-61a6cb7.json",
				"fmm-a2326", "src/main/resources/schemes/bud-bases/section12/fmm-lille_5x5x12_r204_a2326.json");
		String[] perms = { "ABC->ABC", "ABC->ACB", "ABC->BAC", "ABC->BCA", "ABC->CAB", "ABC->CBA" };
		for (Map.Entry<String, String> e : bases.entrySet()) {
			NonCubicBilinearAlgorithm raw = SchemeIO.read(new File(e.getValue()));
			for (String perm : perms) {
				NonCubicBilinearAlgorithm alg = SymmetryTransforms.permuteAxes(raw, perm);
				if (alg == null || alg.n != 5 || alg.m != 12 || alg.p != 5) continue;
				long cost = SerendipitousBudProduct.serendipitousCost(alg, lookup, 4, 2, 5);
				log.info("{} {} → ⟨5,12,5⟩ : serendipitousCost ⊗⟨4,2,5⟩ = {} (plain Kron 6528)",
						e.getKey(), perm, cost);
				int[][] classes = SerendipitousBudProduct.independentClassSizes(alg);
				String[] ax = { "U", "V", "W" };
				for (int i = 0; i < 3; i++) {
					Map<Integer, Integer> h = new TreeMap<>();
					for (int k : classes[i]) if (k >= 2) h.merge(k, 1, Integer::sum);
					log.info("   direction classes {} (size≥2): {}", ax[i], h);
				}
				for (var order : SerendipitousBudProduct.ALL_ORDERINGS) {
					BudDecomposition dec = SerendipitousBudProduct.findBuds(alg, order);
					Map<String, Integer> hist = new TreeMap<>();
					for (Bud b : dec.buds()) {
						hist.merge(b.type() + "x" + b.terms().length, 1, Integer::sum);
					}
					log.info("   order {} : {} (+{} trivial) → costOf={}",
							Arrays.toString(order), hist, dec.trivial().length,
							SerendipitousBudProduct.costOf(dec, lookup, 4, 2, 5));
				}
			}
		}
	}
}
