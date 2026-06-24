package eu.solven.matmul.catalog;

import eu.solven.matmul.recombination.Recombination;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.util.Map;
import java.util.function.Function;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Lineage-level bud inference (the detect-tier promotion of buds): propagating
 * atom buds up the lineage must agree with expanding the composite and reading
 * its buds directly — for the cancellation-free / relabelling ops — and must
 * honestly report a weaker certainty otherwise.
 */
public class TestLineageBudInference {

	private static final String STRASSEN = "src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18_b0.json";
	private static final String LADERMAN = "src/main/resources/schemes/known/section3/laderman_1976-3x3x3_m23_a98_b0.json";

	private static NonCubicBilinearAlgorithm read(String path) throws Exception {
		return SchemeIO.read(eu.solven.matmul.catalog.SchemeResolver.byHint(path));
	}

	/** Atom resolver keyed by a label we control in each test. */
	private static Function<String, NonCubicBilinearAlgorithm> resolver(
			Map<String, NonCubicBilinearAlgorithm> byRef) {
		return byRef::get;
	}

	@Test
	public void kron_buds_inferred_match_expanded_exactly() throws Exception {
		NonCubicBilinearAlgorithm strassen = read(STRASSEN);   // ⟨2,2,2⟩=7
		NonCubicBilinearAlgorithm laderman = read(LADERMAN);   // ⟨3,3,3⟩=23
		NonCubicBilinearAlgorithm actualKron = Compose.kroneckerGeneral(strassen, laderman); // ⟨6,6,6⟩=161

		Lineage.Node lineage = new Lineage.KronProduct(
				new Lineage.Atom("strassen"), new Lineage.Atom("laderman"));
		LineageBudInference.Profile inferred = LineageBudInference.infer(lineage,
				resolver(Map.of("strassen", strassen, "laderman", laderman)));

		LineageBudInference.Profile actual = LineageBudInference.fromExpanded(actualKron);

		assertThat(inferred.certainty()).isEqualTo(LineageBudInference.Certainty.EXACT);
		assertThat(inferred.rank()).isEqualTo(161).isEqualTo(actual.rank());
		// Inferred (Cartesian of atom class sizes, no expansion) == expanded ground truth.
		assertThat(inferred.uClasses()).containsExactly(actual.uClasses());
		assertThat(inferred.vClasses()).containsExactly(actual.vClasses());
		assertThat(inferred.wClasses()).containsExactly(actual.wClasses());
		assertThat(inferred.summary()).isEqualTo(actual.summary());
	}

	@Test
	public void axisflip_is_passthrough_exact() throws Exception {
		NonCubicBilinearAlgorithm laderman = read(LADERMAN);
		LineageBudInference.Profile base = LineageBudInference.fromExpanded(laderman);

		Lineage.Node flipped = new Lineage.AxisFlip(new Lineage.Atom("l"), 5);
		LineageBudInference.Profile inferred =
				LineageBudInference.infer(flipped, resolver(Map.of("l", laderman)));

		assertThat(inferred.certainty()).isEqualTo(LineageBudInference.Certainty.EXACT);
		assertThat(inferred.uClasses()).containsExactly(base.uClasses());
		assertThat(inferred.vClasses()).containsExactly(base.vClasses());
		assertThat(inferred.wClasses()).containsExactly(base.wClasses());
	}

	@Test
	public void transpose_relabels_partitions() throws Exception {
		NonCubicBilinearAlgorithm laderman = read(LADERMAN);
		LineageBudInference.Profile base = LineageBudInference.fromExpanded(laderman);

		// "ABC->BCA": new A=old B, new B=old C, new C=old A.
		Lineage.Node t = new Lineage.Transpose(new Lineage.Atom("l"), "ABC->BCA");
		LineageBudInference.Profile inferred =
				LineageBudInference.infer(t, resolver(Map.of("l", laderman)));

		assertThat(inferred.certainty()).isEqualTo(LineageBudInference.Certainty.EXACT);
		assertThat(inferred.uClasses()).containsExactly(base.vClasses());
		assertThat(inferred.vClasses()).containsExactly(base.wClasses());
		assertThat(inferred.wClasses()).containsExactly(base.uClasses());
	}

	@Test
	public void concat_is_structural_estimate() throws Exception {
		NonCubicBilinearAlgorithm laderman = read(LADERMAN);
		Lineage.Node concat = new Lineage.ConcatCols(
				new Lineage.Atom("l"), new Lineage.Atom("l"));
		LineageBudInference.Profile inferred =
				LineageBudInference.infer(concat, resolver(Map.of("l", laderman)));

		assertThat(inferred.certainty())
				.isEqualTo(LineageBudInference.Certainty.STRUCTURAL_ESTIMATE);
		assertThat(inferred.rank()).isEqualTo(46);  // 23 + 23
	}

	@Test
	public void recombination_is_unknown_without_leaf_resolver() throws Exception {
		// The 2-arg infer has no leaf-by-shape resolver → recombination UNKNOWN.
		NonCubicBilinearAlgorithm laderman = read(LADERMAN);
		Lineage.Node recomb = new Lineage.RecombinationN(
				new Lineage.Atom("l"), new int[] { 1 }, new int[] { 1 }, new int[] { 1 },
				java.util.List.of(new Lineage.Atom("l")));
		LineageBudInference.Profile inferred =
				LineageBudInference.infer(recomb, resolver(Map.of("l", laderman)));

		assertThat(inferred.known()).isFalse();
		assertThat(inferred.certainty()).isEqualTo(LineageBudInference.Certainty.UNKNOWN);
	}

	/** Naive ⟨2,2,2⟩=8 base (4 U-buds of size 2) recombined with Strassen ⟨2,2,2⟩=7
	 *  leaves at a uniform [2,2]³ allocation → ⟨4,4,4⟩=56. Inferring buds from the
	 *  base + allocation (re-deriving sub-shapes, lifting by leaf classes) must match
	 *  expanding the actually-constructed scheme — the recombination indexing proof. */
	@Test
	public void recombination_estimate_matches_constructed_scheme() throws Exception {
		NonCubicBilinearAlgorithm naiveBase = naive222();
		NonCubicBilinearAlgorithm strassen = read(STRASSEN);  // ⟨2,2,2⟩=7
		int[] a = { 2, 2 }, b = { 2, 2 }, c = { 2, 2 };

		Recombination.AlgorithmLookup lookup =
				(n, m, p) -> java.util.Optional.of(strassen);
		Recombination.SotaResolver sota = (n, m, p) -> 7;
		NonCubicBilinearAlgorithm constructed =
				Recombination.constructWithAllocation(naiveBase, lookup, sota, a, b, c);
		assertThat(constructed.r).isEqualTo(56);  // 8 base terms × 7

		LineageBudInference.Profile inferred = LineageBudInference.inferRecombination(
				naiveBase, a, b, c, sz -> LineageBudInference.fromExpanded(strassen));
		LineageBudInference.Profile actual = LineageBudInference.fromExpanded(constructed);

		assertThat(inferred.certainty())
				.isEqualTo(LineageBudInference.Certainty.STRUCTURAL_ESTIMATE);
		assertThat(inferred.rank()).isEqualTo(56);
		// The base's 4 U-buds (size 2) each lift Strassen's 7 U-singletons →
		// 28 result U-classes of size 2. The inference reproduces the expanded truth.
		assertThat(inferred.uClasses()).containsExactly(actual.uClasses());
		assertThat(inferred.vClasses()).containsExactly(actual.vClasses());
		assertThat(inferred.wClasses()).containsExactly(actual.wClasses());
	}

	private static NonCubicBilinearAlgorithm naive222() {
		int r = 8;
		double[][] U = new double[4][r], V = new double[4][r], W = new double[4][r];
		int k = 0;
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				for (int l = 0; l < 2; l++) {
					U[i * 2 + j][k] = 1;
					V[j * 2 + l][k] = 1;
					W[i * 2 + l][k] = 1;
					k++;
				}
		return new NonCubicBilinearAlgorithm(2, 2, 2, U, V, W);
	}

	@Test
	public void unknown_leaf_poisons_kron_to_unknown() throws Exception {
		NonCubicBilinearAlgorithm laderman = read(LADERMAN);
		// Right operand is itself a Recombination → UNKNOWN; the Kron must inherit it.
		Lineage.Node lineage = new Lineage.KronProduct(
				new Lineage.Atom("l"),
				new Lineage.RecombinationN(new Lineage.Atom("l"),
						new int[] { 1 }, new int[] { 1 }, new int[] { 1 },
						java.util.List.of(new Lineage.Atom("l"))));
		LineageBudInference.Profile inferred =
				LineageBudInference.infer(lineage, resolver(Map.of("l", laderman)));

		assertThat(inferred.known()).isFalse();
	}
}
