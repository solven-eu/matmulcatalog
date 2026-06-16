package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.Arrays;
import java.util.List;
import java.util.Optional;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.algebra.Algebra;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.KnownAlgorithm;
import eu.solven.matmul.catalog.KnownAlgorithmCatalog;
import eu.solven.matmul.search.RecursiveComposition;
import eu.solven.matmul.search.RecursiveComposition.Factor;
import eu.solven.matmul.search.RecursiveComposition.Result;

/**
 * Validates the catalog's per-format best-known entries and the recursive
 * composition formulas highlighted in {@code SMALL_MATMUL_CATALOG.md} §4.5.
 */
public class TestRecursiveComposition {

	@Test
	public void strassen_222_squared_gives_49_over_reals() {
		Optional<Result> r = RecursiveComposition.evaluatePower(
				Algebra.nonCommutative(Field.R), Factor.cube(2), 2);
		assertThat(r).isPresent();
		assertThat(r.get().rank).isEqualTo(49);
	}

	@Test
	public void strassen_222_to_the_fifth_gives_16807_over_reals() {
		// Pure Strassen for ⟨32,32,32⟩.
		Optional<Result> r = RecursiveComposition.evaluatePower(
				Algebra.nonCommutative(Field.R), Factor.cube(2), 5);
		assertThat(r).isPresent();
		assertThat(r.get().rank).isEqualTo(16807L);
	}

	@Test
	public void alphatensor_444_over_z2_gives_47() {
		Optional<KnownAlgorithm> best =
				KnownAlgorithmCatalog.bestKnown(4, 4, 4, Algebra.nonCommutative(Field.F2));
		assertThat(best).isPresent();
		assertThat(best.get().rank).isEqualTo(47);
		assertThat(best.get().source).isEqualTo("AlphaTensor");
	}

	@Test
	public void alphaevolve_444_over_complex_gives_48() {
		Optional<KnownAlgorithm> best =
				KnownAlgorithmCatalog.bestKnown(4, 4, 4, Algebra.nonCommutative(Field.C));
		assertThat(best).isPresent();
		assertThat(best.get().rank).isEqualTo(48);
		assertThat(best.get().source).isEqualTo("AlphaEvolve");
		assertThat(best.get().year).isEqualTo(2025);
	}

	@Test
	public void mixed_composition_for_32_over_z2_via_alphatensor_plus_strassen() {
		// ⟨32,32,32⟩ = ⟨4,4,4⟩ · ⟨4,4,4⟩ · ⟨2,2,2⟩ over Z/2.
		// = 47 · 47 · 7 = 15,463 (improves Strassen⁵ = 16,807 by ~8%).
		List<Factor> factors = Arrays.asList(
				Factor.cube(4), Factor.cube(4), Factor.cube(2));
		Optional<Result> r = RecursiveComposition.evaluate(Algebra.nonCommutative(Field.F2), factors);
		assertThat(r).isPresent();
		assertThat(r.get().rank).isEqualTo(15_463L);
	}

	@Test
	public void mixed_composition_for_32_over_complex_via_alphaevolve_plus_strassen() {
		// ⟨32,32,32⟩ via ⟨4,4,4⟩² · ⟨2,2,2⟩ over C: 48 · 48 · 7 = 16,128.
		List<Factor> factors = Arrays.asList(
				Factor.cube(4), Factor.cube(4), Factor.cube(2));
		Optional<Result> r = RecursiveComposition.evaluate(Algebra.nonCommutative(Field.C), factors);
		assertThat(r).isPresent();
		assertThat(r.get().rank).isEqualTo(16_128L);
	}

	@Test
	public void real_ring_composition_falls_back_to_strassen_squared_for_4x4() {
		// Over R/Q/Z we don't have a 47-mult ⟨4,4,4⟩ — best is Strassen² = 49.
		Optional<KnownAlgorithm> best =
				KnownAlgorithmCatalog.bestKnown(4, 4, 4, Algebra.nonCommutative(Field.R));
		assertThat(best).isPresent();
		assertThat(best.get().rank).isEqualTo(49);
	}

	@Test
	public void laderman_333_listed_at_23() {
		Optional<KnownAlgorithm> best =
				KnownAlgorithmCatalog.bestKnown(3, 3, 3, Algebra.nonCommutative(Field.R));
		assertThat(best).isPresent();
		assertThat(best.get().rank).isEqualTo(23);
		assertThat(best.get().source).isEqualTo("Laderman");
		assertThat(best.get().year).isEqualTo(1976);
	}
}
