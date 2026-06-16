package eu.solven.matmul.algebra;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

import org.junit.jupiter.api.Test;

/**
 * Coverage for the type-safe algebra/problem types — fields,
 * commutativity, shape — that should be used throughout the project
 * in place of string-tagged field arguments.
 */
public class TestAlgebraTypes {

	@Test
	public void field_fallback_chain() {
		assertThat(Field.F2.fallbackChain()).containsExactly(Field.F2);
		assertThat(Field.F3.fallbackChain()).containsExactly(Field.F3);
		assertThat(Field.Z.fallbackChain()).containsExactly(Field.Z);
		assertThat(Field.Q.fallbackChain()).containsExactly(Field.Q, Field.Z);
		assertThat(Field.R.fallbackChain()).containsExactly(Field.R, Field.Q, Field.Z);
		assertThat(Field.C.fallbackChain()).containsExactly(Field.C, Field.R, Field.Q, Field.Z);
	}

	@Test
	public void field_fromTag_accepts_historical_abbreviations() {
		assertThat(Field.fromTag("F2")).isEqualTo(Field.F2);
		assertThat(Field.fromTag("Z2")).isEqualTo(Field.F2);
		assertThat(Field.fromTag("Z")).isEqualTo(Field.Z);
		assertThat(Field.fromTag("Q")).isEqualTo(Field.Q);
		assertThat(Field.fromTag("R")).isEqualTo(Field.R);
		assertThat(Field.fromTag("R/Q/Z")).isEqualTo(Field.R); // historical lumped → widest
		assertThat(Field.fromTag("ZT")).isEqualTo(Field.Z);    // ZT = ternary {-1,0,1} sub-class of Z (NOT Q, NOT F2)
		assertThat(Field.fromTag("C")).isEqualTo(Field.C);
		assertThatThrownBy(() -> Field.fromTag("foo"))
				.isInstanceOf(IllegalArgumentException.class);
	}

	@Test
	public void algebra_accepts_within_chain() {
		Algebra nc_R = Algebra.nonCommutative(Field.R);
		Algebra nc_C = Algebra.nonCommutative(Field.C);
		Algebra nc_F2 = Algebra.nonCommutative(Field.F2);
		Algebra cmt_R = Algebra.commutative(Field.R);

		// An NC R-scheme can be used as NC over C (R ⊂ C-friendly).
		assertThat(nc_C.accepts(nc_R)).isTrue();
		// An NC C-scheme can NOT be used as NC over R (C is not in R's chain).
		assertThat(nc_R.accepts(nc_C)).isFalse();
		// Cross-characteristic forbidden.
		assertThat(nc_R.accepts(nc_F2)).isFalse();
		assertThat(nc_F2.accepts(nc_R)).isFalse();
		// Commutative-only can't substitute for NC-required.
		assertThat(nc_R.accepts(cmt_R)).isFalse();
		// NC can substitute for commutative.
		assertThat(cmt_R.accepts(nc_R)).isTrue();
	}

	@Test
	public void problem_shape_and_sortedShape() {
		MatmulProblem p = MatmulProblem.nc(Field.R, 5, 3, 7);
		assertThat(p.shape()).containsExactly(5, 3, 7);
		assertThat(p.sortedShape()).containsExactly(3, 5, 7);
		assertThat(p.maxDim()).isEqualTo(7);
		assertThat(p.isCubic()).isFalse();
		assertThat(p.formatTag()).isEqualTo("⟨5,3,7⟩");
	}

	@Test
	public void problem_isCubic_when_all_axes_equal() {
		assertThat(MatmulProblem.nc(Field.R, 3, 3, 3).isCubic()).isTrue();
		assertThat(MatmulProblem.nc(Field.R, 2, 2, 3).isCubic()).isFalse();
	}

	@Test
	public void problem_rejects_invalid_shape() {
		assertThatThrownBy(() -> MatmulProblem.nc(Field.R, 0, 2, 2))
				.isInstanceOf(IllegalArgumentException.class);
		assertThatThrownBy(() -> MatmulProblem.nc(Field.R, -1, 2, 2))
				.isInstanceOf(IllegalArgumentException.class);
	}

	@Test
	public void problem_toString_carries_algebra_and_shape() {
		MatmulProblem p = MatmulProblem.cmt(Field.F2, 3, 3, 3);
		assertThat(p.toString()).contains("⟨3,3,3⟩").contains("F2").contains("cmt");
	}
}
