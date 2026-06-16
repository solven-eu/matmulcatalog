package eu.solven.matmul;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

import org.junit.jupiter.api.Test;

public class TestShape {

	@Test
	public void basics() {
		Shape s = Shape.of(3, 5, 7);
		assertThat(s.n()).isEqualTo(3);
		assertThat(s.maxDim()).isEqualTo(7);
		assertThat(s.minDim()).isEqualTo(3);
		assertThat(s.volume()).isEqualTo(105L);
		assertThat(s.isCubic()).isFalse();
		assertThat(Shape.of(4, 4, 4).isCubic()).isTrue();
		assertThat(s.toString()).isEqualTo("⟨3,5,7⟩");
		assertThat(s.toFilenameToken()).isEqualTo("3x5x7");
	}

	@Test
	public void order_is_significant_but_equality_is_structural() {
		assertThat(Shape.of(3, 5, 7)).isEqualTo(Shape.of(3, 5, 7));
		assertThat(Shape.of(3, 5, 7)).isNotEqualTo(Shape.of(7, 5, 3));
	}

	@Test
	public void canonical_sorts_and_is_a_distinct_type() {
		CanonicalShape c = Shape.of(7, 3, 5).canonical();
		assertThat(c).isEqualTo(Shape.of(3, 7, 5).canonical());
		assertThat(c.shape()).isEqualTo(Shape.of(3, 5, 7));
		assertThat(c.maxDim()).isEqualTo(7);
		assertThatThrownBy(() -> new CanonicalShape(5, 3, 7)).isInstanceOf(IllegalArgumentException.class);
	}

	@Test
	public void divides_and_cofactor() {
		Shape base = Shape.of(2, 3, 5);
		Shape target = Shape.of(6, 6, 15);
		assertThat(base.divides(target)).isTrue();
		assertThat(base.cofactor(target)).isEqualTo(Shape.of(3, 2, 3));
		assertThat(Shape.of(4, 1, 1).divides(target)).isFalse();
		assertThatThrownBy(() -> Shape.of(4, 1, 1).cofactor(target)).isInstanceOf(IllegalArgumentException.class);
	}

	@Test
	public void parse_and_array_roundtrip() {
		assertThat(Shape.parse("3x5x7")).isEqualTo(Shape.of(3, 5, 7));
		assertThat(Shape.parse("⟨2, 4, 6⟩")).isEqualTo(Shape.of(2, 4, 6));
		assertThat(Shape.parse("derived_recursive-2x3x17_m85")).isEqualTo(Shape.of(2, 3, 17));
		assertThat(Shape.ofArray(Shape.of(9, 8, 7).toArray())).isEqualTo(Shape.of(9, 8, 7));
	}

	@Test
	public void rejects_degenerate_dims() {
		assertThatThrownBy(() -> Shape.of(0, 2, 2)).isInstanceOf(IllegalArgumentException.class);
	}

	@Test
	public void compareTo_is_maxdim_then_volume_then_lex() {
		assertThat(Shape.of(2, 2, 2)).isLessThan(Shape.of(2, 2, 3));   // maxDim 2 < 3
		assertThat(Shape.of(2, 3, 6)).isLessThan(Shape.of(4, 5, 6));   // same maxDim 6, smaller volume
	}
}
