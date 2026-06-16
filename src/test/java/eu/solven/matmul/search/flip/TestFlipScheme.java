package eu.solven.matmul.search.flip;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

import java.util.List;
import java.util.Random;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;

/**
 * Flip-move soundness: every move (flip variants, reductions, splits) must
 * leave the represented tensor EXACTLY invariant — a silently-wrong move would
 * poison every downstream walk, which is precisely the class of silent
 * high-level bug CLAUDE.md requires a guard for. Invariance is asserted with
 * {@link FlipScheme#isExactIntTensor()} (exact integer arithmetic), not the
 * double-based {@code Verifier}, so coefficient growth cannot fake a failure.
 */
public class TestFlipScheme {

	private static NonCubicBilinearAlgorithm strassen;

	@BeforeAll
	static void setUp() {
		strassen = new FieldAwareLookup(Field.Z).find(2, 2, 2).orElseThrow();
	}

	@Test
	public void round_trip_preserves_scheme() {
		FlipScheme s = FlipScheme.of(strassen);
		assertThat(s.rank()).isEqualTo(strassen.r);
		assertThat(s.isExactIntTensor()).isTrue();
		assertThat(Verifier.isExactNonCubic(s.toAlgorithm())).isTrue();
	}

	@Test
	public void non_integer_seed_rejected() {
		double[][] u = { { 0.5 } };
		double[][] v = { { 1 } };
		double[][] w = { { 2 } };
		NonCubicBilinearAlgorithm half = new NonCubicBilinearAlgorithm(1, 1, 1, u, v, w);
		assertThatThrownBy(() -> FlipScheme.of(half))
				.isInstanceOf(IllegalArgumentException.class)
				.hasMessageContaining("non-integer");
	}

	/** One flip per (slot, variant), exact-checked in isolation — pinpoints a
	 *  broken move identity immediately instead of deep inside a random walk. */
	@ParameterizedTest(name = "slot {0} variant{1}")
	@CsvSource({ "U, false", "U, true", "V, false", "V, true", "W, false", "W, true" })
	public void single_flip_is_exact(String slotName, boolean variantB) {
		FlipScheme s = FlipScheme.of(NonCubicBilinearAlgorithm.naive(3, 3, 3));
		FlipScheme.Slot slot = FlipScheme.Slot.valueOf(slotName);
		List<int[]> classes = s.signClasses(slot);
		assertThat(classes).as("naive ⟨3,3,3⟩ must offer %s pivots", slot).isNotEmpty();
		int[] cls = classes.get(0);
		int sign = FlipScheme.signRatio(s.vec(slot, cls[0]), s.vec(slot, cls[1]));
		assertThat(sign).isNotZero();
		s.flip(slot, cls[0], cls[1], sign, variantB);
		assertThat(s.isExactIntTensor())
				.as("flip(%s, variantB=%s) broke the tensor", slot, variantB).isTrue();
	}

	@Test
	public void random_flips_preserve_exactness_unbounded() {
		assertFlipsStayExact(0, 150);
	}

	@Test
	public void random_flips_preserve_exactness_ternary_cap() {
		assertFlipsStayExact(1, 50);
	}

	private static void assertFlipsStayExact(int cap, int minApplied) {
		// naive ⟨3,3,3⟩: rich sign-sharing (every elementary vector is shared),
		// so the walk genuinely exercises all slots even under the ternary cap.
		FlipScheme s = FlipScheme.of(NonCubicBilinearAlgorithm.naive(3, 3, 3));
		Random rng = new Random(42);
		int applied = 0;
		for (int step = 0; step < 600; step++) {
			FlipScheme.Slot slot = FlipScheme.Slot.values()[rng.nextInt(3)];
			List<int[]> classes = s.signClasses(slot);
			if (classes.isEmpty()) {
				continue;
			}
			int[] cls = classes.get(rng.nextInt(classes.size()));
			int i = cls[rng.nextInt(cls.length)];
			int j = cls[rng.nextInt(cls.length)];
			if (i == j) {
				continue;
			}
			int sign = FlipScheme.signRatio(s.vec(slot, i), s.vec(slot, j));
			if (s.flipWithinCap(slot, i, j, sign, rng.nextBoolean(), cap)) {
				applied++;
				if (applied % 25 == 0) {
					assertThat(s.isExactIntTensor())
							.as("tensor must be invariant after %d flips (cap=%d)", applied, cap)
							.isTrue();
				}
			}
		}
		assertThat(applied)
				.as("the loop must actually exercise flips (cap=%d)", cap)
				.isGreaterThanOrEqualTo(minApplied);
		assertThat(s.isExactIntTensor()).isTrue();
		if (cap > 0) {
			assertThat(s.maxAbsCoefficient()).isLessThanOrEqualTo(cap);
		}
	}

	@Test
	public void split_then_reduce_round_trips() {
		// Strassen's vectors have multi-entry supports — splittable; naive's are
		// all single-entry — NOT splittable. Both directions are asserted.
		FlipScheme s = FlipScheme.of(strassen);
		int r0 = s.rank();
		assertThat(s.split(new Random(7))).isTrue();
		assertThat(s.rank()).isEqualTo(r0 + 1);
		assertThat(s.isExactIntTensor())
				.as("split (plus transition) must preserve the tensor").isTrue();
		// The two halves still share the OTHER two slots, so reduce() re-merges.
		s.reduce();
		assertThat(s.rank()).isLessThanOrEqualTo(r0);
		assertThat(s.isExactIntTensor()).isTrue();

		FlipScheme naive = FlipScheme.of(NonCubicBilinearAlgorithm.naive(2, 2, 2));
		assertThat(naive.split(new Random(7)))
				.as("naive vectors are single-entry — nothing to split").isFalse();
	}

	@Test
	public void reduce_is_noop_on_irreducible_seed() {
		// Strassen r=7 is optimal for ⟨2,2,2⟩ — reduce() must not "find" a merge.
		FlipScheme s = FlipScheme.of(strassen);
		int r0 = s.rank();
		s.reduce();
		assertThat(s.rank()).isEqualTo(r0);
		assertThat(s.isExactIntTensor()).isTrue();
	}
}
