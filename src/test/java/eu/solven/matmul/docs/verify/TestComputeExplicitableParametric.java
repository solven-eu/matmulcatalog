package eu.solven.matmul.docs.verify;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.search.LineageReplayer;

/**
 * Guards the 2026-07-07 classifier gap: {@code ComputeExplicitable} treated the
 * parametric-constructor refs ({@code TA_lita(n=N)}, {@code DIS09Lemma4(n=N)}) as
 * unknown → cited-bound, although {@code LineageReplayer.resolveParametric} rebuilds
 * them deterministically from the formula with no catalog resolution. One TA-cube
 * ancestor (⟨28,28,28⟩=10535 = {@code TA_lita(n=28)}) then poisoned EVERY downstream
 * composite to {@code explicitable:false}. The recognizer is shared
 * ({@link LineageReplayer#isParametricRef}) so the two stay in sync by construction.
 */
public class TestComputeExplicitableParametric {

	@Test
	public void parametric_constructor_refs_are_recognised() {
		assertThat(LineageReplayer.isParametricRef("TA_lita(n=28)")).isTrue();
		assertThat(LineageReplayer.isParametricRef("DIS09Lemma4(n=5)")).isTrue();
	}

	@Test
	public void non_parametric_refs_are_not() {
		assertThat(LineageReplayer.isParametricRef("28x28x28")).isFalse();
		assertThat(LineageReplayer.isParametricRef("TA_lita(28)")).isFalse();
		assertThat(LineageReplayer.isParametricRef("ext[3x4x7-r63.json]")).isFalse();
		assertThat(LineageReplayer.isParametricRef("3x4x7@ac0e1ad")).isFalse();
	}
}
