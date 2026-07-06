package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Optional;
import java.util.stream.Stream;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.Compose;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.recombination.Recombination.AlgorithmLookup;

/**
 * Validates {@link Compose#blockSplitCubic} — the constructive procedure
 * for building cubic ⟨n,n,n⟩ algorithms via 2×2×2 block decomposition
 * (naïve, no Sedoglavic algebraic saving yet).
 *
 * For ⟨7,7,7⟩ with split 7 = 4+3, expect rank
 * {@code R(⟨4,4,4⟩) + 3·R(⟨4,4,3⟩) + 3·R(⟨4,3,3⟩) + R(⟨3,3,3⟩) = 49 + 87 + 114 + 23 = 273}.
 */
public class TestBlockSplit {

	// Content-driven catalog-best resolver over Q (char-0 default field; F2-only
	// schemes excluded via fields[], not the long-dead filename filter — post the
	// 2026-06 content-driven migration an F2-only file no longer says "F2" in its
	// name and could silently win the min-rank pick here).
	private static final FieldAwareLookup LOOKUP = new FieldAwareLookup("Q");

	private static AlgorithmLookup catalogLookup() {
		return LOOKUP::find;
	}

	/** Rank of the catalog-best leaf the lookup will hand blockSplitCubic. */
	private static int leaf(int n, int m, int p) {
		return catalogLookup().find(n, m, p).orElseThrow().r;
	}

	/**
	 * Block-split UNIFORM cases (u = v) — sub-products are all ⟨n/2,n/2,n/2⟩
	 * which match exactly in our catalog (no axis permutation needed).
	 * Rank is suboptimal (worse than pure Kronecker), but this validates
	 * that the constructive block-split procedure works end-to-end.
	 */
	@Test
	public void uniform_block_split_444_2_2_constructs_and_verifies_at_56() throws IOException {
		// ⟨4,4,4⟩ with u=v=2: 8 sub-products of ⟨2,2,2⟩ at rank 7 each.
		// Total = 56 (worse than Strassen²=49 but constructive proof block-split works).
		NonCubicBilinearAlgorithm alg = Compose.blockSplitCubic(4, 2, 2, catalogLookup());
		assertThat(alg.n).isEqualTo(4);
		assertThat(alg.r).isEqualTo(8 * leaf(2, 2, 2));
		assertThat(Verifier.isExactNonCubic(alg)).isTrue();
	}

	@Test
	public void uniform_block_split_666_3_3_constructs_and_verifies_at_184() throws IOException {
		// ⟨6,6,6⟩ with u=v=3: 8 sub-products of ⟨3,3,3⟩=23 (Laderman).
		// Total = 184 (worse than fmm-lille's 153 best, but again constructive validation).
		NonCubicBilinearAlgorithm alg = Compose.blockSplitCubic(6, 3, 3, catalogLookup());
		assertThat(alg.n).isEqualTo(6);
		assertThat(alg.r).isEqualTo(8 * leaf(3, 3, 3));
		assertThat(Verifier.isExactNonCubic(alg)).isTrue();
	}

	/**
	 * Mixed split — exercises axis-permutation via
	 * {@link NonCubicBilinearAlgorithm#orientAs}. For ⟨7,7,7⟩ with u=4, v=3,
	 * formula: R(⟨4,4,4⟩=49) + 3·R(⟨4,4,3⟩=29) + 3·R(⟨4,3,3⟩=38) + R(⟨3,3,3⟩=23)
	 * but using canonical sorted ranks this becomes
	 * 49 + 3·29 + 3·38 + 23 = 273 — the naïve block decomp.
	 */
	@Test
	public void mixed_block_split_777_4_3_constructs_and_verifies_at_273() throws IOException {
		NonCubicBilinearAlgorithm alg = Compose.blockSplitCubic(7, 4, 3, catalogLookup());
		assertThat(alg.n).isEqualTo(7);
		// SELF-CONSISTENT expectation: the block decomposition must cost exactly the
		// sum of the leaves the lookup resolves — an invariant of the construction,
		// not of catalog state. (The old hand-pinned constant broke every time the
		// catalog gained a better mixed-shape scheme: 273 → 272 → 271 → …)
		int expected = leaf(4, 4, 4) + 3 * leaf(4, 4, 3) + 3 * leaf(4, 3, 3) + leaf(3, 3, 3);
		assertThat(alg.r).isEqualTo(expected);
		assertThat(alg.r).isLessThanOrEqualTo(273);   // never worse than the 2017-era leaves
		// Full residual at ⟨7,7,7⟩ ≈ 32M ops — acceptable to assert in a test.
		assertThat(Verifier.isExactNonCubic(alg)).isTrue();
	}
}
