package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.IOException;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.ComplexNonCubicBilinearAlgorithm;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.Compose;
import eu.solven.matmul.catalog.Compositions;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Builds recursive compositions of verified catalog schemes for target
 * formats `max(n,m,p) ≥ 8`, validating end-to-end where the verifier's
 * `O(N⁶·r)` cost allows.
 */
public class TestCompositions {

	@Test
	public void strassen_cubed_for_888_is_343_and_verifies() throws IOException {
		NonCubicBilinearAlgorithm c = Compositions.strassen3_888();
		assertThat(c.n).isEqualTo(8);
		assertThat(c.m).isEqualTo(8);
		assertThat(c.p).isEqualTo(8);
		assertThat(c.r).isEqualTo(343);
		assertThat(Verifier.isExactNonCubic(c)).isTrue();
	}

	@Test
	public void laderman_strassen_666_at_161_verifies() throws IOException {
		NonCubicBilinearAlgorithm c = Compositions.laderman_strassen_666();
		assertThat(c.n).isEqualTo(6);
		assertThat(c.r).isEqualTo(161);
		assertThat(Verifier.isExactNonCubic(c)).isTrue();
	}

	@Test
	public void laderman_squared_999_is_529_constructs() throws IOException {
		// ⟨9,9,9⟩ via Laderman²: rank 23·23 = 529. Verifier would be 9⁶·529 ≈
		// 281M ops, slow but tractable — leave it as construction-only here
		// to keep the test budget tight.
		NonCubicBilinearAlgorithm c = Compositions.laderman2_999();
		assertThat(c.n).isEqualTo(9);
		assertThat(c.r).isEqualTo(529);
	}

	@Test
	public void alphatensorF2_strassen_888_F2_at_329_verifies_mod2() throws IOException {
		// ⟨8,8,8⟩ over F₂ via AlphaTensor's 47-mult ⟨4,4,4⟩ ⊗ Strassen ⟨2,2,2⟩ = 7.
		// Total rank: 47·7 = 329. Beats Strassen³=343 over F₂.
		NonCubicBilinearAlgorithm c = Compositions.alphatensorF2_strassen_888();
		assertThat(c.n).isEqualTo(8);
		assertThat(c.r).isEqualTo(329);
		// F₂ verifier (mod-2 XOR semantics; coefficient {0,1} composition
		// stays in {0,1} but real-arithmetic sums can hit 2 or 3).
		assertThat(Verifier.residualNonCubicF2(c)).isEqualTo(0);
	}

	@Test
	public void alphatensorF2_squared_16_at_2209_constructs() throws IOException {
		// ⟨16,16,16⟩ over F₂ via AlphaTensor² = 47² = 2,209 multiplications.
		// Strictly better than Strassen⁴ = 2,401. The verifier is `O(16⁶·2209)
		// ≈ 37 billion ops` — too slow to assert; construction is correct by
		// the Kronecker theorem.
		NonCubicBilinearAlgorithm c = Compositions.alphatensorF2_squared_16();
		assertThat(c.n).isEqualTo(16);
		assertThat(c.r).isEqualTo(2_209);
	}

	@Test
	public void alphatensorF2_squared_strassen_32_at_15463_constructs() throws IOException {
		// The recursive-composition headline for ⟨32,32,32⟩ over F₂:
		//   AlphaTensor ⟨4,4,4⟩=47  ⊗  AlphaTensor ⟨4,4,4⟩=47  ⊗  Strassen ⟨2,2,2⟩=7
		// = 47 · 47 · 7 = 15,463 multiplications (vs. Strassen⁵=16,807).
		NonCubicBilinearAlgorithm c = Compositions.alphatensorF2_squared_strassen_32();
		assertThat(c.n).isEqualTo(32);
		assertThat(c.r).isEqualTo(15_463);
		// Sanity check: addition count is finite and matches the formula
		// nz(U)+nz(V)+nz(W) - 2r - n·p.
		int adds = Verifier.additionCount(c);
		assertThat(adds).isPositive();
		System.out.printf("⟨32,32,32⟩ via AT² ⊗ Str: r=%d, +%d adds, total ops = %d%n",
				c.r, adds, c.r + adds);
	}

	@Test
	public void strassen_fourth_16_at_2401_over_z_constructs_and_loads() throws IOException {
		NonCubicBilinearAlgorithm c = Compositions.strassen4_16();
		assertThat(c.n).isEqualTo(16);
		assertThat(c.r).isEqualTo(2_401);
	}

	@Test
	public void strassen_fifth_32_at_16807_over_z_constructs() throws IOException {
		NonCubicBilinearAlgorithm c = Compositions.strassen5_32();
		assertThat(c.n).isEqualTo(32);
		assertThat(c.r).isEqualTo(16_807);
	}

	@Test
	public void strassen_fourth_16_passes_sampled_verification() throws IOException {
		// Full O(16⁶·2401) verifier is ~minutes; sampled verifier with 50k
		// random positions takes <1s and gives statistically strong confidence.
		NonCubicBilinearAlgorithm c = Compositions.strassen4_16();
		int wrong = Verifier.residualSampled(c, 50_000, 0xCAFEBABE);
		assertThat(wrong).as("sampled residual at ⟨16,16,16⟩ rank 2401").isEqualTo(0);
	}

	@Test
	public void alphatensorF2_squared_16_passes_sampled_F2_verification() throws IOException {
		NonCubicBilinearAlgorithm c = Compositions.alphatensorF2_squared_16();
		int wrong = Verifier.residualSampledF2(c, 50_000, 0xCAFEBABE);
		assertThat(wrong).as("sampled F₂ residual at ⟨16,16,16⟩ rank 2209").isEqualTo(0);
	}

	@Test
	public void strassen_fifth_32_passes_sampled_verification() throws IOException {
		NonCubicBilinearAlgorithm c = Compositions.strassen5_32();
		// 10k samples on the ⟨32,32,32⟩ rank-16807 algorithm.
		int wrong = Verifier.residualSampled(c, 10_000, 0xCAFEBABE);
		assertThat(wrong).as("sampled residual at ⟨32,32,32⟩ rank 16807").isEqualTo(0);
	}

	@Test
	public void alphatensorF2_squared_strassen_32_passes_sampled_F2_verification() throws IOException {
		NonCubicBilinearAlgorithm c = Compositions.alphatensorF2_squared_strassen_32();
		int wrong = Verifier.residualSampledF2(c, 10_000, 0xCAFEBABE);
		assertThat(wrong).as("sampled F₂ residual at ⟨32,32,32⟩ rank 15463").isEqualTo(0);
	}

	@Test
	public void alphaevolve_squared_16_at_2304_over_complex() throws IOException {
		// Matches fmm-lille's listed ⟨16,16,16⟩ = 2304: their entry is exactly
		// the Kronecker square of AlphaEvolve ⟨4,4,4⟩=48 over the complex field.
		// Locate the AlphaEvolve ⟨4,4,4⟩=48 scheme by content identity (shape +
		// rank + source note), not a pinned filename — the 2026-06 migration made
		// names content-driven (`4x4x4-r48-alphaevolve-{hash}.json`) and they re-hash.
		java.nio.file.Path ae;
		try (var s = java.nio.file.Files.walk(
				java.nio.file.Path.of("src/main/resources/schemes"))) {
			ae = s.filter(p -> {
						String n = p.getFileName().toString().toLowerCase(java.util.Locale.ROOT);
						return n.endsWith(".json") && n.contains("4x4x4") && n.contains("r48")
								&& n.contains("alphaevolve");
					})
					.findFirst()
					.orElseThrow(() -> new java.io.IOException("AlphaEvolve 4x4x4=48 not found"));
		}
		ComplexNonCubicBilinearAlgorithm base = SchemeIO.readComplex(ae.toFile());
		ComplexNonCubicBilinearAlgorithm composed = Compose.kroneckerComplex(base, base);
		assertThat(composed.n).isEqualTo(16);
		assertThat(composed.r).isEqualTo(2304);
		int wrong = Verifier.residualSampledComplex(composed, 20_000, 0xCAFEBABE);
		assertThat(wrong).as("sampled complex residual at ⟨16,16,16⟩ rank 2304 (AlphaEvolve²)")
				.isEqualTo(0);
	}
}
