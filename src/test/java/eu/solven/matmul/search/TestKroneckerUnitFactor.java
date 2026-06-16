package eu.solven.matmul.search;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.Optional;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Recombination.AlgorithmLookup;
import eu.solven.matmul.catalog.Recombination.SotaResolver;

/**
 * Unit Kronecker factors: ⟨3,3,18⟩ = ⟨1,1,3⟩⊗⟨3,3,6⟩ = 3·40 = 120. Before the
 * fix the ≥2-only factor enumeration returned empty (n=m=3 prime) and the
 * decomposition was unreachable; now it must be found AND materialise to a
 * Verifier-passing scheme (the degenerate ⟨1,1,3⟩ factor synthesised as naive).
 */
public class TestKroneckerUnitFactor {

	private static SotaResolver sota(FieldAwareLookup lk) {
		return (p, q, r) -> {
			if (p == 0 || q == 0 || r == 0) return 0;
			if (p == 1) return q * r;
			if (q == 1) return p * r;
			if (r == 1) return p * q;
			int v = lk.findRank(p, q, r);
			return v >= Integer.MAX_VALUE / 100 ? p * q * r : v;
		};
	}

	@Test
	public void findsUnitFactorDecomposition() {
		FieldAwareLookup lk = new FieldAwareLookup("R");
		Optional<KroneckerSplitSearch.KroneckerSplit> best =
				KroneckerSplitSearch.findBest(3, 3, 18, sota(lk));
		assertThat(best).as("⟨3,3,18⟩ must admit a Kronecker split now").isPresent();
		// ⟨3,3,6⟩=40 on disk → 3·40 = 120 (matches FMM-Lille).
		assertThat(best.get().totalRank()).isEqualTo(120L);
	}

	@Test
	public void materialisesAndVerifies() throws Exception {
		FieldAwareLookup lk = new FieldAwareLookup("R");
		// AlgorithmLookup over R (find returns the actual on-disk scheme).
		AlgorithmLookup atoms = lk;
		Optional<KroneckerSplitSearch.KroneckerSplit> best =
				KroneckerSplitSearch.findBest(3, 3, 18, sota(lk));
		assertThat(best).isPresent();
		NonCubicBilinearAlgorithm alg = KroneckerSplitSearch.materialise(best.get(), atoms);
		assertThat(alg.n).isEqualTo(3);
		assertThat(alg.m).isEqualTo(3);
		assertThat(alg.p).isEqualTo(18);
		assertThat(alg.r).isEqualTo(120);
		assertThat(Verifier.isExactNonCubic(alg)).as("⟨3,3,18⟩=120 must verify").isTrue();
	}

	@Test
	public void naiveNonCubicVerifies() {
		// the degenerate factor builder itself must be a valid ⟨1,1,3⟩ scheme.
		NonCubicBilinearAlgorithm s = eu.solven.matmul.NaiveMatMul.ofNonCubic(1, 1, 3);
		assertThat(s.r).isEqualTo(3);
		assertThat(Verifier.isExactNonCubic(s)).isTrue();
		// and a genuinely rectangular one (⟨2,3,4⟩=24 naive).
		NonCubicBilinearAlgorithm r = eu.solven.matmul.NaiveMatMul.ofNonCubic(2, 3, 4);
		assertThat(r.r).isEqualTo(24);
		assertThat(Verifier.isExactNonCubic(r)).isTrue();
	}
}
