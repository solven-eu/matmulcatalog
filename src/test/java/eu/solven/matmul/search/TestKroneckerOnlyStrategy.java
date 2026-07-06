package eu.solven.matmul.search;

import static org.assertj.core.api.Assertions.assertThat;

import java.nio.file.Path;
import java.util.List;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.verifiers.Verifier;

/**
 * Guards the {@link RecursiveMaterialiser#STRAT_KRONECKER} strategy: a PLAIN
 * Kronecker-only operator (no block-split search, no bud fusion), selected via
 * {@code --strategies=kron|kronecker}. Historically "kron" silently aliased the
 * full recombination strategy — a restricted sweep the user asked for did far
 * more (and far costlier) work than requested.
 */
public class TestKroneckerOnlyStrategy {

	private static final FieldAwareLookup LOOKUP = new FieldAwareLookup("Q");

	private RecursiveMaterialiser kronOnly() {
		List<eu.solven.matmul.recombination.BlockSplitSearch.NamedBase> pool = List.of();
		RecursiveClosureSota sota = new RecursiveClosureSota(LOOKUP, pool, true, true);
		RecursiveMaterialiser m = new RecursiveMaterialiser(LOOKUP, pool, sota,
				Path.of("target/test-kron-schemes"), false, false);
		m.setStrategies(java.util.Set.of(RecursiveMaterialiser.STRAT_KRONECKER));
		return m;
	}

	@Test
	public void kron_only_builds_strassen_squared_for_4x4x4() {
		// ⟨4,4,4⟩ = ⟨2,2,2⟩⊗⟨2,2,2⟩ = 7·7 — the canonical plain-Kron product. An
		// unbounded upper forces the build; the result must be a verified true
		// Kronecker composition (rank exactly 49 over Q — bud fusion could do no
		// better here and recombination candidates must not leak in).
		RecursiveMaterialiser.Result r = kronOnly().tryKronecker(4, 4, 4, Long.MAX_VALUE);
		assertThat(r).isNotNull();
		assertThat(r.alg().n).isEqualTo(4);
		assertThat(r.alg().m).isEqualTo(4);
		assertThat(r.alg().p).isEqualTo(4);
		assertThat(r.alg().r).isEqualTo(49);
		assertThat(r.lineage()).isInstanceOf(Lineage.KronProduct.class);
		assertThat(Verifier.passesRandomMatmulSpotCheck(r.alg())).isTrue();
	}

	@Test
	public void kron_only_respects_the_upper_bound() {
		// upper == the plain-Kron optimum (49) → no strict improvement possible →
		// the strategy must decline rather than return a tie.
		assertThat(kronOnly().tryKronecker(4, 4, 4, 49L)).isNull();
	}

	@Test
	public void kron_only_handles_prime_axes_via_unit_factors() {
		// ⟨3,3,6⟩: n=m=3 are prime, so the only Kronecker routes go through unit
		// factor pairs (e.g. ⟨1,1,2⟩⊗⟨3,3,3⟩ = 2·23 = 46). SOTA-or-better, not
		// equality: a catalog improvement must never break this guard.
		RecursiveMaterialiser.Result r = kronOnly().tryKronecker(3, 3, 6, Long.MAX_VALUE);
		assertThat(r).isNotNull();
		assertThat(r.alg().r).isLessThanOrEqualTo(46);
		assertThat(Verifier.passesRandomMatmulSpotCheck(r.alg())).isTrue();
	}
}
