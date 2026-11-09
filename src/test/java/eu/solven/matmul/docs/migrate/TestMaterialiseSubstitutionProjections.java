package eu.solven.matmul.docs.migrate;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.LinkedHashMap;
import java.util.Optional;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.search.LineageReplayer;
import eu.solven.matmul.verifiers.Verifier;

/**
 * The substitution-projection closure against the real catalog (SOTA-or-better guards, per
 * CLAUDE.md): the FMM {@code ProjOfAlongMode} rows of issue #23 must stay reachable from the
 * parents on disk, whatever the representative those parents are stored as. A regression here
 * means a parent stopped replaying or the drop stopped killing its guaranteed products.
 */
public class TestMaterialiseSubstitutionProjections {

	private static final FieldAwareLookup LOOKUP = new FieldAwareLookup("Q");
	private static final LineageReplayer REPLAYER = LineageReplayer.withDefaultPool(LOOKUP);

	private static Optional<MaterialiseSubstitutionProjections.Win> best(int n, int m, int p) {
		return MaterialiseSubstitutionProjections.bestFor(new int[] { n, m, p }, Long.MAX_VALUE, LOOKUP, REPLAYER,
				new LinkedHashMap<>());
	}

	/** ⟨13,14,14⟩ from the ⟨14,14,14⟩ cube: FMM 1579 (= 1593 − 14); the drop reached 1578 on 2026-11-09. */
	@Test
	public void cube_14_drops_to_13x14x14_at_or_below_fmm() {
		Optional<MaterialiseSubstitutionProjections.Win> w = best(13, 14, 14);
		assertThat(w).isPresent();
		assertThat(w.get().alg().r).isLessThanOrEqualTo(1579);
		assertThat(Verifier.isExactNonCubic(w.get().alg())).isTrue();
		assertThat(w.get().lineage()).isInstanceOf(Lineage.SubstProject.class);
		// the stub replays to the very scheme the driver built
		assertThat(REPLAYER.replay(w.get().lineage()).r).isEqualTo(w.get().alg().r);
	}

	/** ⟨9,17,17⟩ from the ⟨9,18,18⟩ Kronecker stub: FMM 1574 (two drops, m then p); ours 1570. */
	@Test
	public void kronecker_9x18x18_drops_to_9x17x17_at_or_below_fmm() {
		Optional<MaterialiseSubstitutionProjections.Win> w = best(9, 17, 17);
		assertThat(w).isPresent();
		assertThat(w.get().alg().r).isLessThanOrEqualTo(1574);
		assertThat(Verifier.isExactNonCubic(w.get().alg())).isTrue();
	}
}
