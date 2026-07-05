package eu.solven.matmul.search;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;

/**
 * The resolveLeaf shape-extraction fallback: source-prefixed / canonical-key /
 * named leaf refs that the sweep surfaced as un-replayable must now resolve by
 * shape (so projection parents stop being skipped). Each resolves to a valid
 * scheme of the embedded shape.
 *
 * <p>Pinned {@code shape@hash} refs are the OPPOSITE contract: a resolvable pin
 * replays the exact content; a dangling pin THROWS — never a silent shape-best
 * substitution ("phantom replay", task #91 ⟨17,22,29⟩ 6129→6138; policy in
 * references/PURGE_REFCOUNT_POLICY.md).</p>
 */
public class TestReplayLeafFallback {

	private final FieldAwareLookup lookup = new FieldAwareLookup(Field.Q);
	private final LineageReplayer replayer = LineageReplayer.withDefaultPool(lookup);

	private NonCubicBilinearAlgorithm replay(String ref) {
		return replayer.replay(new Lineage.Atom(ref));
	}

	@Test
	public void source_prefixed_ref_resolves_by_shape() {
		NonCubicBilinearAlgorithm a = replay("alphatensor_Z-2x3x3_m15_a58");
		assertThat(new int[] { a.n, a.m, a.p }).containsExactly(2, 3, 3);
	}

	@Test
	public void perminov_canonical_key_ref_resolves() {
		NonCubicBilinearAlgorithm a = replay("perminov_Z-2x5x15_m118_a2913");
		assertThat(new int[] { a.n, a.m, a.p }).containsExactly(2, 5, 15);
	}

	@Test
	public void angle_named_ref_resolves() {
		NonCubicBilinearAlgorithm a = replay("Strassen<2,2,2>=7");
		assertThat(new int[] { a.n, a.m, a.p }).containsExactly(2, 2, 2);
	}

	private static final Pattern FILE_HASH7 = Pattern.compile("-([0-9a-f]{7})\\.json$");

	@Test
	public void valid_pinned_ref_resolves() {
		// Derive a real pin from a catalog filename ({...}-hash7.json) so the test
		// survives re-hashing: the filename hash7 is a prefix of the content hash.
		String hash7 = lookup.findFiles(2, 2, 2).stream()
				.map(p -> FILE_HASH7.matcher(p.getFileName().toString()))
				.filter(Matcher::find)
				.map(m -> m.group(1))
				.findFirst()
				.orElseThrow();
		NonCubicBilinearAlgorithm a = replay("2x2x2@" + hash7);
		assertThat(new int[] { a.n, a.m, a.p, a.r }).containsExactly(2, 2, 2, 7);
	}

	@Test
	public void dangling_pinned_ref_throws_instead_of_shape_best_substitution() {
		assertThatThrownBy(() -> replay("2x2x2@0000000000000000"))
				.isInstanceOf(IllegalStateException.class)
				.hasMessageContaining("2x2x2@0000000000000000")
				.hasMessageContaining("RepinDanglingLineageRefs");
	}
}
