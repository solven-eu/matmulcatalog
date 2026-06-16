package eu.solven.matmul.docs;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;

/**
 * Regression test for the serendipitous sweep's win-criterion (#159).
 *
 * <p><strong>The bug this guards (silent regression):</strong> the SOTA oracle
 * — the bar a serendipitous {@code r_s} must beat — was built by a filename
 * regex {@code (\d+)x(\d+)x(\d+)_m(\d+)}. The 2026-06 {@code -r{rank}-} rename
 * matched 0 of 14,556 files → empty index → {@code trueSota} returned -1 for
 * every shape → every candidate was skipped → the sweep could never report a
 * win, silently. The prior test exercised the win <em>logic</em> on synthetic
 * maps and so never touched the broken index build. These tests run the real
 * content-driven oracle ({@link FieldAwareLookup#findRank}).</p>
 *
 * <p>Assertions are SOTA-or-better ({@code >0}, {@code ≤}), never equality, so
 * a genuine catalog improvement never breaks the test — only a regression does.</p>
 */
public class TestSerendipitousSweep {

	static FieldAwareLookup lk;

	@BeforeAll
	static void setup() {
		lk = new FieldAwareLookup(Field.Q);
	}

	@Test
	public void sota_oracle_is_content_driven() {
		// The regression guard: known composite shapes MUST resolve to a real
		// SOTA. If the oracle were dead (the filename-regex bug), these were -1.
		assertThat(SerendipitousSweep.trueSota(lk, 9, 9, 9))
				.as("⟨9,9,9⟩ SOTA must be a real value, not the dead-oracle -1").isPositive();
		assertThat(SerendipitousSweep.trueSota(lk, 18, 18, 18)).isPositive();
		assertThat(SerendipitousSweep.catalogRank(lk, 6, 8, 9)).isPositive();
	}

	@Test
	public void catalog_rank_reads_known_shapes() {
		// Content sanity (SOTA-or-better): Laderman ⟨3,3,3⟩=23 is the worst case.
		assertThat(SerendipitousSweep.catalogRank(lk, 3, 3, 3)).isBetween(1L, 23L);
		// ⟨2,2,2⟩ Strassen = 7.
		assertThat(SerendipitousSweep.catalogRank(lk, 2, 2, 2)).isBetween(1L, 7L);
	}

	@Test
	public void trueSota_is_the_min_of_catalog_and_plain_kron() {
		// ⟨18,18,18⟩ factorises (e.g. ⟨3,3,6⟩⊗⟨6,6,3⟩), so plain Kron is buildable
		// and trueSota never exceeds either handle. This is the ⟨18,18,18⟩=3200
		// false-positive guard: a serendipitous r_s only wins BELOW this min.
		long cat = SerendipitousSweep.catalogRank(lk, 18, 18, 18);
		long kron = SerendipitousSweep.kronBest(lk, 18, 18, 18);
		long sota = SerendipitousSweep.trueSota(lk, 18, 18, 18);
		assertThat(kron).as("⟨18,18,18⟩ must have a plain-Kron factorisation").isPositive();
		assertThat(sota).isEqualTo(Math.min(cat < 0 ? Long.MAX_VALUE : cat, kron));
		assertThat(sota).isLessThanOrEqualTo(kron);
	}

	@Test
	public void absent_shape_has_no_sota() {
		// Beyond the ≤32 catalog and unfactorable into known pieces → genuinely -1.
		assertThat(SerendipitousSweep.catalogRank(lk, 41, 41, 41)).isEqualTo(-1L);
		assertThat(SerendipitousSweep.trueSota(lk, 41, 41, 41)).isEqualTo(-1L);
	}
}
