package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.nio.file.Files;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;

/**
 * Guards the born-stamped stub fix: a matrix-less stub written by the
 * materialiser must carry {@code fields[]} (else {@link FieldAwareLookup}
 * treats it as ABSENT — the silent-invisibility bug that hid the freshly
 * materialised ⟨5,32,32⟩=3320 and several derived cube stubs). Also pins the
 * canonical field-inclusion expansion so the {@code StampFields} /
 * materialiser / lookup copies cannot drift apart again.
 */
public class TestBornStampedStub {

	/** The single canonical expander — drift here previously left ~8k integer
	 *  schemes invisible under the F2/F3 selectors. */
	@Test
	public void inclusion_field_names_are_canonical() {
		assertThat(FieldAwareLookup.inclusionFieldNames(Field.Z))
				.containsExactly("F2", "F3", "Z", "Q", "R", "C");
		assertThat(FieldAwareLookup.inclusionFieldNames(Field.Q)).containsExactly("Q", "R", "C");
		assertThat(FieldAwareLookup.inclusionFieldNames(Field.R)).containsExactly("R", "C");
		assertThat(FieldAwareLookup.inclusionFieldNames(Field.C)).containsExactly("C");
		assertThat(FieldAwareLookup.inclusionFieldNames(Field.F2)).containsExactly("F2");
		assertThat(FieldAwareLookup.inclusionFieldNames(Field.F3)).containsExactly("F3");
	}

	@Test
	public void writeStub_with_fields_emits_fields_array(@TempDir File dir) throws Exception {
		NonCubicBilinearAlgorithm one = new NonCubicBilinearAlgorithm(
				1, 1, 1, new double[][] { { 1 } }, new double[][] { { 1 } }, new double[][] { { 1 } });
		File f = new File(dir, "stub-with-fields.json");
		SchemeIO.writeStub(one, f, new Lineage.Atom("test-leaf"),
				List.of("F2", "F3", "Z", "Q", "R", "C"));
		String json = Files.readString(f.toPath());
		assertThat(json).as("stub must carry fields[] (else invisible to the lookup)")
				.contains("\"fields\"").contains("\"Z\"").contains("\"F2\"");
		// fields[] must be readable back as field tags, not dropped.
		assertThat(SchemeIO.fieldTags(SchemeIO.parseJson(f)))
				.contains("F2", "F3", "Z", "Q", "R", "C");
	}

	/**
	 * Field-set inference must be the set-INTERSECTION of leaves' CONTENT fields[]
	 * — never a single-field collapse or a [Z] fallback. The canonical bug:
	 * ⟨5,32,32⟩'s recombination uses the rational ⟨3,8,8⟩=145 ([F3,Q,R,C]) leaf, so
	 * the composed field is [F3,Q,R,C] — NOT [F2,F3,Z,Q,R,C]. Guards against Z/F2
	 * over-claim on rational-leaf schemes.
	 */
	@Test
	public void field_set_is_leaf_intersection_no_Z_overclaim() {
		FieldAwareLookup lk = new FieldAwareLookup(Field.Q);
		// rational leaf: no Z, no F2.
		assertThat(lk.bestFieldsAtShape(3, 8, 8))
				.contains("F3", "Q", "R", "C").doesNotContain("Z", "F2");
		// integer HK leaf: full chain.
		assertThat(lk.bestFieldsAtShape(2, 4, 4))
				.containsExactly("F2", "F3", "Z", "Q", "R", "C");
		// composition with a rational leaf → intersection drops Z/F2.
		assertThat(lk.fieldNamesFromLineage(
				new Lineage.ConcatCols(new Lineage.Atom("2x4x4"), new Lineage.Atom("3x8x8"))))
				.containsExactly("F3", "Q", "R", "C");
		// all-integer composition → full chain.
		assertThat(lk.fieldNamesFromLineage(
				new Lineage.ConcatCols(new Lineage.Atom("2x4x4"), new Lineage.Atom("2x8x8"))))
				.containsExactly("F2", "F3", "Z", "Q", "R", "C");
		// unresolvable leaf → FLOOR to the lookup field's inclusion chain (user
		// 2026-06-13: "always at least shrink into the requested field") — a Q-build
		// resolved that leaf to a Q-valid scheme, so the composition is at least
		// Q-valid. The load-bearing part of this guard is unchanged: NO over-claim —
		// never F2/F3/Z from an unresolved leaf (possible under-claim only).
		assertThat(lk.fieldNamesFromLineage(
				new Lineage.ConcatCols(new Lineage.Atom("2x4x4"), new Lineage.Atom("99x99x99"))))
				.containsExactly("Q", "R", "C");
	}

	@Test
	public void writeStub_without_fields_stays_backcompat(@TempDir File dir) throws Exception {
		NonCubicBilinearAlgorithm one = new NonCubicBilinearAlgorithm(
				1, 1, 1, new double[][] { { 1 } }, new double[][] { { 1 } }, new double[][] { { 1 } });
		File f = new File(dir, "stub-no-fields.json");
		SchemeIO.writeStub(one, f, new Lineage.Atom("test-leaf"));
		assertThat(Files.readString(f.toPath())).doesNotContain("\"fields\"");
	}
}
