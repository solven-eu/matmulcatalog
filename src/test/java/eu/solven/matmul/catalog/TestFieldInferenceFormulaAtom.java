package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

/**
 * Regression guard: {@link FieldAwareLookup#fieldNamesFromLineage} must resolve a
 * PARAMETRIC FORMULA leaf — a ref carrying no {@code NxMxP} shape, e.g.
 * {@code DIS09Lemma4(n=20)} (the Pan/Islam trilinear-aggregation cube). Before the
 * fix the formula atom contributed no shape → the leaf set was empty → the
 * write-time field stamp threw "field inference returned EMPTY" and the shape
 * (⟨19,20,20⟩ = Project(DIS09Lemma4(n=20))) was dropped.
 *
 * <p>The cube is Q-strict (rational, ÷(n+1); valid over Q/R/C, NOT F₂/Z), so its
 * intrinsic field is {@code [Q,R,C]} — taken from the formula identity, NOT a
 * ⟨n,n,n⟩ catalog lookup (which could over-claim integer for a rational leaf).
 */
public class TestFieldInferenceFormulaAtom {

	@Test
	public void dis09Lemma4_projection_stamps_QRC_not_empty() {
		FieldAwareLookup lk = new FieldAwareLookup("Q");
		// ⟨19,20,20⟩ = projection of the ⟨20,20,20⟩ Pan-TA cube (drop one row).
		int[] keepN = range(19), keepM = range(20), keepP = range(20);
		Lineage.Node lineage = new Lineage.Project(
				new Lineage.Atom("DIS09Lemma4(n=20)"), keepN, keepM, keepP);

		assertThat(lk.fieldNamesFromLineage(lineage))
				.containsExactly("Q", "R", "C");
	}

	@Test
	public void bare_dis09Lemma4_atom_is_QRC() {
		FieldAwareLookup lk = new FieldAwareLookup("Q");
		assertThat(lk.fieldNamesFromLineage(new Lineage.Atom("DIS09Lemma4(n=26)")))
				.containsExactly("Q", "R", "C");
	}

	@Test
	public void naive_formula_atom_is_integer_all_fields() {
		FieldAwareLookup lk = new FieldAwareLookup("Q");
		// A naive elementary-product leaf is integer → the full field chain.
		assertThat(lk.fieldNamesFromLineage(new Lineage.Atom("naive-1x1x1")))
				.containsExactly("F2", "F3", "Z", "Q", "R", "C");
	}

	private static int[] range(int k) {
		int[] a = new int[k];
		for (int i = 0; i < k; i++) a[i] = i;
		return a;
	}
}
