package eu.solven.matmul.commutative;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.List;
import java.util.Optional;

import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.recombination.Recombination;
import org.junit.jupiter.api.Test;

/**
 * Sanity guard on the commutative axis: a commutative bound must never be
 * WORSE than a non-commutative rank we already own.
 *
 * <p>Every NC scheme is valid over a commutative ring — commutativity is an
 * extra freedom, never a restriction — so {@code R_c ≤ R} is a theorem. A
 * "best known commutative rank" that exceeds our best known NC rank is
 * therefore a reporting bug, not a fact about the algebra.</p>
 *
 * <p>{@link CommutativeBounds} used to consult only the commutative formulas
 * (Rosowski 2019, DIS09 Table 4) and so under-reported on 209 of the 5456
 * shapes we hold over {@code Q} — smallest ⟨9,18,18⟩ (1624 from the formulas
 * vs 1600 already on disk), widest ⟨32,32,32⟩ (17392 vs 15079). Since this
 * class is the SOTA resolver behind commutative recombination
 * ({@code GenerateDerivedBounds.commutativeRecombineEntries}), the inflated
 * baseline made recombination "improvements" look bigger than they were.</p>
 */
public class TestCommutativeNotWorseThanNc {

	/** Field Q per CLAUDE.md's "sweeps default to --field=Q" (admits Z + Q ingredients). */
	private final FieldAwareLookup nc = new FieldAwareLookup(Field.Q);

	private final CommutativeBounds bounds = new CommutativeBounds();

	/** Shapes chosen to exercise distinct mechanisms: the classic tight case, the
	 *  cubic band where Rosowski wins, and the band where the NC catalog wins. */
	private static final List<int[]> SHAPES = List.of(
			new int[] { 2, 2, 2 }, new int[] { 3, 3, 3 }, new int[] { 4, 4, 4 },
			new int[] { 5, 5, 5 }, new int[] { 8, 8, 8 }, new int[] { 12, 12, 12 },
			new int[] { 16, 16, 16 }, new int[] { 26, 26, 26 }, new int[] { 32, 32, 32 },
			new int[] { 3, 4, 5 }, new int[] { 6, 8, 10 }, new int[] { 9, 18, 18 });

	@Test
	public void commutative_bound_never_exceeds_a_known_nc_rank() {
		for (int[] s : SHAPES) {
			int ncRank = nc.findRank(s[0], s[1], s[2]);
			if (ncRank >= Recombination.SotaResolver.UNKNOWN_RANK) continue;
			Optional<Long> cmt = bounds.bestRank(s[0], s[1], s[2]);
			assertThat(cmt).as("commutative bound known for ⟨%d,%d,%d⟩", s[0], s[1], s[2]).isPresent();
			assertThat(cmt.get())
					.as("⟨%d,%d,%d⟩: commutative bound must not exceed the NC rank %d "
							+ "(every NC scheme is valid commutatively)", s[0], s[1], s[2], ncRank)
					.isLessThanOrEqualTo((long) ncRank);
		}
	}

	@Test
	public void catalog_floor_closes_the_formula_gap() {
		// Documents the defect this floor closes: at these shapes the commutative
		// formulas alone are strictly worse than schemes we already own. Stated as
		// inequalities, so a future formula/catalog improvement can't break them.
		CommutativeBounds formulasOnly = new CommutativeBounds(false);
		for (int[] s : List.of(new int[] { 9, 18, 18 }, new int[] { 32, 32, 32 })) {
			long formula = formulasOnly.bestRank(s[0], s[1], s[2]).orElseThrow();
			long floored = bounds.bestRank(s[0], s[1], s[2]).orElseThrow();
			assertThat(floored).as("⟨%d,%d,%d⟩ floored below formula-only", s[0], s[1], s[2])
					.isLessThan(formula);
			assertThat(floored).isLessThanOrEqualTo(nc.findRank(s[0], s[1], s[2]));
		}
	}

	@Test
	public void commutativity_buys_nothing_at_2x2x2() {
		// Winograd 1971: R_c(⟨2,2,2⟩) = R(⟨2,2,2⟩) = 7. The catalog claimed 6 here
		// until 2026-08 — see TestBoundsVsLowerBounds for the guard on that.
		assertThat(bounds.bestRank(2, 2, 2)).contains(7L);
	}
}
