package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.SerendipitousBudProduct.BudDecomposition;
import eu.solven.matmul.catalog.SerendipitousBudProduct.GridBud;
import eu.solven.matmul.catalog.SerendipitousBudProduct.InnerResolver;
import eu.solven.matmul.verifiers.Verifier;

/**
 * Combined (two-axis) buds in the serendipitous product — Perminov §2.6.4.
 *
 * <p>The silent gap this guards (2026-09-30): the engine fused single-type buds only,
 * so every published {@code s1 ⊗ˢ s2} rank whose base carries a {@code ⟨1,2,2⟩} /
 * {@code ⟨2,1,2⟩} / {@code ⟨2,2,1⟩} grid was unreachable — 41 of the 45 formats where
 * the catalog still sat above Perminov's serendipitous 17–32 band after the importer
 * fixes. No error anywhere: the search just priced ⟨9,9,18⟩ at 920 instead of 913.</p>
 */
public class TestSerendipitousGridBud {

	private static FieldAwareLookup lookup;

	@BeforeAll
	static void buildLookup() {
		lookup = new FieldAwareLookup("Q");
	}

	/** Naive for everything except the one fusion target under test, taken from the catalog. */
	private static InnerResolver naiveExcept(int n, int m, int p) {
		return (a, b, c) -> (a == n && b == m && c == p)
				? lookup.findWithSource(a, b, c).map(FieldAwareLookup.WithSource::alg)
				: Optional.of(NonCubicBilinearAlgorithm.naive(a, b, c));
	}

	/**
	 * The naive ⟨1,2,2⟩ / ⟨2,1,2⟩ / ⟨2,2,1⟩ scheme IS one combined bud of four terms.
	 * Against the unit-axis inner that completes it to ⟨2,2,2⟩ the four rank-2 copies
	 * (8 products; 8 as two single-type buds too) fuse into Strassen's 7 — once per
	 * grid type, each through a different branch of the block construction.
	 */
	@Test
	public void each_grid_type_fuses_into_strassen() {
		int[][] cases = { { 1, 2, 2, 2, 1, 1 }, { 2, 1, 2, 1, 2, 1 }, { 2, 2, 1, 1, 1, 2 } };
		for (int[] c : cases) {
			NonCubicBilinearAlgorithm base = NonCubicBilinearAlgorithm.naive(c[0], c[1], c[2]);
			List<GridBud> grids = SerendipitousBudProduct.candidateDecompositions(base).stream()
					.flatMap(d -> d.grids().stream()).toList();
			assertThat(grids).as("naive ⟨%d,%d,%d⟩ is a single combined bud", c[0], c[1], c[2])
					.anySatisfy(g -> {
						assertThat(new int[] { g.a(), g.b(), g.c() }).containsExactly(c[0], c[1], c[2]);
						assertThat(g.terms()).hasSize(4);
					});

			NonCubicBilinearAlgorithm out = SerendipitousBudProduct.productViaBudsBest(
					base, naiveExcept(2, 2, 2), c[3], c[4], c[5]);

			assertThat(new int[] { out.n, out.m, out.p }).containsExactly(2, 2, 2);
			assertThat(out.r).as("⟨%d,%d,%d⟩ ⊗ˢ ⟨%d,%d,%d⟩", c[0], c[1], c[2], c[3], c[4], c[5]).isEqualTo(7);
			assertThat(Verifier.isExactNonCubic(out)).isTrue();
		}
	}

	/** A ⟨1,2,3⟩ grid (six terms) with a non-unit inner: ⟨1,2,3⟩ ⊗ˢ ⟨2,1,1⟩ → one ⟨2,2,3⟩. */
	@Test
	public void a_rectangular_grid_fuses_into_one_block() {
		NonCubicBilinearAlgorithm base = NonCubicBilinearAlgorithm.naive(1, 2, 3);
		int fused = lookup.findRank(2, 2, 3);
		assertThat(fused).as("the catalog's ⟨2,2,3⟩ beats 6 naive ⟨2,1,1⟩ copies").isLessThan(12);

		NonCubicBilinearAlgorithm out = SerendipitousBudProduct.productViaBudsBest(base, naiveExcept(2, 2, 3), 2, 1, 1);

		assertThat(out.r).isEqualTo(fused);
		assertThat(Verifier.isExactNonCubic(out)).isTrue();
	}

	/**
	 * Proportional — not equal — class members: rescale one term's U by 2 and its V by
	 * 1/2 (the base stays exact). The free factor must absorb the class scales.
	 */
	@Test
	public void class_scales_are_absorbed_by_the_free_factor() {
		NonCubicBilinearAlgorithm naive = NonCubicBilinearAlgorithm.naive(1, 2, 2);
		double[][] u = copy(naive.denseU()), v = copy(naive.denseV()), w = copy(naive.denseW());
		for (double[] row : u) row[3] *= 2;
		for (double[] row : v) row[3] /= 2;
		for (double[] row : u) row[1] *= -1;
		for (double[] row : w) row[1] *= -1;
		NonCubicBilinearAlgorithm base = new NonCubicBilinearAlgorithm(1, 2, 2, u, v, w);
		assertThat(Verifier.isExactNonCubic(base)).isTrue();

		NonCubicBilinearAlgorithm out = SerendipitousBudProduct.productViaBudsBest(base, naiveExcept(2, 2, 2), 2, 1, 1);

		assertThat(out.r).isEqualTo(7);
		assertThat(Verifier.isExactNonCubic(out)).isTrue();
	}

	/**
	 * A base with no grid: the family is the six type orderings (the pre-2026-09-30
	 * family, first and unchanged — so nothing that used to be reachable is lost) plus
	 * their six fewest-options-first variants.
	 */
	@Test
	public void without_a_grid_the_family_starts_with_the_six_orderings() {
		NonCubicBilinearAlgorithm strassen = lookup.findWithSource(2, 2, 2).orElseThrow().alg();
		List<BudDecomposition> family = SerendipitousBudProduct.candidateDecompositions(strassen);
		int orderings = SerendipitousBudProduct.ALL_ORDERINGS.length;
		assertThat(family).hasSize(2 * orderings);
		assertThat(family).allSatisfy(d -> assertThat(d.grids()).isEmpty());
		for (int i = 0; i < orderings; i++) {
			BudDecomposition legacy = SerendipitousBudProduct.findBuds(strassen, SerendipitousBudProduct.ALL_ORDERINGS[i]);
			assertThat(family.get(i).trivial()).containsExactly(legacy.trivial());
			assertThat(family.get(i).buds()).hasSameSizeAs(legacy.buds());
		}
	}

	/**
	 * The catalog-level regression: ⟨9,9,18⟩ from a ⟨3,3,3⟩:23 whose structure is
	 * {@code ⟨1,2,2⟩ + 2⟨2,1,1⟩ + 15⟨1,1,1⟩}, against ⟨6,3,3⟩:40 —
	 * {@code 15·40 + 2·80 + ⟨6,6,6⟩:153 = 913} (Perminov, arXiv:2606.02480). Single-type
	 * buds alone give 920. Assert ≤ (a better ⟨6,6,6⟩ must not break this).
	 */
	@Test
	public void compute_pipeline_reaches_9x9x18_913_via_a_combined_bud() throws Exception {
		List<NonCubicBilinearAlgorithm> bases = new ArrayList<>();
		for (var path : lookup.findFiles(3, 3, 3)) {
			try {
				NonCubicBilinearAlgorithm a = SchemeIO.read(path.toFile());
				if (a.r == 23) bases.add(a);
			} catch (Exception ignore) {
				// stub / non-bilinear — not a reaction base
			}
		}
		assertThat(bases).as("⟨3,3,3⟩:23 representatives on disk").isNotEmpty();

		Optional<SerendipitousSearch.Hit> hit = Optional.empty();
		for (int[] inner : new int[][] { { 3, 3, 6 }, { 3, 6, 3 }, { 6, 3, 3 } }) {
			Optional<SerendipitousSearch.Hit> h = SerendipitousSearch.bestFor(
					3 * inner[0], 3 * inner[1], 3 * inner[2], bases, lookup, 920);
			if (h.isPresent() && (hit.isEmpty() || h.get().rank() < hit.get().rank())) hit = h;
		}

		assertThat(hit).as("a combined bud must beat the single-type 920").isPresent();
		assertThat(hit.get().rank()).isLessThanOrEqualTo(913);
		assertThat(Verifier.isExactNonCubic(hit.get().scheme())).as("⟨9,9,18⟩ result must verify").isTrue();
	}

	/**
	 * One published recipe per mechanism of the structural family, priced (not built)
	 * from the bases on disk — each sat 2–8 above Perminov's rank until its mechanism
	 * existed, with every type ordering and the largest-first grids already in place:
	 * <ul>
	 *   <li>⟨16,15,27⟩ ≤ 3744 — fewest-options-first: a U-triple gives up its third term
	 *       to that term's only W-partner ({@code 18⟨1,1,2⟩ + 3⟨1,2,1⟩}, not 17 + 1 triple + 2);</li>
	 *   <li>⟨6,15,27⟩ ≤ 1475 — same, where the partner must be served before the triple;</li>
	 *   <li>⟨9,30,28⟩ ≤ 4244 — {@code KEEP_U}: 69 grids and four whole U-triples beat 71 grids;</li>
	 *   <li>⟨4,25,26⟩ ≤ 1657 — {@code CORE}: the 2×2 core of a 2×5 grid, over the two
	 *       columns that break no U-triple.</li>
	 * </ul>
	 */
	@Test
	public void the_structural_family_prices_perminovs_recipes() throws Exception {
		int[][] rows = {
				// base n,m,p, base rank, inner n2,m2,p2, published rank
				{ 4, 5, 9, 132, 4, 3, 3, 3744 },
				{ 2, 5, 9, 73, 3, 3, 3, 1475 },
				{ 3, 10, 14, 312, 3, 3, 2, 4244 },
				{ 2, 5, 13, 104, 2, 5, 2, 1657 } };
		for (int[] r : rows) {
			long best = Long.MAX_VALUE;
			for (var path : lookup.findFiles(r[0], r[1], r[2])) {
				NonCubicBilinearAlgorithm base;
				try {
					base = SchemeIO.read(path.toFile());
				} catch (Exception e) {
					continue; // stub — not a reaction base
				}
				if (base.r != r[3] || base.n != r[0] || base.m != r[1] || base.p != r[2]) continue;
				best = Math.min(best, SerendipitousBudProduct.serendipitousCost(base, lookup, r[4], r[5], r[6]));
			}
			assertThat(best).as("⟨%d,%d,%d⟩:%d ⊗ˢ ⟨%d,%d,%d⟩", r[0], r[1], r[2], r[3], r[4], r[5], r[6])
					.isLessThanOrEqualTo(r[7]);
		}
	}

	private static double[][] copy(double[][] m) {
		double[][] out = new double[m.length][];
		for (int i = 0; i < m.length; i++) out[i] = m[i].clone();
		return out;
	}
}
