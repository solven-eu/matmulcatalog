package eu.solven.matmul.papers.ta;

import java.util.List;
import java.util.Optional;
import java.util.OptionalLong;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.papers.dis2009.PanTrilinearAggregation;
import eu.solven.matmul.papers.khoruzhii2026.LitaTaConstruction;
import eu.solven.matmul.papers.khoruzhii2026.LitaTrilinearAggregation;
import eu.solven.matmul.papers.schwartzzwecher2025.TaNew25Construction;

/**
 * The trilinear-aggregation (TA) family registry: every {@link TrilinearAggregation}
 * member, plus helpers to take the BEST (min) rank across members at a cubic
 * {@code ⟨n,n,n⟩} target — the "alternative ranks by accepting multiple muls"
 * abstraction.
 *
 * <p>Members (each a thin adapter over its paper's existing implementation):</p>
 * <ul>
 *   <li>{@link #PAN} — Pan 1980 (SIAM), even only, bound-only.</li>
 *   <li>{@link #DIS} — Islam 2009 / DIS09 §3, even+odd, BUILDABLE.</li>
 *   <li>{@link #HS}  — Pan / Hadas–Schwartz 1982, even ({@code n≠16}), bound-only.</li>
 *   <li>{@link #SZ}  — Schwartz–Zwecher 2025, even ({@code n≠16}), BUILDABLE.</li>
 *   <li>{@link #LITA}— Khoruzhii–Gelß–Pokutta 2026, even+odd ({@code n≥19}), BUILDABLE.</li>
 * </ul>
 *
 * <p>{@code DIS} is the member historically (mis)named "Pan" — its formulas are
 * Islam/DIS09, not Pan's own (whose genuine even-only bound is {@code PAN}). See
 * {@link PanTrilinearAggregation} for the underlying formulas + attribution.</p>
 */
public enum TrilinearAggregations implements TrilinearAggregation {

	/** Pan 1980 (SIAM): {@code (2n³+9n²−6n)/4}, even only. Bound-only. */
	PAN("TA_pan", "Pan 1980 SIAM J. Comput. 9(2)") {
		@Override
		public OptionalLong cubicRank(int n) {
			return wrap(n < 2 ? -1L : PanTrilinearAggregation.panSiam1980Bound(n));
		}
	},

	/** Islam 2009 / DIS09 §3: even {@code (n³+12n²+11n)/3}, odd {@code (n³+15n²+14n−6)/3}. Buildable. */
	DIS("TA_dis", "Islam 2009 MSc / DIS09 §3") {
		@Override
		public OptionalLong cubicRank(int n) {
			return n < 2 ? OptionalLong.empty() : OptionalLong.of(PanTrilinearAggregation.cubicBound(n));
		}

		@Override
		public boolean canBuild(int n) {
			return n >= 2;
		}

		@Override
		public Optional<NonCubicBilinearAlgorithm> build(int n) {
			return n < 2 ? Optional.empty() : Optional.of(PanTrilinearAggregation.build(n));
		}
	},

	/** Pan / Hadas–Schwartz 1982: {@code (4n³+45n²+128n+108)/12}, even, {@code n≠16}. Bound-only. */
	HS("TA_hs", "Pan / Hadas–Schwartz 1982") {
		@Override
		public OptionalLong cubicRank(int n) {
			return wrap(n < 2 ? -1L : PanTrilinearAggregation.panHadasSchwartz1982Bound(n));
		}
	},

	/** Schwartz–Zwecher 2025 (arXiv:2508.01748 Thm 3.4): even, {@code n≠16}. BUILDABLE
	 *  via {@link TaNew25Construction} (Appendix B port; n0=4/6/8 exact-verified, larger
	 *  even n0 spot-checked; rank == the bound). */
	SZ("TA_sz", "Schwartz–Zwecher 2025 (arXiv:2508.01748)") {
		@Override
		public OptionalLong cubicRank(int n) {
			return wrap(n < 2 ? -1L : PanTrilinearAggregation.schwartzZwecher2025Bound(n));
		}

		@Override
		public boolean canBuild(int n) {
			return n >= 2 && n % 2 == 0 && n != 16;
		}

		@Override
		public Optional<NonCubicBilinearAlgorithm> build(int n) {
			return canBuild(n) ? Optional.of(TaNew25Construction.build(n)) : Optional.empty();
		}
	},

	/** Khoruzhii–Gelß–Pokutta 2026 (LITA): even+odd, {@code n≥19}. Buildable (when the
	 *  construction port lands; bound-only until then). */
	LITA("TA_lita", "Khoruzhii–Gelß–Pokutta 2026 (github.com/khoruzhii/lita)") {
		@Override
		public OptionalLong cubicRank(int n) {
			return n < LitaTrilinearAggregation.MIN_N
					? OptionalLong.empty()
					: OptionalLong.of(LitaTrilinearAggregation.cubicRank(n));
		}

		@Override
		public boolean canBuild(int n) {
			return n >= LitaTrilinearAggregation.MIN_N;
		}

		@Override
		public Optional<NonCubicBilinearAlgorithm> build(int n) {
			// Faithful Maple→Java port (LitaTaConstruction). Odd N exact-verify
			// (sparse); even N are dense and verify via the repo's spot-check path.
			return n < LitaTrilinearAggregation.MIN_N
					? Optional.empty()
					: Optional.of(LitaTaConstruction.build(n));
		}
	};

	private final String tag;
	private final String paperRef;

	TrilinearAggregations(String tag, String paperRef) {
		this.tag = tag;
		this.paperRef = paperRef;
	}

	@Override
	public String tag() {
		return tag;
	}

	@Override
	public String paperRef() {
		return paperRef;
	}

	/** {@code -1} sentinel (the paper-bound convention for "not applicable") → empty. */
	private static OptionalLong wrap(long boundOrMinusOne) {
		return boundOrMinusOne < 0 ? OptionalLong.empty() : OptionalLong.of(boundOrMinusOne);
	}

	/** All members, in attribution order. */
	public static List<TrilinearAggregation> all() {
		return List.of(values());
	}

	/** The genuine Pan-FAMILY members (excludes LITA) — backs {@code bestPanTaBound}. */
	private static final List<TrilinearAggregation> PAN_FAMILY = List.of(PAN, DIS, HS, SZ);

	/** A member's rank and identity at a target N. */
	public record Ranked(TrilinearAggregation method, long rank) {}

	/** Best (min) rank over the given members at {@code ⟨n,n,n⟩}, with the winning
	 *  member, or empty when none applies. Ties resolve to earliest in {@code members}. */
	public static Optional<Ranked> best(int n, List<TrilinearAggregation> members) {
		Ranked best = null;
		for (TrilinearAggregation m : members) {
			OptionalLong r = m.cubicRank(n);
			if (r.isEmpty()) {
				continue;
			}
			if (best == null || r.getAsLong() < best.rank()) {
				best = new Ranked(m, r.getAsLong());
			}
		}
		return Optional.ofNullable(best);
	}

	/** Best rank over ALL members (incl. LITA) at {@code ⟨n,n,n⟩}. */
	public static Optional<Ranked> best(int n) {
		return best(n, all());
	}

	/** Best (min) rank over all members, or empty when none applies. */
	public static OptionalLong bestRank(int n) {
		return best(n).map(r -> OptionalLong.of(r.rank())).orElse(OptionalLong.empty());
	}

	/** Best rank over the genuine Pan family (PAN/DIS/HS/SZ; excludes LITA) — the
	 *  semantics of the historical {@code bestPanTaBound}. */
	public static OptionalLong bestPanFamilyRank(int n) {
		return best(n, PAN_FAMILY).map(r -> OptionalLong.of(r.rank())).orElse(OptionalLong.empty());
	}

	/** Best member that can also {@link TrilinearAggregation#build} at this N (DIS, or LITA
	 *  once ported), with its rank — for the materialiser (never offer an unbuildable rank). */
	public static Optional<Ranked> bestBuildable(int n) {
		Ranked best = null;
		for (TrilinearAggregation m : all()) {
			if (!m.canBuild(n)) {
				continue;
			}
			OptionalLong r = m.cubicRank(n);
			if (r.isPresent() && (best == null || r.getAsLong() < best.rank())) {
				best = new Ranked(m, r.getAsLong());
			}
		}
		return Optional.ofNullable(best);
	}
}
