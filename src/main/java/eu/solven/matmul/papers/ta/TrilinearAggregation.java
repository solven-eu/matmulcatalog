package eu.solven.matmul.papers.ta;

import java.util.Optional;
import java.util.OptionalLong;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * One member of the <b>trilinear-aggregation (TA) family</b> for the cubic
 * matmul tensor {@code ⟨N,N,N⟩}. Each member gives an <em>alternative</em> rank
 * (multiplication count) for the SAME shape, valid over its own
 * method-specific domain. The known members (see {@link TrilinearAggregations}):
 *
 * <ul>
 *   <li>{@code TA_pan}  — Pan 1980 (SIAM), even only.</li>
 *   <li>{@code TA_dis}  — Islam 2009 / DIS09 §3, even <b>and</b> odd. Buildable.</li>
 *   <li>{@code TA_hs}   — Pan / Hadas–Schwartz 1982, even, {@code n≠16}.</li>
 *   <li>{@code TA_sz}   — Schwartz–Zwecher 2025, even, {@code n≠16}.</li>
 *   <li>{@code TA_lita} — Khoruzhii–Gelß–Pokutta 2026, even and odd, {@code n≥19}. Buildable.</li>
 * </ul>
 *
 * <p>All are non-commutative and Q-rational (char-0; not F₂/F₃). The point of
 * the interface is to let callers gather every member's rank at a target N and
 * pick the best — "alternative ranks by accepting multiple muls" — without
 * hard-coding which formula wins where (the crossover differs by N and parity;
 * the n=16 pole and the n≥19 floor are encoded as empty results).</p>
 */
public interface TrilinearAggregation {

	/** Short tag, e.g. {@code "TA_dis"}, {@code "TA_lita"}. */
	String tag();

	/** Bibliographic reference. */
	String paperRef();

	/**
	 * Rank of {@code ⟨n,n,n⟩} via this method, or {@link OptionalLong#empty()}
	 * when the method does not apply at {@code n} (wrong parity, the {@code n=16}
	 * pole, below a domain floor, …). Never throws for {@code n ≥ 1}.
	 */
	OptionalLong cubicRank(int n);

	/**
	 * Whether {@link #build(int)} can materialise {@code ⟨n,n,n⟩} — a CHEAP
	 * capability check (no construction). Bound-only members return {@code false}.
	 */
	default boolean canBuild(int n) {
		return false;
	}

	/**
	 * Explicit factor-matrix construction for {@code ⟨n,n,n⟩}, or
	 * {@link Optional#empty()} for bound-only members (no constructor available).
	 * When present, its rank equals {@link #cubicRank(int)}. Prefer
	 * {@link #canBuild(int)} when you only need to test availability — building can
	 * be expensive (dense even-N LITA schemes).
	 */
	default Optional<NonCubicBilinearAlgorithm> build(int n) {
		return Optional.empty();
	}
}
