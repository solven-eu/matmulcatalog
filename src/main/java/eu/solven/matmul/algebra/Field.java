package eu.solven.matmul.algebra;

import java.util.List;

/**
 * The arithmetic field a scheme's factor matrices live in.
 *
 * <p>Splitting Z / Q / R is mathematically correct — they're distinct
 * — but in practice most non-commutative matmul algorithms are
 * Z-coefficient and trivially valid over Q and R. The
 * {@link #fallbackChain()} encodes that an algorithm valid over a
 * smaller field is also valid over any larger one (downward
 * transfer).</p>
 *
 * <ul>
 *   <li>{@link #Z} — integers only (typically {±1, 0}, sometimes
 *       {-2..2}). Lifts to Q, R, C and reduces mod p to F_p.</li>
 *   <li>{@link #Q} — rationals; may include {@code 1/2}, {@code 1/3}.
 *       Lifts to R, C. Does NOT generally lift to F_p (e.g., {@code 1/2}
 *       is undefined in F₂).</li>
 *   <li>{@link #R} — real numbers; may include irrationals
 *       ({@code √2}, {@code π}). Lifts to C only.</li>
 *   <li>{@link #C} — complex numbers; the broadest characteristic-0
 *       option (e.g., AlphaEvolve's ⟨4,4,4⟩=48).</li>
 *   <li>{@link #F2} — GF(2), characteristic 2.</li>
 *   <li>{@link #F3} — GF(3), characteristic 3.</li>
 * </ul>
 *
 * <p>Combine with the commutativity flag via {@link Algebra}.</p>
 */
public enum Field {

	/** Integers. Lifts to Q, R, C; reduces mod 2 to F₂. */
	Z("Z"),

	/** Rationals (may use {@code 1/2}, {@code 1/3}). Lifts to R, C; does NOT lift to F_p. */
	Q("Q"),

	/** Reals (may include irrationals). Lifts to C only. */
	R("R"),

	/** Complex numbers. The broadest characteristic-0 field tracked here. */
	C("C"),

	/** GF(2). */
	F2("F2"),

	/** GF(3). */
	F3("F3");

	private final String tag;

	Field(String tag) {
		this.tag = tag;
	}

	/** Canonical display tag (e.g. {@code "Z"}, {@code "F2"}). */
	public String tag() {
		return tag;
	}

	/**
	 * Returns the chain of fields a scheme can come from to be valid
	 * over THIS field, in priority order (most-restrictive / native first).
	 *
	 * <ul>
	 *   <li>{@link #Z} → only {@link #Z}</li>
	 *   <li>{@link #Q} → {@link #Q}, {@link #Z}</li>
	 *   <li>{@link #R} → {@link #R}, {@link #Q}, {@link #Z}</li>
	 *   <li>{@link #C} → {@link #C}, {@link #R}, {@link #Q}, {@link #Z}</li>
	 *   <li>{@link #F2} → only {@link #F2}</li>
	 *   <li>{@link #F3} → only {@link #F3}</li>
	 * </ul>
	 *
	 * <p>Note: this captures DIRECT containment. Mod-p reduction
	 * (e.g. Z → F2) is a different transfer and not encoded here —
	 * it would require coefficient-level validation that the scheme's
	 * matrices don't use coefficients that vanish mod p.</p>
	 */
	public List<Field> fallbackChain() {
		return switch (this) {
			case Z -> List.of(Z);
			case Q -> List.of(Q, Z);
			case R -> List.of(R, Q, Z);
			case C -> List.of(C, R, Q, Z);
			case F2 -> List.of(F2);
			case F3 -> List.of(F3);
		};
	}

	/**
	 * Parse a tag string ({@code "Z"}, {@code "Q"}, {@code "R"},
	 * {@code "C"}, {@code "F2"}, {@code "Z2"}, {@code "F3"},
	 * {@code "Z3"}, {@code "ZT"}, {@code "R/Q/Z"}) back to a Field.
	 * {@code "ZT"} (Perminov's "ternary" Z-target — integer coefficients
	 * restricted to {-1,0,1}) parses to {@link #Z}: ZT is a SUB-CLASS of Z,
	 * not a field of its own, and definitely NOT F₂/Z₂. Historical
	 * {@code "R/Q/Z"} (single-bucket) parses to {@link #R} (widest of the three).
	 */
	public static Field fromTag(String tag) {
		if (tag == null) throw new IllegalArgumentException("null tag");
		switch (tag) {
			case "F2":
			case "Z2":
				return F2;
			case "F3":
			case "Z3":
				return F3;
			case "Z":
			case "ZT":  // ZT = integer scheme with ternary {-1,0,1} coefficients ⊂ Z
				return Z;
			case "Q":
				return Q;
			case "R":
			case "R/Q/Z":  // historical single-bucket → widest characteristic-0 non-C
				return R;
			case "C":
				return C;
			default:
				throw new IllegalArgumentException("unknown field tag: " + tag);
		}
	}
}
