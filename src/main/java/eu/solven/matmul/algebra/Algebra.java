package eu.solven.matmul.algebra;

/**
 * The arithmetic setting an algorithm operates in: a {@link Field}
 * plus a flag for whether the algorithm exploits commutativity of the
 * underlying scalars.
 *
 * <p>Commutativity is <strong>orthogonal</strong> to the field choice:</p>
 * <ul>
 *   <li>{@code (Z_Q_R, commutative=true)} — commutative-only scheme over
 *       integers/rationals/reals (e.g. Waksman 1970). Does NOT lift to
 *       recursive matmul over matrix entries.</li>
 *   <li>{@code (Z_Q_R, commutative=false)} — standard non-commutative
 *       bilinear scheme (e.g. Strassen 1969). Lifts to recursive matmul
 *       and is what produces practical fast matmul implementations.</li>
 *   <li>{@code (F2, commutative=true)} — commutative-only over F₂.
 *       F₂ scalars commute, but the SCHEME still exploits that for
 *       extra savings.</li>
 *   <li>{@code (F2, commutative=false)} — standard bilinear over F₂
 *       (e.g. AlphaTensor's {@code ⟨4,4,4⟩=47}).</li>
 * </ul>
 *
 * <p>Every rank claim in this catalog must specify both axes — see
 * {@code RANK_KNOWLEDGE.md} §1.2 for why.</p>
 */
public record Algebra(Field field, boolean commutative) {

	public Algebra {
		if (field == null) throw new IllegalArgumentException("field is required");
	}

	/** Shorthand for the standard non-commutative setting in a given field. */
	public static Algebra nonCommutative(Field field) {
		return new Algebra(field, false);
	}

	/** Shorthand for the commutative-only setting in a given field. */
	public static Algebra commutative(Field field) {
		return new Algebra(field, true);
	}

	/**
	 * Whether an algorithm valid under {@code other} is also valid under
	 * {@code this}. True iff (a) the source field is in our fallback
	 * chain, AND (b) the source isn't commutative-only when we want
	 * non-commutative-friendly behaviour.
	 *
	 * <p>Common case: a non-commutative {@code Z_Q_R} scheme is valid
	 * under any non-commutative request in {@code C, Z_Q_R} but NOT
	 * under {@code F2} (different characteristic).</p>
	 */
	public boolean accepts(Algebra other) {
		if (!field.fallbackChain().contains(other.field)) return false;
		// A commutative-only scheme can't be used where non-commutative is required.
		if (!this.commutative && other.commutative) return false;
		return true;
	}

	@Override
	public String toString() {
		return field.tag() + (commutative ? " (cmt)" : " (NC)");
	}
}
