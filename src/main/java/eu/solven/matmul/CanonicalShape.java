package eu.solven.matmul;

/**
 * A {@link Shape} guaranteed to be in canonical (sorted) order {@code n ≤ m ≤ p}.
 * Distinct type so APIs that key by shape-CLASS — where rank is invariant under
 * the S₃ axis action — can demand canonicality at compile time rather than
 * re-sorting defensively. Obtain one via {@link Shape#canonical()}.
 *
 * <p>Use this for catalog "best rank at this shape" maps and cross-orientation
 * dedup; use a plain {@link Shape} wherever orientation matters (composition,
 * lineage, on-disk factor layout).</p>
 */
public record CanonicalShape(int n, int m, int p) {

	public CanonicalShape {
		if (n < 1 || n > m || m > p) {
			throw new IllegalArgumentException(
					"CanonicalShape requires 1 ≤ n ≤ m ≤ p, got ⟨" + n + "," + m + "," + p + "⟩");
		}
	}

	/** This canonical key as a (sorted) {@link Shape}. */
	public Shape shape() {
		return Shape.of(n, m, p);
	}

	public int maxDim() {
		return p;  // sorted ⇒ p is the max
	}

	public long volume() {
		return (long) n * m * p;
	}

	public boolean isCubic() {
		return n == p;  // sorted ⇒ n==p ⟹ all equal
	}

	public String toFilenameToken() {
		return n + "x" + m + "x" + p;
	}

	@Override
	public String toString() {
		return "⟨" + n + "," + m + "," + p + "⟩";
	}
}
