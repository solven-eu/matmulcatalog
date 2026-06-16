package eu.solven.matmul;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * A matrix-multiplication shape {@code ⟨n,m,p⟩}: multiply an {@code n×m} matrix
 * by an {@code m×p} matrix. Immutable value type — the single replacement for the
 * {@code (int n, int m, int p)} triples and {@code int[]} shapes scattered across
 * the codebase.
 *
 * <p><b>Order is significant.</b> {@code ⟨n,m,p⟩} is NOT {@code ⟨p,m,n⟩}: the axis
 * order is the orientation, which composition (Kronecker/concat) and lineage
 * depend on. The matmul <em>tensor</em> is S₃-symmetric, so for shape-class keys
 * (rank is orientation-invariant) use {@link #canonical()} / {@link CanonicalShape},
 * not the raw shape.</p>
 *
 * <p>Algebra context (field, commutativity) is tracked separately — a {@code Shape}
 * is purely the three dimensions (see the field-discipline rule in CLAUDE.md).</p>
 */
public record Shape(int n, int m, int p) implements Comparable<Shape> {

	public Shape {
		if (n < 1 || m < 1 || p < 1) {
			throw new IllegalArgumentException("Shape dims must be ≥ 1, got ⟨" + n + "," + m + "," + p + "⟩");
		}
	}

	public static Shape of(int n, int m, int p) {
		return new Shape(n, m, p);
	}

	/** Largest axis — the band a scheme is filed under. */
	public int maxDim() {
		return Math.max(n, Math.max(m, p));
	}

	/** Smallest axis. */
	public int minDim() {
		return Math.min(n, Math.min(m, p));
	}

	/** {@code n·m·p} (long to avoid overflow at the 32³ end). */
	public long volume() {
		return (long) n * m * p;
	}

	public boolean isCubic() {
		return n == m && m == p;
	}

	/** True if any axis is 1 (a trivial / naïve-degenerate shape). */
	public boolean hasWidthOneAxis() {
		return n == 1 || m == 1 || p == 1;
	}

	/** The sorted ({@code n≤m≤p}) shape-class key — rank is invariant under the S₃
	 *  axis action, so this is the right key for "best rank at this shape". */
	public CanonicalShape canonical() {
		int[] s = { n, m, p };
		java.util.Arrays.sort(s);
		return new CanonicalShape(s[0], s[1], s[2]);
	}

	/** True if {@code this} divides {@code target} on every axis (Kronecker cofactor exists). */
	public boolean divides(Shape target) {
		return target.n % n == 0 && target.m % m == 0 && target.p % p == 0;
	}

	/** The Kronecker cofactor {@code target / this} (requires {@link #divides}). */
	public Shape cofactor(Shape target) {
		if (!divides(target)) {
			throw new IllegalArgumentException(this + " does not divide " + target);
		}
		return new Shape(target.n / n, target.m / m, target.p / p);
	}

	public int[] toArray() {
		return new int[] { n, m, p };
	}

	public static Shape ofArray(int[] a) {
		if (a == null || a.length != 3) {
			throw new IllegalArgumentException("expected int[3], got " + java.util.Arrays.toString(a));
		}
		return new Shape(a[0], a[1], a[2]);
	}

	private static final Pattern PARSE =
			Pattern.compile("\\D*?(\\d+)\\s*[x,]\\s*(\\d+)\\s*[x,]\\s*(\\d+)");

	/** Parse {@code "3x5x7"}, {@code "3,5,7"}, or {@code "⟨3,5,7⟩"} (first triple found). */
	public static Shape parse(String s) {
		Matcher mm = PARSE.matcher(s);
		if (!mm.find()) {
			throw new IllegalArgumentException("no ⟨n,m,p⟩ in '" + s + "'");
		}
		return new Shape(Integer.parseInt(mm.group(1)), Integer.parseInt(mm.group(2)), Integer.parseInt(mm.group(3)));
	}

	/** Compact {@code "3x5x7"} form — the filename / lineage token convention. */
	public String toFilenameToken() {
		return n + "x" + m + "x" + p;
	}

	@Override
	public String toString() {
		return "⟨" + n + "," + m + "," + p + "⟩";
	}

	/** Order by max axis, then by volume, then lexicographically — the
	 *  "rank by format" sort used in the catalog. */
	@Override
	public int compareTo(Shape o) {
		int c = Integer.compare(maxDim(), o.maxDim());
		if (c != 0) return c;
		c = Long.compare(volume(), o.volume());
		if (c != 0) return c;
		c = Integer.compare(n, o.n);
		if (c != 0) return c;
		c = Integer.compare(m, o.m);
		if (c != 0) return c;
		return Integer.compare(p, o.p);
	}
}
