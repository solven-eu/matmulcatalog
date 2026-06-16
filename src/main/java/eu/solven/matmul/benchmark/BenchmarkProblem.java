package eu.solven.matmul.benchmark;

/**
 * A single decomposition problem: target tensor + rank to try. {@code id} is
 * used to name the per-problem artifact directory; keep it filesystem-safe
 * (alphanumerics, underscore, hyphen).
 *
 * <p>Includes optional per-slot zero-forcing masks ({@code forceZeroSlots})
 * for restricted-target problems (e.g. upper-triangular, diagonal-plus-one).
 * For unrestricted dense targets, pass {@code null}.</p>
 */
public final class BenchmarkProblem {

	public final String id;
	public final Alphabet alphabet;
	/** Matrix-mult format: {@code A} is {@code n×m}, {@code B} is {@code m×p}, {@code C} is {@code n×p}. */
	public final int n;
	public final int m;
	public final int p;
	/** Short label used in the MD roll-up ("dense", "upper-tri", "diag+1", ...). */
	public final String variantLabel;
	/** Shape {@code [n·m][m·p][n·p]}. */
	public final int[][][] target;
	/**
	 * Per-slot zero-forcing masks: {@code [0] → U}, {@code [1] → V}, {@code [2] → W}.
	 * Any entry may be {@code null} for no restriction on that slot; the outer
	 * array may also be {@code null} for no restrictions anywhere.
	 */
	public final boolean[][] forceZeroSlots;
	public final int rank;

	public BenchmarkProblem(String id, Alphabet alphabet, int n, int m, int p,
			String variantLabel, int[][][] target, boolean[][] forceZeroSlots, int rank) {
		this.id = id;
		this.alphabet = alphabet;
		this.n = n;
		this.m = m;
		this.p = p;
		this.variantLabel = variantLabel;
		this.target = target;
		this.forceZeroSlots = forceZeroSlots;
		this.rank = rank;
	}

	public int dimU() { return n * m; }
	public int dimV() { return m * p; }
	public int dimW() { return n * p; }

	public String formatLabel() {
		return String.format("⟨%d,%d,%d⟩", n, m, p);
	}
}
