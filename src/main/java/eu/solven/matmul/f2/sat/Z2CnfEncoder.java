package eu.solven.matmul.f2.sat;

import java.util.ArrayList;
import java.util.List;

import eu.solven.matmul.BilinearAlgorithm;

/**
 * CNF encoder for Z/2 bilinear matmul search.
 *
 * Searches for a rank-r decomposition {@code T = ⊕_k u_k ⊗ v_k ⊗ w_k} over GF(2)
 * (i.e., addition is XOR, multiplication is AND), where the target tensor T is
 * arbitrary (caller supplies it; typically dense or triangular matmul over Z/2).
 *
 * <p><b>Supports non-cubic targets</b>: the three factor matrices may have
 * different numbers of rows ({@code dimU, dimV, dimW}). For matmul
 * {@code ⟨n, m, p⟩}: {@code dimU = n·m, dimV = m·p, dimW = n·p}.</p>
 *
 * <p><b>Encoding scheme:</b></p>
 * <ul>
 *   <li><b>Scalars</b>: one boolean per {@code U[i, k]}, {@code V[i, k]},
 *       {@code W[i, k]}. Total {@code (dimU + dimV + dimW)·r} variables.</li>
 *   <li><b>Cubic products</b>: one boolean {@code P[a, b, c, k] =
 *       U[a, k] ∧ V[b, k] ∧ W[c, k]}. Tseitin clauses: {@code 4} per product.</li>
 *   <li><b>XOR-sum constraints</b>: per {@code (a, b, c)}, the XOR over k of
 *       {@code P[a, b, c, k]} must equal {@code target[a][b][c]}. Encoded as a
 *       chain of partial-sum auxiliary variables {@code S[a, b, c, k]} with
 *       {@code 4} XOR clauses per chain link plus a final unit clause.</li>
 * </ul>
 *
 * <p>Total clause count for {@code n=2, r=7}: ~3400 clauses, ~980 variables.
 * Tractable for SAT4J in milliseconds.</p>
 */
public class Z2CnfEncoder {

	private final int dimU;
	private final int dimV;
	private final int dimW;
	private final int r;
	private final int[][][] target;
	private final boolean columnLexOrdering;
	/** Per-slot zero-forcing masks ({@code [0] → U, [1] → V, [2] → W}). Null entries → no restriction for that slot. */
	private final boolean[][] forceZeroPositions;

	/** {@code scalarVar[slot][i][k]} → DIMACS var (1-indexed). {@code scalarVar[0]} has dimU rows, {@code [1]} has dimV, {@code [2]} has dimW. */
	private int[][][] scalarVar;
	private int[][][][] productVar; // [a][b][c][k]
	private int[][][][] sumVar; // [a][b][c][k] — partial XOR sum at index k

	private int nextVar = 1;
	private final List<int[]> clauses = new ArrayList<>();

	/** Cubic default: lex-ordering on columns enabled (necessary for SAT4J on hard ranks). */
	public Z2CnfEncoder(int n, int r, int[][][] target) {
		this(n * n, n * n, n * n, r, target, true, null);
	}

	public Z2CnfEncoder(int n, int r, int[][][] target, boolean columnLexOrdering) {
		this(n * n, n * n, n * n, r, target, columnLexOrdering, null);
	}

	/**
	 * Cubic with single per-position mask shared across U/V/W (legacy API).
	 *
	 * @param forceZeroPositions length-{@code n²} mask; positions where true have their
	 *                           U/V/W rows forced to 0 via unit clauses. Used to
	 *                           model restricted-input/output matmul (e.g.
	 *                           diagonal-only inputs, triangular, etc.). Pass
	 *                           {@code null} for no restriction.
	 */
	public Z2CnfEncoder(int n, int r, int[][][] target, boolean columnLexOrdering,
			boolean[] forceZeroPositions) {
		this(n * n, n * n, n * n, r, target, columnLexOrdering,
				forceZeroPositions == null ? null
						: new boolean[][] { forceZeroPositions, forceZeroPositions, forceZeroPositions });
	}

	/**
	 * Full non-cubic constructor. Slot dimensions must match the target tensor's
	 * shape: {@code target.length == dimU}, {@code target[0].length == dimV},
	 * {@code target[0][0].length == dimW}.
	 *
	 * @param forceZeroPositions per-slot zero-forcing masks ({@code [0]} → U,
	 *                           {@code [1]} → V, {@code [2]} → W); any element
	 *                           may be {@code null} to skip restrictions for
	 *                           that slot, and the outer array may also be
	 *                           {@code null} for no restrictions anywhere.
	 */
	public Z2CnfEncoder(int dimU, int dimV, int dimW, int r, int[][][] target,
			boolean columnLexOrdering, boolean[][] forceZeroPositions) {
		if (target.length != dimU || target[0].length != dimV || target[0][0].length != dimW) {
			throw new IllegalArgumentException(String.format(
					"target shape [%d][%d][%d] doesn't match dims [%d][%d][%d]",
					target.length, target[0].length, target[0][0].length, dimU, dimV, dimW));
		}
		this.dimU = dimU;
		this.dimV = dimV;
		this.dimW = dimW;
		this.r = r;
		this.target = target;
		this.columnLexOrdering = columnLexOrdering;
		this.forceZeroPositions = forceZeroPositions;
		encode();
	}

	private int dimOf(int slot) {
		switch (slot) {
		case 0:
			return dimU;
		case 1:
			return dimV;
		case 2:
			return dimW;
		default:
			throw new IllegalArgumentException("slot=" + slot);
		}
	}

	private int newVar() {
		return nextVar++;
	}

	private void addClause(int... literals) {
		clauses.add(literals);
	}

	private void encode() {
		allocateVars();
		encodeProductTseitin();
		encodeXorSums();
		if (forceZeroPositions != null) {
			encodeForcedZeros();
		}
		if (columnLexOrdering) {
			encodeColumnLexOrdering();
		}
	}

	/** Unit clauses forcing U/V/W rows for the masked positions to 0. */
	private void encodeForcedZeros() {
		for (int slot = 0; slot < 3; slot++) {
			if (forceZeroPositions[slot] == null) continue;
			boolean[] mask = forceZeroPositions[slot];
			int dim = dimOf(slot);
			for (int i = 0; i < dim && i < mask.length; i++) {
				if (!mask[i]) continue;
				for (int k = 0; k < r; k++) {
					addClause(-scalarVar[slot][i][k]);
				}
			}
		}
	}

	/**
	 * Symmetry-breaking: require columns of the {@code [U; V; W]} concatenation
	 * (length {@code dimU + dimV + dimW} per column) to be strictly lex-increasing
	 * in k. Breaks {@code S_r} permutation of rank-1 terms (factor {@code r!}
	 * reduction) and forbids two identical columns (which would represent a
	 * smaller-rank decomposition in disguise).
	 *
	 * Standard lex-leader encoding via per-bit "equal-so-far" auxiliary variables:
	 *   eq[0]   = TRUE
	 *   eq[i+1] ↔ eq[i] ∧ (a_i ↔ b_i)
	 *   forbid    eq[i] ∧ a_i ∧ ¬b_i  (would make a > b at position i)
	 *   force     ¬eq[N]              (strict inequality: a ≠ b)
	 */
	private void encodeColumnLexOrdering() {
		int colWidth = dimU + dimV + dimW;
		for (int k = 0; k < r - 1; k++) {
			int[] aBits = columnBits(k);
			int[] bBits = columnBits(k + 1);
			int[] eq = new int[colWidth + 1];
			eq[0] = newVar();
			addClause(eq[0]); // eq[0] = TRUE

			for (int i = 0; i < colWidth; i++) {
				eq[i + 1] = newVar();
				int eqI = eq[i];
				int eqNext = eq[i + 1];
				int aI = aBits[i];
				int bI = bBits[i];

				// eq[i+1] ↔ eq[i] ∧ (a_i ↔ b_i)
				// → direction (three implications collapsed):
				addClause(-eqNext, eqI);            // eq[i+1] → eq[i]
				addClause(-eqNext, -aI, bI);        // eq[i+1] → (a_i → b_i)
				addClause(-eqNext, aI, -bI);        // eq[i+1] → (b_i → a_i)
				// ← direction:
				addClause(-eqI, -aI, -bI, eqNext);  // eq[i] ∧ a_i ∧ b_i  → eq[i+1]
				addClause(-eqI, aI, bI, eqNext);    // eq[i] ∧ ¬a_i ∧ ¬b_i → eq[i+1]

				// a ≤ b: forbid (eq[i] ∧ a_i ∧ ¬b_i)
				addClause(-eqI, -aI, bI);
			}
			// strict: a ≠ b → ¬eq[N]
			addClause(-eq[colWidth]);
		}
	}

	private int[] columnBits(int k) {
		int[] bits = new int[dimU + dimV + dimW];
		int idx = 0;
		for (int slot = 0; slot < 3; slot++) {
			int dim = dimOf(slot);
			for (int i = 0; i < dim; i++) {
				bits[idx++] = scalarVar[slot][i][k];
			}
		}
		return bits;
	}

	private void allocateVars() {
		scalarVar = new int[3][][];
		scalarVar[0] = new int[dimU][r];
		scalarVar[1] = new int[dimV][r];
		scalarVar[2] = new int[dimW][r];
		for (int slot = 0; slot < 3; slot++) {
			int dim = dimOf(slot);
			for (int i = 0; i < dim; i++) {
				for (int k = 0; k < r; k++) {
					scalarVar[slot][i][k] = newVar();
				}
			}
		}
		productVar = new int[dimU][dimV][dimW][r];
		for (int a = 0; a < dimU; a++) {
			for (int b = 0; b < dimV; b++) {
				for (int c = 0; c < dimW; c++) {
					for (int k = 0; k < r; k++) {
						productVar[a][b][c][k] = newVar();
					}
				}
			}
		}
		sumVar = new int[dimU][dimV][dimW][r];
		for (int a = 0; a < dimU; a++) {
			for (int b = 0; b < dimV; b++) {
				for (int c = 0; c < dimW; c++) {
					for (int k = 0; k < r; k++) {
						sumVar[a][b][c][k] = newVar();
					}
				}
			}
		}
	}

	/**
	 * P[a, b, c, k] ↔ U[a, k] ∧ V[b, k] ∧ W[c, k]
	 *   forward:  P → U, P → V, P → W  (3 clauses)
	 *   backward: U ∧ V ∧ W → P        (1 clause)
	 */
	private void encodeProductTseitin() {
		for (int a = 0; a < dimU; a++) {
			for (int b = 0; b < dimV; b++) {
				for (int c = 0; c < dimW; c++) {
					for (int k = 0; k < r; k++) {
						int p = productVar[a][b][c][k];
						int u = scalarVar[0][a][k];
						int v = scalarVar[1][b][k];
						int w = scalarVar[2][c][k];
						addClause(-p, u);
						addClause(-p, v);
						addClause(-p, w);
						addClause(-u, -v, -w, p);
					}
				}
			}
		}
	}

	/**
	 * For each (a, b, c): T[a, b, c] = ⊕_k P[a, b, c, k].
	 *
	 * Chained encoding via auxiliary partial sums:
	 *   S[0]  = P[0]            (alias via 2 equivalence clauses)
	 *   S[k]  = S[k-1] ⊕ P[k]   (4 XOR clauses per step)
	 *   S[r-1] = target          (1 unit clause)
	 */
	private void encodeXorSums() {
		for (int a = 0; a < dimU; a++) {
			for (int b = 0; b < dimV; b++) {
				for (int c = 0; c < dimW; c++) {
					int s0 = sumVar[a][b][c][0];
					int p0 = productVar[a][b][c][0];
					// S[0] = P[0]
					addClause(-s0, p0);
					addClause(s0, -p0);

					// S[k] = S[k-1] XOR P[k]
					for (int k = 1; k < r; k++) {
						int sPrev = sumVar[a][b][c][k - 1];
						int pK = productVar[a][b][c][k];
						int sK = sumVar[a][b][c][k];
						encodeXor(sPrev, pK, sK);
					}

					// Final: S[r-1] = target
					int sFinal = sumVar[a][b][c][r - 1];
					if (target[a][b][c] == 1) {
						addClause(sFinal);
					} else {
						addClause(-sFinal);
					}
				}
			}
		}
	}

	/**
	 * x ↔ a ⊕ b, four CNF clauses (each forbids one row of the XOR truth table
	 * where the prediction doesn't match):
	 *   (0,0,1), (0,1,0), (1,0,0), (1,1,1) forbidden.
	 */
	private void encodeXor(int a, int b, int x) {
		addClause(a, b, -x);
		addClause(a, -b, x);
		addClause(-a, b, x);
		addClause(-a, -b, -x);
	}

	public List<int[]> getClauses() {
		return clauses;
	}

	public int getVarCount() {
		return nextVar - 1;
	}

	public int getDimU() {
		return dimU;
	}

	public int getDimV() {
		return dimV;
	}

	public int getDimW() {
		return dimW;
	}

	public int getRank() {
		return r;
	}

	/**
	 * Decode a SAT assignment back into a {@link BilinearAlgorithm}. Only valid
	 * when the target is cubic and {@code dimU = dimV = dimW = n²} for some
	 * integer n. For non-cubic targets call {@link #decodeRaw} instead.
	 */
	public BilinearAlgorithm decode(boolean[] assignment) {
		if (dimU != dimV || dimV != dimW) {
			throw new IllegalStateException(
					"Non-cubic encoder; use decodeRaw() and inspect U/V/W directly.");
		}
		int n = (int) Math.round(Math.sqrt(dimU));
		if (n * n != dimU) {
			throw new IllegalStateException("Cubic dim must be a perfect square, got " + dimU);
		}
		double[][][] uvw = decodeRaw(assignment);
		return new BilinearAlgorithm(n, uvw[0], uvw[1], uvw[2]);
	}

	/**
	 * Non-cubic-aware decode. Returns {@code {U, V, W}} as a 3-element array of
	 * {@code double[][]} with shapes {@code [dimU][r], [dimV][r], [dimW][r]}.
	 * Entries are 0.0 or 1.0 (the GF(2) values).
	 */
	public double[][][] decodeRaw(boolean[] assignment) {
		double[][] U = new double[dimU][r];
		double[][] V = new double[dimV][r];
		double[][] W = new double[dimW][r];
		for (int i = 0; i < dimU; i++) {
			for (int k = 0; k < r; k++) {
				U[i][k] = assignment[scalarVar[0][i][k] - 1] ? 1.0 : 0.0;
			}
		}
		for (int i = 0; i < dimV; i++) {
			for (int k = 0; k < r; k++) {
				V[i][k] = assignment[scalarVar[1][i][k] - 1] ? 1.0 : 0.0;
			}
		}
		for (int i = 0; i < dimW; i++) {
			for (int k = 0; k < r; k++) {
				W[i][k] = assignment[scalarVar[2][i][k] - 1] ? 1.0 : 0.0;
			}
		}
		return new double[][][] { U, V, W };
	}
}
