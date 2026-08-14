package eu.solven.matmul.papers.schwartzzwecher2025;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.fraction.BigFraction;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Faithful port of Schwartz–Zwecher 2025 (arXiv:2508.01748) TA-New25 cubic
 * matmul, per Appendix B (explicit trilinear form) — see
 * {@code references/SZ_TA_NEW25_PORT_SPEC.md}. Even {@code n0 != 16}.
 *
 * <p>Rank {@code tNew = n0³/3 + 15/4 n0² + 61/6 n0 + 8}. The construction is the
 * φ-transform + two aggregation families (Ŝ, Ṡ\Ṡ1) + off-diagonal cancellation
 * (S̃\S̃1 via a fixed ⟨2,2,2;7⟩) + the united diagonal {@code R(i)}.</p>
 *
 * <p><b>VERIFIED (2026-08-06): build(n0) reconstructs the matmul tensor EXACTLY.</b>
 * n0=4/6/8 pass {@code Verifier.isExactNonCubic}; n0=10..32 pass the random
 * spot-check; rank == {@code tNew} at every tested n0 (n0=28/30/32 reproduce the
 * held dense SZ imports 10550/12688/15096). Wired into {@code TrilinearAggregations.SZ}.</p>
 *
 * <p><b>Construction</b> (all three matrices share the SAME transform):
 * {@code A*=B*=C*=φ(X)=(I2⊗L)·X·(I2⊗R)}, with {@code L=[I_half; −uᵀ]},
 * {@code R=[I_half−(1/d)J | −(1/d)u]}, {@code d=n0/2+1}. {@code Tr(φA·φB·φC)=Tr(ABC)}
 * because {@code R·L=I}. Four product families sum to {@code Tr(A*B*C*)}:
 * (a) Ŝ aggregation-symmetric, (b) Ṡ\Ṡ1 aggregation-barred, (c) off-diagonal
 * cancellation ({@code −d·}Σ of a 2×2 block trace via a fixed Strassen ⟨2,2,2;7⟩),
 * (d) the united diagonal {@code R(i)} (7 products/i). Factor recovery is the
 * φ-pullback shared by U/V/W: coeff of C[x][y] in a functional {@code Σγ_pq C*_pq}
 * is {@code Σ_pq γ_pq·PL[p][x]·PR[y][q]}, landing on output {@code (AB)_{yx}}
 * (trace dual). Coefficients are exact {@code BigFraction} (γ=1−9/d rational).</p>
 *
 * <p>The two bugs that cost the debug run (both found by diffing the arXiv LaTeX
 * e-print, not the garbled PDF): an asymmetric {@code B←L·B·Lᵀ} transform (must be
 * uniform φ), and a transpose in family (b)'s third C-term ({@code C*_{q,s}} not
 * {@code C*_{s,q}}). Harness: {@code docs.explore.ProbeTaNew25} (per-output-cell
 * residual). Spec: {@code references/SZ_TA_NEW25_PORT_SPEC.md}. Value is
 * provenance/coverage: LITA strictly beats SZ on rank.</p>
 */
public final class TaNew25Construction {
	private TaNew25Construction() {}

	public static long tNew(int n0) {
		// n0³/3 + 15/4 n0² + 61/6 n0 + 8, integer for even n0.
		BigFraction t = new BigFraction(n0).pow(3).divide(3)
				.add(new BigFraction(15, 4).multiply(new BigFraction(n0).pow(2)))
				.add(new BigFraction(61, 6).multiply(new BigFraction(n0)))
				.add(new BigFraction(8));
		return t.longValue();
	}

	public static NonCubicBilinearAlgorithm build(int n0) {
		if (n0 < 2 || (n0 & 1) == 1 || n0 == 16) {
			throw new IllegalArgumentException("TA-New25 needs even n0 != 16, got " + n0);
		}
		return new Builder(n0, false).build();
	}

	/**
	 * The construction's SYMBOLIC (transformed-space) forms, one per product:
	 * {@code out[product][axis]} with axis 0/1/2 = A*, B*, C*, each a length-{@code dim²}
	 * coefficient vector over the transformed variable {@code X*_{p,q}} (index
	 * {@code p·dim+q}, {@code dim=n0+2}). This is the space where SZ/LITA-style
	 * unification lives (kin = proportional on two axes) — the delivered factor
	 * matrices (see {@link #build}) hide it, because the φ-pullback densifies the
	 * forms into all-distinct directions. Same family emission as {@link #build};
	 * only Astar/Bstar are swapped from the φ-pullback to symbolic units and the
	 * pullback recovery is skipped. For structural analysis only ({@code ProbeTaKinGraph}).
	 */
	public static double[][][] buildSymbolicForms(int n0) {
		if (n0 < 2 || (n0 & 1) == 1 || n0 == 16) {
			throw new IllegalArgumentException("TA-New25 needs even n0 != 16, got " + n0);
		}
		Builder b = new Builder(n0, true);
		b.emitAll();
		int r = b.u.size();
		double[][][] out = new double[r][3][];
		for (int l = 0; l < r; l++) {
			out[l][0] = toDoubles(b.u.get(l));
			out[l][1] = toDoubles(b.v.get(l));
			out[l][2] = toDoubles(b.w.get(l));
		}
		return out;
	}

	private static double[] toDoubles(BigFraction[] f) {
		double[] d = new double[f.length];
		for (int i = 0; i < f.length; i++) d[i] = f[i].doubleValue();
		return d;
	}

	/** Symbolic transformed-space forms PLUS per-product class tags (0=generic,1=boundary,
	 *  2=correction,3=diagonal). For the boundary support-isolation analysis. */
	public record SymbolicTagged(double[][][] forms, int[] productClass) {}

	public static SymbolicTagged buildSymbolicTagged(int n0) {
		if (n0 < 2 || (n0 & 1) == 1 || n0 == 16) {
			throw new IllegalArgumentException("TA-New25 needs even n0 != 16, got " + n0);
		}
		Builder b = new Builder(n0, true);
		b.emitAll();
		int r = b.u.size();
		double[][][] out = new double[r][3][];
		int[] classes = new int[r];
		for (int l = 0; l < r; l++) {
			out[l][0] = toDoubles(b.u.get(l));
			out[l][1] = toDoubles(b.v.get(l));
			out[l][2] = toDoubles(b.w.get(l));
			classes[l] = b.cls.get(l);
		}
		return new SymbolicTagged(out, classes);
	}

	/** The built scheme plus each product's structural class: 0=generic aggregation,
	 *  1=boundary aggregation (2-distinct triples — the whole {@code 2·N²} tail),
	 *  2=off-diagonal correction, 3=united diagonal R(i). Aligned with product columns. */
	public record Tagged(NonCubicBilinearAlgorithm alg, int[] productClass) {}

	public static Tagged buildTagged(int n0) {
		if (n0 < 2 || (n0 & 1) == 1 || n0 == 16) {
			throw new IllegalArgumentException("TA-New25 needs even n0 != 16, got " + n0);
		}
		Builder b = new Builder(n0, false);
		NonCubicBilinearAlgorithm alg = b.build();
		int[] classes = new int[b.cls.size()];
		for (int i = 0; i < classes.length; i++) classes[i] = b.cls.get(i);
		return new Tagged(alg, classes);
	}

	private static final class Builder {
		final int n0, half, d, dim, n2, abLen;
		final BigFraction gamma, dF;
		// Astar/Bstar[p][q] = form (length n2) over A/B entries. Cstar[p][q] = padded-C
		// UNIT (length dim²) at p·dim+q. Output recovered in toAlg via C_out=(I2⊗Lᵀ)C*(I2⊗L).
		final BigFraction[][][] Astar, Bstar, Cstar;
		BigFraction[][] fPL, fPLt, fPR;   // I2⊗L, I2⊗Lᵀ, I2⊗R — for the C-recovery
		final List<BigFraction[]> u = new ArrayList<>(), v = new ArrayList<>(), w = new ArrayList<>();
		// Per-product structural class (for boundary-reducibility analysis):
		// 0=generic aggregation (3-distinct), 1=boundary aggregation (2-distinct, the N² tail),
		// 2=off-diagonal correction (c), 3=united diagonal R(i). Tracked only when tagging.
		int curClass = 0;
		final List<Integer> cls = new ArrayList<>();

		Builder(int n0, boolean symbolic) {
			this.n0 = n0;
			this.half = n0 / 2;
			this.d = half + 1;
			this.dim = 2 * d;         // = n0 + 2
			this.n2 = n0 * n0;
			this.abLen = symbolic ? dim * dim : n2;   // A/B form length: symbolic units vs φ-pullback
			this.dF = new BigFraction(d);
			this.gamma = BigFraction.ONE.subtract(new BigFraction(9, d));  // 1 − 9/d
			// φ(X) = (I2⊗L)·X·(I2⊗R). PL: dim×n0, PR: n0×dim.
			BigFraction[][] L = new BigFraction[d][half];      // [I_{half}; −uᵀ]
			for (int r = 0; r < d; r++) for (int c = 0; c < half; c++) {
				L[r][c] = (r < half) ? (r == c ? BigFraction.ONE : BigFraction.ZERO)
						: BigFraction.MINUS_ONE;
			}
			BigFraction[][] R = new BigFraction[half][d];      // [I_{half} − (1/d)J | −(1/d)u]
			BigFraction invD = new BigFraction(1, d);
			for (int a = 0; a < half; a++) for (int c = 0; c < d; c++) {
				R[a][c] = (c < half) ? (a == c ? BigFraction.ONE : BigFraction.ZERO).subtract(invD)
						: invD.negate();
			}
			BigFraction[][] PL = new BigFraction[dim][n0];
			BigFraction[][] PR = new BigFraction[n0][dim];      // I2⊗R  (for A's right)
			BigFraction[][] PLt = new BigFraction[n0][dim];     // I2⊗Lᵀ (for B's right)
			for (BigFraction[] row : PL) java.util.Arrays.fill(row, BigFraction.ZERO);
			for (BigFraction[] row : PR) java.util.Arrays.fill(row, BigFraction.ZERO);
			for (BigFraction[] row : PLt) java.util.Arrays.fill(row, BigFraction.ZERO);
			for (int b = 0; b < 2; b++) {
				for (int r = 0; r < d; r++) for (int c = 0; c < half; c++) PL[b * d + r][b * half + c] = L[r][c];
				for (int a = 0; a < half; a++) for (int c = 0; c < d; c++) PR[b * half + a][b * d + c] = R[a][c];
				for (int a = 0; a < half; a++) for (int c = 0; c < d; c++) PLt[b * half + a][b * d + c] = L[c][a]; // Lᵀ
			}
			// A ← L·A·R ; B ← L·B·Lᵀ.  Store PL/PR for the C-recovery (C*=φ(C)).
			fPL = PL; fPLt = PLt; fPR = PR;
			// UNIFORM φ for all three (paper §"Explicit Description"): A*=B*=C*=φ(X)=(I2⊗L)X(I2⊗R).
			// Tr(φA·φB·φC)=Tr(ABC) since R·L=I. (An asymmetric B←L·B·Lᵀ is WRONG here.)
			// Symbolic mode keeps A*/B* as transformed-space UNITS (no pullback) for structural analysis.
			Astar = symbolic ? unitStar() : star(PL, PR);
			Bstar = symbolic ? unitStar() : star(PL, PR);
			// C* is the PADDED output variable: Cstar[p][q] = unit at (p·dim+q) in the dim²
			// padded space. Output C_out=(I2⊗Lᵀ)·C*·(I2⊗L) is applied in toAlg.
			Cstar = unitStar();
		}

		/** star[p][q] = unit at (p·dim+q) in the dim²-dim transformed space (X*_{p,q} itself). */
		BigFraction[][][] unitStar() {
			BigFraction[][][] s = new BigFraction[dim][dim][];
			for (int p = 0; p < dim; p++) for (int q = 0; q < dim; q++) {
				s[p][q] = new BigFraction[dim * dim];
				java.util.Arrays.fill(s[p][q], BigFraction.ZERO);
				s[p][q][p * dim + q] = BigFraction.ONE;
			}
			return s;
		}

		/** star[p][q][a*n0+b] = PL[p][a]·PR[b][q]. */
		BigFraction[][][] star(BigFraction[][] PL, BigFraction[][] PR) {
			BigFraction[][][] s = new BigFraction[dim][dim][];
			for (int p = 0; p < dim; p++) for (int q = 0; q < dim; q++) {
				BigFraction[] form = new BigFraction[n2];
				java.util.Arrays.fill(form, BigFraction.ZERO);
				for (int a = 0; a < n0; a++) {
					if (PL[p][a].equals(BigFraction.ZERO)) continue;
					for (int b = 0; b < n0; b++) {
						form[a * n0 + b] = PL[p][a].multiply(PR[b][q]);
					}
				}
				s[p][q] = form;
			}
			return s;
		}

		int bar(int i) { return (i + d) % dim; }   // ī = i + n0/2+1 (mod n0+2)

		/** Accumulator for one linear form (A/B: length n2 over inputs; C: length dim² padded). */
		final class Form {
			final BigFraction[] c;
			Form(int len) { c = new BigFraction[len]; java.util.Arrays.fill(c, BigFraction.ZERO); }
			Form add(BigFraction[][][] star, int p, int q, BigFraction s) {
				BigFraction[] f = star[p][q];
				for (int t = 0; t < c.length; t++) if (!f[t].equals(BigFraction.ZERO)) c[t] = c[t].add(s.multiply(f[t]));
				return this;
			}
			Form add(BigFraction[][][] star, int p, int q) { return add(star, p, q, BigFraction.ONE); }
		}

		void emit(Form a, Form b, Form c) { u.add(a.c); v.add(b.c); w.add(c.c); cls.add(curClass); }

		int distinct(int i, int j, int k) {
			int c = 1;
			if (j != i) c++;
			if (k != i && k != j) c++;
			return c;
		}

		/** (a) aggregation-symmetric product for member (p,q,s). */
		void familyA(int p, int q, int s) {
			emit(new Form(abLen).add(Astar, p, q).add(Astar, q, s).add(Astar, s, p),
				new Form(abLen).add(Bstar, q, s).add(Bstar, s, p).add(Bstar, p, q),
				new Form(dim * dim).add(Cstar, s, p).add(Cstar, p, q).add(Cstar, q, s));
		}

		/** (b) aggregation-barred product for member (p,q,s); bar() is an involution. */
		void familyB(int p, int q, int s) {
			BigFraction N1 = BigFraction.MINUS_ONE;
			emit(new Form(abLen).add(Astar, p, q, N1).add(Astar, bar(q), s).add(Astar, s, bar(p)),
				new Form(abLen).add(Bstar, q, bar(s)).add(Bstar, s, p).add(Bstar, bar(p), q),
				new Form(dim * dim).add(Cstar, bar(s), p, N1).add(Cstar, p, bar(q)).add(Cstar, q, s));
		}

		NonCubicBilinearAlgorithm build() {
			emitAll();
			return toAlg();
		}

		/** Emit all four families into u/v/w (works in both normal and symbolic mode). */
		void emitAll() {
			BigFraction G = gamma, ONE = BigFraction.ONE, N1 = BigFraction.MINUS_ONE;
			BigFraction Ginv = ONE.divide(G), G2inv = Ginv.multiply(Ginv);
			// (a) Ŝ = {(i,j,k),(ī,j̄,k̄) : i≤j<k or k<j≤i} — BOTH members per qualifying base.
			for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) for (int k = 0; k < d; k++) {
				if (!((i <= j && j < k) || (k < j && j <= i))) continue;
				curClass = distinct(i, j, k) == 2 ? 1 : 0;   // 2-distinct → boundary (N² tail)
				familyA(i, j, k);
				familyA(bar(i), bar(j), bar(k));
			}
			// (b) Ṡ\Ṡ1: Ṡ = {(i,j,k),(ī,j̄,k̄)}, Ṡ1 = {(i,i,i)} unbarred only.
			for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) for (int k = 0; k < d; k++) {
				curClass = distinct(i, j, k) == 2 ? 1 : 0;
				if (!(i == j && j == k)) familyB(i, j, k);         // unbarred, minus the 3 diagonals
				familyB(bar(i), bar(j), bar(k));                   // barred partner, always
			}
			// (c) off-diagonal cancellation: −d·Σ_{i≠j} Trace of the 2×2 block via a fixed ⟨2,2,2;7⟩ (Strassen).
			curClass = 2;
			for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) {
				if (i == j) continue;
				strassenTrace(i, j, dF.negate());
			}
			// (d) united diagonal R(i)
			curClass = 3;
			for (int i = 0; i < d; i++) {
				int ib = bar(i);
				BigFraction dGam = dF.multiply(G), dGinv = dF.multiply(Ginv), dG2inv = dF.multiply(G2inv);
				// 1
				emit(new Form(abLen).add(Astar, ib, i).add(Astar, i, ib).add(Astar, i, i, N1),
					new Form(abLen).add(Bstar, ib, i).add(Bstar, i, ib).add(Bstar, i, i),
					new Form(dim * dim).add(Cstar, ib, ib, dF.multiply(ONE.subtract(G)).divide(G))
						.add(Cstar, ib, i, G.subtract(dF).divide(G).negate())
						.add(Cstar, i, ib, G.negate().add(dF).divide(G).negate())
						.add(Cstar, i, i, ONE.subtract(dF)));
				// 2
				emit(new Form(abLen).add(Astar, i, ib),
					new Form(abLen).add(Bstar, ib, ib, G.negate().subtract(ONE).divide(G))
						.add(Bstar, ib, i, Ginv.negate())
						.add(Bstar, i, ib, ONE.subtract(G2inv))
						.add(Bstar, i, i, G.subtract(ONE).divide(G)),
					new Form(dim * dim).add(Cstar, ib, ib, dF).add(Cstar, ib, i, dF).add(Cstar, i, ib, dGinv).add(Cstar, i, i, dF));
				// 3
				emit(new Form(abLen).add(Astar, i, ib).add(Astar, i, i, G),
					new Form(abLen).add(Bstar, ib, ib, G.add(ONE).divide(G)).add(Bstar, ib, i, G.add(ONE).divide(G))
						.add(Bstar, i, ib, G2inv).add(Bstar, i, i, Ginv),
					new Form(dim * dim).add(Cstar, i, ib, dGinv).add(Cstar, i, i, dF));
				// 4
				emit(new Form(abLen).add(Astar, ib, i).add(Astar, i, i, G.negate().subtract(ONE)),
					new Form(abLen).add(Bstar, ib, ib).add(Bstar, ib, i).add(Bstar, i, ib, Ginv).add(Bstar, i, i),
					new Form(dim * dim).add(Cstar, i, ib, dG2inv).add(Cstar, i, i, dF.add(dGinv)));
				// 5
				emit(new Form(abLen).add(Astar, ib, ib).add(Astar, ib, i).add(Astar, i, ib, Ginv.negate()).add(Astar, i, i, N1),
					new Form(abLen).add(Bstar, ib, ib, G.negate().subtract(ONE)).add(Bstar, i, ib, Ginv.negate()),
					new Form(dim * dim).add(Cstar, ib, ib, dF.multiply(G.subtract(ONE)).divide(G)).add(Cstar, ib, i, dGinv.negate()));
				// 6
				emit(new Form(abLen).add(Astar, ib, i).add(Astar, i, i, N1),
					new Form(abLen).add(Bstar, ib, ib, G.negate().subtract(ONE)).add(Bstar, ib, i, N1)
						.add(Bstar, i, ib, G.negate().subtract(ONE).divide(G)).add(Bstar, i, i, N1),
					new Form(dim * dim).add(Cstar, ib, ib, dF.multiply(ONE.subtract(G)).divide(G)).add(Cstar, ib, i, dGinv)
						.add(Cstar, i, ib, dF.multiply(G.subtract(ONE)).divide(G.multiply(G)).negate()).add(Cstar, i, i, dGinv));
				// 7
				emit(new Form(abLen).add(Astar, ib, ib).add(Astar, i, ib, G.negate().subtract(ONE).divide(G)),
					new Form(abLen).add(Bstar, ib, ib, N1).add(Bstar, i, ib, G.subtract(ONE).divide(G)),
					new Form(dim * dim).add(Cstar, ib, ib, dGinv).add(Cstar, ib, i, dF.negate().subtract(dGinv).negate()));
			}
		}

		/**
		 * {@code sc·Trace(Ablk·Bblk·Cblk)} for the 2×2 star blocks (page 19), via
		 * Strassen ⟨2,2,2;7⟩. Blocks: a11=(i,j) a12=(ī,j) a21=(i,j̄) a22=(ī,j̄)
		 * (col↔bar i, row↔bar j); same for B. Cblk=[[Cij,−Cīj],[−Cij̄,Cīj̄]].
		 * {@code Trace((AB)·Cblk)=Σ (AB)ₓᵧ·Cblkᵧₓ} ⟹ per-Mₖ W-form derived below.
		 */
		void strassenTrace(int i, int j, BigFraction sc) {
			int ib = bar(i), jb = bar(j);
			BigFraction N = BigFraction.MINUS_ONE;
			// M1=(a11+a22)(b11+b22), W=Cij+Cīj̄
			emit(new Form(abLen).add(Astar, i, j).add(Astar, ib, jb), new Form(abLen).add(Bstar, i, j).add(Bstar, ib, jb),
				new Form(dim * dim).add(Cstar, i, j, sc).add(Cstar, ib, jb, sc));
			// M2=(a21+a22)b11, W=−Cīj−Cīj̄
			emit(new Form(abLen).add(Astar, i, jb).add(Astar, ib, jb), new Form(abLen).add(Bstar, i, j),
				new Form(dim * dim).add(Cstar, ib, j, sc.negate()).add(Cstar, ib, jb, sc.negate()));
			// M3=a11(b12−b22), W=−Cij̄+Cīj̄
			emit(new Form(abLen).add(Astar, i, j), new Form(abLen).add(Bstar, ib, j).add(Bstar, ib, jb, N),
				new Form(dim * dim).add(Cstar, i, jb, sc.negate()).add(Cstar, ib, jb, sc));
			// M4=a22(b21−b11), W=Cij−Cīj
			emit(new Form(abLen).add(Astar, ib, jb), new Form(abLen).add(Bstar, i, jb).add(Bstar, i, j, N),
				new Form(dim * dim).add(Cstar, i, j, sc).add(Cstar, ib, j, sc.negate()));
			// M5=(a11+a12)b22, W=−Cij−Cij̄
			emit(new Form(abLen).add(Astar, i, j).add(Astar, ib, j), new Form(abLen).add(Bstar, ib, jb),
				new Form(dim * dim).add(Cstar, i, j, sc.negate()).add(Cstar, i, jb, sc.negate()));
			// M6=(a21−a11)(b11+b12), W=Cīj̄
			emit(new Form(abLen).add(Astar, i, jb).add(Astar, i, j, N), new Form(abLen).add(Bstar, i, j).add(Bstar, ib, j),
				new Form(dim * dim).add(Cstar, ib, jb, sc));
			// M7=(a12−a22)(b21+b22), W=Cij
			emit(new Form(abLen).add(Astar, ib, j).add(Astar, ib, jb, N), new Form(abLen).add(Bstar, i, jb).add(Bstar, ib, jb),
				new Form(dim * dim).add(Cstar, i, j, sc));
		}

		NonCubicBilinearAlgorithm toAlg() {
			int r = u.size();
			double[][] U = new double[n2][r], V = new double[n2][r], W = new double[n2][r];
			for (int col = 0; col < r; col++) {
				for (int t = 0; t < n2; t++) {
					U[t][col] = u.get(col)[t].doubleValue();
					V[t][col] = v.get(col)[t].doubleValue();
				}
				// Consistent φ-pullback (same as U/V's `star`): coeff of C[x][y] in a
				// C*-functional Σ γ_pq C*_pq is Σ_pq γ_pq·PL[p][x]·PR[y][q]; that contributes
				// to output (AB)_{yx} (trace dual), so W[y·n0+x] = Σ γ_pq PL[p][x] PR[y][q].
				BigFraction[] wPad = w.get(col);
				for (int y = 0; y < n0; y++) for (int x = 0; x < n0; x++) {
					BigFraction acc = BigFraction.ZERO;
					for (int p = 0; p < dim; p++) {
						if (fPL[p][x].equals(BigFraction.ZERO)) continue;
						for (int q = 0; q < dim; q++) {
							BigFraction wv = wPad[p * dim + q];
							if (wv.equals(BigFraction.ZERO)) continue;
							acc = acc.add(fPL[p][x].multiply(fPR[y][q]).multiply(wv));
						}
					}
					W[y * n0 + x][col] = acc.doubleValue();
				}
			}
			return new NonCubicBilinearAlgorithm(n0, n0, n0, U, V, W);
		}
	}
}
