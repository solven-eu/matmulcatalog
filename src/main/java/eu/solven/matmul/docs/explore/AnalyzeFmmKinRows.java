package eu.solven.matmul.docs.explore;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.io.MapleSchemeParser;
import lombok.extern.slf4j.Slf4j;

/**
 * Test the "kin-row unification" hypothesis for FMM-Lille's ⟨27,28,28⟩.
 *
 * <p>FMM publishes a {@code _raw.mpl} scheme with 10442 products — the SAME count
 * as our coordinate projection of the ⟨28,28,28⟩=10556 cube — yet claims rank
 * 10413. The gap is not a stronger projection: it is products that became
 * <em>proportional</em> after projection (same rank-1 input form {@code (u,v)} up
 * to a scalar), which can be merged ("kin-row unification"). This probe parses the
 * raw scheme and counts DISTINCT rank-1 input forms; if that equals 10413, the
 * mechanism is product-merging, not linear projection.</p>
 *
 * <pre>mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.explore.AnalyzeFmmKinRows</pre>
 */
@Slf4j
public final class AnalyzeFmmKinRows {
	private AnalyzeFmmKinRows() {}

	private static final double EPS = 1e-9;

	public static void main(String[] args) throws Exception {
		File raw = new File("references/fmm-lille/27x28x28/27x28x28_raw.mpl");
		NonCubicBilinearAlgorithm a = MapleSchemeParser.parseRawFmmLille(raw, 27, 28, 28);
		int r = a.r, nm = 27 * 28, mp = 28 * 28;
		log.info("parsed FMM ⟨27,28,28⟩ raw: products(r)={}", r);

		// Distinct rank-1 INPUT forms (u,v) up to scalar: two products with the same
		// (u-direction, v-direction) are scalar multiples → the same multiplication →
		// kin, mergeable into one.
		// Two products merge into one (their rank-1 sum is rank-1) iff they share two
		// of three factor directions: (u,v) [proportional product], (u,w) [same left
		// factor + same output column], or (v,w) [same right factor + same column].
		Map<String, Integer> uv = new HashMap<>(), uw = new HashMap<>(), vw = new HashMap<>();
		for (int k = 0; k < r; k++) {
			String uk = dirKey(a.denseU(), nm, k), vk = dirKey(a.denseV(), mp, k), wk = dirKey(a.denseW(), nm, k);
			uv.merge(uk + "|" + vk, 1, Integer::sum);
			uw.merge(uk + "|" + wk, 1, Integer::sum);
			vw.merge(vk + "|" + wk, 1, Integer::sum);
		}
		log.info("merge by (u,v) [proportional]      : distinct={}  → {} merges", uv.size(), r - uv.size());
		log.info("merge by (u,w) [left + output kin] : distinct={}  → {} merges", uw.size(), r - uw.size());
		log.info("merge by (v,w) [right + output kin]: distinct={}  → {} merges", vw.size(), r - vw.size());
		// Group products by shared LEFT factor; each group's slice is u⊗S with
		// S = Σ_{k∈group} v_k⊗w_k. The slice needs only rank(S) products, so unifying
		// it saves (groupSize − rank(S)). Sum over all groups = the kin-row saving.
		Map<String, java.util.List<Integer>> leftGroups = new HashMap<>();
		for (int k = 0; k < r; k++)
			leftGroups.computeIfAbsent(dirKey(a.denseU(), nm, k), x -> new java.util.ArrayList<>()).add(k);
		Map<Integer, Integer> sizeHisto = new java.util.TreeMap<>();
		int totalSaved = 0;
		Map<Integer, Integer> savedBySize = new java.util.TreeMap<>();
		for (java.util.List<Integer> grp : leftGroups.values()) {
			sizeHisto.merge(grp.size(), 1, Integer::sum);
			if (grp.size() < 2) continue;
			int rk = sliceRank(a.denseV(), mp, a.denseW(), nm, grp);
			int saved = grp.size() - rk;
			totalSaved += saved;
			if (saved > 0) savedBySize.merge(grp.size(), saved, Integer::sum);
		}
		log.info("left-factor group sizes: {}", sizeHisto + " (size→#groups)");
		log.info("slice-rank savings by group size: {}", savedBySize + " (size→products saved)");
		log.info("SUMMARY: raw products={}  kin-row unify saves {}  →  {}  (FMM headline=10413, residual={})",
				r, totalSaved, r - totalSaved, (r - totalSaved) - 10413);
	}

	/** rank(Σ_{k∈grp} v_k w_kᵀ): orthonormalise the v_k (basis Qv, rv≤|grp|), then the
	 *  rank equals the #independent rows of (Qvᵀ V)·Wᵀ — i.e. of the rv combinations
	 *  Σ_k (q_b·v_k) w_k. All small (|grp| is 2–3 here). */
	private static int sliceRank(double[][] V, int mp, double[][] W, int nm, java.util.List<Integer> grp) {
		int g = grp.size();
		double[][] vs = new double[g][mp], ws = new double[g][nm];
		for (int t = 0; t < g; t++) {
			int k = grp.get(t);
			for (int x = 0; x < mp; x++) vs[t][x] = V[x][k];
			for (int x = 0; x < nm; x++) ws[t][x] = W[x][k];
		}
		// Orthonormal basis Qv of the v_t (Gram-Schmidt); keep the basis vectors.
		java.util.List<double[]> qv = gramSchmidt(vs);
		// Rows b: r_b = Σ_t (qv_b · v_t) w_t  (length nm). Their rank = rank(S).
		java.util.List<double[]> rows = new java.util.ArrayList<>();
		for (double[] q : qv) {
			double[] row = new double[nm];
			for (int t = 0; t < g; t++) {
				double dot = 0; for (int x = 0; x < mp; x++) dot += q[x] * vs[t][x];
				if (Math.abs(dot) > EPS) for (int x = 0; x < nm; x++) row[x] += dot * ws[t][x];
			}
			rows.add(row);
		}
		return gramSchmidt(rows.toArray(new double[0][])).size();
	}

	/** Return an orthonormal basis (list) of the given row-vectors; size = rank. */
	private static java.util.List<double[]> gramSchmidt(double[][] vecs) {
		java.util.List<double[]> basis = new java.util.ArrayList<>();
		for (double[] v0 : vecs) {
			double[] v = v0.clone();
			for (double[] b : basis) {
				double dot = 0; for (int x = 0; x < v.length; x++) dot += v[x] * b[x];
				for (int x = 0; x < v.length; x++) v[x] -= dot * b[x];
			}
			double nrm = 0; for (double x : v) nrm += x * x; nrm = Math.sqrt(nrm);
			if (nrm > 1e-6) { for (int x = 0; x < v.length; x++) v[x] /= nrm; basis.add(v); }
		}
		return basis;
	}

	/** Canonical direction key for column k of factor F (rows × r): divide by the
	 *  signed largest-magnitude entry, round, and stringify the support+ratios. */
	private static String dirKey(double[][] F, int rows, int k) {
		int piv = -1; double pivVal = 0;
		for (int i = 0; i < rows; i++) {
			double v = F[i][k];
			if (Math.abs(v) > Math.abs(pivVal)) { pivVal = v; piv = i; }
		}
		if (piv < 0) return "0"; // zero form (shouldn't happen)
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < rows; i++) {
			double v = F[i][k] / pivVal;
			if (Math.abs(v) > EPS) sb.append(i).append(':').append(Math.round(v * 1e6)).append(';');
		}
		return sb.toString();
	}
}
