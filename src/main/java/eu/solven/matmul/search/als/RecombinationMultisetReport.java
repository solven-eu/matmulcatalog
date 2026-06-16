package eu.solven.matmul.search.als;

import eu.solven.matmul.catalog.Recombination;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.search.RecombinationMultisetOrbit;

/**
 * Parameterised report: <em>given a target shape and a base shape</em>, tabulate
 * every recombination multiset the base can realise (the exact GL-orbit, via
 * {@link RecombinationMultisetOrbit}) and the number of multiplications each
 * yields when recombined to the target.
 *
 * <p><b>The allocation is independent per axis (n / m / p).</b> Each target axis
 * is split into base-many parts; a multiset's rank is minimised over the
 * <em>joint</em> per-axis split space (not the symmetric diagonal). The report
 * shows that {@code best} plus the winning {@code (n / m / p)} split, and — when
 * the target+base are cubic — the symmetric splits {@code a+b} as reference
 * columns (the per-axis optimum can beat them).</p>
 *
 * <p>For a ⟨2,2,2⟩ base the 6 <b>ternary</b> multisets are flagged + cited by
 * curated/known base hash (constructable); the rest need a rational change-of-basis
 * (predictions only). Run: {@code mvn -o exec:java -Dexec.mainClass=…RecombinationMultisetReport
 * -Dexec.args="--target=7x7x7 --base=2x2x2"}.</p>
 */
public final class RecombinationMultisetReport {
	private RecombinationMultisetReport() {}

	private static final int[][] PERMS3 = {
			{ 0, 1, 2 }, { 0, 2, 1 }, { 1, 0, 2 }, { 1, 2, 0 }, { 2, 0, 1 }, { 2, 1, 0 } };
	private static final String KNOWN2 = "src/main/resources/schemes/known/section2/";
	private static final String CURATED2 = "src/main/resources/schemes/curated/section2/";

	public static void main(String[] args) throws Exception {
		int[] target = parseShape(arg(args, "--target", "7x7x7"));
		int[] baseShape = parseShape(arg(args, "--base", "2x2x2"));
		FieldAwareLookup lk = new FieldAwareLookup(Field.Q);
		NonCubicBilinearAlgorithm base = resolveIntegerBase(baseShape, lk);
		if (base == null) {
			System.out.printf("No integer base scheme found at ⟨%d,%d,%d⟩%n", baseShape[0], baseShape[1], baseShape[2]);
			return;
		}
		String md = report(target, base, lk);
		Path out = Path.of("references", "recomb-multisets-" + show(target) + "-base-" + show(baseShape) + ".md");
		Files.createDirectories(out.getParent());
		Files.writeString(out, md);
		System.out.println("Wrote " + out);
	}

	/** The report markdown for {@code target} recombined from {@code base}. */
	public static String report(int[] target, NonCubicBilinearAlgorithm base, FieldAwareLookup lk) throws Exception {
		int[] baseShape = { base.n, base.m, base.p };
		RecombinationMultisetOrbit.Result orbit = RecombinationMultisetOrbit.enumerate(base, 3);

		// Per-axis split lists — INDEPENDENT per axis. Descending parts so block-index
		// 0 = larger block. The recombination rank is minimised over the joint product.
		List<int[]> splitsN = descCompositions(target[0], baseShape[0]);
		List<int[]> splitsM = descCompositions(target[1], baseShape[1]);
		List<int[]> splitsP = descCompositions(target[2], baseShape[2]);
		splitsN.sort(Comparator.comparingInt(c -> c[0]));
		splitsM.sort(Comparator.comparingInt(c -> c[0]));
		splitsP.sort(Comparator.comparingInt(c -> c[0]));

		// Symmetric reference columns (same split on all axes) — only when cubic.
		boolean cubic = target[0] == target[1] && target[1] == target[2]
				&& baseShape[0] == baseShape[1] && baseShape[1] == baseShape[2];
		List<int[]> symSplits = cubic ? splitsN : List.of();

		boolean is222 = baseShape[0] == 2 && baseShape[1] == 2 && baseShape[2] == 2;
		Set<String> ternaryKeys = is222 ? ternaryKeys() : Set.of();
		Map<String, String[]> covered = is222 ? coveredBases() : Map.of();

		record Row(String pattern, long[] sym, long best, String bestAlloc, boolean ternary, String base) {}
		List<Row> rows = new ArrayList<>();
		for (int[][] idx : orbit.representativeShapes.values()) {
			long[] sym = new long[symSplits.size()];
			for (int si = 0; si < symSplits.size(); si++) {
				int[] c = symSplits.get(si);
				sym[si] = rankAt(idx, c, c, c, lk);
			}
			long best = Long.MAX_VALUE;
			int[] bn = null, bm = null, bp = null;
			for (int[] sN : splitsN) for (int[] sM : splitsM) for (int[] sP : splitsP) {
				long r = rankAt(idx, sN, sM, sP, lk);
				if (r > 0 && r < best) { best = r; bn = sN; bm = sM; bp = sP; }
			}
			String bestAlloc = best == Long.MAX_VALUE ? "?" : join(bn) + " / " + join(bm) + " / " + join(bp);
			boolean tern = false; String baseRef = "—";
			if (is222) {
				String key = canon(toSizes(idx, 9, 8));
				tern = ternaryKeys.contains(key);
				if (tern && covered.containsKey(key)) baseRef = covered.get(key)[0] + " (@" + covered.get(key)[1] + ")";
			}
			rows.add(new Row(patternOf(idx, is222), sym, best == Long.MAX_VALUE ? -1 : best, bestAlloc, tern, baseRef));
		}
		rows.sort(Comparator.comparingLong(Row::best));

		long[] symMin = new long[symSplits.size()];
		java.util.Arrays.fill(symMin, Long.MAX_VALUE);
		long bestMin = Long.MAX_VALUE;
		for (Row r : rows) {
			for (int c = 0; c < r.sym().length; c++) if (r.sym()[c] > 0) symMin[c] = Math.min(symMin[c], r.sym()[c]);
			if (r.best() > 0) bestMin = Math.min(bestMin, r.best());
		}

		// ── table: # | multiset | <symmetric a+b …> | best | best alloc | [ternary | base] ──
		List<String[]> table = new ArrayList<>();
		List<String> header = new ArrayList<>(List.of("#", "multiset"));
		for (int[] c : symSplits) header.add(join(c));
		header.add("best");
		header.add("best alloc (n / m / p)");
		if (is222) { header.add("ternary"); header.add("base (hash)"); }
		table.add(header.toArray(String[]::new));
		int i = 1;
		for (Row r : rows) {
			List<String> cells = new ArrayList<>();
			cells.add(Integer.toString(i++));
			cells.add(r.pattern());
			for (int c = 0; c < r.sym().length; c++) {
				String v = r.sym()[c] < 0 ? "?" : Long.toString(r.sym()[c]);
				if (r.sym()[c] == symMin[c]) v = "**" + v + "**";
				cells.add(v);
			}
			String bv = r.best() < 0 ? "?" : Long.toString(r.best());
			if (r.best() == bestMin) bv = "**" + bv + "**";
			cells.add(bv);
			cells.add(r.bestAlloc());
			if (is222) { cells.add(r.ternary() ? "✓" : "—"); cells.add(r.ternary() ? r.base() : "—"); }
			table.add(cells.toArray(String[]::new));
		}

		long ternCount = rows.stream().filter(Row::ternary).count();
		StringBuilder md = new StringBuilder();
		md.append("# ⟨").append(join3(target)).append("⟩ recombination multisets — base ⟨")
				.append(join3(baseShape)).append("⟩\n\n");
		md.append("Every recombination multiset the ⟨").append(join3(baseShape))
				.append("⟩ base can realise (exact GL-orbit via `RecombinationMultisetOrbit`), and the number of ")
				.append("multiplications each yields when recombined to ⟨").append(join3(target)).append("⟩.\n\n");
		md.append("**The allocation is independent per axis (n, m, p).** Each axis is split into ")
				.append(baseShape[0]).append("/").append(baseShape[1]).append("/").append(baseShape[2])
				.append(" parts; the `best` column is the minimum over the *joint* per-axis split space and ")
				.append("`best alloc` is the winning `(n / m / p)` split.");
		if (cubic) md.append(" The `a+b` columns are the *symmetric* splits (same on all three axes) for reference — ")
				.append("the per-axis optimum can beat them.");
		md.append("\n\n");
		md.append("- **").append(rows.size()).append(" multisets** total over ℚ");
		if (is222) md.append(", **").append(ternCount).append(" ternary** (constructable, cited by base hash; the other ")
				.append(rows.size() - ternCount).append(" need a rational change-of-basis — predictions only)");
		md.append(".\n");
		md.append("- `multiset` is allocation-independent");
		if (is222) md.append(" (B = larger block ⌈⌉, S = smaller block ⌊⌋; the per-product triple is axis-sorted)");
		md.append("; **bold** = column-best.\n\n");
		md.append(renderTable(table));

		// Per-base structural metadata: split each big block into plain (B) vs mixed
		// (M = B⊕S). Same cost, but explains why multisets with equal B/S counts differ.
		if (is222 && !covered.isEmpty()) {
			md.append("\n## Block structure (per ternary base)\n\n");
			md.append("The cost-multiset above uses only block *size* (B/S). But a `B` reaches the big size ")
					.append("two ways: a **plain** big block (touches only the big block) or a **mixed** `B⊕S` ")
					.append("sum (touches both). Same cost — yet the mixed block carries small-block support that ")
					.append("blocks the further shaving a plain `S` would permit, which is *why* two bases can ")
					.append("share a B/S cost-multiset and still be distinct. Legend: **B** plain big, ")
					.append("**M** mixed big (B⊕S), **S** small.\n\n");
			List<String[]> st = new ArrayList<>();
			st.add(new String[] { "base (hash)", "fine multiset (B / M / S)" });
			for (String key : ternaryKeys) {
				String[] cov = covered.get(key);
				if (cov == null) continue;
				NonCubicBilinearAlgorithm b = readBase(cov[0]);
				if (b != null) st.add(new String[] { cov[0] + " (@" + cov[1] + ")", finePattern(seedFrom(b)) });
			}
			md.append(renderTable(st));
		}

		md.append("\n_Generated by `RecombinationMultisetReport` (target ⟨").append(join3(target))
				.append("⟩, base ⟨").append(join3(baseShape)).append("⟩)._\n");
		return md.toString();
	}

	/** Σ catalog rank of the 7 sub-shapes at per-axis splits {@code sN/sM/sP}; -1 if any unknown. */
	private static long rankAt(int[][] idx, int[] sN, int[] sM, int[] sP, FieldAwareLookup lk) {
		long sum = 0;
		for (int[] tr : idx) {
			int A = sN[tr[0]], B = sM[tr[1]], C = sP[tr[2]];
			// A sub-shape with a unit axis is trivial: rank = product (no catalog entry).
			long r = (A == 1 || B == 1 || C == 1) ? (long) A * B * C : lk.findRank(A, B, C);
			if (r <= 0 || r >= Recombination.SotaResolver.UNKNOWN_RANK) return -1;
			sum += r;
		}
		return sum;
	}

	// ── multiset rendering ──────────────────────────────────────────────────
	/** Allocation-independent pattern: counts of each product's sorted block-index triple. */
	private static String patternOf(int[][] idx, boolean is222) {
		TreeMap<String, Integer> ct = new TreeMap<>();
		for (int[] tr : idx) {
			int[] s = tr.clone();
			java.util.Arrays.sort(s);
			StringBuilder sb = new StringBuilder();
			for (int x : s) sb.append(is222 ? (x == 0 ? "B" : "S") : Integer.toString(x));
			ct.merge(sb.toString(), 1, Integer::sum);
		}
		List<String> parts = new ArrayList<>();
		ct.forEach((k, v) -> parts.add(v + "×" + k));
		return String.join(" + ", parts);
	}

	/**
	 * Per-product <em>fine</em> pattern that distinguishes a PLAIN big block from a
	 * MIXED one. The cost-multiset uses only the block <em>size</em> (B/S), but a
	 * product can reach the big size two ways on an axis: by touching only the big
	 * block ("plain B") or by touching <em>both</em> blocks ("mixed", a B⊕S sum).
	 * Both cost the same, but the mixed block carries small-block support that
	 * blocks the further shaving a plain small would allow — so two multisets with
	 * equal B/S counts can be structurally different. Legend: B plain big, M mixed
	 * big (B⊕S), S small. Computed from the actual base scheme.
	 */
	private static String finePattern(int[][][] seed) {
		TreeMap<String, Integer> ct = new TreeMap<>();
		for (int k = 0; k < 7; k++) {
			int[] u = seed[k][0], v = seed[k][1], w = seed[k][2];
			int uRowN = 0, uColM = 0, vRowM = 0, vColP = 0, wRowN = 0, wColP = 0;
			for (int a = 0; a < 4; a++) if (u[a] != 0) { uRowN |= 1 << (a >> 1); uColM |= 1 << (a & 1); }
			for (int b = 0; b < 4; b++) if (v[b] != 0) { vRowM |= 1 << (b >> 1); vColP |= 1 << (b & 1); }
			for (int c = 0; c < 4; c++) if (w[c] != 0) { wRowN |= 1 << (c >> 1); wColP |= 1 << (c & 1); }
			char[] t = { cls(uRowN, wRowN), cls(uColM, vRowM), cls(vColP, wColP) };
			java.util.Arrays.sort(t);            // axis-symmetric grouping (B < M < S)
			ct.merge(new String(t), 1, Integer::sum);
		}
		List<String> parts = new ArrayList<>();
		ct.forEach((k, v) -> parts.add(v + "×" + k));
		return String.join(" + ", parts);
	}

	/** Effective block on an axis from its two factor views: B plain / M mixed (B⊕S) / S small. */
	private static char cls(int view1, int view2) {
		boolean bothBig = (view1 & 1) != 0 && (view2 & 1) != 0;   // sub-dim = big iff both views touch block0
		if (!bothBig) return 'S';
		boolean mixed = (view1 & 2) != 0 || (view2 & 2) != 0;     // a view also touches the small block
		return mixed ? 'M' : 'B';
	}

	private static NonCubicBilinearAlgorithm readBase(String filename) {
		for (String dir : new String[] { KNOWN2, CURATED2 }) {
			File f = new File(dir, filename);
			if (f.exists()) {
				try { return SchemeIO.readBilinear(f); } catch (Exception e) { return null; }
			}
		}
		return null;
	}

	private static int[][] toSizes(int[][] idx, int big, int small) {
		int[][] sz = new int[idx.length][3];
		for (int k = 0; k < idx.length; k++) for (int ax = 0; ax < 3; ax++) sz[k][ax] = idx[k][ax] == 0 ? big : small;
		return sz;
	}

	private static String canon(int[][] sh) {
		String best = null;
		for (int[] perm : PERMS3) {
			String[] s = new String[sh.length];
			for (int k = 0; k < sh.length; k++) s[k] = sh[k][perm[0]] + "," + sh[k][perm[1]] + "," + sh[k][perm[2]];
			java.util.Arrays.sort(s);
			String key = String.join("|", s);
			if (best == null || key.compareTo(best) < 0) best = key;
		}
		return best;
	}

	// ── ⟨2,2,2⟩ ternary classification ──────────────────────────────────────
	private static int[][][] seedFrom(NonCubicBilinearAlgorithm a) {
		double[][] U = a.denseU(), V = a.denseV(), W = a.denseW();
		int[][][] seed = new int[7][3][4];
		for (int k = 0; k < 7; k++) for (int i = 0; i < 4; i++) {
			seed[k][0][i] = (int) Math.round(U[i][k]);
			seed[k][1][i] = (int) Math.round(V[i][k]);
			seed[k][2][i] = (int) Math.round(W[i][k]);
		}
		return seed;
	}

	private static Set<String> ternaryKeys() throws Exception {
		List<int[][][]> seeds = new ArrayList<>();
		File[] kf = new File(KNOWN2).listFiles((d, n) -> n.endsWith(".json"));
		if (kf != null) for (File f : kf) {
			NonCubicBilinearAlgorithm a;
			try { a = SchemeIO.readBilinear(f); } catch (Exception e) { continue; }
			if (a.n == 2 && a.m == 2 && a.p == 2 && a.r == 7 && Verifier.isExactNonCubic(a)) seeds.add(seedFrom(a));
		}
		return Ternary2x2x2Orbit.sweep(seeds, -1, 1).representatives.keySet();
	}

	private static Map<String, String[]> coveredBases() throws Exception {
		Map<String, String[]> covered = new LinkedHashMap<>();
		for (String dir : new String[] { KNOWN2, CURATED2 }) {
			File[] files = new File(dir).listFiles((d, n) -> n.endsWith(".json"));
			if (files == null) continue;
			for (File f : files) {
				NonCubicBilinearAlgorithm a;
				try { a = SchemeIO.readBilinear(f); } catch (Exception e) { continue; }
				if (a.n != 2 || a.m != 2 || a.p != 2 || a.r != 7 || !Verifier.isExactNonCubic(a)) continue;
				covered.putIfAbsent(Ternary2x2x2Orbit.canonicalMultisetKey(seedFrom(a)),
						new String[] { f.getName(), SchemeIO.contentHash(a).substring(0, 7) });
			}
		}
		return covered;
	}

	// ── base resolution + allocations ───────────────────────────────────────
	private static NonCubicBilinearAlgorithm resolveIntegerBase(int[] s, FieldAwareLookup lk) {
		for (Path p : lk.findFiles(s[0], s[1], s[2])) {
			NonCubicBilinearAlgorithm a;
			try { a = SchemeIO.readBilinear(p.toFile()); } catch (Exception e) { continue; }
			if (a.n != s[0] || a.m != s[1] || a.p != s[2] || !Verifier.isExactNonCubic(a)) continue;
			if (isInteger(a)) return a;
		}
		return null;
	}

	private static boolean isInteger(NonCubicBilinearAlgorithm a) {
		for (double[][] M : new double[][][] { a.denseU(), a.denseV(), a.denseW() })
			for (double[] row : M) for (double x : row) if (Math.abs(x - Math.rint(x)) > 1e-9) return false;
		return true;
	}

	/** Descending compositions of {@code n} into exactly {@code k} positive parts. */
	private static List<int[]> descCompositions(int n, int k) {
		List<int[]> out = new ArrayList<>();
		descComp(n, k, n, new int[k], 0, out);
		return out;
	}

	private static void descComp(int rem, int k, int maxP, int[] cur, int idx, List<int[]> out) {
		if (k == 1) { if (rem >= 1 && rem <= maxP) { cur[idx] = rem; out.add(cur.clone()); } return; }
		for (int p = Math.min(maxP, rem - (k - 1)); p >= 1; p--) {
			cur[idx] = p; descComp(rem - p, k - 1, p, cur, idx + 1, out);
		}
	}

	// ── markdown table padding ──────────────────────────────────────────────
	private static String renderTable(List<String[]> rows) {
		int cols = rows.get(0).length;
		int[] w = new int[cols];
		for (String[] r : rows) for (int c = 0; c < cols; c++) w[c] = Math.max(w[c], r[c].length());
		StringBuilder sb = new StringBuilder();
		for (int ri = 0; ri < rows.size(); ri++) {
			String[] r = rows.get(ri);
			sb.append("|");
			for (int c = 0; c < cols; c++) sb.append(' ').append(pad(r[c], w[c])).append(" |");
			sb.append('\n');
			if (ri == 0) {
				sb.append("|");
				for (int c = 0; c < cols; c++) sb.append(' ').append("-".repeat(Math.max(3, w[c]))).append(" |");
				sb.append('\n');
			}
		}
		return sb.toString();
	}

	private static String pad(String s, int w) {
		StringBuilder b = new StringBuilder(s);
		while (b.length() < w) b.append(' ');
		return b.toString();
	}

	// ── small parse/format helpers ──────────────────────────────────────────
	private static String arg(String[] a, String key, String def) {
		for (String s : a) if (s.startsWith(key + "=")) return s.substring(key.length() + 1);
		return def;
	}

	private static int[] parseShape(String s) {
		String[] t = s.split("x");
		return new int[] { Integer.parseInt(t[0]), Integer.parseInt(t[1]), Integer.parseInt(t[2]) };
	}

	private static String show(int[] s) { return s[0] + "x" + s[1] + "x" + s[2]; }

	private static String join3(int[] s) { return s[0] + "," + s[1] + "," + s[2]; }

	private static String join(int[] c) {
		StringBuilder b = new StringBuilder();
		for (int i = 0; i < c.length; i++) { if (i > 0) b.append('+'); b.append(c[i]); }
		return b.toString();
	}
}
