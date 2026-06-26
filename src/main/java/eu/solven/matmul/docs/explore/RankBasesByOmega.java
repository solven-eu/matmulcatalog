package eu.solven.matmul.docs.explore;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import lombok.extern.slf4j.Slf4j;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.ObjectMapper;
import tools.jackson.databind.json.JsonMapper;

/**
 * Rank every catalog scheme by its implied {@code ω} value, so we can
 * pick outer bases by ω first and shape second. Output is a Markdown
 * table sorted by ω ascending, optionally filtered to non-commutative
 * leaf schemes (the natural candidate set for the BlockSplitSearch
 * pool).
 *
 * <p>For a scheme of shape {@code ⟨n,m,p⟩} with rank {@code r}, the
 * implied exponent is
 * <pre>
 *   ω = 3·log(r) / log(n·m·p)
 * </pre>
 * — the value of the matrix-multiplication exponent that {@code r}
 * recursive uses of this scheme would yield asymptotically. Lower is
 * better.
 *
 * <p>This driver doesn't change anything in the pool itself — it just
 * surfaces the data. See task linking the result to
 * {@link eu.solven.matmul.search.RecombinationPoolConfig}.
 */
@Slf4j
public final class RankBasesByOmega {

	private static final Path CATALOG = Path.of("docs/catalog.json");

	private RankBasesByOmega() {}

	public static void main(String[] args) throws IOException {
		int maxDim = 12;
		Path outPath = Path.of("docs/bases-by-omega.md");
		boolean leafOnly = true;
		boolean ncOnly = true;
		int targetDim = 17;     // reference cubic target ⟨targetDim,targetDim,targetDim⟩
		int maxImbalance = 3;   // for the imb-bounded column
		for (String a : args) {
			if (a.startsWith("--max-dim=")) maxDim = Integer.parseInt(a.substring("--max-dim=".length()));
			else if (a.startsWith("--out=")) outPath = Path.of(a.substring("--out=".length()));
			else if (a.startsWith("--target=")) targetDim = Integer.parseInt(a.substring("--target=".length()));
			else if (a.startsWith("--max-imbalance=")) maxImbalance = Integer.parseInt(a.substring("--max-imbalance=".length()));
			else if ("--include-derived".equals(a)) leafOnly = false;
			else if ("--include-commutative".equals(a)) ncOnly = false;
			else if ("--help".equals(a) || "-h".equals(a)) {
				System.out.println(
						"Usage: RankBasesByOmega [--max-dim=N] [--out=PATH] "
								+ "[--target=K] [--max-imbalance=M] "
								+ "[--include-derived] [--include-commutative]");
				return;
			} else throw new IllegalArgumentException("Unknown arg: " + a);
		}
		final int TARGET_DIM = targetDim;
		final int MAX_IMB = maxImbalance;

		ObjectMapper mapper = JsonMapper.builder().build();
		JsonNode catalog = mapper.readTree(CATALOG.toFile());
		JsonNode schemes = catalog.get("schemes");

		List<Row> rows = new ArrayList<>();
		for (JsonNode s : schemes) {
			JsonNode fmt = s.get("format");
			if (fmt == null || !fmt.isArray() || fmt.size() != 3) continue;
			int n = fmt.get(0).asInt(), m = fmt.get(1).asInt(), p = fmt.get(2).asInt();
			int max = Math.max(n, Math.max(m, p));
			if (max > maxDim) continue;
			if (n < 2 || m < 2 || p < 2) continue;
			int rank = s.has("rank") ? s.get("rank").asInt() : -1;
			if (rank < 1) continue;
			String field = s.has("field") ? s.get("field").asText("") : "";
			String source = s.has("source") ? s.get("source").asText("") : "";
			String src = source.toLowerCase();
			String file = s.has("file") ? s.get("file").asText("") : "";
			// catalog.json doesn't carry the per-scheme "commutative" flag.
			// Apply the same source/filename heuristic FrontierClosure uses
			// (commutative authors + _commutative filename suffix).
			boolean commutative = s.has("commutative") && s.get("commutative").asBoolean(false);
			if (!commutative) {
				if (file.contains("_commutative")) commutative = true;
				else if (src.startsWith("waksman") || src.startsWith("rosowski")
						|| src.startsWith("makarov") || src.startsWith("islam")
						|| src.startsWith("smith ") || src.startsWith("probert")) commutative = true;
			}
			if (ncOnly && commutative) continue;
			if (leafOnly && src.startsWith("derived-")) continue;
			boolean cubic = (n == m) && (m == p);
			double omega = 3.0 * Math.log(rank) / Math.log((long) n * m * p);

			// ── Allocation counts at the reference target ⟨T,T,T⟩. ─
			// Only meaningful when every axis fits (n,m,p ≤ T). For
			// larger bases write 0.
			long allocsRaw = 0, allocsImb = 0, allocsFlip = 0, allocsPerm = 0;
			if (n <= TARGET_DIM && m <= TARGET_DIM && p <= TARGET_DIM) {
				List<int[]> allocsAxisN = compositionsPositive(n, TARGET_DIM);
				List<int[]> allocsAxisM = compositionsPositive(m, TARGET_DIM);
				List<int[]> allocsAxisP = compositionsPositive(p, TARGET_DIM);
				List<int[]> imbN = filterByImbalance(allocsAxisN, MAX_IMB);
				List<int[]> imbM = filterByImbalance(allocsAxisM, MAX_IMB);
				List<int[]> imbP = filterByImbalance(allocsAxisP, MAX_IMB);
				allocsRaw = (long) allocsAxisN.size() * allocsAxisM.size() * allocsAxisP.size();
				allocsImb = (long) imbN.size() * imbM.size() * imbP.size();
				// Axis-flip canonicalisation — per axis, keep allocations
				// that are lex-min of the {alloc, rev(alloc)} pair.
				long flipN = imbN.stream().filter(a -> !lexGreaterThanReverse(a)).count();
				long flipM = imbM.stream().filter(a -> !lexGreaterThanReverse(a)).count();
				long flipP = imbP.stream().filter(a -> !lexGreaterThanReverse(a)).count();
				allocsFlip = flipN * flipM * flipP;
				// Axis-permutation canonicalisation — applies only at
				// cubic targets with cubic bases. For each (allocA,
				// allocB, allocC), keep only the lex-min of its S₃
				// orbit. Enumerated with a cap to avoid blow-up on
				// large bases; we emit -1 ("--") when the enumeration
				// would exceed CAP.
				if (cubic) {
					final long CAP = 5_000_000L;
					if (allocsFlip > CAP) {
						allocsPerm = -1;
					} else {
						List<int[]> flipCanonN = imbN.stream()
								.filter(a -> !lexGreaterThanReverse(a)).toList();
						long permCount = 0;
						for (int[] aA : flipCanonN) {
							for (int[] aB : flipCanonN) {
								for (int[] aC : flipCanonN) {
									if (isS3Canonical(aA, aB, aC)) permCount++;
								}
							}
						}
						allocsPerm = permCount;
					}
				} else {
					// Non-cubic base — S₃ on (allocA, allocB, allocC)
					// doesn't apply directly. Future work: per-base
					// permutation symmetry detection.
					allocsPerm = allocsFlip;  // no further reduction
				}
			}

			rows.add(new Row(n, m, p, rank, max, cubic, field, commutative, source, omega,
					allocsRaw, allocsImb, allocsFlip, allocsPerm));
		}

		// Sort by ω ascending (primary), then by max-dim ascending (secondary)
		rows.sort(Comparator.<Row>comparingDouble(r -> r.omega)
				.thenComparingInt(r -> r.maxDim)
				.thenComparingInt(r -> r.rank));

		Files.createDirectories(outPath.getParent());
		try (PrintWriter pw = new PrintWriter(outPath.toFile())) {
			pw.printf("# Candidate outer bases ranked by implied ω%n%n");
			pw.printf("Filter: leafOnly=%s, ncOnly=%s, maxDim≤%d.  ω = 3·log(r)/log(n·m·p).%n%n",
					leafOnly, ncOnly, maxDim);
			pw.printf("Allocation columns are evaluated at the reference cubic target ⟨%d,%d,%d⟩.%n",
					TARGET_DIM, TARGET_DIM, TARGET_DIM);
			pw.printf("  - **raw**:  product of per-axis compositions, no constraints.%n");
			pw.printf("  - **imb≤%d**:  filtered by per-axis max−min ≤ %d.%n", MAX_IMB, MAX_IMB);
			pw.printf("  - **flip-canon**:  + per-axis allocation kept only if lex-smaller than its reverse "
					+ "(safe iff axis-flip orbit of base is covered by mask sweep or pool expansion).%n");
			pw.printf("  - **+perm-canon**:  + cubic axis-permutation orbit canonicalised on (aA,aB,aC); "
					+ "applies to cubic bases on the cubic target only. `--` means enumeration capped.%n%n");
			pw.printf("Source: regenerate via `eu.solven.matmul.docs.explore.RankBasesByOmega`.%n%n");
			pw.printf("| ω | Shape | Rank | Cubic? | Field | Source "
					+ "| raw | imb≤%d | flip-canon | +perm-canon |%n", MAX_IMB);
			pw.printf("|---|---|---|---|---|---|---:|---:|---:|---:|%n");
			for (Row r : rows) {
				pw.printf("| %.4f | ⟨%d,%d,%d⟩ | %d | %s | %s | %s "
								+ "| %s | %s | %s | %s |%n",
						r.omega, r.n, r.m, r.p, r.rank,
						r.cubic ? "yes" : "no", r.field, r.source,
						fmt(r.allocsRaw), fmt(r.allocsImb),
						fmt(r.allocsFlip), fmt(r.allocsPerm));
			}
		}
		log.info("Wrote {} rows to {}", rows.size(), outPath);
		log.info("Top 15 by ω:");
		for (int i = 0; i < Math.min(15, rows.size()); i++) {
			Row r = rows.get(i);
			log.info("  ω={} ⟨{},{},{}⟩={} ({})",
					String.format("%.4f", r.omega), r.n, r.m, r.p, r.rank, r.source);
		}
	}

	private record Row(int n, int m, int p, int rank, int maxDim, boolean cubic,
			String field, boolean commutative, String source, double omega,
			long allocsRaw, long allocsImb, long allocsFlip, long allocsPerm) {}

	/** All positive-part compositions of {@code budget} into {@code parts} parts. */
	private static List<int[]> compositionsPositive(int parts, int budget) {
		List<int[]> out = new ArrayList<>();
		int[] buf = new int[parts];
		emitCompositions(parts, budget, 0, buf, out);
		return out;
	}

	private static void emitCompositions(int parts, int remaining, int idx,
			int[] buf, List<int[]> out) {
		if (idx == parts - 1) {
			if (remaining >= 1) {
				buf[idx] = remaining;
				out.add(buf.clone());
			}
			return;
		}
		int max = remaining - (parts - 1 - idx);
		for (int v = 1; v <= max; v++) {
			buf[idx] = v;
			emitCompositions(parts, remaining - v, idx + 1, buf, out);
		}
	}

	private static List<int[]> filterByImbalance(List<int[]> allocs, int maxImb) {
		List<int[]> out = new ArrayList<>();
		for (int[] a : allocs) {
			int mn = a[0], mx = a[0];
			for (int v : a) { if (v < mn) mn = v; if (v > mx) mx = v; }
			if (mx - mn <= maxImb) out.add(a);
		}
		return out;
	}

	/**
	 * True iff {@code alloc} is lexicographically strictly greater than
	 * its reverse — i.e.\ it's the "larger" representative of its
	 * reverse-orbit. Palindromes (alloc == rev(alloc)) return false (kept).
	 */
	private static boolean lexGreaterThanReverse(int[] alloc) {
		int n = alloc.length;
		for (int i = 0; i < n / 2; i++) {
			int j = n - 1 - i;
			if (alloc[i] != alloc[j]) return alloc[i] > alloc[j];
		}
		return false;
	}

	/**
	 * Check whether the triple {@code (aA, aB, aC)} is the
	 * lexicographically smallest among its 6 S₃ permutations. Used to
	 * canonicalise the axis-permutation orbit for cubic-base/cubic-target
	 * search.
	 */
	private static boolean isS3Canonical(int[] aA, int[] aB, int[] aC) {
		// Compare against the 5 other permutations
		if (compareTriple(aA, aB, aC, aA, aC, aB) > 0) return false;
		if (compareTriple(aA, aB, aC, aB, aA, aC) > 0) return false;
		if (compareTriple(aA, aB, aC, aB, aC, aA) > 0) return false;
		if (compareTriple(aA, aB, aC, aC, aA, aB) > 0) return false;
		if (compareTriple(aA, aB, aC, aC, aB, aA) > 0) return false;
		return true;
	}

	private static int compareTriple(int[] a1, int[] a2, int[] a3,
			int[] b1, int[] b2, int[] b3) {
		int c = compareInts(a1, b1);
		if (c != 0) return c;
		c = compareInts(a2, b2);
		if (c != 0) return c;
		return compareInts(a3, b3);
	}

	private static int compareInts(int[] a, int[] b) {
		int n = Math.min(a.length, b.length);
		for (int i = 0; i < n; i++) {
			if (a[i] != b[i]) return Integer.compare(a[i], b[i]);
		}
		return Integer.compare(a.length, b.length);
	}

	private static String fmt(long v) {
		if (v < 0) return "--";
		if (v == 0) return "0";
		if (v >= 1_000_000_000L) return String.format("%.1fG", v / 1e9);
		if (v >= 1_000_000L) return String.format("%.1fM", v / 1e6);
		if (v >= 1_000L) return String.format("%.1fk", v / 1e3);
		return String.valueOf(v);
	}
}
