package eu.solven.matmul.docs;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import eu.solven.matmul.ComplexNonCubicBilinearAlgorithm;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.KnownAlgorithm;
import eu.solven.matmul.catalog.KnownAlgorithmCatalog;
import eu.solven.matmul.catalog.SchemeIO;
import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;

/**
 * One-shot generator: walks {@code src/main/resources/schemes/} and
 * {@link KnownAlgorithmCatalog}, then writes {@code COVERAGE.md} at the repo
 * root showing what's verified vs. what's only known from the literature.
 *
 * <p>Re-run via:
 * <pre>
 *   mvn -q test-compile
 *   java -cp target/classes:target/test-classes:&dollar;CLASSPATH \
 *        eu.solven.matmul.catalog.CoverageMatrixGenerator
 * </pre>
 *
 * <p>Output rows are grouped by <b>max-dimension</b> ({@code max(n,m,p)}) and
 * formats within each group are sorted lexicographically with components in
 * ascending order (so {@code ⟨2,2,3⟩}, {@code ⟨2,3,3⟩}, {@code ⟨3,3,3⟩} all
 * appear in section "max-dim 3"). Per-format rows show the per-field best
 * known rank, with `✓` for verified schemes (factor matrices on disk) and `📄`
 * for literature-only entries.</p>
 */
public class CoverageMatrixGenerator {

	private static final File SCHEMES_DIR = new File("src/main/resources/schemes");
	private static final File OUTPUT = new File("generated/COVERAGE.md");

	private static final Pattern SCHEME_FILE = Pattern.compile(
			"(?<source>.*)_(?<n>\\d+)x(?<m>\\d+)x(?<p>\\d+)_(?:r|m)(?<rank>\\d+)[^/]*\\.json");

	public static void main(String[] args) throws IOException {
		List<Entry> entries = new ArrayList<>();
		entries.addAll(scanSchemes());
		entries.addAll(scanCatalog());

		// Group ALL catalog entries per canonical format for the "history" column.
		Map<CanonicalFormat, List<Entry>> historyByFormat = new LinkedHashMap<>();
		for (Entry e : entries) {
			if (e.year == 0) continue; // skip scheme-file entries (no year)
			historyByFormat.computeIfAbsent(e.canonical(), x -> new ArrayList<>()).add(e);
		}
		for (List<Entry> hist : historyByFormat.values()) {
			hist.sort(Comparator.comparingInt((Entry x) -> x.year).thenComparingInt(x -> x.rank));
		}

		// Best per (canonical format, field): lower rank wins, ties broken by verified > literature.
		Map<Key, Entry> best = new LinkedHashMap<>();
		for (Entry e : entries) {
			Key k = new Key(e.canonical(), e.field);
			Entry prev = best.get(k);
			if (prev == null || e.rank < prev.rank
					|| (e.rank == prev.rank && e.verified && !prev.verified)) {
				best.put(k, e);
			}
		}

		// Group canonical formats by max-dim.
		Map<Integer, List<CanonicalFormat>> formatsByMaxDim = new TreeMap<>();
		for (Key k : best.keySet()) {
			formatsByMaxDim.computeIfAbsent(k.format.maxDim(), x -> new ArrayList<>())
					.add(k.format);
		}
		for (List<CanonicalFormat> list : formatsByMaxDim.values()) {
			List<CanonicalFormat> sorted = list.stream().distinct()
					.sorted(Comparator.comparingInt((CanonicalFormat f) -> f.a)
							.thenComparingInt(f -> f.b)
							.thenComparingInt(f -> f.c))
					.toList();
			list.clear();
			list.addAll(sorted);
		}

		writeMd(formatsByMaxDim, best, historyByFormat);
		System.out.println("wrote " + OUTPUT.getAbsolutePath());
	}

	private static List<Entry> scanSchemes() throws IOException {
		List<Entry> out = new ArrayList<>();
		if (!SCHEMES_DIR.isDirectory()) return out;
		File[] files;
		try (var s = Files.walk(SCHEMES_DIR.toPath())) {
			files = s.filter(p -> p.toString().endsWith(".json"))
					.map(Path::toFile)
					.sorted()
					.toArray(File[]::new);
		}
		for (File f : files) {
			Matcher m = SCHEME_FILE.matcher(f.getName());
			if (!m.matches()) continue;
			int n = Integer.parseInt(m.group("n"));
			int mm = Integer.parseInt(m.group("m"));
			int p = Integer.parseInt(m.group("p"));
			int rank = Integer.parseInt(m.group("rank"));
			String source = m.group("source");

			String body = Files.readString(f.toPath());
			JsonNode root;
			try {
				root = JsonMapper.builder().build().readTree(body);
			} catch (Exception e) {
				System.err.println("skip " + f.getName() + ": parse error " + e.getMessage());
				continue;
			}
			boolean complex = root.path("complex").asBoolean(false);
			boolean z2 = SchemeIO.isZ2(root);
			String fieldStr = root.path("field").asText(null);
			FieldBucket fb = bucketize(fieldStr, z2, complex);

			Integer additions = null;
			boolean integerCoefs = false;
			try {
				if (complex) {
					ComplexNonCubicBilinearAlgorithm alg = SchemeIO.readComplex(f);
					additions = Verifier.additionCount(alg);
				} else {
					NonCubicBilinearAlgorithm alg = SchemeIO.read(f);
					additions = Verifier.additionCount(alg);
					integerCoefs = hasOnlyIntegerCoefs(alg);
				}
			} catch (Exception e) {
				// Stays null — coverage row still gets the verified ✓.
			}
			out.add(new Entry(new Format(n, mm, p), fb, rank, capitalize(source),
					0, true, additions, f.getName()));

			// Lift integer-coefficient R-algorithms into the F₂ column: any scheme
			// with entries in {-1, 0, 1, …} ⊂ Z reduces mod 2 to a valid F₂ scheme
			// at the same rank. Half-integer / complex schemes don't lift.
			if (fb == FieldBucket.REAL_RING && integerCoefs) {
				out.add(new Entry(new Format(n, mm, p), FieldBucket.F2, rank,
						capitalize(source) + " [Z→F₂]", 0, true, additions, f.getName()));
			}
		}
		return out;
	}

	private static boolean hasOnlyIntegerCoefs(NonCubicBilinearAlgorithm alg) {
		double[][] srcU = alg.denseU();
		double[][] srcV = alg.denseV();
		double[][] srcW = alg.denseW();
		return isAllIntegers(srcU) && isAllIntegers(srcV) && isAllIntegers(srcW);
	}

	private static boolean isAllIntegers(double[][] m) {
		for (double[] row : m) {
			for (double x : row) {
				if (x != Math.rint(x)) return false;
			}
		}
		return true;
	}

	private static List<Entry> scanCatalog() {
		List<Entry> out = new ArrayList<>();
		for (KnownAlgorithm a : KnownAlgorithmCatalog.all()) {
			FieldBucket fb;
			if (a.algebra.commutative()) {
				fb = FieldBucket.COMMUTATIVE;
			} else {
				fb = switch (a.algebra.field()) {
					case F2 -> FieldBucket.F2;
					case F3 -> FieldBucket.OTHER;
					case Z, Q, R -> FieldBucket.REAL_RING;
					case C -> FieldBucket.COMPLEX;
				};
			}
			out.add(new Entry(new Format(a.n, a.m, a.p), fb, a.rank, a.source,
					a.year, false, null, null));
		}
		return out;
	}

	private static FieldBucket bucketize(String fieldStr, boolean z2, boolean complex) {
		if (z2) return FieldBucket.F2;
		if (complex) return FieldBucket.COMPLEX;
		if (fieldStr == null) return FieldBucket.REAL_RING;
		String s = fieldStr.toLowerCase().replace(" ", "");
		if (s.contains("c")) return FieldBucket.COMPLEX;
		if (s.equals("z/2") || s.equals("f_2") || s.equals("f2") || s.equals("gf(2)"))
			return FieldBucket.F2;
		return FieldBucket.REAL_RING;
	}

	private static String capitalize(String s) {
		if (s == null || s.isEmpty()) return s;
		return Character.toUpperCase(s.charAt(0)) + s.substring(1);
	}

	private static void writeMd(Map<Integer, List<CanonicalFormat>> formatsByMaxDim,
			Map<Key, Entry> best,
			Map<CanonicalFormat, List<Entry>> historyByFormat) throws IOException {
		try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(OUTPUT)))) {
			pw.println("# Matrix-Multiplication Algorithm Coverage");
			pw.println();
			pw.println("Auto-generated from `src/main/resources/schemes/` + `KnownAlgorithmCatalog`.");
			pw.println("Re-run via `java … CoverageMatrixGenerator`.");
			pw.println();
			pw.println("Cell legend:");
			pw.println();
			pw.println("- **✓ r=N** — verified scheme on disk at `src/main/resources/schemes/`. N is the multiplication count.");
			pw.println("- **📄 r=N** — known from the literature, no scheme file yet.");
			pw.println("- **—** — not tracked at all.");
			pw.println();
			pw.println("`+N adds` is the addition count (`nz(U)+nz(V)+nz(W) − 2r − n·p`), shown for verified schemes only.");
			pw.println();
			pw.println("Formats are grouped by `max(n,m,p)` (the **section number**), then sorted with");
			pw.println("components in ascending order. Permutations of `⟨n,m,p⟩` are merged (the rank is");
			pw.println("invariant under axis permutation).");
			pw.println();

			FieldBucket[] cols = { FieldBucket.F2, FieldBucket.REAL_RING,
					FieldBucket.COMPLEX, FieldBucket.COMMUTATIVE };

			for (Map.Entry<Integer, List<CanonicalFormat>> e : formatsByMaxDim.entrySet()) {
				int maxDim = e.getKey();
				pw.printf("## Section %d — max-dimension = %d%n%n", maxDim, maxDim);
				pw.print("| format |");
				for (FieldBucket fb : cols) pw.printf(" %s |", fb.label);
				pw.println(" best (any field) | adds | history |");
				pw.print("|---|");
				for (int i = 0; i < cols.length; i++) pw.print("---|");
				pw.println("---|---|---|");

				for (CanonicalFormat cf : e.getValue()) {
					pw.printf("| `%s` |", cf);
					Integer additions = null;
					int bestRank = Integer.MAX_VALUE;
					List<String> bestFields = new ArrayList<>();
					for (FieldBucket fb : cols) {
						Entry entry = best.get(new Key(cf, fb));
						pw.printf(" %s |", formatCell(entry));
						if (entry != null && entry.additions != null && additions == null) {
							additions = entry.additions;
						}
						if (entry != null) {
							if (entry.rank < bestRank) {
								bestRank = entry.rank;
								bestFields.clear();
								bestFields.add(fb.label);
							} else if (entry.rank == bestRank) {
								bestFields.add(fb.label);
							}
						}
					}
					if (bestRank == Integer.MAX_VALUE) {
						pw.print(" — |");
					} else {
						pw.printf(" **%d** in %s |", bestRank, String.join(", ", bestFields));
					}
					pw.printf(" %s |", additions == null ? "—" : ("+" + additions));
					pw.printf(" %s |%n", formatHistory(historyByFormat.get(cf)));
				}
				pw.println();
			}

			pw.println("---");
			pw.println();
			pw.println("## How to fill gaps");
			pw.println();
			pw.println("- A `📄` cell becomes `✓` once a scheme file is committed to");
			pw.println("  `src/main/resources/schemes/` and round-trips through `SchemeIO` + `Verifier`.");
			pw.println("- A `—` cell becomes `📄` once an entry is added to");
			pw.println("  `eu.solven.matmul.catalog.KnownAlgorithmCatalog`.");
			pw.println("- Important format-field cells we know exist but don't have schemes for:");
			pw.println("  AlphaTensor's full F₂ result set (~50 schemes), Smirnov's `{-1,0,+1}` catalog");
			pw.println("  for `⟨6,6,6⟩` / `⟨7,7,7⟩` / etc.");
		}
	}

	private static String formatCell(Entry e) {
		if (e == null) return "—";
		String tag = e.verified ? "✓" : "📄";
		return String.format("**%s %d** %s", tag, e.rank, e.source);
	}

	/**
	 * Render the chronological catalog history for a format as
	 * {@code year src=rank/field} bullets joined with {@code <br>} so they
	 * render as a small bulleted list inside one Markdown table cell. Returns
	 * {@code —} if no catalog entries exist (only scheme files).
	 */
	private static String formatHistory(List<Entry> hist) {
		if (hist == null || hist.isEmpty()) return "—";
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < hist.size(); i++) {
			Entry e = hist.get(i);
			if (i > 0) sb.append("<br>");
			sb.append("• ")
					.append(e.year).append(' ')
					.append(e.source).append(" r=").append(e.rank)
					.append(" (").append(e.field.label).append(')');
		}
		return sb.toString();
	}

	// ───────────────────────────────────────────────────────────────────────────

	enum FieldBucket {
		F2("F₂"), REAL_RING("R/Q/Z"), COMPLEX("C"), COMMUTATIVE("commutative"), OTHER("other");
		final String label;
		FieldBucket(String l) { this.label = l; }
	}

	record Format(int n, int m, int p) {
		CanonicalFormat canonical() {
			int[] a = { n, m, p };
			Arrays.sort(a);
			return new CanonicalFormat(a[0], a[1], a[2]);
		}
	}

	record CanonicalFormat(int a, int b, int c) {
		int maxDim() { return c; }
		@Override public String toString() { return String.format("⟨%d,%d,%d⟩", a, b, c); }
	}

	record Key(CanonicalFormat format, FieldBucket field) {}

	record Entry(Format raw, FieldBucket field, int rank, String source, int year,
			boolean verified, Integer additions, String fileName) {
		CanonicalFormat canonical() { return raw.canonical(); }
	}
}
