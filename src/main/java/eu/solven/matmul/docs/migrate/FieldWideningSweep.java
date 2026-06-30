package eu.solven.matmul.docs.migrate;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

import eu.solven.matmul.ComplexNonCubicBilinearAlgorithm;
import eu.solven.matmul.NonBilinearAlgorithm;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.SchemeIO;
import lombok.extern.slf4j.Slf4j;
import tools.jackson.databind.JsonNode;

/**
 * Audit every catalog scheme to detect whether it is valid in a STRICTER
 * field than the one currently tagged.
 *
 * <p>Our pipeline tags schemes by <i>source</i> (e.g. {@code "field": "Q"}
 * because Perminov's Q sub-catalog is the importer) rather than by the
 * actual coefficient profile. So a {@code "field": "Q"} scheme whose
 * U/V/W coefficients are all integers is also Z-valid; an
 * {@code "field": "C"} scheme whose imaginary parts are all zero is
 * R-valid; if those reals are integers, it is Z-valid. The catalog
 * lookup misses those schemes on stricter-field queries.</p>
 *
 * <h3>Field lattice (strict → permissive)</h3>
 * <pre>
 *   F2_pm1  -- subset of Z with coefficients in {-1, 0, +1} that is also
 *             F2-compatible (no negative entries — i.e. all in {0, 1})
 *   PM1     -- coefficients in {-1, 0, +1} (Z, but extra-clean — these are
 *             the schemes most likely to survive symbolic transforms)
 *   Z       -- coefficients in Z (within tolerance)
 *   Q       -- coefficients in Q (rational with bounded denominator)
 *   R       -- coefficients in R (all imaginary parts within tolerance of 0)
 *   C       -- coefficients in C
 * </pre>
 *
 * <p>The tool proposes the most restrictive valid tag based on the actual
 * U/V/W entries. It does NOT modify scheme files — it produces a markdown
 * proposal report at {@code references/field-widening-proposals.md}.</p>
 *
 * <h3>Verification of proposals</h3>
 * <p>For schemes where a narrower field is proposed, the tool re-runs
 * {@link Verifier#isExactNonCubic} (or the F2 / complex / non-bilinear
 * variant) to confirm the scheme still verifies — this just sanity-checks
 * that the parsed data round-trips, since narrowing the algebraic class
 * does not change the bilinear identity. The pass count is reported.</p>
 *
 * <p>Run with:</p>
 * <pre>
 *   mvn -q exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.FieldWideningSweep
 * </pre>
 */
@Slf4j
public final class FieldWideningSweep {

	private static final Path SCHEMES_ROOT = Path.of("src/main/resources/schemes");
	private static final Path REPORT_PATH =
			Path.of("references/field-widening-proposals.md");

	/** Tolerance for treating a double as the closest integer / zero. */
	static final double EPS = 1e-9;

	/** Max denominator we consider acceptable for a "rational" classification. */
	static final int Q_MAX_DENOM = 4096;

	/** Matches {@code _{n}x{m}x{p}_r{rank}} in scheme filenames. */
	private static final Pattern SHAPE = Pattern.compile("[_-](\\d+)x(\\d+)x(\\d+)_(?:r|m)(\\d+)");

	private FieldWideningSweep() {}

	public static void main(String[] args) throws Exception {
		List<Path> files;
		try (Stream<Path> walk = Files.walk(SCHEMES_ROOT)) {
			files = walk.filter(p -> p.toString().endsWith(".json"))
					.sorted()
					.toList();
		}
		log.info("Scanning {} scheme files under {}", files.size(), SCHEMES_ROOT);

		List<Proposal> proposals = new ArrayList<>();
		AtomicInteger processed = new AtomicInteger();
		AtomicInteger unreadable = new AtomicInteger();
		long t0 = System.nanoTime();

		for (Path p : files) {
			try {
				Proposal pr = analyze(p);
				if (pr != null) proposals.add(pr);
			} catch (Exception e) {
				unreadable.incrementAndGet();
				log.warn("UNREADABLE {}: {}: {}", SCHEMES_ROOT.relativize(p),
						e.getClass().getSimpleName(), e.getMessage());
			}
			int done = processed.incrementAndGet();
			if (done % 500 == 0 || done == files.size()) {
				long elapsedMs = (System.nanoTime() - t0) / 1_000_000L;
				log.info("[progress] {}/{} files scanned, {} proposals so far, {}ms",
						done, files.size(), proposals.size(), elapsedMs);
			}
		}

		// Re-verify schemes with proposed changes to count how many round-trip.
		int verifyPass = 0;
		int verifyFail = 0;
		for (Proposal pr : proposals) {
			if (!pr.isWidening()) continue;
			try {
				if (reverifies(pr.absPath)) verifyPass++;
				else verifyFail++;
			} catch (Exception e) {
				verifyFail++;
			}
		}
		log.info("Re-verification of widening proposals: {} pass, {} fail",
				verifyPass, verifyFail);

		writeReport(files.size(), unreadable.get(), proposals, verifyPass, verifyFail);
		log.info("Report written to {}", REPORT_PATH);
	}

	// --------------------------------------------------------------------
	// Analysis
	// --------------------------------------------------------------------

	/**
	 * Analyze one scheme file. Returns a {@link Proposal} unless the file is
	 * already at the strictest applicable class (in which case the proposal
	 * is still emitted but tagged as {@code no-op} so callers can audit
	 * coverage).
	 */
	static Proposal analyze(Path p) throws IOException {
		// Resolve by content (tolerant of the 2026-06 filename rename) rather than
		// the literal path — keeps the sweep robust to relabelled schemes.
		File f = eu.solven.matmul.catalog.SchemeResolver.byHint(p.toString());
		JsonNode root = SchemeIO.parseJson(f);
		String currentField = readCurrentField(root, p);

		CoefficientProfile profile;
		if (SchemeIO.isComplex(root)) {
			ComplexNonCubicBilinearAlgorithm cx = SchemeIO.readComplex(root);
			profile = profileComplex(cx);
		} else if (SchemeIO.isNonBilinear(root)) {
			NonBilinearAlgorithm nb = SchemeIO.readNonBilinear(root);
			profile = profileNonBilinear(nb);
		} else {
			NonCubicBilinearAlgorithm alg = SchemeIO.isReduced(root)
					? SchemeIO.readReduced(root) : SchemeIO.read(root);
			profile = profileBilinear(alg);
		}

		FieldTag detected = profile.toFieldTag();
		boolean f2Tagged = SchemeIO.isZ2(root);
		boolean complexTagged = SchemeIO.isComplex(root);
		FieldTag current = FieldTag.fromText(currentField, f2Tagged, complexTagged);

		Path rel = relativizeSafe(p);
		return new Proposal(rel, p, currentField, current, detected, profile);
	}

	/**
	 * {@code Path.relativize} throws when the argument doesn't share the
	 * same root (common for temp-file fixtures in tests). Fall back to the
	 * filename in that case.
	 */
	private static Path relativizeSafe(Path p) {
		try {
			Path abs = p.toAbsolutePath().normalize();
			Path root = SCHEMES_ROOT.toAbsolutePath().normalize();
			if (abs.startsWith(root)) return root.relativize(abs);
		} catch (Exception ignored) {}
		return p.getFileName();
	}

	private static String readCurrentField(JsonNode root, Path p) {
		// Unified fields[] (task #174): the strictest characteristic-0 tag is the
		// "current" field (Z ⊂ Q ⊂ R ⊂ C); a prime-field-only scheme reports F2/F3.
		java.util.List<String> tags = SchemeIO.fieldTags(root);
		if (!tags.isEmpty()) {
			for (String t : new String[] { "Z", "Q", "R", "C", "F3", "F2" }) {
				if (tags.contains(t)) return t;
			}
		}
		if (root.has("field")) {
			JsonNode v = root.get("field");
			if (v != null && v.isTextual()) return v.asString();
		}
		// Infer from filename suffix
		String name = p.getFileName().toString();
		if (name.endsWith("_F2.json") || name.contains("_F2_")) return "F_2 (filename)";
		if (name.endsWith("_Z.json") || name.contains("_Z_")) return "Z (filename)";
		if (name.endsWith("_Q.json") || name.contains("_Q_")) return "Q (filename)";
		if (name.endsWith("_R.json") || name.contains("_R_")) return "R (filename)";
		if (name.endsWith("_C.json") || name.contains("_C_") || name.contains("xC")) return "C (filename)";
		if (name.endsWith("_ZT.json") || name.contains("_ZT_")) return "ZT (filename)";
		if (SchemeIO.isZ2(root)) return "F_2 (z2:true)";
		if (SchemeIO.isComplex(root)) return "C (complex:true)";
		return "(untagged)";
	}

	// --------------------------------------------------------------------
	// Coefficient profiling
	// --------------------------------------------------------------------

	/**
	 * Aggregated stats over all U/V/W entries: imaginary-part max, integer
	 * deviation, denominator bound, whether sign-set fits {-1,0,+1} etc.
	 */
	static final class CoefficientProfile {
		boolean hasComplex;            // any |imag| > EPS
		double maxImagAbs;             // largest |imag| seen
		double maxIntegerDeviation;    // largest |x - round(x)| over reals
		double maxAbs;                 // largest |x| seen (helps spot scale)
		int maxDenom;                  // smallest denom q such that q·x ∈ Z for all
		                              // x in coefficients (-1 if no finite q works)
		boolean allInPM1;              // all real parts ∈ {-1, 0, +1} (within EPS)
		boolean f2Compatible;          // all real parts ∈ {0, +1} (within EPS)
		boolean allRealInZ;            // all real parts integer (within EPS)
		long entriesScanned;

		CoefficientProfile() {
			maxImagAbs = 0;
			maxIntegerDeviation = 0;
			maxAbs = 0;
			maxDenom = 1;
			allInPM1 = true;
			f2Compatible = true;
			allRealInZ = true;
			entriesScanned = 0;
			hasComplex = false;
		}

		void accumulateReal(double x) {
			entriesScanned++;
			double ax = Math.abs(x);
			if (ax > maxAbs) maxAbs = ax;
			double dev = Math.abs(x - Math.round(x));
			if (dev > maxIntegerDeviation) maxIntegerDeviation = dev;
			if (dev > EPS) {
				allRealInZ = false;
				// Try to find a finite denominator q ≤ Q_MAX_DENOM s.t. q·x is integer.
				int q = smallestDenominator(x);
				if (q < 0) {
					maxDenom = -1;
				} else if (maxDenom > 0 && q > maxDenom) {
					maxDenom = q;
				}
			}
			// PM1 / F2 check: nearest integer must be in {-1, 0, +1} / {0, 1}.
			long rounded = Math.round(x);
			if (dev > EPS || rounded < -1 || rounded > 1) {
				allInPM1 = false;
				f2Compatible = false;
			} else {
				if (rounded < 0) f2Compatible = false;
			}
		}

		void accumulateComplex(double re, double im) {
			accumulateReal(re);
			double aim = Math.abs(im);
			if (aim > EPS) hasComplex = true;
			if (aim > maxImagAbs) maxImagAbs = aim;
		}

		FieldTag toFieldTag() {
			if (hasComplex) return FieldTag.C;
			if (allInPM1) {
				return f2Compatible ? FieldTag.F2_PM1 : FieldTag.PM1;
			}
			if (allRealInZ) return FieldTag.Z;
			if (maxDenom > 0) return FieldTag.Q;
			return FieldTag.R;
		}
	}

	/**
	 * Return the smallest positive integer q ≤ {@link #Q_MAX_DENOM} such
	 * that q·x is within EPS of an integer; -1 if no such q exists in range
	 * (treated as "irrational").
	 */
	static int smallestDenominator(double x) {
		for (int q = 1; q <= Q_MAX_DENOM; q++) {
			double qx = q * x;
			if (Math.abs(qx - Math.round(qx)) < EPS) return q;
		}
		return -1;
	}

	private static CoefficientProfile profileBilinear(NonCubicBilinearAlgorithm alg) {
		CoefficientProfile pr = new CoefficientProfile();
		accumulateMatrix(pr, alg.denseU());
		accumulateMatrix(pr, alg.denseV());
		accumulateMatrix(pr, alg.denseW());
		return pr;
	}

	private static CoefficientProfile profileComplex(ComplexNonCubicBilinearAlgorithm cx) {
		CoefficientProfile pr = new CoefficientProfile();
		accumulateComplexMatrix(pr, cx.uRe, cx.uIm);
		accumulateComplexMatrix(pr, cx.vRe, cx.vIm);
		accumulateComplexMatrix(pr, cx.wRe, cx.wIm);
		return pr;
	}

	private static CoefficientProfile profileNonBilinear(NonBilinearAlgorithm nb) {
		CoefficientProfile pr = new CoefficientProfile();
		accumulateMatrix(pr, nb.Ua);
		accumulateMatrix(pr, nb.Ub);
		accumulateMatrix(pr, nb.Va);
		accumulateMatrix(pr, nb.Vb);
		accumulateMatrix(pr, nb.W);
		return pr;
	}

	private static void accumulateMatrix(CoefficientProfile pr, double[][] m) {
		for (double[] row : m) {
			for (double x : row) pr.accumulateReal(x);
		}
	}

	private static void accumulateComplexMatrix(CoefficientProfile pr,
			double[][] re, double[][] im) {
		for (int i = 0; i < re.length; i++) {
			for (int j = 0; j < re[i].length; j++) {
				pr.accumulateComplex(re[i][j], im[i][j]);
			}
		}
	}

	// --------------------------------------------------------------------
	// Field-tag lattice
	// --------------------------------------------------------------------

	enum FieldTag {
		F2_PM1(0, "F_2 (subset of Z, {0,1})"),
		PM1(1, "Z (only {-1,0,+1})"),
		F2(1, "F_2"),
		Z(2, "Z"),
		Q(3, "Q"),
		R(4, "R"),
		C(5, "C"),
		UNKNOWN(99, "unknown");

		final int permissiveness;
		final String label;

		FieldTag(int permissiveness, String label) {
			this.permissiveness = permissiveness;
			this.label = label;
		}

		/**
		 * Parse a free-form field tag string (from JSON {@code "field"} or
		 * filename suffix) into the most natural enum value.
		 */
		static FieldTag fromText(String s, boolean z2Tagged, boolean complexTagged) {
			if (z2Tagged) return F2;
			if (complexTagged) return C;
			if (s == null) return UNKNOWN;
			String t = s.trim().toLowerCase(java.util.Locale.ROOT);
			if (t.startsWith("f_2") || t.startsWith("f2") || t.contains("z_2")) return F2;
			if (t.startsWith("0.5*c") || t.endsWith("c") || t.equals("c")
					|| t.contains("complex")) {
				return C;
			}
			if (t.startsWith("0.5*z") || t.startsWith("z")
					|| t.startsWith("zt")) return Z;
			if (t.startsWith("q")) return Q;
			if (t.startsWith("r")) return R;
			return UNKNOWN;
		}
	}

	// --------------------------------------------------------------------
	// Proposals
	// --------------------------------------------------------------------

	static final class Proposal {
		final Path relPath;
		final Path absPath;
		final String currentText;
		final FieldTag current;
		final FieldTag detected;
		final CoefficientProfile profile;

		Proposal(Path relPath, Path absPath, String currentText,
				FieldTag current, FieldTag detected, CoefficientProfile profile) {
			this.relPath = relPath;
			this.absPath = absPath;
			this.currentText = currentText;
			this.current = current;
			this.detected = detected;
			this.profile = profile;
		}

		boolean isWidening() {
			// A widening (= narrower-field proposal) is when the detected tag
			// is strictly more restrictive than the current. Exception: don't
			// flag F2-compatible PM1 as a "widening" over an existing F2 tag
			// — the F2 tag is already the strictest.
			if (current == FieldTag.UNKNOWN) return true;
			if (current == FieldTag.F2 && detected == FieldTag.F2_PM1) return false;
			return detected.permissiveness < current.permissiveness;
		}

		String transitionKey() {
			return current.name() + " -> " + detected.name();
		}
	}

	// --------------------------------------------------------------------
	// Re-verification (after a hypothetical narrowing)
	// --------------------------------------------------------------------

	static boolean reverifies(Path p) throws IOException {
		File f = p.toFile();
		JsonNode root = SchemeIO.parseJson(f);
		if (SchemeIO.isComplex(root)) {
			ComplexNonCubicBilinearAlgorithm cx = SchemeIO.readComplex(root);
			return Verifier.isExactComplex(cx);
		}
		if (SchemeIO.isNonBilinear(root)) {
			NonBilinearAlgorithm nb = SchemeIO.readNonBilinear(root);
			return Verifier.isExactNonBilinear(nb);
		}
		NonCubicBilinearAlgorithm alg = SchemeIO.isReduced(root)
				? SchemeIO.readReduced(root) : SchemeIO.read(root);
		if (SchemeIO.isZ2(root)) return Verifier.isExactNonCubicF2(alg);
		return Verifier.isExactNonCubic(alg);
	}

	// --------------------------------------------------------------------
	// Report writing
	// --------------------------------------------------------------------

	private static void writeReport(int totalScanned, int unreadable,
			List<Proposal> proposals, int verifyPass, int verifyFail) throws IOException {
		Files.createDirectories(REPORT_PATH.getParent());

		// Summary by transition.
		Map<String, Integer> transitionCounts = new TreeMap<>();
		int wideningCount = 0;
		for (Proposal p : proposals) {
			if (!p.isWidening()) continue;
			wideningCount++;
			transitionCounts.merge(p.transitionKey(), 1, Integer::sum);
		}

		List<Proposal> widenings = new ArrayList<>();
		for (Proposal p : proposals) if (p.isWidening()) widenings.add(p);
		widenings.sort(Comparator.comparing((Proposal x) -> x.transitionKey())
				.thenComparing(x -> x.relPath.toString()));

		try (BufferedWriter bw = new BufferedWriter(new FileWriter(REPORT_PATH.toFile()));
				PrintWriter out = new PrintWriter(bw)) {
			out.println("# Field-widening proposals");
			out.println();
			out.println("Auto-generated by `FieldWideningSweep`. Do NOT hand-edit —");
			out.println("re-run the driver after the catalog changes.");
			out.println();
			out.println("Each row proposes a stricter field tag than the one currently");
			out.println("stored in the JSON. The detection is purely coefficient-based:");
			out.println("we never apply a wider tag than what the U/V/W entries warrant.");
			out.println();
			out.println("## Summary");
			out.println();
			out.printf("- Scanned: **%d** scheme files%n", totalScanned);
			out.printf("- Unreadable: **%d**%n", unreadable);
			out.printf("- Widening proposals: **%d**%n", wideningCount);
			out.printf("- Re-verification of widening proposals: **%d** pass, **%d** fail%n",
					verifyPass, verifyFail);
			out.println();
			out.println("### Transitions (current -> proposed)");
			out.println();
			out.println("| Transition | Count |");
			out.println("| --- | ---: |");
			for (Map.Entry<String, Integer> e : transitionCounts.entrySet()) {
				out.printf("| %s | %d |%n", e.getKey(), e.getValue());
			}
			out.println();

			out.println("## Detected widenings");
			out.println();
			out.println("| File | Current | Proposed | maxAbs | maxIntDev | maxDenom | maxImag |");
			out.println("| --- | --- | --- | ---: | ---: | ---: | ---: |");
			for (Proposal p : widenings) {
				out.printf("| `%s` | %s | %s | %s | %s | %s | %s |%n",
						p.relPath.toString().replace('\\', '/'),
						escape(p.current.name() + " (" + p.currentText + ")"),
						escape(p.detected.name()),
						fmt(p.profile.maxAbs),
						fmt(p.profile.maxIntegerDeviation),
						p.profile.maxDenom < 0 ? "irr" : Integer.toString(p.profile.maxDenom),
						fmt(p.profile.maxImagAbs));
			}
		}
	}

	private static String escape(String s) {
		return s.replace("|", "\\|");
	}

	private static String fmt(double x) {
		if (x == 0.0) return "0";
		if (x < 1e-3) return String.format(java.util.Locale.ROOT, "%.2e", x);
		return String.format(java.util.Locale.ROOT, "%.4g", x);
	}
}
