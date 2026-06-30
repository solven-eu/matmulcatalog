package eu.solven.matmul.docs.verify;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.EnumSet;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.Semaphore;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

import eu.solven.matmul.ComplexNonCubicBilinearAlgorithm;
import eu.solven.matmul.NonBilinearAlgorithm;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.verifiers.LineageVerifier;
import eu.solven.matmul.util.ProgressMonitor;
import lombok.extern.slf4j.Slf4j;
import tools.jackson.databind.JsonNode;

/**
 * Generic catalog-integrity tool: walk every {@code *.json} under
 * {@code src/main/resources/schemes/}, load it via {@link SchemeIO},
 * verify the matmul tensor identity via {@link Verifier}, and report
 * a structured summary.
 *
 * <p>Exit code: {@code 0} if every scheme verifies, {@code 1}
 * otherwise — so CI can fail the workflow on any broken scheme.</p>
 *
 * <p>Concurrency: a fixed pool of {@code availableProcessors - 1}
 * worker threads. Progress is logged every 250 files with a
 * <strong>shape-weighted ETA</strong>: heavier shapes ({@code ⟨24,24,24⟩
 * = r·(n·m + m·p + n·p)}) count proportionally more than tiny ones, so
 * the remaining-time estimate doesn't collapse when the easy files clear
 * fast and the hard ones drag.</p>
 *
 * <p>Run locally:</p>
 * <pre>
 *   mvn -q exec:java -Dexec.mainClass=eu.solven.matmul.docs.verify.VerifyAllSchemes
 * </pre>
 *
 * <p>Or from CI: see {@code .github/workflows/verify-catalog.yml}.</p>
 */
@Slf4j
public final class VerifyAllSchemes {

	private static final Path SCHEMES_ROOT = Path.of("src/main/resources/schemes");

	/** Weight (= r·(nm+mp+np)) above which a scheme densifies to roughly &ge;30 MB of
	 *  factor matrices and is throttled to bounded concurrency to cap peak heap. */
	private static final long HEAVY_WEIGHT = 10_000_000L;

	/** Compositional verifier for lineage-only STUB schemes (maxDim &gt; materialise
	 *  cap): they carry no explicit U/V/W matrices, so they're verified by exact-
	 *  checking their primitive leaves + trusting the operators — never expanded.
	 *  Shared (thread-safe: concurrent atom cache + read-only lookup). */
	private static final LineageVerifier STUB_VERIFIER =
			new LineageVerifier(new FieldAwareLookup(Field.Q));

	/** Matches {@code _{n}x{m}x{p}_r{rank}} in scheme filenames. */
	private static final Pattern SHAPE = Pattern.compile("[_-](\\d+)x(\\d+)x(\\d+)_(?:r|m)(\\d+)");

	private VerifyAllSchemes() {}

	public static void main(String[] args) throws Exception {
		List<Path> files;
		try (Stream<Path> walk = Files.walk(SCHEMES_ROOT)) {
			files = walk.filter(p -> p.toString().endsWith(".json"))
					.sorted()
					.toList();
		}

		// Per-file weight from filename: r·(n·m + m·p + n·p) ≈ cost of one
		// random-matmul spot check on that shape. Drives the ETA so it
		// stays accurate when the easy half clears quickly and the
		// remaining files are 30× heavier.
		long[] weights = new long[files.size()];
		long totalWeight = 0L;
		for (int i = 0; i < files.size(); i++) {
			weights[i] = weightOf(files.get(i));
			totalWeight += weights[i];
		}
		log.info("Verifying {} scheme files ({} weight units) under {}",
				files.size(), totalWeight, SCHEMES_ROOT);

		int nThreads = Math.max(1, Runtime.getRuntime().availableProcessors() - 1);
		ExecutorService pool = Executors.newFixedThreadPool(nThreads);
		log.info("Parallel verify on {} threads", nThreads);

		AtomicInteger verified = new AtomicInteger();
		AtomicInteger broken = new AtomicInteger();
		AtomicInteger unreadable = new AtomicInteger();
		List<String> brokenList = new ArrayList<>();
		long t0 = System.nanoTime();
		int total = files.size();

		// Shape-weighted ETA + heap stats via the shared ProgressMonitor.
		long[] weightsRef = weights;
		ProgressMonitor monitor = ProgressMonitor.builder()
				.total(total)
				.weight(i -> weightsRef[(int) i])
				.logEveryMillis(2_000)
				.logger(log)
				.label("verify-catalog")
				.build();
		monitor.start();

		// Bound peak memory. A big dense scheme (SZ n≥24 cubes, FMM cubes) densifies
		// to hundreds of MB of double[][]; with all 15 worker threads loading one at
		// once a 5 GB heap GC-thrashes then OOMs. So THROTTLE the heavy schemes to a
		// few concurrent (HEAVY_WEIGHT ≈ ≥~30 MB of factor matrices) while the small
		// ones — the vast majority — run unthrottled. Stubs cost ~nothing (compositional)
		// and never trip the gate.
		int heavyConcurrency = Math.max(1, nThreads / 4);
		Semaphore heavyGate = new Semaphore(heavyConcurrency);
		log.info("Heavy-scheme throttle: {} concurrent above weight {}", heavyConcurrency, HEAVY_WEIGHT);

		List<Future<?>> futures = new ArrayList<>(total);
		for (int i = 0; i < files.size(); i++) {
			final Path p = files.get(i);
			final int idx = i;
			final boolean heavy = weights[i] > HEAVY_WEIGHT;
			futures.add(pool.submit(() -> {
				try {
					if (heavy) heavyGate.acquireUninterruptibly();
					try {
						verifyOne(p, brokenList, verified, broken, unreadable);
					} finally {
						if (heavy) heavyGate.release();
					}
				} finally {
					monitor.tick(idx);
				}
			}));
		}
		for (Future<?> f : futures) {
			try { f.get(); } catch (Exception ignored) {}
		}
		pool.shutdown();
		monitor.done();

		long elapsedMs = (System.nanoTime() - t0) / 1_000_000L;
		brokenList.sort(Comparator.naturalOrder());
		for (String line : brokenList) log.warn("  {}", line);

		log.info("");
		log.info("Catalog verification complete: {} verified, {} broken, {} unreadable, {}ms",
				verified.get(), broken.get(), unreadable.get(), elapsedMs);

		if (broken.get() > 0 || unreadable.get() > 0) {
			System.exit(1);
		}
	}

	private static void verifyOne(Path p, List<String> brokenList,
			AtomicInteger verified, AtomicInteger broken, AtomicInteger unreadable) {
		File f = p.toFile();
		// Load the file ONCE into a JsonNode; the kind-checks and readers all
		// accept JsonNode so dispatch costs zero re-reads. Previous version
		// re-read the file 2–3× per scheme (isNonBilinear, isComplex,
		// isReduced, then the matching read*), which dominated wall-clock
		// on large catalogs.
		JsonNode root;
		try {
			root = SchemeIO.parseJson(f);
		} catch (Exception e) {
			unreadable.incrementAndGet();
			synchronized (brokenList) {
				brokenList.add(String.format("UNREADABLE  %s  (parse: %s: %s)",
						SCHEMES_ROOT.relativize(p),
						e.getClass().getSimpleName(),
						truncate(String.valueOf(e.getMessage()), 80)));
			}
			return;
		}
		try {
			// Lineage-only STUB (maxDim > materialise cap): no explicit matrices to
			// read — verify compositionally (exact leaves + trusted operators).
			if (SchemeIO.isStub(root)) {
				java.util.Optional<Lineage.Node> ln = SchemeIO.readLineage(root);
				LineageVerifier.Result vr = ln.isPresent()
						? STUB_VERIFIER.verify(ln.get())
						: LineageVerifier.Result.fail("stub missing lineage");
				record(vr.certified(), verified, broken, brokenList,
						String.format("BROKEN  %s  (stub: %s)",
								SCHEMES_ROOT.relativize(p), vr.detail()));
				return;
			}
			if (SchemeIO.isNonBilinear(root)) {
				NonBilinearAlgorithm nb = SchemeIO.readNonBilinear(root);
				boolean ok = Verifier.passesRandomMatmulSpotCheckNB(nb);
				record(ok, verified, broken, brokenList,
						String.format("BROKEN  %s  ⟨%d,%d,%d⟩=r%d  (non-bilinear)",
								SCHEMES_ROOT.relativize(p), nb.n, nb.m, nb.p, nb.r));
				return;
			}
			if (SchemeIO.isComplex(root)) {
				ComplexNonCubicBilinearAlgorithm cx = SchemeIO.readComplex(root);
				boolean ok = Verifier.isExactComplex(cx);
				record(ok, verified, broken, brokenList,
						String.format("BROKEN  %s  ⟨%d,%d,%d⟩=r%d  (complex)",
								SCHEMES_ROOT.relativize(p), cx.n, cx.m, cx.p, cx.r));
				return;
			}
			NonCubicBilinearAlgorithm alg = SchemeIO.isReduced(root)
					? SchemeIO.readReduced(root) : SchemeIO.read(root);
			// Verify the matmul-tensor identity over EVERY algebra the scheme
			// actually claims (see CheckField / fieldsToCheck). The scheme is
			// verified iff it holds over all of them. This replaces an ad-hoc
			// "F2 ? modular : real" fork that silently mis-routed F3-only schemes
			// to the real check (→ spurious BROKEN).
			EnumSet<CheckField> toCheck = fieldsToCheck(root);
			boolean ok = toCheck.stream().allMatch(cf -> cf.holds(alg));
			record(ok, verified, broken, brokenList,
					String.format("BROKEN  %s  ⟨%d,%d,%d⟩=r%d  (%s)",
							SCHEMES_ROOT.relativize(p), alg.n, alg.m, alg.p, alg.r,
							fieldLabel(toCheck)));
		} catch (Exception e) {
			unreadable.incrementAndGet();
			synchronized (brokenList) {
				brokenList.add(String.format("UNREADABLE  %s  (%s: %s)",
						SCHEMES_ROOT.relativize(p),
						e.getClass().getSimpleName(),
						truncate(String.valueOf(e.getMessage()), 80)));
			}
		}
	}

	/**
	 * An algebra over which a bilinear scheme's matmul-tensor identity is checked.
	 * These are distinct <em>characteristics</em>, not the catalog's field tags: an
	 * integer/rational (characteristic-0) identity reduces mod p (the reduction is
	 * a ring homomorphism), so a single {@link #CHAR0} check certifies every char-0
	 * field <em>and</em> its finite reductions; only a finite-field-ONLY scheme
	 * needs its own modular check.
	 *
	 * <p>{@link #CHAR0} uses a random real-valued spot check (cheap, and a wrong
	 * integer identity fails it with probability ≈ 1); {@link #F2}/{@link #F3} use
	 * exact modular arithmetic (a real-valued check would fail spuriously since,
	 * e.g., {@code 1+1=0} in F₂ but {@code ≠0} in ℝ).</p>
	 */
	private enum CheckField {
		CHAR0("char-0") {
			@Override boolean holds(NonCubicBilinearAlgorithm a) {
				return Verifier.passesRandomMatmulSpotCheck(a);
			}
		},
		F2("F2") {
			@Override boolean holds(NonCubicBilinearAlgorithm a) {
				return Verifier.residualNonCubicF2(a) == 0;
			}
		},
		F3("F3") {
			@Override boolean holds(NonCubicBilinearAlgorithm a) {
				return Verifier.residualNonCubicF3(a) == 0;
			}
		};

		private final String label;
		CheckField(String label) { this.label = label; }

		abstract boolean holds(NonCubicBilinearAlgorithm alg);
	}

	/**
	 * The minimal set of algebras to actually verify for a bilinear scheme, from
	 * its declared {@code fields[]}. Any characteristic-0 tag (Z/Q/R/C) collapses
	 * to a single {@link CheckField#CHAR0} check (it subsumes the super-fields and
	 * the finite reductions); a finite-field-only scheme is checked in each finite
	 * field it claims. Falls back to the filename F₂ signal when {@code fields[]}
	 * is absent (legacy / synthetic names).
	 */
	private static EnumSet<CheckField> fieldsToCheck(JsonNode root) {
		List<String> tags = SchemeIO.fieldTags(root);
		boolean char0 = tags.contains("Z") || tags.contains("Q")
				|| tags.contains("R") || tags.contains("C");
		if (char0) {
			return EnumSet.of(CheckField.CHAR0);
		}
		EnumSet<CheckField> out = EnumSet.noneOf(CheckField.class);
		if (tags.contains("F2")) out.add(CheckField.F2);
		if (tags.contains("F3")) out.add(CheckField.F3);
		if (out.isEmpty()) {
			out.add(SchemeIO.isZ2(root) ? CheckField.F2 : CheckField.CHAR0);
		}
		return out;
	}

	/** Human-readable label of the algebras checked, for the BROKEN report line. */
	private static String fieldLabel(EnumSet<CheckField> fs) {
		return fs.stream().map(c -> c.label).collect(java.util.stream.Collectors.joining(","));
	}

	private static void record(boolean ok, AtomicInteger verified, AtomicInteger broken,
			List<String> brokenList, String brokenMsg) {
		if (ok) {
			verified.incrementAndGet();
		} else {
			broken.incrementAndGet();
			synchronized (brokenList) {
				brokenList.add(brokenMsg);
			}
		}
	}

	private static String truncate(String s, int n) {
		if (s == null) return "";
		return s.length() <= n ? s : s.substring(0, n) + "…";
	}

	/**
	 * Per-file work-weight from the filename: {@code r·(n·m + m·p + n·p)},
	 * matching the dominant cost of {@code Verifier.passesRandomMatmulSpotCheck}
	 * (three matrix scans). Files without a parseable shape fall back to 1.
	 */
	private static long weightOf(Path p) {
		Matcher m = SHAPE.matcher(p.getFileName().toString());
		if (!m.find()) return 1L;
		long n = Long.parseLong(m.group(1));
		long mm = Long.parseLong(m.group(2));
		long pp = Long.parseLong(m.group(3));
		long r = Long.parseLong(m.group(4));
		return r * (n * mm + mm * pp + n * pp);
	}
}
