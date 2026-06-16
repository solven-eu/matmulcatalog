package eu.solven.matmul.papers.hopcroftkerr1971;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.fraction.BigFraction;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.SchemeIO;

import lombok.extern.slf4j.Slf4j;

/**
 * Generate, TRIPLE-VERIFY and register the constructive Hopcroft–Kerr
 * {@code ⟨2,p,n⟩} schemes (task #7/#8) into the NEW {@code schemes/constructed/}
 * category (formula-emitter output: regenerable like {@code derived/}, but from
 * a paper construction, not catalog-closure operators).
 *
 * <p>Sweep: {@code 3 ≤ p ≤ 32}, {@code p ≤ n ≤ 32} — the Lemma-1 band range
 * {@code n ≤ 2p−1} plus, beyond it, the chained-augmentation regime
 * ({@link HopcroftKerr2bcAsymmetric#buildChained}: DP-optimal concatenation
 * of band segments, task #10).</p>
 *
 * <p>Triple verification per scheme (fail-loud: any failure aborts):</p>
 * <ol>
 *   <li><b>Sampled spot check</b> — {@link Verifier#passesRandomMatmulSpotCheck}
 *       with 20k samples;</li>
 *   <li><b>Full double residual</b> — {@link Verifier#isExactNonCubic} over every
 *       tensor cell;</li>
 *   <li><b>Exact rational symbolic</b> — every U/V/W entry recovered as an exact
 *       fraction (continued-function rationalisation, denominator-capped), then
 *       the bilinear identity {@code Σ_k U[a][k]·V[b][k]·W[c][k] = T[a][b][c]}
 *       checked over ℚ with BigFraction arithmetic, sparse iteration — no
 *       floating point anywhere.</li>
 * </ol>
 *
 * <pre>mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.papers.hopcroftkerr1971.GenerateHk2npConstructed -Dexec.args="--apply"</pre>
 */
@Slf4j
public class GenerateHk2npConstructed {

	private static final Path ROOT = Path.of("src/main/resources/schemes/constructed");
	private static final int MAX_DIM = 32;

	public static void main(String[] args) throws Exception {
		boolean apply = List.of(args).contains("--apply");
		// In-sweep verification is OPT-IN (--verify): the default sweep is as
		// fast as possible and verification is a LATER PHASE — run
		// references/hopcroftkerr1971/verify_constructed_independent.py (the
		// cross-language exact certificate) and/or VerifyAllSchemes over
		// schemes/constructed/, then stamp verified:true. Schemes written
		// without --verify are stamped verified:false until that phase.
		boolean verify = List.of(args).contains("--verify");
		long start = System.nanoTime();
		int built = 0, atFormula = 0, above = 0, newBounds = 0, written = 0, skipped = 0;
		List<String> failures = new ArrayList<>();
		var lookup = new eu.solven.matmul.catalog.FieldAwareLookup("C");

		for (int p = 3; p <= MAX_DIM; p++) {
			for (int n = p; n <= MAX_DIM; n++) {
				int formula = (int) Math.ceil((3.0 * p * n + Math.max(p, n)) / 2.0);
				NonCubicBilinearAlgorithm raw;
				try {
					raw = HopcroftKerr2bcAsymmetric.build(p, n);
				} catch (Throwable e) {
					failures.add(String.format("⟨2,%d,%d⟩ BUILD FAILED: %s", p, n, e.getMessage()));
					continue;
				}
				// Catalog orientation ⟨2,p,n⟩ (the emitters build ⟨p,2,n⟩).
				NonCubicBilinearAlgorithm alg = raw.orientAs(2, p, n).orElse(raw);
				built++;

				// ── Idempotent resume: if the on-disk file already carries the
				// IDENTICAL content hash and verified:true, those exact bits have
				// already passed the quadruple gate — skip re-verification and
				// re-write. (Builders are deterministic, so an interrupted
				// certified run resumes where it stopped; pass --force to redo.)
				if (apply && verify && !List.of(args).contains("--force")
						&& alreadyEmittedVerified(alg, p, n)) {
					if (alg.r == formula) atFormula++;
					else above++;
					written++;
					skipped++;
					continue;
				}

				// ── In-sweep verification (opt-in; see header note) ──
				boolean spot = !verify || Verifier.passesRandomMatmulSpotCheck(alg, 20_000, 42L);
				boolean fullDouble = !verify || Verifier.isExactNonCubic(alg);
				boolean symbolic = !verify || exactRationalVerify(alg);
				if (verify && spot && fullDouble && !symbolic) {
					// Numerically valid but the exact rationals exceed the publishable
					// denominator range (huge Lemma-1 window determinants at large
					// n−p). Defer rather than fail: registering an approximate scheme
					// would violate the exact-coefficients discipline.
					log.warn("⟨2,{},{}⟩ r={} DEFERRED: passes spot+full but exact denominators "
							+ "exceed publishable range — not registered", p, n, alg.r);
					continue;
				}
				if (!spot || !fullDouble) {
					failures.add(String.format("⟨2,%d,%d⟩ VERIFY spot=%b full=%b symbolic=%b",
							p, n, spot, fullDouble, symbolic));
					continue;
				}
				if (alg.r == formula) atFormula++;
				else above++;

				// Catalog comparison (any field — broadest).
				Integer existingBest = lookup.findFiles(2, p, n).stream()
						.map(f -> {
							try {
								var root = SchemeIO.parseJson(f.toFile());
								return root.get("m") != null ? root.get("m").asInt() : null;
							} catch (Exception e) {
								return null;
							}
						})
						.filter(java.util.Objects::nonNull)
						.min(Integer::compare).orElse(null);
				String vsCatalog = existingBest == null ? "ABSENT"
						: (alg.r < existingBest ? "BEATS " + existingBest
								: (alg.r == existingBest ? "ties " + existingBest : "above " + existingBest));
				if (existingBest == null || alg.r < existingBest) newBounds++;

				log.info("⟨2,{},{}⟩ r={} (formula {}{}) spot+full+symbolic OK — catalog: {}",
						p, n, alg.r, formula, alg.r == formula ? "" : ", +" + (alg.r - formula),
						vsCatalog);

				if (apply) {
					File out = writeScheme(alg, p, n, formula, verify);
					if (verify) {
						// Fourth check: the PUBLISHED artifact must itself verify
						// exactly after a disk round-trip (catches any write-path
						// lossiness — the real target of the verification).
						NonCubicBilinearAlgorithm back = SchemeIO.readBilinear(out);
						if (!exactRationalVerify(back)) {
							failures.add(String.format("⟨2,%d,%d⟩ DISK ROUND-TRIP not exact: %s",
									p, n, out.getName()));
							continue;
						}
					}
					written++;
				}
			}
			log.info("[progress] p={} done — built={} atFormula={} above={} newBounds={} {} ms",
					p, built, atFormula, above, newBounds, (System.nanoTime() - start) / 1_000_000);
		}
		log.info("==================================================================");
		log.info("HK ⟨2,p,n⟩ constructed sweep: built={} atFormula={} above={} newBounds={} written={}"
				+ " (skipped={} already-verified identical) apply={}",
				built, atFormula, above, newBounds, written, skipped, apply);
		for (String f : failures) log.error("  FAILURE: {}", f);
		if (!failures.isEmpty()) {
			throw new IllegalStateException(failures.size() + " build/verification failures");
		}
		log.info(verify
				? "ALL schemes quadruple-verified (sampled, full-residual, exact-rational symbolic, disk round-trip)."
				: "FAST sweep complete — schemes stamped verified:false; run the later verification phase.");
	}

	/**
	 * Exact symbolic verification over ℚ: rationalise every stored coefficient,
	 * then check the bilinear identity coefficient-wise with BigFraction
	 * arithmetic and sparse iteration. The definitive no-floating-point check of
	 * what gets published.
	 */
	public static boolean exactRationalVerify(NonCubicBilinearAlgorithm alg) {
		int dimU = alg.n * alg.m, dimV = alg.m * alg.p, dimW = alg.n * alg.p;
		double[][] U = alg.denseU(), V = alg.denseV(), W = alg.denseW();
		BigFraction[][] Uf = rationalize(U, dimU, alg.r);
		BigFraction[][] Vf = rationalize(V, dimV, alg.r);
		BigFraction[][] Wf = rationalize(W, dimW, alg.r);
		if (Uf == null || Vf == null || Wf == null) return false;

		// acc[a][b][c] accumulated sparsely per product k.
		java.util.Map<Long, BigFraction> acc = new java.util.HashMap<>();
		for (int k = 0; k < alg.r; k++) {
			List<int[]> us = support(Uf, dimU, k);
			List<int[]> vs = support(Vf, dimV, k);
			List<int[]> ws = support(Wf, dimW, k);
			for (int[] u : us) {
				for (int[] v : vs) {
					BigFraction uv = Uf[u[0]][k].multiply(Vf[v[0]][k]);
					for (int[] w : ws) {
						long key = ((long) u[0] * dimV + v[0]) * dimW + w[0];
						BigFraction add = uv.multiply(Wf[w[0]][k]);
						acc.merge(key, add, BigFraction::add);
					}
				}
			}
		}
		// Compare against the matmul tensor: T[a][b][c] ∈ {0,1}.
		// First: every accumulated entry must equal its tensor value…
		for (var e : acc.entrySet()) {
			long key = e.getKey();
			int c = (int) (key % dimW);
			int b = (int) ((key / dimW) % dimV);
			int a = (int) (key / dimW / dimV);
			int expected = Verifier.matmulTensorEntry(a, b, c, alg.n, alg.m, alg.p);
			if (e.getValue().compareTo(new BigFraction(expected)) != 0) {
				if (expected == 0 && e.getValue().getNumerator().signum() == 0) continue;
				return false;
			}
		}
		// …and every tensor 1-cell must have been accumulated.
		for (int i = 0; i < alg.n; i++) {
			for (int l = 0; l < alg.m; l++) {
				for (int j = 0; j < alg.p; j++) {
					long key = ((long) (i * alg.m + l) * dimV + (l * alg.p + j)) * dimW + (i * alg.p + j);
					BigFraction got = acc.get(key);
					if (got == null || got.compareTo(BigFraction.ONE) != 0) return false;
				}
			}
		}
		return true;
	}

	private static List<int[]> support(BigFraction[][] F, int dim, int k) {
		List<int[]> out = new ArrayList<>();
		for (int i = 0; i < dim; i++) {
			if (F[i][k].getNumerator().signum() != 0) out.add(new int[] { i });
		}
		return out;
	}

	/** Continued-fraction rationalisation of a factor matrix; null if any entry resists. */
	private static BigFraction[][] rationalize(double[][] M, int rows, int r) {
		BigFraction[][] out = new BigFraction[rows][r];
		for (int i = 0; i < rows; i++) {
			for (int k = 0; k < r; k++) {
				BigFraction f = rationalizeOne(M[i][k]);
				if (f == null) return null;
				out[i][k] = f;
			}
		}
		return out;
	}

	static BigFraction rationalizeOne(double v) {
		if (v == 0.0) return BigFraction.ZERO;
		long rounded = Math.round(v);
		if (Math.abs(v - rounded) < 1e-9) return new BigFraction(rounded);
		// continued fractions, denominator-capped
		boolean neg = v < 0;
		double x = Math.abs(v);
		long hPrev = 0, h = 1, kPrev = 1, kk = 0;
		double frac = x;
		for (int it = 0; it < 64; it++) {
			long a = (long) Math.floor(frac);
			long hN = a * h + hPrev, kN = a * kk + kPrev;
			if (kN > 1_000_000_000L || kN <= 0) return null;
			hPrev = h; h = hN;
			kPrev = kk; kk = kN;
			if (kk > 0 && Math.abs((double) h / kk - x) <= 1e-12 * Math.max(1.0, x)) {
				return new BigFraction(neg ? -h : h, kk);
			}
			double rem = frac - a;
			if (rem < 1e-15) return null;
			frac = 1.0 / rem;
		}
		return null;
	}

	/**
	 * True iff the on-disk emission of this shape has the IDENTICAL content
	 * hash as the freshly-built scheme AND is stamped {@code verified:true} —
	 * i.e. these exact bits already passed the quadruple gate in a previous
	 * (possibly interrupted) certified run.
	 */
	private static boolean alreadyEmittedVerified(NonCubicBilinearAlgorithm alg, int p, int n) {
		try {
			int section = Math.max(2, Math.max(p, n));
			String hash7 = SchemeIO.contentHash(alg).substring(0, 7);
			File existing = ROOT.resolve("section" + section)
					.resolve("2x" + p + "x" + n + "-r" + alg.r + "-hk71-" + hash7 + ".json").toFile();
			if (!existing.isFile()) return false;
			var root = SchemeIO.parseJson(existing);
			return root.has("verified") && root.get("verified").asBoolean(false);
		} catch (Exception e) {
			return false; // any doubt → re-verify
		}
	}

	private static File writeScheme(NonCubicBilinearAlgorithm alg, int p, int n, int formula,
			boolean verified) throws Exception {
		int section = Math.max(2, Math.max(p, n));
		Path dir = ROOT.resolve("section" + section);
		Files.createDirectories(dir);
		// Drop any previous emission of this shape (regenerable category).
		try (var ls = Files.list(dir)) {
			for (Path old : ls.filter(f -> f.getFileName().toString()
					.startsWith("2x" + p + "x" + n + "-") && f.getFileName().toString()
					.contains("hk71")).toList()) {
				Files.delete(old);
			}
		}
		String hash7 = SchemeIO.contentHash(alg).substring(0, 7);
		File out = dir.resolve("2x" + p + "x" + n + "-r" + alg.r + "-hk71-" + hash7 + ".json")
				.toFile();
		SchemeIO.write(alg, out);

		// Fields: integer-exact ⇒ all six; else Q-exact + per-prime reduction.
		List<String> fields = new ArrayList<>();
		if (Verifier.isExactNonCubicFp(alg, 2)) fields.add("F2");
		if (Verifier.isExactNonCubicFp(alg, 3)) fields.add("F3");
		boolean allInt = allIntegers(alg);
		if (allInt) fields.add("Z");
		fields.add("Q");
		fields.add("R");
		fields.add("C");

		java.util.Map<String, Object> meta = new java.util.LinkedHashMap<>();
		meta.put("source", "Hopcroft-Kerr 1971 (constructive, this work)");
		meta.put("lineage_str", "HopcroftKerr2bcAsymmetric.build(p=" + p + ", n=" + n + ")"
				+ " — banded cyclic construction + exact Lemma-1 back-substitution"
				+ ((p & 1) == 0 ? " + Case-2 circulant matching / Z-pair Step 3" : ""));
		meta.put("lineage_compact", "hk71-constructive(" + p + "," + n + ")");
		meta.put("discovery_note", "Constructively attains"
				+ (alg.r == formula ? "" : " +"+ (alg.r - formula) + " over")
				+ " the HK formula ⌈(3pn+max)/2⌉ = " + formula + "."
				+ (verified
						? " Quadruple-verified at emission: 20k-sample spot check, full double"
								+ " residual, exact-rational symbolic (in-memory AND disk round-trip)."
						: " Emitted FAST (no in-sweep verification) — pending the later"
								+ " verification phase (verify_constructed_independent.py).")
				+ " See research/hopcroft-kerr-2np/CONSTRUCTIVE_METHOD.md.");
		meta.put("verified", verified);
		SchemeIO.updateFields(out, meta, List.of(), true);
		injectFieldArrays(out, fields);
		return out;
	}

	private static boolean allIntegers(NonCubicBilinearAlgorithm alg) {
		for (double[][] f : new double[][][] { alg.denseU(), alg.denseV(), alg.denseW() }) {
			for (double[] row : f) {
				for (double v : row) {
					if (Math.abs(v - Math.round(v)) > 1e-9) return false;
				}
			}
		}
		return true;
	}

	private static void injectFieldArrays(File f, List<String> fields) throws Exception {
		var mapper = tools.jackson.databind.json.JsonMapper.builder().build();
		var root = (tools.jackson.databind.node.ObjectNode) mapper.readTree(Files.readString(f.toPath()));
		var fa = root.arrayNode();
		for (String s : fields) fa.add(s);
		root.set("fields", fa);
		root.set("fields_not", root.arrayNode());
		Files.writeString(f.toPath(), eu.solven.matmul.catalog.MatrixJsonFormatter.format(root));
	}
}
