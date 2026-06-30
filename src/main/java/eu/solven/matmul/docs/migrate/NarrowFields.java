package eu.solven.matmul.docs.migrate;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.SchemeIO;
import lombok.extern.slf4j.Slf4j;
import tools.jackson.databind.JsonNode;

/**
 * Expands each bilinear scheme's {@code fields[]} from the single-element tag the
 * structural #174 migration left, to the FULL set of fields the scheme is valid
 * over — so the Z ⊂ Q ⊂ R ⊂ C fallback is explicit rather than inferred (#175).
 *
 * <p>Per the user's 2026-06-03 insight, R/C are free given Z or Q (containment):
 * we never run 6 separate verifications. Instead:</p>
 * <ul>
 *   <li><b>integer</b> coefficients → exact-verify over Z once. Z-exact ⇒ valid
 *       in EVERY field (the integer identity reduces mod any prime), so
 *       {@code [F2,F3,Z,Q,R,C]}. If NOT Z-exact (AlphaTensor-F2 style: {0,1}
 *       coefficients whose identity only holds mod 2), test F2/F3 and tag those.</li>
 *   <li><b>rational</b> coefficients → verify over Q once; Q-valid ⇒
 *       {@code [Q,R,C]} (½ excludes Z/F2/F3).</li>
 * </ul>
 *
 * <p>Exact verification is O(n⁶·r); to stay tractable we only verify schemes with
 * {@code max(n,m,p) ≤ MAX_VERIFY_DIM}. Larger schemes trust their already-corrected
 * base tag and expand mechanically along the containment chain.</p>
 *
 * <p>Rewrites a file only when {@code fields[]} actually changes (minimal diff).</p>
 *
 * <pre>mvn -q exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.NarrowFields -Dexec.args="--apply"</pre>
 */
@Slf4j
public final class NarrowFields {

	private static final Path SCHEMES_ROOT = Path.of("src/main/resources/schemes");
	/** Exact verification cap — above this, expand mechanically from the base tag. */
	private static final int MAX_VERIFY_DIM = 10;

	private NarrowFields() {}

	public static void main(String[] args) throws Exception {
		boolean apply = List.of(args).contains("--apply");
		List<Path> files;
		try (Stream<Path> walk = Files.walk(SCHEMES_ROOT)) {
			files = walk.filter(p -> p.toString().endsWith(".json")).toList();
		}
		AtomicInteger changed = new AtomicInteger();
		AtomicInteger verified = new AtomicInteger();
		AtomicInteger mechanical = new AtomicInteger();
		AtomicInteger errors = new AtomicInteger();

		files.parallelStream().forEach(file -> {
			try {
				JsonNode root = SchemeIO.parseJson(file.toFile());
				if (SchemeIO.isNonBilinear(root) || SchemeIO.isComplex(root)) return;
				List<String> cur = SchemeIO.fieldTags(root);
				if (cur.isEmpty()) return;
				// Stubs have no matrices to verify, but the mechanical containment
				// chain still applies to their existing base tag (Z ⇒ F2/F3/Q/R/C is
				// a theorem given the Z claim) — without this, the ~7.8k derived
				// stub under-claimers stay invisible to the SPA's F2/F3 selectors.
				if (SchemeIO.isStub(root) || (!root.has("u_sparse") && !root.has("u"))) {
					List<String> expanded = expandFromBase(cur);
					if (expanded.isEmpty()) { errors.incrementAndGet(); return; }
					mechanical.incrementAndGet();
					if (!sameSet(cur, expanded)) {
						if (apply) rewriteFields(file, expanded);
						changed.incrementAndGet();
					}
					return;
				}

				NonCubicBilinearAlgorithm alg = SchemeIO.read(root);
				int maxDim = Math.max(alg.n, Math.max(alg.m, alg.p));
				List<String> fields;
				if (maxDim <= MAX_VERIFY_DIM) {
					fields = verifyFields(alg);
					verified.incrementAndGet();
				} else {
					fields = expandFromBase(cur);
					mechanical.incrementAndGet();
				}
				if (fields.isEmpty()) { errors.incrementAndGet(); return; }

				if (!sameSet(cur, fields) && apply) {
					rewriteFields(file, fields);
					changed.incrementAndGet();
				} else if (!sameSet(cur, fields)) {
					changed.incrementAndGet();
				}
			} catch (Exception e) {
				errors.incrementAndGet();
				log.warn("failed on {}: {}", file.getFileName(), e.getMessage());
			}
		});
		log.info("NarrowFields: changed={} (verified={}, mechanical={}) errors={} apply={}",
				changed.get(), verified.get(), mechanical.get(), errors.get(), apply);
	}

	/** Exact field set via the cheap Z/Q-base + mechanical-containment rule. */
	static List<String> verifyFields(NonCubicBilinearAlgorithm alg) {
		if (allIntegers(alg)) {
			if (Verifier.isExactNonCubic(alg)) {
				// Integer-exact ⇒ valid over every field (reduces mod any prime).
				return List.of("F2", "F3", "Z", "Q", "R", "C");
			}
			// Not Z-exact but integer: prime-field-only (e.g. AlphaTensor F₂).
			java.util.List<String> out = new java.util.ArrayList<>();
			if (Verifier.isExactNonCubicF2(alg)) out.add("F2");
			if (Verifier.isExactNonCubicF3(alg)) out.add("F3");
			return out;
		}
		// Rational coefficients: verify once, then R/C are free. A Q-exact scheme
		// ALSO reduces mod a prime p when every denominator is coprime to p
		// (num·den⁻¹ mod p): 1/2 is representable in F₃ (2⁻¹ ≡ 2) but not F₂;
		// 1/3 in F₂ (3 ≡ 1) but not F₃. Verifier.residualNonCubicFp checks
		// representability + exactness per prime.
		if (!Verifier.isExactNonCubic(alg)) return List.of();
		java.util.List<String> out = new java.util.ArrayList<>();
		if (Verifier.isExactNonCubicFp(alg, 2)) out.add("F2");
		if (Verifier.isExactNonCubicFp(alg, 3)) out.add("F3");
		out.addAll(List.of("Q", "R", "C"));
		return out;
	}

	/** Mechanical expansion from the structural base tag (no verification). */
	static List<String> expandFromBase(List<String> cur) {
		if (cur.contains("Z")) return List.of("F2", "F3", "Z", "Q", "R", "C");
		if (cur.contains("Q")) return List.of("Q", "R", "C");
		if (cur.contains("R")) return List.of("R", "C");
		if (cur.contains("C")) return List.of("C");
		if (cur.contains("F2")) return List.of("F2");
		if (cur.contains("F3")) return List.of("F3");
		return cur;
	}

	private static boolean allIntegers(NonCubicBilinearAlgorithm alg) {
		for (double[][] mat : new double[][][] { alg.denseU(), alg.denseV(), alg.denseW() }) {
			for (double[] row : mat) {
				for (double v : row) {
					if (v != Math.rint(v)) return false;
				}
			}
		}
		return true;
	}

	private static boolean sameSet(List<String> a, List<String> b) {
		return new java.util.HashSet<>(a).equals(new java.util.HashSet<>(b));
	}

	private static final List<String> ALL = List.of("F2", "F3", "Z", "Q", "R", "C");
	private static final Pattern FIELDS = Pattern.compile("\"fields\"\\s*:\\s*\\[[^\\]]*\\]");
	private static final Pattern FIELDS_NOT = Pattern.compile("\"fields_not\"\\s*:\\s*\\[[^\\]]*\\]");

	/** Textual patch of fields[] / fields_not[] — minimal diff, preserves formatting. */
	private static void rewriteFields(Path file, List<String> fields) throws java.io.IOException {
		String s = Files.readString(file);
		String fieldsArr = "\"fields\": [" + quoteJoin(fields) + "]";
		java.util.List<String> not = new java.util.ArrayList<>();
		for (String f : ALL) if (!fields.contains(f)) not.add(f);
		String notArr = "\"fields_not\": [" + quoteJoin(not) + "]";

		Matcher mf = FIELDS.matcher(s);
		if (mf.find()) s = s.substring(0, mf.start()) + fieldsArr + s.substring(mf.end());
		Matcher mn = FIELDS_NOT.matcher(s);
		if (mn.find()) s = s.substring(0, mn.start()) + notArr + s.substring(mn.end());
		// Route through the single shared formatter so the on-disk file stays
		// canonical (and is validated as JSON before being written).
		Files.writeString(file, eu.solven.matmul.catalog.MatrixJsonFormatter.format(s));
	}

	private static String quoteJoin(List<String> xs) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < xs.size(); i++) {
			if (i > 0) sb.append(", ");
			sb.append('"').append(xs.get(i)).append('"');
		}
		return sb.toString();
	}
}
