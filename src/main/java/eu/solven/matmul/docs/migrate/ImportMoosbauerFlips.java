package eu.solven.matmul.docs.migrate;

import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.net.http.HttpClient;
import java.net.http.HttpRequest;
import java.net.http.HttpResponse;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.Duration;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.SchemeIO;

import lombok.extern.slf4j.Slf4j;
import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;

/**
 * Import the <b>original flip-graph</b> matrix-multiplication schemes from
 * <a href="https://github.com/jakobmoosbauer/flips">jakobmoosbauer/flips</a>
 * (Kauers–Moosbauer, "Flip Graphs for Matrix Multiplication", arXiv:2212.01175).
 *
 * <p>Why a <em>direct</em> importer rather than relying on Perminov's re-encoding:
 * these schemes currently reach the catalog laundered through Perminov's
 * {@code known/} subtree and get their attribution reconstructed from folder
 * names. Importing from the origin gives true provenance ({@code source_scheme_url}
 * points at the author's own file), the authors' own integer coefficients (not a
 * ternary-ised re-encoding), and <b>both</b> field variants per shape. We keep
 * Perminov's add-reduced versions alongside — they may still win on the additions
 * axis at equal rank.</p>
 *
 * <p>Each {@code solutions/*.exp} file is one rank-one term per line in the
 * trilinear convention {@code ⟨n,m,p⟩ = Σ a_{ij} ⊗ b_{jk} ⊗ c_{ki}}:
 * {@code (A-form)*(B-form)*(C-form)}. Unlike the {@code ±1}-only Kauers-2026
 * "structured" files, the flip-graph {@code mod0} files carry <b>integer
 * coefficients</b> (e.g. {@code -2*a11+2*a12+a21-a22}, {@code c21+2*c22}); the
 * term regex captures the magnitude, not just the sign. The {@code .m} (Maple)
 * siblings are duplicates of the {@code .exp} files and are ignored.</p>
 *
 * <p>Filename {@code {nmp}-{rank}-mod{0|0a|2}.exp}: {@code nmp} is the
 * <em>sorted</em> shape, so the true orientation is read from the actual index
 * ranges (as in {@code ImportKauers2026Structured}). {@code mod0}/{@code mod0a}
 * are characteristic-0 (integer ⇒ exact over ℚ, reduces to every field);
 * {@code mod2} is an F₂-only scheme that need not lift to characteristic 0
 * (e.g. {@code 555-95-mod2} = F₂⟨5,5,5⟩=95 vs {@code 555-97-mod0} = ℤ⟨5,5,5⟩=97).</p>
 *
 * <p>Run (defaults: maxDim 32, no limit; fetches over HTTP from GitHub):</p>
 * <pre>MAVEN_OPTS="-Xmx2g" mvn -q -ntp exec:java \
 *   -Dexec.mainClass=eu.solven.matmul.docs.migrate.ImportMoosbauerFlips \
 *   -Dexec.args="--max-dim=16"</pre>
 *
 * <p>After importing, stamp + regenerate: {@code StampFields}, {@code StampAdditions},
 * {@code GenerateCatalogManifest} (the latter computes {@code zt} from coefficients).
 * Discovery-vs-rediscovery tagging (which of these ranks the flip graph
 * <em>first</em> established vs merely re-found) is a follow-up audit per the
 * {@code research/DISCOVERIES_PENDING_ANALYSIS.md} workflow — {@code source} here
 * records where we obtained the scheme, matching how the Perminov-routed copies
 * are already labelled.</p>
 */
@Slf4j
public final class ImportMoosbauerFlips {

	private ImportMoosbauerFlips() {}

	private static final String CONTENTS_URL =
			"https://api.github.com/repos/jakobmoosbauer/flips/contents/solutions?ref=main";
	private static final String RAW_BASE =
			"https://raw.githubusercontent.com/jakobmoosbauer/flips/main/solutions/";
	private static final String BLOB_BASE =
			"https://github.com/jakobmoosbauer/flips/blob/main/solutions/";
	private static final String USER_AGENT = "solven-matmul-catalog/flips-import";
	private static final Path KNOWN = Path.of("src/main/resources/schemes/known");

	/** Origin attribution — the original flip-graph paper (arXiv:2212.01175);
	 *  matches the {@code jakobmoosbauer_flips} case in {@code PerminovKnownAttribution}. */
	private static final String SOURCE = "Kauers-Moosbauer 2023";
	private static final String PAPER_URL = "https://arxiv.org/abs/2212.01175";

	/** {@code {nmp}-{rank}-mod{0|0a|2}.exp} — nmp is the SORTED shape (label only). */
	private static final Pattern NAME =
			Pattern.compile("^(\\d+)-(\\d+)-(mod0a|mod0|mod2)\\.exp$");
	/** One signed, integer-coefficient term: optional sign, optional {@code k*},
	 *  then {@code [abc]} + two single-digit indices (max dim is 6). */
	private static final Pattern TERM =
			Pattern.compile("([+-]?)\\s*(\\d*)\\s*\\*?\\s*([abc])(\\d)(\\d)");

	private static final JsonMapper MAPPER = JsonMapper.builder().build();

	public static void main(String[] args) throws Exception {
		int maxDim = intArg(args, "--max-dim", 32);
		int limit = intArg(args, "--limit", 0);

		log.info("listing flips solutions/ …");
		JsonNode contents = MAPPER.readTree(fetch(CONTENTS_URL));
		List<String> names = new ArrayList<>();
		for (JsonNode e : contents) {
			String name = e.path("name").asString();
			if (name.endsWith(".exp")) {
				names.add(name);
			}
		}
		names.sort(null);
		log.info("found {} .exp solution files", names.size());

		int wrote = 0, skipExisting = 0, skipRange = 0, fail = 0;
		long t0 = System.nanoTime();
		for (String name : names) {
			Matcher fm = NAME.matcher(name);
			if (!fm.matches()) {
				log.warn("[skip] unexpected name {}", name);
				continue;
			}
			String mod = fm.group(3);
			boolean f2 = mod.equals("mod2");

			try {
				String body = fetch(RAW_BASE + name);
				NonCubicBilinearAlgorithm alg = parse(body);
				if (alg == null) {
					log.warn("[skip] {} — no parseable terms", name);
					fail++;
					continue;
				}
				int maxd = Math.max(alg.n, Math.max(alg.m, alg.p));
				if (maxd > maxDim) {
					skipRange++;
					continue;
				}
				boolean ok = f2 ? Verifier.isExactNonCubicF2(alg) : Verifier.isExactNonCubic(alg);
				if (!ok) {
					// A mod0 file that fails ℚ-exactness, or a mod2 that fails F₂, means
					// the parse or the source is wrong — fail loud, never write it.
					log.warn("[FAIL] {} (⟨{},{},{}⟩ r={}) did NOT verify over {}", name,
							alg.n, alg.m, alg.p, alg.r, f2 ? "F2" : "Q");
					fail++;
					continue;
				}

				String hash7 = SchemeIO.contentHash(alg).substring(0, 7);
				String note = "flips_" + mod;
				Path dir = KNOWN.resolve("section" + maxd);
				Files.createDirectories(dir);
				File out = dir.resolve(alg.n + "x" + alg.m + "x" + alg.p
						+ "-r" + alg.r + "-" + note + "-" + hash7 + ".json").toFile();
				if (out.exists()) {
					skipExisting++;
					continue;
				}
				SchemeIO.write(alg, out);

				Map<String, Object> meta = new LinkedHashMap<>();
				meta.put("source", SOURCE);
				meta.put("imported_via", "jakobmoosbauer/flips");
				meta.put("original_source_path", "solutions/" + name);
				meta.put("source_scheme_url", BLOB_BASE + name);
				meta.put("source_paper_url", PAPER_URL);
				meta.put("commutative", false);
				meta.put("verified", true);
				// mod0/mod0a: integer coefficients ⇒ exact over ℤ, reduces to every
				// field (theorem). mod2: F₂-only scheme that need not lift.
				meta.put("fields", f2
						? List.of("F2")
						: List.of("F2", "F3", "Z", "Q", "R", "C"));
				SchemeIO.addFields(out, meta, /* apply */ true);

				int adds = Verifier.additionCount(alg);
				log.info("⟨{},{},{}⟩ r={} a={} {} → {}", alg.n, alg.m, alg.p, alg.r, adds,
						f2 ? "F2" : "Z", out.getName());
				wrote++;
				if (limit > 0 && wrote >= limit) {
					break;
				}
			} catch (RuntimeException | IOException | InterruptedException e) {
				log.warn("[ERR] {}: {}", name, e.toString());
				fail++;
			}
		}
		long ms = (System.nanoTime() - t0) / 1_000_000L;
		log.info("Done: {} imported, {} skipped-existing, {} out-of-range, {} failed ({}ms). "
				+ "Next: StampFields / StampAdditions / GenerateCatalogManifest.",
				wrote, skipExisting, skipRange, fail, ms);
	}

	/**
	 * Parse a flip-graph {@code .exp} body (one rank-one term per line) into a
	 * {@link NonCubicBilinearAlgorithm}. Returns {@code null} if no term parses.
	 * The true orientation ⟨n,m,p⟩ is recovered from the index ranges, not the
	 * (sorted) filename.
	 */
	static NonCubicBilinearAlgorithm parse(String body) {
		List<String[]> terms = new ArrayList<>();  // [Aform, Bform, Cform]
		List<Integer> signs = new ArrayList<>();   // overall ±1 per term
		for (String raw : body.split("\\R")) {
			String line = raw.trim();
			if (line.isEmpty() || !line.contains("*")) {
				continue;
			}
			// Peel an overall term sign and a paren pair wrapping the WHOLE product
			// (e.g. `-(a25*(…)*c12)` negates the whole rank-one term). Loops because a
			// term can carry both (`-( … )`). The sign is folded into one factor below;
			// negating the product ≡ negating any single factor.
			int sign = 1;
			boolean changed = true;
			while (changed) {
				changed = false;
				if (line.startsWith("-")) { sign = -sign; line = line.substring(1).trim(); changed = true; }
				else if (line.startsWith("+")) { line = line.substring(1).trim(); changed = true; }
				if (line.startsWith("(") && wrapsWhole(line)) {
					line = line.substring(1, line.length() - 1).trim();
					changed = true;
				}
			}
			// Factors are separated by top-level '*' (depth 0). A '*' inside a
			// coefficient (e.g. 2*a11) always sits inside parentheses, so depth-aware
			// splitting keeps coefficients intact AND keeps bare single-variable
			// factors (e.g. the un-parenthesised `c12` in `(-a21)*(-b11+b13)*c12`),
			// which a parens-only extractor would silently drop.
			String[] g = new String[3];
			for (String factor : splitTopLevel(line)) {
				char which = firstLetter(factor);
				if (which == '?') {
					continue;
				}
				g[which == 'a' ? 0 : which == 'b' ? 1 : 2] = factor;
			}
			if (g[0] == null || g[1] == null || g[2] == null) {
				continue;
			}
			terms.add(g);
			signs.add(sign);
		}
		int r = terms.size();
		if (r == 0) {
			return null;
		}
		// pass 1: recover the true (unsorted) shape from index ranges.
		int n = 0, m = 0, p = 0;
		for (String[] t : terms) {
			for (int[] e : parseForm(t[0])) { n = Math.max(n, e[0]); m = Math.max(m, e[1]); }  // a_{ij}
			for (int[] e : parseForm(t[1])) { m = Math.max(m, e[0]); p = Math.max(p, e[1]); }  // b_{jk}
			for (int[] e : parseForm(t[2])) { p = Math.max(p, e[0]); n = Math.max(n, e[1]); }  // c_{ki}
		}
		// pass 2: fill U/V/W (dim×rank).
		double[][] U = new double[n * m][r], V = new double[m * p][r], W = new double[n * p][r];
		for (int k = 0; k < r; k++) {
			String[] t = terms.get(k);
			int sign = signs.get(k);  // fold the overall term sign into the a-form
			for (int[] e : parseForm(t[0])) U[(e[0] - 1) * m + (e[1] - 1)][k] += sign * e[2];  // a_{ij}
			for (int[] e : parseForm(t[1])) V[(e[0] - 1) * p + (e[1] - 1)][k] += e[2];  // b_{jk}
			for (int[] e : parseForm(t[2])) W[(e[1] - 1) * p + (e[0] - 1)][k] += e[2];  // c_{ki}→C_{ik}
		}
		return new NonCubicBilinearAlgorithm(n, m, p, U, V, W);
	}

	/** Does the leading {@code (} of {@code s} (which must start with {@code (}) close
	 *  only at the final character — i.e. one paren pair wraps the whole expression? */
	private static boolean wrapsWhole(String s) {
		int depth = 0;
		for (int i = 0; i < s.length(); i++) {
			char c = s.charAt(i);
			if (c == '(') {
				depth++;
			} else if (c == ')') {
				depth--;
				if (depth == 0) {
					return i == s.length() - 1;
				}
			}
		}
		return false;
	}

	/** Split a {@code FACTOR*FACTOR*FACTOR} line on {@code *} at parenthesis depth 0,
	 *  so coefficient stars ({@code 2*a11}, always inside parens) stay within a factor. */
	private static List<String> splitTopLevel(String line) {
		List<String> parts = new ArrayList<>();
		int depth = 0, start = 0;
		for (int i = 0; i < line.length(); i++) {
			char c = line.charAt(i);
			if (c == '(') {
				depth++;
			} else if (c == ')') {
				depth--;
			} else if (c == '*' && depth == 0) {
				parts.add(line.substring(start, i));
				start = i + 1;
			}
		}
		parts.add(line.substring(start));
		return parts;
	}

	private static char firstLetter(String s) {
		for (int i = 0; i < s.length(); i++) {
			char c = s.charAt(i);
			if (c == 'a' || c == 'b' || c == 'c') {
				return c;
			}
		}
		return '?';
	}

	/** Parse a linear form into {@code [idx1, idx2, coeff]} triples. A bare term
	 *  ({@code a11}) has coefficient ±1; {@code k*a11} has coefficient ±k. */
	private static List<int[]> parseForm(String form) {
		List<int[]> out = new ArrayList<>();
		Matcher m = TERM.matcher(form);
		while (m.find()) {
			int sign = "-".equals(m.group(1)) ? -1 : 1;
			String mag = m.group(2);
			int coeff = sign * (mag.isEmpty() ? 1 : Integer.parseInt(mag));
			out.add(new int[] { Integer.parseInt(m.group(4)), Integer.parseInt(m.group(5)), coeff });
		}
		return out;
	}

	private static String fetch(String url) throws IOException, InterruptedException {
		HttpClient client = HttpClient.newBuilder()
				.followRedirects(HttpClient.Redirect.NORMAL)
				.connectTimeout(Duration.ofSeconds(30))
				.build();
		HttpRequest req = HttpRequest.newBuilder(URI.create(url))
				.header("User-Agent", USER_AGENT)
				.timeout(Duration.ofSeconds(60))
				.GET()
				.build();
		IOException last = null;
		for (int attempt = 0; attempt <= 2; attempt++) {
			try {
				HttpResponse<String> resp = client.send(req, HttpResponse.BodyHandlers.ofString());
				if (resp.statusCode() == 200) {
					return resp.body();
				}
				last = new IOException("HTTP " + resp.statusCode() + " for " + url);
			} catch (IOException e) {
				last = e;
			}
			Thread.sleep(500L * (attempt + 1));
		}
		throw last;
	}

	private static int intArg(String[] args, String key, int dflt) {
		for (String a : args) {
			if (a.startsWith(key + "=")) {
				return Integer.parseInt(a.substring(key.length() + 1));
			}
		}
		return dflt;
	}
}
