package eu.solven.matmul.papers.kauers2026;

import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.SchemeIO;
import lombok.extern.slf4j.Slf4j;

/**
 * Import the Kauers–Moosbauer–Wood 2026 "structured" decompositions
 * (arXiv:2602.11041, REFERENCES.md [83]) from their public repo
 * {@code github.com/mkauers/matrix-multiplication/structured}, pre-downloaded to
 * {@code target/kmw-exp/*.exp}.
 *
 * <p>Each {@code .exp} line is one rank-one term in the trilinear convention
 * {@code ⟨n,m,p⟩ = Σ a_{ij} ⊗ b_{jk} ⊗ c_{ki}}:
 * {@code (A-form)*(B-form)*(C-form)} with ±1 coefficients. The directory name is
 * the <em>sorted</em> shape; the true orientation is read from the index ranges
 * ({@code n}=max i in a / max i in c, {@code m}=max j in a / max j in b,
 * {@code p}=max k in b / max k in c). The C-form's {@code c_{ki}} maps to output
 * {@code C_{ik}} → {@code W} row {@code (i-1)·p+(k-1)}.</p>
 */
@Slf4j
public final class ImportKauers2026Structured {
	private ImportKauers2026Structured() {}

	private static final Pattern TERM = Pattern.compile("([+-]?)\\s*([abc])([1-9])([1-9])");
	private static final Pattern GROUP = Pattern.compile("\\(([^)]*)\\)");

	public static void main(String[] args) throws Exception {
		File dir = new File("target/kmw-exp");
		File[] exps = dir.listFiles((d, name) -> name.endsWith(".exp"));
		if (exps == null) {
			throw new IllegalStateException("no .exp files in " + dir);
		}
		java.util.Arrays.sort(exps);
		int ok = 0, f2 = 0, bad = 0;
		for (File exp : exps) {
			List<String> lines = Files.readAllLines(exp.toPath());
			List<String[]> terms = new ArrayList<>();  // [Aform, Bform, Cform]
			for (String line : lines) {
				if (line.isBlank() || !line.contains("*")) {
					continue;
				}
				String[] g = new String[3];
				Matcher gm = GROUP.matcher(line);
				while (gm.find()) {
					String body = gm.group(1);
					char which = firstLetter(body);
					int idx = which == 'a' ? 0 : which == 'b' ? 1 : 2;
					g[idx] = body;
				}
				if (g[0] == null || g[1] == null || g[2] == null) {
					continue;
				}
				terms.add(g);
			}
			int r = terms.size();
			if (r == 0) {
				continue;
			}
			// pass 1: dims
			int n = 0, m = 0, p = 0;
			for (String[] t : terms) {
				for (int[] e : parse(t[0])) { n = Math.max(n, e[0]); m = Math.max(m, e[1]); }
				for (int[] e : parse(t[1])) { m = Math.max(m, e[0]); p = Math.max(p, e[1]); }
				for (int[] e : parse(t[2])) { p = Math.max(p, e[0]); n = Math.max(n, e[1]); }
			}
			// pass 2: fill U/V/W
			double[][] U = new double[n * m][r], V = new double[m * p][r], W = new double[n * p][r];
			for (int k = 0; k < r; k++) {
				String[] t = terms.get(k);
				for (int[] e : parse(t[0])) U[(e[0] - 1) * m + (e[1] - 1)][k] += e[2];  // a_{ij}
				for (int[] e : parse(t[1])) V[(e[0] - 1) * p + (e[1] - 1)][k] += e[2];  // b_{jk}
				for (int[] e : parse(t[2])) W[(e[1] - 1) * p + (e[0] - 1)][k] += e[2];  // c_{ki}→C_{ik}
			}
			NonCubicBilinearAlgorithm alg = new NonCubicBilinearAlgorithm(n, m, p, U, V, W);
			boolean exactQ = Verifier.isExactNonCubic(alg);
			boolean exactF2 = !exactQ && Verifier.isExactNonCubicF2(alg);
			if (!exactQ && !exactF2) {
				log.warn("⟨{},{},{}⟩ r={} from {} — does NOT verify (Q or F2); skip", n, m, p, r, exp.getName());
				bad++;
				continue;
			}
			int adds = Verifier.additionCount(alg);
			File outDir = new File("src/main/resources/schemes/known/section" + Math.max(n, Math.max(m, p)));
			outDir.mkdirs();
			// Field (Q vs F2) lives in fields[] content, not the filename; the content
			// hash distinguishes the two variants (different U/V/W → different name).
			File out = new File(outDir, SchemeIO.canonicalName(alg, "kauers_2026"));
			SchemeIO.write(alg, out);
			if (exactQ) ok++; else f2++;
			log.info("⟨{},{},{}⟩ r={} a={} {} → {}", n, m, p, r, adds, exactQ ? "Q" : "F2", out.getName());
		}
		log.info("=== imported {} (Q) + {} (F2); {} unverifiable ===", ok, f2, bad);
	}

	private static char firstLetter(String s) {
		for (int i = 0; i < s.length(); i++) {
			char c = s.charAt(i);
			if (c == 'a' || c == 'b' || c == 'c') return c;
		}
		return '?';
	}

	/** Parse a form into [idx1, idx2, ±1] triples. */
	private static List<int[]> parse(String form) {
		List<int[]> out = new ArrayList<>();
		Matcher m = TERM.matcher(form);
		while (m.find()) {
			int sign = "-".equals(m.group(1)) ? -1 : 1;
			out.add(new int[] { Integer.parseInt(m.group(3)), Integer.parseInt(m.group(4)), sign });
		}
		return out;
	}
}
