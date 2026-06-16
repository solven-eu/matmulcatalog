package eu.solven.matmul.docs.migrate;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.Compose;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.MatrixJsonFormatter;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.search.LineageReplayer;
import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;
import tools.jackson.databind.node.ArrayNode;
import tools.jackson.databind.node.ObjectNode;

/**
 * Audit + fix every {@code fmm_lille} file. Most are <em>plain derivations</em>
 * (Kronecker / concat of catalog atoms); a few are genuine atoms imported under
 * the wrong attribution. For each file:
 *
 * <ol>
 *   <li><b>Derivation</b> — search for a Kron or 2-part concat decomposition into
 *       catalog schemes whose lineage <em>replays to a Verifier-passing scheme of
 *       the exact ⟨n,m,p⟩ and rank</em>. If found, rewrite as a lineage-only stub
 *       (source {@code Derived_FMM}) and relocate to {@code derived/}. Rank is
 *       preserved; only the redundant FMM coefficients are dropped.</li>
 *   <li><b>Genuine atom</b> — fix attribution where unambiguous: ⟨3,3,6⟩=40 →
 *       Smirnov 2013; the inner-dimension-2 family (a {@code 2} axis, rank = the
 *       Hopcroft-Kerr closed form) → Hopcroft-Kerr 1971. Otherwise leave the file
 *       and log it for human review.</li>
 * </ol>
 *
 * <p>Conservative by construction: never deletes unreproducible data, never lowers
 * a rank, every derivation gated on replay+verify. Dry run by default; {@code --apply}
 * writes.</p>
 */
public final class FixFmmLille {
	private FixFmmLille() {}

	private static final JsonMapper MAPPER = JsonMapper.builder().build();
	private static final String ROOT = "src/main/resources/schemes";

	public static void main(String[] args) throws Exception {
		boolean apply = List.of(args).contains("--apply");
		FieldAwareLookup lookup = new FieldAwareLookup(Field.Q);
		LineageReplayer replayer = LineageReplayer.withDefaultPool(lookup);

		List<Path> files;
		try (var s = Files.walk(Path.of(ROOT))) {
			files = s.filter(p -> p.getFileName().toString().matches(".*fmm_lille.*\\.json"))
					.sorted().toList();
		}
		int derived = 0, attributed = 0, genuine = 0, failed = 0;
		List<String> reviewList = new ArrayList<>();
		for (Path p : files) {
			File f = p.toFile();
			JsonNode root;
			try { root = SchemeIO.parseJson(f); } catch (Exception e) { failed++; continue; }
			JsonNode nArr = root.get("n");
			if (nArr == null || !nArr.isArray() || nArr.size() != 3) { failed++; continue; }
			int n = nArr.get(0).asInt(), m = nArr.get(1).asInt(), pp = nArr.get(2).asInt();
			int r = root.get("m").asInt();

			Lineage.Node lineage = findDerivation(replayer, n, m, pp, r);
			if (lineage != null) {
				derived++;
				NonCubicBilinearAlgorithm alg = replayer.replay(lineage);  // already verified inside
				List<String> fields = fieldsOf(root);
				int maxDim = Math.max(n, Math.max(m, pp));
				String fname = String.format("%dx%dx%d-r%d-derived_fmm-%s.json",
						n, m, pp, r, SchemeIO.shortHash(alg));
				Path dst = Path.of(ROOT, "derived", "section" + maxDim, fname);
				System.out.printf("DERIVED  ⟨%d,%d,%d⟩=%d  %s  → %s%n", n, m, pp, r,
						Lineage.prettyCompact(lineage), dst.getFileName());
				if (apply) {
					Files.createDirectories(dst.getParent());
					SchemeIO.writeStub(alg, dst.toFile(), lineage);
					stamp(dst, fields, "Derived_FMM", null);
					Files.deleteIfExists(p);  // remove the redundant materialised FMM import
				}
				continue;
			}
			// genuine atom — attribution
			String attribSource = attribution(n, m, pp, r);
			if (attribSource != null) {
				attributed++;
				System.out.printf("ATTRIB   ⟨%d,%d,%d⟩=%d  fmm-lille → %s%n", n, m, pp, r, attribSource);
				if (apply) reattribute(f, attribSource);
			} else {
				genuine++;
				reviewList.add(String.format("⟨%d,%d,%d⟩=%d  %s", n, m, pp, r, f.getName()));
			}
		}
		System.out.printf("%n%s: %d derived→derived/, %d re-attributed, %d genuine-kept, %d skipped%n",
				apply ? "APPLIED" : "DRY RUN", derived, attributed, genuine, failed);
		System.out.println("\n=== genuine atoms kept as fmm-lille (need human attribution) ===");
		reviewList.forEach(s -> System.out.println("  " + s));
		if (!apply) System.out.println("\nRe-run with --apply to write.");
	}

	/**
	 * Find a Kron or concat decomposition whose lineage replays to a verifying
	 * ⟨n,m,p⟩=r scheme. Returns the lineage, or null if none validates.
	 */
	private static Lineage.Node findDerivation(LineageReplayer rep, int n, int m, int p, int r) {
		// Kron: every per-axis factor split × the catalog rank of each factor.
		for (int[] nf : factorPairs(n)) for (int[] mf : factorPairs(m)) for (int[] pf : factorPairs(p)) {
			int n1 = nf[0], n2 = nf[1], m1 = mf[0], m2 = mf[1], p1 = pf[0], p2 = pf[1];
			if (n1 * m1 * p1 == 1 || n2 * m2 * p2 == 1) continue;  // trivial ⟨1,1,1⟩ factor
			Lineage.Node lin = new Lineage.KronProduct(
					new Lineage.Atom(ref(n1, m1, p1)), new Lineage.Atom(ref(n2, m2, p2)));
			if (replaysTo(rep, lin, n, m, p, r)) return lin;
		}
		// Concat: split ONE axis into two parts, both catalog schemes.
		for (int a = 1; a < n; a++) { Lineage.Node l = concat(rep, ref(a, m, p), ref(n - a, m, p),
				Axis.N, n, m, p, r); if (l != null) return l; }
		for (int b = 1; b < m; b++) { Lineage.Node l = concat(rep, ref(n, b, p), ref(n, m - b, p),
				Axis.M, n, m, p, r); if (l != null) return l; }
		for (int c = 1; c < p; c++) { Lineage.Node l = concat(rep, ref(n, m, c), ref(n, m, p - c),
				Axis.P, n, m, p, r); if (l != null) return l; }
		return null;
	}

	private enum Axis { N, M, P }

	private static Lineage.Node concat(LineageReplayer rep, String left, String right, Axis ax,
			int n, int m, int p, int r) {
		Lineage.Node lin = switch (ax) {
			case N -> new Lineage.ConcatRows(new Lineage.Atom(left), new Lineage.Atom(right));
			case M -> new Lineage.SumInner(new Lineage.Atom(left), new Lineage.Atom(right));
			case P -> new Lineage.ConcatCols(new Lineage.Atom(left), new Lineage.Atom(right));
		};
		return replaysTo(rep, lin, n, m, p, r) ? lin : null;
	}

	private static boolean replaysTo(LineageReplayer rep, Lineage.Node lin, int n, int m, int p, int r) {
		try {
			NonCubicBilinearAlgorithm a = rep.replay(lin);
			return a.n == n && a.m == m && a.p == p && a.r == r
					&& Verifier.passesRandomMatmulSpotCheck(a);
		} catch (RuntimeException e) {
			return false;
		}
	}

	/** Unambiguous attribution for a genuine atom, or null to leave as fmm-lille. */
	private static String attribution(int n, int m, int p, int r) {
		int[] s = { n, m, p };
		java.util.Arrays.sort(s);
		if (s[0] == 3 && s[1] == 3 && s[2] == 6 && r == 40) return "Smirnov 2013";
		// Hopcroft-Kerr inner-dimension-2 family: a 2 axis, rank = ⌈(3·x·y+max(x,y))/2⌉.
		if (s[0] == 2) {
			int x = s[1], y = s[2];
			int hk = (int) Math.ceil((3.0 * x * y + Math.max(x, y)) / 2.0);
			if (r == hk) return "Hopcroft-Kerr 1971";
		}
		return null;
	}

	private static void reattribute(File f, String source) throws java.io.IOException {
		String existing = Files.readString(f.toPath());
		ObjectNode root = (ObjectNode) MAPPER.readTree(existing);
		root.put("source", source);
		root.put("attribution", source);
		root.put("prior_art", "imported via FMM-Lille; discovery attributed to " + source);
		Files.writeString(f.toPath(), MatrixJsonFormatter.format(root));
		FieldAwareLookup.onSchemeWritten(f);
	}

	private static void stamp(Path file, List<String> fields, String source, String attribution)
			throws java.io.IOException {
		ObjectNode root = (ObjectNode) MAPPER.readTree(Files.readString(file));
		ArrayNode arr = root.arrayNode();
		fields.forEach(arr::add);
		root.set("fields", arr);
		root.put("source", source);
		root.put("commutative", false);
		root.put("verified", true);
		if (attribution != null) root.put("attribution", attribution);
		Files.writeString(file, MatrixJsonFormatter.format(root));
		FieldAwareLookup.onSchemeWritten(file.toFile());
	}

	private static List<String> fieldsOf(JsonNode root) {
		Set<String> out = new LinkedHashSet<>(SchemeIO.fieldTags(root));
		if (out.isEmpty()) out.addAll(List.of("Q", "R", "C"));
		return new ArrayList<>(out);
	}

	private static String ref(int a, int b, int c) { return a + "x" + b + "x" + c; }

	private static List<int[]> factorPairs(int n) {
		List<int[]> out = new ArrayList<>();
		out.add(new int[] { 1, n });
		out.add(new int[] { n, 1 });
		for (int a = 2; a * a <= n; a++) {
			if (n % a == 0) {
				out.add(new int[] { a, n / a });
				if (a != n / a) out.add(new int[] { n / a, a });
			}
		}
		return out;
	}
}
