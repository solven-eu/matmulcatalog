package eu.solven.matmul.docs.migrate;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.CatalogPolicy;
import eu.solven.matmul.catalog.Compose;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.ProjectionSearch;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.catalog.VectorProjection;
import eu.solven.matmul.search.LineageReplayer;
import eu.solven.matmul.verifiers.Verifier;

import lombok.extern.slf4j.Slf4j;

/**
 * Coordinate-project Kronecker-product parents built from <em>every representative</em> of
 * their factors — the stored files and their {@link VectorProjection} rebases — and persist a
 * strict improvement as a {@code Project(KronProduct(…))} stub pinned on the exact factor files.
 *
 * <p>The closure's projection scatter ({@code ProjectFmmGaps}) resolves ONE parent per shape (the
 * rank-best file) and prunes it by its stored coordinate margin, so a Kronecker stub whose
 * factor representatives have margin 0 is never projected. Here the representative is the
 * variable: for each factor shape every best-rank file is offered, plus the GL-rebased
 * representatives that turn a common rank-1 direction into a coordinate. A rebased
 * representative that a win stands on is written — only then — as an explicit scheme under
 * {@code schemes/margin-bases/sectionN/} (rank-tied, margin-rich: the class
 * {@code references/PROJECTION_MARGIN_TRADEOFF.md} reserves that folder for) through the
 * contributed-scheme gate, {@code discovery:false}, attribution = the original's. 2026-11-09:
 * ⟨12,23,24⟩ 3790 → 3776 (FMM 3788) via a rebased ⟨3,3,6⟩ ⊗ ⟨4,8,4⟩.</p>
 *
 * <p>For the FMM {@code ProjOfAlongMode} rows themselves (a drop that saves m / n / p products
 * from any representative) see {@link MaterialiseSubstitutionProjections}; this driver covers
 * what coordinate margins can add on top of that. Tier of every rank written: <b>bound</b>.</p>
 *
 * <pre>
 *   mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.MaterialiseVectorProjections \
 *       -Dexec.args="--shapes-file=generated/upstream-gaps.txt [--execute] [--max-factor-dim=8] [--slack=0.05]"
 * </pre>
 * Dry-run by default (nothing written; a rebased factor it would pin is reported as {@code *}).
 */
@Slf4j
public final class MaterialiseVectorProjections {
	private MaterialiseVectorProjections() {}

	private static final Path SCHEMES_ROOT = Path.of("src/main/resources/schemes");
	private static final String NOTE = "rebased_vector_margin";

	/** A factor representative; {@code file} is null for a rebase that no win has pinned yet. */
	static final class Rep {
		File file;
		final String hash;
		final NonCubicBilinearAlgorithm alg;
		final Rep origin;
		final String originSource;

		Rep(File file, String hash, NonCubicBilinearAlgorithm alg, Rep origin, String originSource) {
			this.file = file;
			this.hash = hash;
			this.alg = alg;
			this.origin = origin;
			this.originSource = originSource;
		}

		boolean rebased() {
			return origin != null;
		}
	}

	record Win(int[] shape, int ours, NonCubicBilinearAlgorithm alg, Lineage.Node lineage, String how, List<Rep> reps) {}

	public static void main(String[] args) throws Exception {
		List<int[]> targets = new ArrayList<>();
		boolean execute = false;
		int maxFactorDim = 8;
		double slack = 0.05;
		for (String a : args) {
			if (a.equals("--execute")) {
				execute = true;
			} else if (a.startsWith("--shapes-file=")) {
				targets.addAll(readShapeFile(a.substring("--shapes-file=".length())));
			} else if (a.startsWith("--shapes=")) {
				for (String tok : a.substring("--shapes=".length()).split("[;,]+")) {
					if (!tok.isBlank()) {
						targets.add(parseShape(tok));
					}
				}
			} else if (a.startsWith("--max-factor-dim=")) {
				maxFactorDim = Integer.parseInt(a.substring("--max-factor-dim=".length()));
			} else if (a.startsWith("--slack=")) {
				slack = Double.parseDouble(a.substring("--slack=".length()));
			} else {
				throw new IllegalArgumentException("unknown arg: " + a);
			}
		}
		if (targets.isEmpty()) {
			throw new IllegalArgumentException("usage: --shapes-file=PATH | --shapes=NxMxP,… [--execute]");
		}
		targets.sort(Comparator.comparingLong(t -> (long) t[0] * t[1] * t[2]));
		log.info("MaterialiseVectorProjections: {} target(s), maxFactorDim={}, slack={}, {}", targets.size(),
				maxFactorDim, slack, execute ? "--execute" : "DRY-RUN");

		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		Map<String, List<Rep>> repCache = new LinkedHashMap<>();
		List<Win> wins = new ArrayList<>();
		long t0 = System.currentTimeMillis();
		int done = 0;
		for (int[] t : targets) {
			done++;
			int ours = lookup.findRank(t[0], t[1], t[2]);
			if (!known(ours)) {
				log.info("⟨{},{},{}⟩: no catalog rank — skipped", t[0], t[1], t[2]);
				continue;
			}
			Optional<Win> w = bestFor(t, ours, lookup, repCache, maxFactorDim, slack);
			if (w.isPresent()) {
				wins.add(w.get());
				log.info("[win] ⟨{},{},{}⟩ {} → {}  via {}", t[0], t[1], t[2], ours, w.get().alg().r, w.get().how());
				if (execute) {
					persist(w.get(), lookup);
				}
			} else {
				log.info("⟨{},{},{}⟩ = {}: no Kronecker parent projects below it", t[0], t[1], t[2], ours);
			}
			log.info("[progress] {}/{} targets, {} win(s), {}ms elapsed", done, targets.size(), wins.size(),
					System.currentTimeMillis() - t0);
		}
		log.info("MaterialiseVectorProjections done: {} win(s) over {} target(s) ({})", wins.size(), targets.size(),
				execute ? "written" : "DRY-RUN, nothing written");
		for (Win w : wins) {
			System.out.printf("%dx%dx%d  %d -> %d  %s%n", w.shape()[0], w.shape()[1], w.shape()[2], w.ours(), w.alg().r,
					w.how());
		}
	}

	/** The best projected Kronecker parent for one target, over all parents / splits / representatives. */
	static Optional<Win> bestFor(int[] t, int ours, FieldAwareLookup lookup, Map<String, List<Rep>> repCache,
			int maxFactorDim, double slack) throws IOException {
		long upper = ours;
		Win best = null;
		for (int dn = 0; dn <= 1; dn++) {
			for (int dm = 0; dm <= 1; dm++) {
				for (int dp = 0; dp <= 1; dp++) {
					if (dn + dm + dp == 0) {
						continue;
					}
					int[] P = { t[0] + dn, t[1] + dm, t[2] + dp };
					for (int[][] split : splits(P, maxFactorDim)) {
						int[] f1 = split[0], f2 = split[1];
						int r1 = lookup.findRank(f1[0], f1[1], f1[2]), r2 = lookup.findRank(f2[0], f2[1], f2[2]);
						if (!known(r1) || !known(r2) || (long) r1 * r2 > ours * (1 + slack)) {
							continue;
						}
						for (Oriented a : orientedReps(f1, lookup, repCache)) {
							for (Oriented b : orientedReps(f2, lookup, repCache)) {
								NonCubicBilinearAlgorithm k = Compose.kroneckerGeneral(a.alg(), b.alg());
								Optional<ProjectionSearch.Hit> hit = ProjectionSearch.bestFor(t[0], t[1], t[2], List.of(k),
										upper, 1);
								if (hit.isEmpty()) {
									continue;
								}
								ProjectionSearch.Hit h = hit.get();
								Lineage.Node tree = new Lineage.Project(new Lineage.KronProduct(a.leaf(), b.leaf()), h.keepN(),
										h.keepM(), h.keepP());
								String how = String.format("Project(%dx%dx%d = %dx%dx%d%s ⊗ %dx%dx%d%s)", P[0], P[1], P[2],
										f1[0], f1[1], f1[2], a.rep().rebased() ? "*" : "", f2[0], f2[1], f2[2],
										b.rep().rebased() ? "*" : "");
								best = new Win(t, ours, h.scheme(), tree, how, List.of(a.rep(), b.rep()));
								upper = h.rank();
							}
						}
					}
				}
			}
		}
		return Optional.ofNullable(best);
	}

	private static boolean known(int r) {
		return r > 0 && r < eu.solven.matmul.recombination.Recombination.SotaResolver.UNKNOWN_RANK;
	}

	/** All Kronecker factorisations {@code P = f1 ⊗ f2} with both factors non-trivial and within the dim cap;
	 *  ordered pairs (the product is not symmetric in the operands' roles). */
	static List<int[][]> splits(int[] P, int maxFactorDim) {
		List<int[][]> out = new ArrayList<>();
		for (int a = 1; a <= P[0]; a++) {
			if (P[0] % a != 0) continue;
			for (int b = 1; b <= P[1]; b++) {
				if (P[1] % b != 0) continue;
				for (int c = 1; c <= P[2]; c++) {
					if (P[2] % c != 0) continue;
					int[] f1 = { a, b, c }, f2 = { P[0] / a, P[1] / b, P[2] / c };
					if (a * b * c == 1 || f2[0] * f2[1] * f2[2] == 1) continue;
					if (Math.max(a, Math.max(b, c)) > maxFactorDim
							|| Math.max(f2[0], Math.max(f2[1], f2[2])) > maxFactorDim) continue;
					out.add(new int[][] { f1, f2 });
				}
			}
		}
		return out;
	}

	/** A representative oriented to the exact factor shape, with its pinned leaf. */
	record Oriented(Rep rep, NonCubicBilinearAlgorithm alg, Lineage.Node leaf) {}

	private static List<Oriented> orientedReps(int[] f, FieldAwareLookup lookup, Map<String, List<Rep>> cache)
			throws IOException {
		int[] s = f.clone();
		Arrays.sort(s);
		String key = s[0] + "x" + s[1] + "x" + s[2];
		List<Rep> reps = cache.get(key);
		if (reps == null) {
			reps = representatives(s, lookup);
			cache.put(key, reps);
		}
		List<Oriented> out = new ArrayList<>();
		for (Rep r : reps) {
			Optional<int[]> perm = r.alg.orientPermFor(f[0], f[1], f[2]);
			if (perm.isEmpty()) {
				continue;
			}
			NonCubicBilinearAlgorithm oriented = r.alg.orientByPerm(perm.get());
			Lineage.Node atom = new Lineage.Atom(r.alg.n + "x" + r.alg.m + "x" + r.alg.p + "@" + r.hash);
			Lineage.Node leaf = Arrays.equals(perm.get(), NonCubicBilinearAlgorithm.ORIENT_PERMS[0]) ? atom
					: new Lineage.OrientAs(atom, f[0], f[1], f[2], perm.get().clone());
			out.add(new Oriented(r, oriented, leaf));
		}
		return out;
	}

	/**
	 * Every explicit best-rank file of a (canonical) factor shape plus the {@link VectorProjection}
	 * variants of each, de-duplicated by content hash. Nothing is written here: a rebase is
	 * persisted by {@link #persist} only when a win pins it.
	 */
	static List<Rep> representatives(int[] s, FieldAwareLookup lookup) throws IOException {
		int best = lookup.findRank(s[0], s[1], s[2]);
		Map<String, Rep> byHash = new LinkedHashMap<>();
		for (Path p : lookup.findFilesNonCommutative(s[0], s[1], s[2])) {
			tools.jackson.databind.JsonNode root;
			NonCubicBilinearAlgorithm alg;
			try {
				root = SchemeIO.parseJson(p.toFile());
				if (SchemeIO.isStub(root)) {
					continue;
				}
				alg = SchemeIO.readBilinear(p.toFile());
			} catch (Exception e) {
				log.debug("skipping {}: {}", p.getFileName(), e.toString());
				continue;
			}
			if (alg.r != best) {
				continue;
			}
			String hash = SchemeIO.readHash(root);
			if (hash == null || hash.isBlank()) {
				hash = SchemeIO.contentHash(alg);
			}
			String source = root.hasNonNull("attribution_for_rank") ? root.get("attribution_for_rank").asText()
					: root.hasNonNull("source") ? root.get("source").asText() : p.getFileName().toString();
			Rep origin = new Rep(p.toFile(), hash, alg, null, source);
			byHash.putIfAbsent(hash, origin);
			for (NonCubicBilinearAlgorithm v : VectorProjection.variants(alg)) {
				String vh = SchemeIO.contentHash(v);
				if (byHash.containsKey(vh)) {
					continue;
				}
				byHash.put(vh, new Rep(null, vh, v, origin, source));
			}
		}
		log.info("⟨{},{},{}⟩ = {}: {} representative(s) incl. rebases — direction margins of the stored ones: {}",
				s[0], s[1], s[2], best, byHash.size(),
				byHash.values().stream().filter(r -> !r.rebased())
						.map(r -> Arrays.toString(VectorProjection.directionMargins(r.alg))).toList());
		return new ArrayList<>(byHash.values());
	}

	/** Persist a rebased representative under {@code margin-bases/sectionN/} through the contributed gate. */
	private static File writeMarginBase(Rep rep, FieldAwareLookup lookup) throws IOException {
		NonCubicBilinearAlgorithm v = rep.alg;
		Rep origin = rep.origin;
		ImportContributedSchemes.Gate g = ImportContributedSchemes.gate(v, lookup);
		if (!g.passed()) {
			throw new IllegalStateException("rebased representative of " + origin.file.getName() + " failed the gate");
		}
		String src = "GL change of basis (vector-projection margin) of " + origin.file.getName();
		ImportContributedSchemes.Options opts = new ImportContributedSchemes.Options(List.of(), NOTE, src, null, null,
				"false", origin.originSource, null,
				"Same rank, GL-equivalent representative with coordinate projection margins "
						+ Arrays.toString(ProjectionSearch.axisMargins(v)) + " (original "
						+ Arrays.toString(ProjectionSearch.axisMargins(origin.alg)) + "); see catalog.VectorProjection",
				"margin-bases", "Q", true);
		File out = ImportContributedSchemes.write(g, opts, SCHEMES_ROOT);
		FieldAwareLookup.onSchemeWritten(out);
		log.info("wrote margin base {} (margins {})", SCHEMES_ROOT.relativize(out.toPath()),
				Arrays.toString(ProjectionSearch.axisMargins(v)));
		return out;
	}

	/** Write the win's rebased factors (if any), then the win as a pinned stub (or dense file below the
	 *  materialise cap), and prove the disk artifact replays to the built scheme. */
	private static void persist(Win w, FieldAwareLookup lookup) throws IOException {
		for (Rep r : w.reps()) {
			if (r.file == null) {
				r.file = writeMarginBase(r, lookup);
			}
		}
		NonCubicBilinearAlgorithm alg = w.alg();
		if (!Verifier.passesRandomMatmulSpotCheck(alg)) {
			throw new IllegalStateException("projected scheme fails the spot check: " + w.how());
		}
		int maxDim = Math.max(alg.n, Math.max(alg.m, alg.p));
		Path dir = SCHEMES_ROOT.resolve("derived").resolve("section" + maxDim);
		Files.createDirectories(dir);
		File out = dir.resolve(SchemeIO.canonicalName(alg, "derived")).toFile();
		if (out.exists()) {
			log.info("{} already on disk", out.getName());
			return;
		}
		List<String> fields = lookup.fieldNamesFromLineage(w.lineage());
		if (maxDim > CatalogPolicy.MATERIALISE_MAX_DIM) {
			SchemeIO.writeStub(alg, out, w.lineage(), fields);
		} else {
			SchemeIO.write(alg, out, w.lineage());
			Map<String, Object> meta = new LinkedHashMap<>();
			meta.put("fields", fields);
			meta.put("source", "derived");
			meta.put("verified", true);
			SchemeIO.updateFields(out, meta, List.of(), true);
		}
		FieldAwareLookup.onSchemeWritten(out);
		NonCubicBilinearAlgorithm back = LineageReplayer.withDefaultPool(lookup).replayFromFile(out);
		if (back.r != alg.r || !SchemeIO.contentHash(back).equals(SchemeIO.contentHash(alg))) {
			Files.delete(out.toPath());
			throw new IllegalStateException("replay of " + out.getName() + " does not reproduce the built scheme (r="
					+ back.r + " vs " + alg.r + ") — file removed; lineage " + Lineage.prettyString(w.lineage()));
		}
		log.info("wrote {} (r={}, fields={}, lineage={})", SCHEMES_ROOT.relativize(out.toPath()), alg.r, fields,
				Lineage.prettyString(w.lineage()));
	}

	/** One shape per line, leading token {@code NxMxP} (the gap list's layout); comments and blanks skipped. */
	static List<int[]> readShapeFile(String path) throws IOException {
		List<int[]> out = new ArrayList<>();
		for (String line : Files.readAllLines(Path.of(path))) {
			String t = line.trim();
			if (t.isEmpty() || t.startsWith("#")) {
				continue;
			}
			out.add(parseShape(t.split("\\s+")[0]));
		}
		return out;
	}

	static int[] parseShape(String tok) {
		String[] p = tok.trim().split("[xX*,]+");
		if (p.length != 3) {
			throw new IllegalArgumentException("bad shape token: '" + tok + "'");
		}
		return new int[] { Integer.parseInt(p[0]), Integer.parseInt(p[1]), Integer.parseInt(p[2]) };
	}
}
