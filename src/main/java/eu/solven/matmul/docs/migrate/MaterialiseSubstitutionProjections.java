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
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.catalog.SubstitutionProjection;
import eu.solven.matmul.search.LineageReplayer;
import eu.solven.matmul.verifiers.Verifier;

import lombok.extern.slf4j.Slf4j;

/**
 * Replay FMM-Lille's {@code ProjOfAlongMode} lemma over the catalog ({@link SubstitutionProjection}):
 * for each target ⟨n,m,p⟩ and each parent one index larger on one, two or three axes that we
 * hold on disk, drop the extra indices by substitution — every axis order, and on each axis the
 * index whose drop leaves the fewest products — and persist a strict improvement as a
 * {@code SubstProject(…(parent@hash))} stub that replays bit-for-bit.
 *
 * <p>This is the operator behind the 61 FMM index moves of 2026-11-04..06 (issue #23):
 * ⟨9,17,18⟩:1591 = ⟨9,18,18⟩:1600 − 9, ⟨9,17,17⟩:1574 = 1591 − 17, ⟨6,21,23⟩:1779 = 1800 − 21,
 * ⟨13,14,14⟩:1579 = 1593 − 14, ⟨16,17,17⟩:2706 = 2723 − 17, … — all from products the catalog
 * already held, none reachable by the coordinate projection the closure runs (0 saved on every
 * axis for those Kronecker stubs). Dropping the axes in another order can beat FMM's figure
 * (⟨9,17,17⟩: p then m = 1573). Tier of every rank written: <b>bound</b> (greedy index per
 * step, one drop per axis).</p>
 *
 * <pre>
 *   mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.MaterialiseSubstitutionProjections \
 *       -Dexec.args="--shapes-file=generated/upstream-gaps.txt [--execute]"
 * </pre>
 * Dry-run by default. Heavy parents (rank ~10k) take a few seconds each; run in the background
 * with the heap-dump flags.
 */
@Slf4j
public final class MaterialiseSubstitutionProjections {
	private MaterialiseSubstitutionProjections() {}

	private static final Path SCHEMES_ROOT = Path.of("src/main/resources/schemes");

	/** A replayable on-disk parent, oriented to the shape the target needs. */
	record Parent(Path file, NonCubicBilinearAlgorithm alg, Lineage.Node leaf) {}

	record Win(int[] shape, int ours, NonCubicBilinearAlgorithm alg, Lineage.Node lineage, String how) {}

	public static void main(String[] args) throws Exception {
		List<int[]> targets = new ArrayList<>();
		boolean execute = false;
		for (String a : args) {
			if (a.equals("--execute")) {
				execute = true;
			} else if (a.startsWith("--shapes-file=")) {
				targets.addAll(MaterialiseVectorProjections.readShapeFile(a.substring("--shapes-file=".length())));
			} else if (a.startsWith("--shapes=")) {
				for (String tok : a.substring("--shapes=".length()).split("[;,]+")) {
					if (!tok.isBlank()) {
						targets.add(MaterialiseVectorProjections.parseShape(tok));
					}
				}
			} else {
				throw new IllegalArgumentException("unknown arg: " + a);
			}
		}
		if (targets.isEmpty()) {
			throw new IllegalArgumentException("usage: --shapes-file=PATH | --shapes=NxMxP,… [--execute]");
		}
		// Smallest first: a freshly written ⟨9,17,18⟩ serves as the parent of ⟨9,17,17⟩ in the same run.
		targets.sort(Comparator.comparingLong(t -> (long) t[0] * t[1] * t[2]));
		log.info("MaterialiseSubstitutionProjections: {} target(s), {}", targets.size(), execute ? "--execute" : "DRY-RUN");

		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		LineageReplayer replayer = LineageReplayer.withDefaultPool(lookup);
		Map<String, Optional<Parent>> parentCache = new LinkedHashMap<>();
		List<Win> wins = new ArrayList<>();
		long t0 = System.currentTimeMillis();
		int done = 0;
		for (int[] t : targets) {
			done++;
			int ours = lookup.findRank(t[0], t[1], t[2]);
			boolean held = ours > 0 && ours < eu.solven.matmul.recombination.Recombination.SotaResolver.UNKNOWN_RANK;
			long upper = held ? ours : Long.MAX_VALUE;
			Optional<Win> w = bestFor(t, upper, lookup, replayer, parentCache);
			if (w.isPresent()) {
				wins.add(w.get());
				log.info("[win] ⟨{},{},{}⟩ {} → {}  via {}", t[0], t[1], t[2], held ? ours : "-", w.get().alg().r,
						w.get().how());
				if (execute) {
					persist(w.get(), lookup, replayer);
				}
			} else {
				log.info("⟨{},{},{}⟩ = {}: no parent drops below it", t[0], t[1], t[2], held ? ours : "-");
			}
			log.info("[progress] {}/{} targets, {} win(s), {}ms elapsed", done, targets.size(), wins.size(),
					System.currentTimeMillis() - t0);
		}
		log.info("MaterialiseSubstitutionProjections done: {} win(s) over {} target(s) ({})", wins.size(),
				targets.size(), execute ? "written" : "DRY-RUN, nothing written");
		for (Win w : wins) {
			System.out.printf("%dx%dx%d  %d -> %d  %s%n", w.shape()[0], w.shape()[1], w.shape()[2], w.ours(), w.alg().r,
					w.how());
		}
	}

	/** Best substitution-projection of any one-larger parent down to {@code t}, strictly below {@code upper}. */
	static Optional<Win> bestFor(int[] t, long upper, FieldAwareLookup lookup, LineageReplayer replayer,
			Map<String, Optional<Parent>> parentCache) {
		Win best = null;
		for (int dn = 0; dn <= 1; dn++) {
			for (int dm = 0; dm <= 1; dm++) {
				for (int dp = 0; dp <= 1; dp++) {
					if (dn + dm + dp == 0) {
						continue;
					}
					int[] P = { t[0] + dn, t[1] + dm, t[2] + dp };
					int pr = lookup.findRank(P[0], P[1], P[2]);
					if (pr <= 0 || pr >= eu.solven.matmul.recombination.Recombination.SotaResolver.UNKNOWN_RANK) {
						continue;
					}
					// The lemma guarantees roughly these savings whatever the order (incidental kills can add
					// a few more): skip a parent that is hopeless even with twice the guaranteed saving.
					long guaranteed = (dn > 0 ? Math.min(P[1], P[2]) : 0) + (dm > 0 ? Math.min(P[0], P[2]) : 0)
							+ (dp > 0 ? Math.min(P[0], P[1]) : 0);
					if (upper != Long.MAX_VALUE && pr - 2L * guaranteed >= upper) {
						continue;
					}
					Optional<Parent> parent = parentCache.computeIfAbsent(P[0] + "x" + P[1] + "x" + P[2],
							k -> resolveParent(P, lookup, replayer));
					if (parent.isEmpty()) {
						continue;
					}
					List<Integer> axes = new ArrayList<>();
					if (dn > 0) axes.add(0);
					if (dm > 0) axes.add(1);
					if (dp > 0) axes.add(2);
					for (List<Integer> order : permutations(axes)) {
						NonCubicBilinearAlgorithm cur = parent.get().alg();
						Lineage.Node node = parent.get().leaf();
						StringBuilder how = new StringBuilder(P[0] + "x" + P[1] + "x" + P[2] + "=" + cur.r);
						for (int axis : order) {
							int idx = SubstitutionProjection.bestIndex(cur, axis);
							cur = SubstitutionProjection.drop(cur, axis, idx);
							node = new Lineage.SubstProject(node, axis, idx);
							how.append(" ⊖").append("nmp".charAt(axis)).append(idx).append("→").append(cur.r);
						}
						if (cur.r < upper && (best == null || cur.r < best.alg().r)) {
							best = new Win(t, upper == Long.MAX_VALUE ? -1 : (int) upper, cur, node, how.toString());
						}
					}
				}
			}
		}
		return Optional.ofNullable(best);
	}

	private static List<List<Integer>> permutations(List<Integer> items) {
		List<List<Integer>> out = new ArrayList<>();
		if (items.size() <= 1) {
			out.add(new ArrayList<>(items));
			return out;
		}
		for (int i = 0; i < items.size(); i++) {
			List<Integer> rest = new ArrayList<>(items);
			Integer head = rest.remove(i);
			for (List<Integer> tail : permutations(rest)) {
				List<Integer> perm = new ArrayList<>();
				perm.add(head);
				perm.addAll(tail);
				out.add(perm);
			}
		}
		return out;
	}

	/**
	 * The cheapest on-disk scheme of shape {@code P} that replays to at most its recorded rank,
	 * oriented to {@code P}, with a leaf pinning the exact file ({@code shape@hash}, wrapped in an
	 * explicit {@code OrientAs} when the file's native orientation differs). The same discipline as
	 * the materialiser's parent resolution: never a bare shape ref (projection is parent-sensitive).
	 */
	static Optional<Parent> resolveParent(int[] P, FieldAwareLookup lookup, LineageReplayer replayer) {
		List<Path> files = new ArrayList<>(lookup.findFilesNonCommutative(P[0], P[1], P[2]));
		Map<Path, Integer> recorded = new LinkedHashMap<>();
		for (Path f : files) {
			try {
				recorded.put(f, SchemeIO.parseJson(f.toFile()).path("m").asInt(Integer.MAX_VALUE));
			} catch (Exception e) {
				recorded.put(f, Integer.MAX_VALUE);
			}
		}
		files.sort(Comparator.comparingInt(recorded::get));
		for (Path f : files) {
			try {
				tools.jackson.databind.JsonNode root = SchemeIO.parseJson(f.toFile());
				NonCubicBilinearAlgorithm alg = replayer.replayFromFile(f.toFile());
				if (alg.r > recorded.get(f)) {
					log.warn("{} replays to r={} above its recorded {} — skipped as a parent", f.getFileName(), alg.r,
							recorded.get(f));
					continue;
				}
				String hash = SchemeIO.readHash(root);
				if (hash == null || hash.isBlank()) {
					hash = SchemeIO.contentHash(alg);
				}
				Lineage.Node atom = new Lineage.Atom(alg.n + "x" + alg.m + "x" + alg.p + "@" + hash);
				Optional<int[]> perm = alg.orientPermFor(P[0], P[1], P[2]);
				if (perm.isEmpty()) {
					continue;
				}
				NonCubicBilinearAlgorithm oriented = alg.orientByPerm(perm.get());
				Lineage.Node leaf = Arrays.equals(perm.get(), NonCubicBilinearAlgorithm.ORIENT_PERMS[0]) ? atom
						: new Lineage.OrientAs(atom, P[0], P[1], P[2], perm.get().clone());
				return Optional.of(new Parent(f, oriented, leaf));
			} catch (Exception e) {
				log.warn("{} does not replay ({}) — skipped as a parent", f.getFileName(), e.toString());
			}
		}
		return Optional.empty();
	}

	/** Write the win (stub above the materialise cap, dense below) with coefficient-derived fields, and
	 *  prove the disk artifact replays to the built scheme. */
	static void persist(Win w, FieldAwareLookup lookup, LineageReplayer replayer) throws IOException {
		NonCubicBilinearAlgorithm alg = w.alg();
		// Verification: the sampled identity check (the documented big-shape path — the symbolic
		// expansion behind Verifier.isExactNonCubic is quadratic in the factor density and the
		// substitution makes V dense: >6 min on one ⟨13,14,14⟩). The lemma is exact by construction
		// and TestSubstitutionProjection proves it symbolically on small instances.
		if (!Verifier.passesRandomMatmulSpotCheck(alg, 8, 1e-6)) {
			throw new IllegalStateException("substitution-projected scheme fails the spot check: " + w.how());
		}
		List<String> fields = fieldsOf(alg);
		int maxDim = Math.max(alg.n, Math.max(alg.m, alg.p));
		Path dir = SCHEMES_ROOT.resolve("derived").resolve("section" + maxDim);
		Files.createDirectories(dir);
		File out = dir.resolve(SchemeIO.canonicalName(alg, "derived")).toFile();
		if (out.exists()) {
			log.info("{} already on disk", out.getName());
			return;
		}
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
		NonCubicBilinearAlgorithm back = replayer.replayFromFile(out);
		if (back.r != alg.r || !SchemeIO.contentHash(back).equals(SchemeIO.contentHash(alg))) {
			Files.delete(out.toPath());
			throw new IllegalStateException("replay of " + out.getName() + " does not reproduce the built scheme (r="
					+ back.r + " vs " + alg.r + ") — file removed; lineage " + Lineage.prettyString(w.lineage()));
		}
		log.info("wrote {} (r={}, fields={}, lineage={})", SCHEMES_ROOT.relativize(out.toPath()), alg.r, fields,
				Lineage.prettyString(w.lineage()));
	}

	/**
	 * Fields from the built coefficients, in the catalog's order: the substitution solves a linear
	 * system, so an integer parent can yield a rational child ({@code Q} but not {@code Z}); validity
	 * mod 2 / mod 3 is the sampled identity with denominators inverted mod p (false when a denominator
	 * is not invertible — exactly the CLAUDE.md "1/2 works in F₃ but not F₂" rule).
	 */
	static List<String> fieldsOf(NonCubicBilinearAlgorithm alg) {
		List<String> fields = new ArrayList<>();
		if (Verifier.passesRandomMatmulSpotCheckFp(alg, 2, 8)) fields.add("F2");
		if (Verifier.passesRandomMatmulSpotCheckFp(alg, 3, 8)) fields.add("F3");
		if (allIntegers(alg)) fields.add("Z");
		fields.add("Q");
		fields.add("R");
		fields.add("C");
		return fields;
	}

	private static boolean allIntegers(NonCubicBilinearAlgorithm alg) {
		for (double[][] f : List.of(alg.denseU(), alg.denseV(), alg.denseW())) {
			for (double[] row : f) {
				for (double x : row) {
					if (Math.abs(x - Math.rint(x)) > 1e-9) {
						return false;
					}
				}
			}
		}
		return true;
	}
}
