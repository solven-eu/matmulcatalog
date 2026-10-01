package eu.solven.matmul.docs.migrate;

import java.io.File;
import java.nio.file.Path;
import java.time.LocalDate;
import java.time.ZoneOffset;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.io.MapleTensorParser;

import lombok.extern.slf4j.Slf4j;

/**
 * Import FMM-Lille <b>tensor artifacts</b> — the explicit scheme behind an index rank —
 * for the formats where the catalog is below FMM and FMM publishes the matrices.
 *
 * <p>When a cross-check row is WORSE, FMM's per-shape artifact
 * ({@code https://fmm.univ-lille.fr/{n}x{m}x{p}_tensor.mpl.bz2}) settles what kind of
 * gap it is (see {@code references/fmm-artifact-audit.md}): the artifact's rank is
 * above ours (the index is not backed — no gap), it is a recipe-only placeholder, or
 * it is a real scheme at the index rank. This driver lands the third kind: parse
 * ({@link MapleTensorParser}), run the contributed-scheme gate
 * ({@link ImportContributedSchemes#gate}: exact symbolic proof, spot-check, fields
 * from the coefficients), write the canonical file under {@code schemes/known/}.</p>
 *
 * <p>Input: already-decompressed {@code *_tensor.mpl} files (the JDK has no bzip2):</p>
 * <pre>
 *   curl -sLO https://fmm.univ-lille.fr/7x11x30_tensor.mpl.bz2 &amp;&amp; bunzip2 7x11x30_tensor.mpl.bz2
 *   mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.ImportFmmTensorArtifacts \
 *       -Dexec.args="7x11x30_tensor.mpl [--execute]"
 * </pre>
 *
 * <p>Dry-run by default. A file is skipped — never written — when it is a placeholder or
 * a parametric family, fails the gate, or does not beat the catalog's rank over Q (an artifact at or above
 * our rank adds nothing; pass {@code --keep-ties} to load a tie as an explicit
 * witness). Attribution: FMM-Lille is an aggregator, so {@code source} is
 * {@code "FMM-Lille"} and {@code discovery} is {@code "TBD"};
 * {@code GenerateCatalogManifest} resolves the originator from the FMM digest's
 * {@code references[]} where the page names one.</p>
 */
@Slf4j
public final class ImportFmmTensorArtifacts {
	private ImportFmmTensorArtifacts() {}

	private static final Path SCHEMES_ROOT = Path.of("src/main/resources/schemes");
	static final String NOTE = "fmm_lille";
	static final String SOURCE = "FMM-Lille";

	public static void main(String[] args) throws Exception {
		boolean execute = List.of(args).contains("--execute");
		boolean keepTies = List.of(args).contains("--keep-ties");
		List<File> inputs = new ArrayList<>();
		for (String a : args) if (!a.startsWith("--")) inputs.add(new File(a));
		if (inputs.isEmpty()) {
			throw new IllegalArgumentException("usage: ImportFmmTensorArtifacts <shape>_tensor.mpl... [--execute] [--keep-ties]");
		}
		FieldAwareLookup lookup = new FieldAwareLookup(Field.fromTag("Q"), SCHEMES_ROOT);
		int written = 0, skipped = 0, failed = 0;
		long t0 = System.nanoTime();
		for (File in : inputs) {
			long s = System.nanoTime();
			NonCubicBilinearAlgorithm alg;
			try {
				alg = MapleTensorParser.parse(in);
			} catch (MapleTensorParser.UnsupportedArtifactException e) {
				log.warn("{}: NOT A SCHEME — {}", in.getName(), e.getMessage());
				skipped++;
				continue;
			}
			ImportContributedSchemes.Gate g = ImportContributedSchemes.gate(alg, lookup);
			log.info("{} → ⟨{},{},{}⟩ r={}  exact={} spot={}  fields={}{}  catalog(Q): {}  ({} ms)", in.getName(), alg.n,
					alg.m, alg.p, alg.r, g.exactSymbolic(), g.spotCheck(), g.fields(), g.ternary() ? " ZT" : "",
					g.vsCatalog(), (System.nanoTime() - s) / 1_000_000L);
			if (!g.passed()) {
				log.error("{}: GATE FAILED — not imported", in.getName());
				failed++;
				continue;
			}
			boolean beats = g.catalogBest() == null || alg.r < g.catalogBest();
			boolean ties = g.catalogBest() != null && alg.r == g.catalogBest();
			if (!beats && !(ties && keepTies)) {
				log.info("{}: artifact rank {} does not beat the catalog's {} — nothing to import", in.getName(), alg.r,
						g.catalogBest());
				skipped++;
				continue;
			}
			if (!execute) {
				log.info("{}: would import (DRY-RUN, pass --execute)", in.getName());
				continue;
			}
			String shape = alg.n + "x" + alg.m + "x" + alg.p;
			File out = ImportContributedSchemes.write(g, options(g, shape));
			Map<String, Object> provenance = new LinkedHashMap<>();
			provenance.put("source_scheme_url", "https://fmm.univ-lille.fr/" + shape + "_tensor.mpl.bz2");
			SchemeIO.addFields(out, provenance, /* apply */ true);
			written++;
			log.info("[progress] {}/{} — wrote {} ({} ms elapsed)", written, inputs.size(),
					SCHEMES_ROOT.relativize(out.toPath()), (System.nanoTime() - t0) / 1_000_000L);
		}
		log.info("ImportFmmTensorArtifacts: {} inputs, {} written, {} skipped, {} failed ({})", inputs.size(), written,
				skipped, failed, execute ? "--execute" : "DRY-RUN");
		if (failed > 0) {
			throw new IllegalStateException(failed + " artifact(s) failed the import gate");
		}
	}

	private static ImportContributedSchemes.Options options(ImportContributedSchemes.Gate g, String shape) {
		String today = LocalDate.now(ZoneOffset.UTC).toString();
		String discoveryNote = "FMM-Lille tensor artifact " + shape + "_tensor.mpl.bz2 (fmm.univ-lille.fr), " + g.alg().r
				+ " triads, downloaded and imported " + today + " by ImportFmmTensorArtifacts (exact BigInteger symbolic"
				+ " proof, random spot-check, mod-2 / mod-3 reductions, canonical-disk round-trip). Against this catalog"
				+ " over Q at import time: " + g.vsCatalog() + ". FMM-Lille is an aggregator: who first reached this rank"
				+ " is on the per-format page when it names a source — `discovery` stays TBD until that is read.";
		String notes = "Re-encoded from FMM-Lille's Maple Triad listing (third matrix p×n, transposed here) into the"
				+ " catalog's canonical format — no coefficient was changed.";
		return new ImportContributedSchemes.Options(List.of(), NOTE, SOURCE, LocalDate.now(ZoneOffset.UTC).getYear(),
				"https://fmm.univ-lille.fr/" + shape + ".html", "TBD", null, discoveryNote, notes, "known", "Q", true);
	}
}
