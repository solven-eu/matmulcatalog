package eu.solven.matmul.catalog;

import java.util.Map;
import java.util.Optional;
import java.util.Set;

/**
 * Re-attribution of schemes that live under Perminov's
 * <a href="https://github.com/dronperminov/FastMatrixMultiplication">FastMatrixMultiplication</a>
 * {@code schemes/known/<sub-source>/} subtree.
 *
 * <p>Perminov files a scheme under {@code known/<sub>} precisely because it
 * <b>originates elsewhere</b> — he merely re-encodes / mirrors it (his README
 * "Analyzed Schemes &amp; Data Sources" names the upstream for each sub-folder).
 * Crediting {@code "Perminov 2023"} for those misattributes the historical
 * record: e.g. {@code ⟨2,7,7⟩=76} sits under {@code known/meta_flip_graph/} and is
 * <b>Kauers &amp; Wood 2025</b>'s meta-flip-graph result (arXiv:2510.19787), NOT
 * Perminov's. Only {@code schemes/results/*} is Perminov's own work.</p>
 *
 * <p>This is the Perminov analogue of the FMM-Lille aggregator exception in
 * {@code GenerateCatalogManifest}: an aggregator's republication is attributed to
 * the true originator, not to the aggregator. The link back to Perminov's file is
 * preserved separately in {@code source_scheme_url} / {@code original_source_path}
 * (and {@code imported_via}), so provenance is not lost — only the {@code source}
 * (= discoverer) is corrected.</p>
 *
 * <p>The mapping is keyed off the {@code original_source_path} (Perminov's own
 * directory layout), NOT the local filename — content/provenance drives metadata
 * (the repo's "read from CONTENT, never the filename" rule).</p>
 */
public final class PerminovKnownAttribution {
	private PerminovKnownAttribution() {}

	/** What to do with a scheme found at a given Perminov path. */
	public enum Disposition {
		/** Originates with an external author — {@code source} becomes that author. */
		EXTERNAL,
		/** Perminov's own discovery ({@code schemes/results/*}) — keep "Perminov 2023". */
		PERMINOV_OWN,
		/**
		 * A redundant mirror we already obtain by other means — {@code tensor}
		 * (FMM, pulled via our own FMM channel) and {@code matmulcatalog} (this
		 * very repo). Importers must NOT freshly pull these; existing on-disk
		 * copies are re-attributed (to their true origin), never silently deleted
		 * — they may be load-bearing recombination bases for derived schemes.
		 */
		SKIP_FRESH_IMPORT
	}

	/** Corrected attribution for a Perminov-repo path. */
	public record Attribution(String source, Disposition disposition) {
		public boolean isPerminovOwn() {
			return disposition == Disposition.PERMINOV_OWN;
		}
	}

	/**
	 * Map a Perminov-repo-relative path to its corrected attribution.
	 *
	 * @param originalSourcePath e.g. {@code "schemes/known/meta_flip_graph/277/k….m"}
	 *        or {@code "schemes/results/Z/2x7x7_m76_….json"}. {@code null}/blank →
	 *        empty (caller keeps the existing source).
	 * @return the corrected attribution, or empty when the path carries no signal.
	 */
	public static Optional<Attribution> forPath(String originalSourcePath) {
		if (originalSourcePath == null || originalSourcePath.isBlank()) {
			return Optional.empty();
		}
		String p = originalSourcePath.replace('\\', '/');
		int k = p.indexOf("known/");
		if (k < 0) {
			int s = p.indexOf("schemes/");
			String thirdParty = THIRD_PARTY_IN_RESULTS.get(s < 0 ? p : p.substring(s));
			if (thirdParty != null) {
				return external(thirdParty);
			}
			// Anything else NOT under known/ (notably schemes/results/*) is Perminov's own.
			return Optional.of(new Attribution("Perminov 2023", Disposition.PERMINOV_OWN));
		}
		String rest = p.substring(k + "known/".length());
		int slash = rest.indexOf('/');
		String sub = slash < 0 ? rest : rest.substring(0, slash);
		String file = rest.substring(rest.lastIndexOf('/') + 1);
		switch (sub) {
		case "meta_flip_graph":
			// arXiv:2510.19787 "Exploring the Meta Flip Graph for Matrix Multiplication"
			// (Oct 2025). NOTE: the *meta* flip graph is 2025; the *original* flip graph
			// is Kauers-Moosbauer 2022/2023 (the jakobmoosbauer_flips folder below).
			return external("Kauers & Wood 2025");
		case "alpha_tensor":
			return external("AlphaTensor 2022");
		case "alpha_evolve":
			return external("AlphaEvolve 2025");
		case "a_60_addition":
			return external("Stapleton 2025 (a=60)");
		case "jakobmoosbauer_flips":
			return external("Kauers-Moosbauer 2023");
		case "jakobmoosbauer_symmetric_flips":
			return external("Moosbauer-Poole 2025");
		case "fmm_add_reduction":
			return external("FMM add-reduction (werekorren)");
		case "classic":
			// classic/<Author>-<dims>-<rank>…  — the filename names the literature author.
			return external(classicAuthor(file));
		case "tensor":
			// FMM tensor-format dump (we ingest FMM by other means) → attribute to the
			// aggregator; GenerateCatalogManifest's FMM-Lille exception further resolves
			// it to the per-format literature citation where the FMM digest has one.
			return Optional.of(new Attribution("FMM-Lille", Disposition.SKIP_FRESH_IMPORT));
		case "matmulcatalog":
			// Our OWN catalog, re-hosted by Perminov — circular; never re-import.
			return Optional.of(new Attribution("MatMulCatalog (this repo)", Disposition.SKIP_FRESH_IMPORT));
		case "MerlijnW70_fmm_schemes":
			// M. S. Witteveen's dataset (doi:10.5281/zenodo.22943995, CC BY 4.0), mirrored by
			// Perminov since 2026-09-29. We pull it from the ORIGIN repo
			// (docs.migrate.ImportWitteveenSchemes, sync-reference-catalogs.yml) — the
			// mirror lags a release behind.
			return Optional.of(new Attribution("Witteveen 2026", Disposition.SKIP_FRESH_IMPORT));
		case "lita":
			// Khoruzhii–Serafin–Gelß–Pokutta LITA cubes (github.com/khoruzhii/lita), 13³–16³
			// re-encoded by Perminov. Ours come from the origin archives
			// (docs.migrate.ImportKhoruzhiiLita, TA_lita_npz stubs).
			return Optional.of(new Attribution("Khoruzhii 2026", Disposition.SKIP_FRESH_IMPORT));
		default:
			// An unseen known/<sub>: still NOT Perminov's. Fall back to the aggregator
			// label so it is at least not miscredited to Perminov.
			return Optional.of(new Attribution("FMM-Lille", Disposition.EXTERNAL));
		}
	}

	/**
	 * Every {@code schemes/known/<sub>} the switch in {@link #forPath} maps explicitly.
	 * A sub-folder outside this set is a NEW third-party source Perminov started
	 * mirroring: {@code forPath} still answers (aggregator fallback, never "Perminov"),
	 * but nobody imports it — {@code ImportPerminovSchemes} only lists
	 * {@code schemes/results/*}. That is how {@code MerlijnW70_fmm_schemes}
	 * (⟨11,13,15⟩=1364, a new best over every ring) and {@code lita} sat upstream
	 * un-imported behind a green sync job (2026-09-29/30). The importer therefore
	 * WARNS on every sub-folder {@link #isMappedSubtree} does not know.
	 */
	public static final Set<String> MAPPED_SUBTREES = Set.of("meta_flip_graph", "alpha_tensor", "alpha_evolve",
			"a_60_addition", "jakobmoosbauer_flips", "jakobmoosbauer_symmetric_flips", "fmm_add_reduction",
			"classic", "tensor", "matmulcatalog", "MerlijnW70_fmm_schemes", "lita");

	public static boolean isMappedSubtree(String sub) {
		return MAPPED_SUBTREES.contains(sub);
	}

	/**
	 * Third-party schemes Perminov filed under his OWN {@code schemes/results/} tree, so
	 * the directory rule would credit them to him. Keyed by repo-relative path.
	 *
	 * <p>Witteveen's five ⟨2,p,n⟩ ternary schemes: published in
	 * {@code MerlijnW70/fmm-schemes} release 1.0 (2026-09-24); Perminov's commit
	 * {@code f9f91b21} "add new schemes from @MerlijnW70" (2026-09-29) put twenty of that
	 * dataset under {@code known/MerlijnW70_fmm_schemes/} and these five under
	 * {@code results/ZT/}. Content hashes match the origin files for m280 / m304 / m324 /
	 * m348 (checked 2026-09-30; the origin's m374 was replaced in release 1.2). Our sync
	 * imported m280 and m348 as "Perminov 2023" within the half-hour — the
	 * misattribution this table repairs (via {@code AttributePerminovKnown}) and prevents.</p>
	 */
	private static final Map<String, String> THIRD_PARTY_IN_RESULTS = Map.of(
			"schemes/results/ZT/2x12x15_m280_ZT.json", "Witteveen 2026",
			"schemes/results/ZT/2x13x15_m304_ZT.json", "Witteveen 2026",
			"schemes/results/ZT/2x13x16_m324_ZT.json", "Witteveen 2026",
			"schemes/results/ZT/2x14x16_m348_ZT.json", "Witteveen 2026",
			"schemes/results/ZT/2x15x16_m374_ZT.json", "Witteveen 2026");

	private static Optional<Attribution> external(String source) {
		return Optional.of(new Attribution(source, Disposition.EXTERNAL));
	}

	/** Resolve {@code classic/<Author>-…} to the literature source (year embedded so
	 *  {@code yearOfSource} can derive the chronology). Unknown author → FMM-Lille,
	 *  which the manifest then resolves per-format. */
	static String classicAuthor(String file) {
		String f = file.toLowerCase();
		if (f.startsWith("strassen")) return "Strassen 1969";
		if (f.startsWith("winograd")) return "Winograd 1971";
		if (f.startsWith("hopcroft")) return "Hopcroft-Kerr 1971";
		if (f.startsWith("laderman")) return "Laderman 1976";
		if (f.startsWith("makarov")) return "Makarov 1986";
		if (f.startsWith("schonhage") || f.startsWith("schönhage")) return "Schönhage 1981";
		if (f.startsWith("smirnov")) return "Smirnov 2013";
		if (f.startsWith("pan")) return "Pan 1978";
		return "FMM-Lille";
	}
}
