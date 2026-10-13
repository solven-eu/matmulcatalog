package eu.solven.matmul.docs.migrate;

import java.util.List;
import java.util.Optional;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * The registry of third-party GitHub repositories that publish explicit schemes as
 * Perminov-style JSON ({@code n}, {@code m}, {@code u}, {@code v}, {@code w}) and that
 * the sync job pulls on a schedule ({@link ImportUpstreamSchemeRepos}). FMM-Lille and
 * Perminov have their own channels ({@code SyncReferenceCatalogs},
 * {@code ImportPerminovSchemes}, {@code ImportFmmTensorArtifacts}); the KGP LITA
 * archives have {@code ImportKhoruzhiiLita --download}.
 *
 * <p>Adding a source = one entry here (plus a REFERENCES.md entry, a
 * {@code docs/catalog.js} source → reference mapping, and — if Perminov mirrors it —
 * a {@code PerminovKnownAttribution} case). The sync job then imports every
 * {@code (shape, rank, class)} the catalog lacks, through the contributed-scheme gate,
 * stamped {@code discovery: "TBD"} for the audit ({@code RestampContributedProvenance}).</p>
 */
public final class UpstreamSchemeRepos {
	private UpstreamSchemeRepos() {}

	/**
	 * One upstream repository.
	 *
	 * @param key      short handle ({@code --repo=key})
	 * @param repo     {@code owner/name}
	 * @param branch   branch to list
	 * @param name     pattern over the repo-relative path: groups 1–3 the shape, 4 the rank,
	 *                 optional group 5 a coefficient-class tag ({@code ZT}/{@code Z}/{@code Q});
	 *                 without a tag the class is unknown and the coarsest ({@code Q}) is
	 *                 assumed for the "already held" test
	 * @param note     filename label of the imported files
	 * @param source   {@code source} stamped on the imported files (the catalog's attribution key)
	 * @param year     {@code year} stamped
	 * @param citation how to credit the dataset (licence included) — goes into {@code discovery_note}
	 */
	public record UpstreamSchemeRepo(String key, String repo, String branch, Pattern name, String note, String source,
			int year, String citation) {

		public String treeUrl() {
			return "https://api.github.com/repos/" + repo + "/git/trees/" + branch + "?recursive=1";
		}

		public String rawBase() {
			return "https://raw.githubusercontent.com/" + repo + "/" + branch + "/";
		}

		public String blobBase() {
			return "https://github.com/" + repo + "/blob/" + branch + "/";
		}

		/** Shape / rank / class of an upstream path, or empty when it is not a scheme file. */
		public Optional<ImportUpstreamSchemeRepos.Candidate> parse(String path) {
			Matcher m = name.matcher(path);
			if (!m.matches()) return Optional.empty();
			String cls = m.groupCount() >= 5 && m.group(5) != null ? KnownSchemeKeys.classOfTag(m.group(5))
					: KnownSchemeKeys.Q;
			return Optional.of(new ImportUpstreamSchemeRepos.Candidate(path, Integer.parseInt(m.group(1)),
					Integer.parseInt(m.group(2)), Integer.parseInt(m.group(3)), Integer.parseInt(m.group(4)), cls));
		}
	}

	/** Merlijn S. Witteveen's dataset — REFERENCES [90]. */
	public static final UpstreamSchemeRepo WITTEVEEN = new UpstreamSchemeRepo("witteveen", "MerlijnW70/fmm-schemes",
			"main", Pattern.compile("^schemes/(\\d+)x(\\d+)x(\\d+)_m(\\d+)_([A-Za-z0-9]+)\\.json$"), "witteveen_2026",
			"Witteveen 2026", 2026,
			"Merlijn S. Witteveen, New fast matrix multiplication schemes, doi:10.5281/zenodo.22943995, "
					+ "github.com/MerlijnW70/fmm-schemes — licensed CC BY 4.0");

	/** The Cloudroom Math Campaign's dataset — REFERENCES [91], issue #12. */
	public static final UpstreamSchemeRepo CLOUDROOM = new UpstreamSchemeRepo("cloudroom", "matyasstoch/cloudroom-math",
			"main", Pattern.compile("^schemes/(\\d+)x(\\d+)x(\\d+)_rank(\\d+)\\.json$"), "cloudroom_2026",
			"Cloudroom Math 2026", 2026,
			"Cloudroom Math Campaign, Cloudroom Math: new matrix multiplication algorithms found by AI agents, "
					+ "github.com/matyasstoch/cloudroom-math — licensed Apache-2.0");

	public static final List<UpstreamSchemeRepo> ALL = List.of(WITTEVEEN, CLOUDROOM);

	public static UpstreamSchemeRepo byKey(String key) {
		return ALL.stream().filter(r -> r.key().equals(key)).findFirst()
				.orElseThrow(() -> new IllegalArgumentException("unknown upstream repo '" + key + "'; known: "
						+ ALL.stream().map(UpstreamSchemeRepo::key).toList()));
	}
}
