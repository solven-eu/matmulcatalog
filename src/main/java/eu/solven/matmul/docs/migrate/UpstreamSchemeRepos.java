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
 *
 * <p>A repository may also publish <em>addition-reduced</em> programs (Perminov's
 * {@code _reduced} layout: {@code u_fresh}/{@code v_fresh}/{@code w_fresh}). Those are
 * keyed by the pattern group {@link UpstreamSchemeRepo#slpGroup}: a file where that
 * group matches is pulled even when its {@code (shape, rank, class)} is already held,
 * because what it carries is the <em>program</em> — attached to the scheme we hold
 * when the matrices are identical, written as a new scheme otherwise.</p>
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
	 * @param slpGroup the pattern group that matches only on an addition-reduced (SLP) file,
	 *                 or 0 when the repository publishes none
	 * @param rankRows the files are in the rank-row sparse layout ({@code io.RankRowSchemeJson}:
	 *                 {@code w} has one row per product), not Perminov's
	 * @param explicitMaxDim largest dimension imported as an explicit catalog file. Above it a
	 *                 scheme is only <em>listed</em> (the rank digest
	 *                 {@code references/catalogs/upstream-repo-ranks.json}, read by the upstream
	 *                 scan): the repositories that set it publish hundreds of rational schemes of
	 *                 rank 2000–12000 whose canonical JSON runs to tens of MB each
	 */
	public record UpstreamSchemeRepo(String key, String repo, String branch, Pattern name, String note, String source,
			int year, String citation, int slpGroup, boolean rankRows, int explicitMaxDim) {

		/** A Perminov-layout repository, imported explicitly at every dimension. */
		public UpstreamSchemeRepo(String key, String repo, String branch, Pattern name, String note, String source,
				int year, String citation, int slpGroup) {
			this(key, repo, branch, name, note, source, year, citation, slpGroup, false, Integer.MAX_VALUE);
		}

		/** A repository without addition-reduced programs. */
		public UpstreamSchemeRepo(String key, String repo, String branch, Pattern name, String note, String source,
				int year, String citation) {
			this(key, repo, branch, name, note, source, year, citation, 0);
		}

		public String treeUrl() {
			return "https://api.github.com/repos/" + repo + "/git/trees/" + branch + "?recursive=1";
		}

		public String rawBase() {
			return "https://raw.githubusercontent.com/" + repo + "/" + branch + "/";
		}

		public String blobBase() {
			return "https://github.com/" + repo + "/blob/" + branch + "/";
		}

		/** Shape / rank / class (and whether it is an SLP file) of an upstream path, or empty when it is not a scheme file. */
		public Optional<ImportUpstreamSchemeRepos.Candidate> parse(String path) {
			Matcher m = name.matcher(path);
			if (!m.matches()) return Optional.empty();
			String cls = m.groupCount() >= 5 && m.group(5) != null ? KnownSchemeKeys.classOfTag(m.group(5))
					: KnownSchemeKeys.Q;
			boolean reduced = slpGroup > 0 && m.groupCount() >= slpGroup && m.group(slpGroup) != null;
			return Optional.of(new ImportUpstreamSchemeRepos.Candidate(path, Integer.parseInt(m.group(1)),
					Integer.parseInt(m.group(2)), Integer.parseInt(m.group(3)), Integer.parseInt(m.group(4)), cls,
					reduced));
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

	/**
	 * mizopl / we*star — REFERENCES [92], issue #20. Two directories under
	 * {@code solutions/fast-matrix-multiplication/}: {@code zt_rank_schemes/} (ternary / integer
	 * rank schemes, Perminov sparse layout, {@code {shape}_m{rank}_{ZT|Z}.json}) and
	 * {@code slp_addition_reduced/} (addition-reduced programs, Perminov {@code _reduced}
	 * layout, {@code {shape}_m{rank}_cr{reduced}_cn{naive}_ZT_reduced.json} — group 6).
	 * Perminov mirrors the rank schemes under {@code known/mizopl/} (not the programs);
	 * we pull from the origin.
	 */
	public static final UpstreamSchemeRepo MIZOPL = new UpstreamSchemeRepo("mizopl", "mizopl/we-star", "main",
			Pattern.compile("^solutions/fast-matrix-multiplication/(?:zt_rank_schemes|slp_addition_reduced)/"
					+ "(\\d+)x(\\d+)x(\\d+)_m(\\d+)(?:_cr\\d+_cn\\d+)?_([A-Za-z0-9]+?)(_reduced)?\\.json$"),
			"mizopl_2026", "mizopl 2026", 2026,
			"mizopl (we*star), Fast matrix multiplication: ZT/Z rank schemes and addition-reduced SLPs, "
					+ "github.com/mizopl/we-star (solutions/fast-matrix-multiplication) — licensed MIT", 6);

	/**
	 * Greg Sidebottom's {@code mat_mul} — REFERENCES [93]. One folder per scheme under
	 * {@code schemes/} (records up to 16, with PLinOpt programs beside the JSON), {@code large/}
	 * (deletions from the LITA squares 17³–32³, each followed by LITA's 2-reduction) and
	 * {@code products/} (shared-factor products and deletions from them); the scheme itself is
	 * {@code {shape}_m{rank}[_{ZT|Z|Q}].json.gz} in the rank-row layout. Perminov mirrors four
	 * under {@code known/gsidebottom/}; we pull from the origin.
	 */
	public static final UpstreamSchemeRepo SIDEBOTTOM = new UpstreamSchemeRepo("sidebottom", "gsidebottom/mat_mul",
			"main",
			Pattern.compile("^(?:schemes|large|products)/[^/]+/(\\d+)x(\\d+)x(\\d+)_m(\\d+)(?:_(ZT|Z|Q))?\\.json\\.gz$"),
			"sidebottom_2026", "Sidebottom 2026", 2026,
			"Greg Sidebottom, mat_mul: matrix multiplication schemes — rank and operation-count records, in PLinOpt "
					+ "format, with Lean 4 proofs, github.com/gsidebottom/mat_mul — licensed Apache-2.0",
			0, true, 16);

	/**
	 * Haotian He's oblique-projection records — REFERENCES [94], announced as Perminov's
	 * issue #18 (2026-11-24). {@code schemes_classic_le16/} (10 schemes) and
	 * {@code schemes_large/} (48, up to 30), {@code {shape}_m{rank}_{dyadic|Q}.json.gz} in the
	 * rank-row layout. The filename tag names the denominators (dyadic = ℤ[1/2]), not one of
	 * our coefficient classes — every file is a Q scheme — so it is not captured.
	 */
	public static final UpstreamSchemeRepo HE_OBLIQUE = new UpstreamSchemeRepo("he-oblique",
			"heht24-cyber/fmm-oblique-records", "main",
			Pattern.compile("^schemes_(?:classic_le16|large)/(\\d+)x(\\d+)x(\\d+)_m(\\d+)_[A-Za-z0-9]+\\.json\\.gz$"),
			"he_2026", "He 2026", 2026,
			"Haotian He, Oblique hyperplane projection deletion for fast matrix multiplication, "
					+ "github.com/heht24-cyber/fmm-oblique-records (release v1.0) — licensed MIT",
			0, true, 16);

	/**
	 * Arhan Canli's Quantavane records ⟨9,10,11⟩=640 and ⟨9,11,11⟩=708 — REFERENCES [95]. The
	 * origin ({@code arhancanli/quantavane-fmm-records}) publishes Maple tensors
	 * ({@code .mpl.bz2}); Perminov's mirror {@code known/arhancanli/} holds the two records as
	 * JSON, which is what this entry pulls. The origin's intermediate schemes (641, 642, 643,
	 * 709) are not mirrored and not imported.
	 */
	public static final UpstreamSchemeRepo CANLI = new UpstreamSchemeRepo("canli",
			"dronperminov/FastMatrixMultiplication", "master",
			Pattern.compile("^schemes/known/arhancanli/(\\d+)x(\\d+)x(\\d+)_m(\\d+)\\.json$"), "canli_2026",
			"Canli 2026", 2026,
			"Arhan Canli, Quantavane FMM records, github.com/arhancanli/quantavane-fmm-records — read from Andrew "
					+ "Perminov's mirror (FastMatrixMultiplication schemes/known/arhancanli)");

	/**
	 * Joshua Stapleton's 81 addition-reduced ternary programs — REFERENCES [96]. Published in
	 * Perminov's repository only ({@code known/stapleton/additions_reduced_ZT/}, his commit of
	 * 2026-11-12), in his {@code _reduced} layout — group 6 marks the program.
	 */
	public static final UpstreamSchemeRepo STAPLETON = new UpstreamSchemeRepo("stapleton",
			"dronperminov/FastMatrixMultiplication", "master",
			Pattern.compile("^schemes/known/stapleton/additions_reduced_ZT/"
					+ "(\\d+)x(\\d+)x(\\d+)_m(\\d+)(?:_cr\\d+_cn\\d+)?_([A-Za-z0-9]+?)(_reduced)?\\.json$"),
			"stapleton_2026", "Stapleton 2026", 2026,
			"Joshua Stapleton, addition-reduced ternary schemes, read from Andrew Perminov's "
					+ "FastMatrixMultiplication (schemes/known/stapleton/additions_reduced_ZT)", 6);

	/**
	 * Chandan Shukla's three integer decompositions — REFERENCES [97], the "Shukla 2026" FMM-Lille
	 * cites for ⟨10,19,24⟩=2723 and ⟨14,18,19⟩=2829 (the third, ⟨10,13,36⟩=2786, is past the
	 * importer's default {@code --max-dim=32}). {@code tensorDecomposition/{shape}_rank{r}.json},
	 * rank rows of {@code [index, value]} under {@code dimensions} / {@code rank}; the
	 * {@code 14x19x18} file is the same scheme in another orientation and reads as already held.
	 * No licence file; the README offers the coefficient data "for free reuse … with attribution".
	 */
	public static final UpstreamSchemeRepo SHUKLA = new UpstreamSchemeRepo("shukla",
			"chandanshukla1989/fmm_decompositions", "main",
			Pattern.compile("^tensorDecomposition/(\\d+)x(\\d+)x(\\d+)_rank(\\d+)\\.json$"), "shukla_2026",
			"Shukla 2026", 2026,
			"Chandan Shukla, Verified bilinear decompositions for three small matrix formats, "
					+ "github.com/chandanshukla1989/fmm_decompositions — coefficient data offered for free reuse with "
					+ "attribution (README)",
			0, true, Integer.MAX_VALUE);

	public static final List<UpstreamSchemeRepo> ALL = List.of(WITTEVEEN, CLOUDROOM, MIZOPL, SIDEBOTTOM, HE_OBLIQUE,
			CANLI, STAPLETON, SHUKLA);

	public static UpstreamSchemeRepo byKey(String key) {
		return ALL.stream().filter(r -> r.key().equals(key)).findFirst()
				.orElseThrow(() -> new IllegalArgumentException("unknown upstream repo '" + key + "'; known: "
						+ ALL.stream().map(UpstreamSchemeRepo::key).toList()));
	}
}
