package eu.solven.matmul.docs.migrate;

import java.util.ArrayList;
import java.util.List;

/**
 * Pull new schemes from Merlijn S. Witteveen's dataset
 * <a href="https://github.com/MerlijnW70/fmm-schemes">MerlijnW70/fmm-schemes</a>
 * (doi:10.5281/zenodo.22943995, <b>CC BY 4.0</b>) — REFERENCES.md [90].
 *
 * <p>Why a dedicated channel existed: the dataset is cited by FMM-Lille (⟨11,13,15⟩,
 * ⟨11,14,14⟩) and mirrored by Perminov under {@code schemes/known/MerlijnW70_fmm_schemes/},
 * but neither route reached the catalog — we import no FMM artifact automatically,
 * and {@link ImportPerminovSchemes} only lists Perminov's own {@code schemes/results/}.
 * Release 1.0 (2026-09-24, with a new best over every ring at ⟨11,13,15⟩ and
 * ⟨11,14,14⟩) therefore sat un-imported for six days behind green sync jobs.</p>
 *
 * <p>Since 2026-10-13 this is the {@code witteveen} entry of the
 * {@link UpstreamSchemeRepos} registry; this class is {@link ImportUpstreamSchemeRepos}
 * restricted to it ({@code --repo=witteveen}), kept for the existing command lines.</p>
 */
public final class ImportWitteveenSchemes {
	private ImportWitteveenSchemes() {}

	public static void main(String[] args) throws Exception {
		List<String> forwarded = new ArrayList<>(List.of(args));
		forwarded.add("--repo=" + UpstreamSchemeRepos.WITTEVEEN.key());
		ImportUpstreamSchemeRepos.main(forwarded.toArray(new String[0]));
	}
}
