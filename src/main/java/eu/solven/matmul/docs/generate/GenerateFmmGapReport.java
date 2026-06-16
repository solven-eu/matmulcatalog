package eu.solven.matmul.docs.generate;

import eu.solven.matmul.catalog.Recombination;

import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;

import eu.solven.matmul.catalog.FieldAwareLookup;

/**
 * Catalog-vs-FMM gap report (2026-06-04). Diffs every format in the FMM-Lille
 * digest ({@code references/catalogs/fmm-lille-catalog.json}, 5426 formats with published
 * non-commutative ranks) against OUR on-disk best
 * ({@link FieldAwareLookup#findRank} over R — the source of truth the search
 * uses). Pure data diff: no search, no predictions, no false positives.
 *
 * <p>Emits two lists: shapes where FMM strictly beats our catalog (import/derive
 * targets), and shapes FMM lists that we have NO scheme for at all.</p>
 */
public final class GenerateFmmGapReport {

	private GenerateFmmGapReport() {}

	record Gap(int n, int m, int p, int fmm, int ours, boolean haveOurs, String refs) {}

	public static void main(String[] args) throws Exception {
		JsonMapper mapper = JsonMapper.builder().build();
		JsonNode root;
		try (FileReader r = new FileReader("references/catalogs/fmm-lille-catalog.json")) {
			root = mapper.readTree(r);
		}
		JsonNode entries = root.get("entries");
		FieldAwareLookup lk = new FieldAwareLookup("R");
		int sentinel = Recombination.SotaResolver.UNKNOWN_RANK;

		List<Gap> worse = new ArrayList<>();   // FMM < ours (we have a worse scheme)
		List<Gap> missing = new ArrayList<>(); // FMM has it, we have nothing
		int compared = 0, tieOrBetter = 0, skippedUnverified = 0;

		for (JsonNode e : entries) {
			JsonNode fmt = e.get("format");
			if (fmt == null || !fmt.isArray() || fmt.size() != 3) continue;
			int n = fmt.get(0).asInt(), m = fmt.get(1).asInt(), p = fmt.get(2).asInt();
			JsonNode rankNode = e.get("rank");
			if (rankNode == null) continue;
			int fmm = rankNode.asInt();
			// Skip degenerate / trivial (any axis 1 → naive, no scheme needed).
			if (n <= 1 || m <= 1 || p <= 1) continue;
			// Skip SYNTHESIZED FMM bounds (empty page / HK-formula, no backing
			// scheme) — they are not real FMM results and create phantom gaps.
			JsonNode ver = e.get("verified");
			if (ver != null && !ver.asBoolean()) { skippedUnverified++; continue; }

			// FMM provenance: a non-empty references list (e.g. "moosbauer:2025",
			// "smirnov:2013") means a published result; empty = FMM-derived
			// composition (concat / Kronecker / recursion), which we should be
			// able to re-derive ourselves.
			StringBuilder refb = new StringBuilder();
			JsonNode refs = e.get("references");
			if (refs != null && refs.isArray()) {
				for (JsonNode rf : refs) { if (refb.length() > 0) refb.append(';'); refb.append(rf.asString()); }
			}
			String refStr = refb.toString();

			int ours = lk.findRank(n, m, p);
			boolean haveOurs = ours < sentinel;
			compared++;
			if (!haveOurs) {
				missing.add(new Gap(n, m, p, fmm, -1, false, refStr));
			} else if (fmm < ours) {
				worse.add(new Gap(n, m, p, fmm, ours, true, refStr));
			} else {
				tieOrBetter++;
			}
		}

		worse.sort(Comparator.<Gap>comparingInt(g -> g.ours() - g.fmm()).reversed());
		missing.sort(Comparator.comparingInt(g -> g.n() * g.m() * g.p()));

		StringBuilder sb = new StringBuilder();
		sb.append("# Catalog vs FMM-Lille digest — rank gaps\n\n");
		sb.append("Pure data diff of `references/catalogs/fmm-lille-catalog.json` (FMM published NC ranks) "
				+ "vs our on-disk best (`FieldAwareLookup(\"R\").findRank`). Non-trivial formats only.\n\n");
		sb.append(String.format("- formats compared (we have a scheme): %d%n", compared - missing.size()));
		sb.append(String.format("- **FMM strictly better than us: %d** (import/derive targets)%n", worse.size()));
		sb.append(String.format("- we tie or beat FMM: %d%n", tieOrBetter));
		sb.append(String.format("- FMM has it, we have NOTHING: %d%n", missing.size()));
		sb.append(String.format("- SKIPPED synthesized/unverified FMM bounds (empty page / HK-formula): %d%n%n",
				skippedUnverified));

		long refd = worse.stream().filter(g -> !g.refs().isEmpty()).count();
		sb.append(String.format("### split of the %d FMM-better gaps by FMM provenance%n", worse.size()));
		sb.append(String.format("- with a published REFERENCE (real import/cite targets): **%d**%n", refd));
		sb.append(String.format("- NO reference = FMM-derived composition (re-derivable by us): **%d**%n%n",
				worse.size() - refd));

		sb.append("## (A) FMM-better WITH a reference — published results to import (sorted by gap)\n\n");
		sb.append("| shape | FMM | ours | gap | reference |\n| --- | ---: | ---: | ---: | --- |\n");
		for (Gap g : worse) {
			if (g.refs().isEmpty()) continue;
			sb.append(String.format("| ⟨%d,%d,%d⟩ | %d | %d | %d | %s |%n",
					g.n(), g.m(), g.p(), g.fmm(), g.ours(), g.ours() - g.fmm(), g.refs()));
		}
		sb.append("\n## (B) FMM-better with NO reference — FMM-derived, re-derivable by us (sorted by gap)\n\n");
		sb.append("| shape | FMM | ours | gap |\n| --- | ---: | ---: | ---: |\n");
		for (Gap g : worse) {
			if (!g.refs().isEmpty()) continue;
			sb.append(String.format("| ⟨%d,%d,%d⟩ | %d | %d | %d |%n",
					g.n(), g.m(), g.p(), g.fmm(), g.ours(), g.ours() - g.fmm()));
		}
		sb.append("\n## FMM lists a shape we have no scheme for (first 200 by size)\n\n");
		sb.append("| shape | FMM |\n| --- | ---: |\n");
		int cap = Math.min(200, missing.size());
		for (int i = 0; i < cap; i++) {
			Gap g = missing.get(i);
			sb.append(String.format("| ⟨%d,%d,%d⟩ | %d |%n", g.n(), g.m(), g.p(), g.fmm()));
		}
		if (missing.size() > cap) sb.append(String.format("%n_(+%d more)_%n", missing.size() - cap));

		// Tracked location (not ephemeral target/) so it is a durable hand-off
		// artifact — e.g. for the metaflip agent hunting high-projection-margin
		// bases to close the (B) FMM-derived gaps.
		File out = new File("docs/comparison/fmm-gap-report.md");
		try (java.io.FileWriter w = new java.io.FileWriter(out)) {
			w.write(sb.toString());
		}
		System.out.println("Wrote " + out);
		System.out.printf("FMM strictly better: %d | missing: %d | tie-or-better: %d%n",
				worse.size(), missing.size(), tieOrBetter);
		System.out.println("--- top FMM-better gaps ---");
		for (int i = 0; i < Math.min(25, worse.size()); i++) {
			Gap g = worse.get(i);
			System.out.printf("  ⟨%d,%d,%d⟩  FMM=%d  ours=%d  gap=%d%n",
					g.n(), g.m(), g.p(), g.fmm(), g.ours(), g.ours() - g.fmm());
		}
	}
}
