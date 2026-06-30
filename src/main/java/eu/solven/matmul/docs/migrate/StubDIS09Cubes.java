package eu.solven.matmul.docs.migrate;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.papers.dis2009.PanTrilinearAggregation;
import eu.solven.matmul.search.LineageReplayer;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import tools.jackson.databind.JsonNode;

/**
 * One-shot migration: strip the fully-materialised {@code *dis09_Q*} cube
 * files ({@code ⟨n,n,n⟩} from the DIS09 Appendix Lemma 4 / Pan trilinear
 * aggregation formula) to lineage-only <strong>stubs</strong>.
 *
 * <p>These ~29 files total ~221&nbsp;MB on disk yet are 100% reproducible
 * from their {@code DIS09Lemma4(n=N)} lineage via
 * {@link PanTrilinearAggregation#build(int)} — which the
 * {@link LineageReplayer} already resolves. Storing the explicit
 * {@code u_sparse/v_sparse/w_sparse} blocks is pure redundancy; the recipe is
 * the scheme. (They also violate the {@code MATERIALISE_MAX_DIM=16} policy:
 * everything above 16 should be a stub.)</p>
 *
 * <p>For each file this driver (a) rebuilds via the formula and re-checks
 * shape + rank + a random spot-check, refusing to touch a file whose lineage
 * no longer reproduces it; then (b) removes the three factor blocks and sets
 * {@code scheme_type:"stub"} via {@link SchemeIO#updateFields}, preserving all
 * curated scalar metadata ({@code fields}, {@code commutative}, {@code source},
 * {@code additions}, {@code lineage}, …).</p>
 *
 * <p>Idempotent: a file already a stub (no {@code u_sparse}) is skipped.
 * Pass {@code --apply} to write; default is a dry run.</p>
 *
 * <pre>
 *   mvn -q -o exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.StubDIS09Cubes -Dexec.args="--apply"
 * </pre>
 */
public final class StubDIS09Cubes {
	private StubDIS09Cubes() {}

	// The factor blocks (the bulk) plus realization-specific, recomputable
	// derived data (the SLP addition schedule) — all reproduced by replay /
	// EnrichSchemeMetrics, so none belongs in a lineage-only stub.
	private static final List<String> STRIP_KEYS = List.of(
			"u_sparse", "v_sparse", "w_sparse", "u", "v", "w",
			"slp", "min_additions", "additions_optimal");

	public static void main(String[] args) throws Exception {
		boolean apply = List.of(args).contains("--apply");
		Path root = Path.of("src/main/resources/schemes");
		LineageReplayer replayer = LineageReplayer.withDefaultPool(new FieldAwareLookup(Field.Q));

		List<Path> files;
		try (var s = Files.walk(root)) {
			files = s.filter(p -> p.getFileName().toString().matches(".*dis09_Q.*\\.json"))
					.sorted().toList();
		}

		long before = 0, after = 0;
		int stubbed = 0, already = 0, refused = 0;
		for (Path p : files) {
			File f = p.toFile();
			long sz = Files.size(p);
			before += sz;
			JsonNode rn = SchemeIO.parseJson(f);
			if (SchemeIO.isStub(rn) || !(rn.has("u_sparse") || rn.has("u"))) {
				already++;
				after += sz;
				continue;
			}
			int declRank = rn.get("m").asInt();
			int n = rn.get("n").get(0).asInt();
			Lineage.Node lin = SchemeIO.readLineage(rn).orElse(null);

			// Safety gate: the lineage must reproduce a verifying scheme of the
			// declared shape & rank, else we'd be deleting irrecoverable data.
			NonCubicBilinearAlgorithm alg;
			try {
				alg = replayer.replay(lin);
			} catch (RuntimeException e) {
				System.out.printf("REFUSE %s — replay failed: %s%n", f.getName(), e.getMessage());
				refused++;
				after += sz;
				continue;
			}
			boolean ok = alg.isCubic() && alg.n == n && alg.r == declRank
					&& Verifier.passesRandomMatmulSpotCheck(alg);
			if (!ok) {
				System.out.printf("REFUSE %s — replay mismatch (got ⟨%d,%d,%d⟩ r=%d)%n",
						f.getName(), alg.n, alg.m, alg.p, alg.r);
				refused++;
				after += sz;
				continue;
			}

			boolean changed = SchemeIO.updateFields(f,
					Map.of("scheme_type", "stub"), STRIP_KEYS, apply);
			long newSz = apply ? Files.size(p) : estimateStubSize(f);
			after += newSz;
			stubbed++;
			System.out.printf("%-6s %-26s n=%-2d r=%-6d  %8d → %5d bytes%n",
					apply ? "STUB" : "(dry)", f.getName(), n, declRank, sz, newSz);
			if (!changed && apply) {
				System.out.printf("   note: %s reported no change%n", f.getName());
			}
		}
		System.out.printf("%n%s: %d stubbed, %d already-stub, %d refused; disk %.1f MB → %.1f MB (saved %.1f MB)%n",
				apply ? "APPLIED" : "DRY RUN", stubbed, already, refused,
				before / 1048576.0, after / 1048576.0, (before - after) / 1048576.0);
		if (!apply) {
			System.out.println("Re-run with --apply to write.");
		}
	}

	/** Rough stub size for the dry-run report (the real write is measured directly). */
	private static long estimateStubSize(File f) {
		return 600;  // a DIS09 stub is ~0.5 KB (n/rank/fields/source/lineage)
	}
}
