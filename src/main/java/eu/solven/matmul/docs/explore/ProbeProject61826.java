package eu.solven.matmul.docs.explore;

import java.io.File;
import java.util.List;
import java.util.Optional;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.ProjectionSearch;
import eu.solven.matmul.search.LineageReplayer;

/**
 * Throwaway probe for ⟨6,18,26⟩ (+2 over master). Master reaches 1726 by PROJECTING its
 * ⟨6,18,27⟩ (dropping P-col 25); the branch used a ConcatCols (1728) and never projected.
 * The branch's ⟨6,18,27⟩ is a DIFFERENT scheme (hash 8e3a6ae vs master's 9223a53). Replay
 * each parent and ask {@link ProjectionSearch#bestFor} for the best ⟨6,18,26⟩ projection:
 *  - branch parent reaches 1726 → SEARCH GAP (branch never tried projecting its own parent).
 *  - branch parent only 1728, master 1726 → WORSE PARENT (need master's high-margin one).
 */
public class ProbeProject61826 {

	public static void main(String[] args) {
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		LineageReplayer rep = LineageReplayer.withDefaultPool(lookup);

		String SC = "/private/tmp/claude-501/-Users-blacelle-workspace4-strassen/"
				+ "2949a555-de05-4f76-9ee6-f1f1ecc4c1cf/scratchpad/";
		File branchFile = new File("src/main/resources/schemes/derived/section27/"
				+ "6x18x27-r1744-derived-8e3a6ae.json");
		File masterFile = new File(SC + "master-6x18x27.json");

		probe("branch ⟨6,18,27⟩ (8e3a6ae)", rep, branchFile);
		probe("master ⟨6,18,27⟩ (9223a53)", rep, masterFile);

		// Does master's CLAIMED ⟨6,18,26⟩=1726 actually replay+verify, or is it a phantom?
		try {
			NonCubicBilinearAlgorithm child = rep.replayFromFile(new File(SC + "master-6x18x26.json"));
			boolean ok = eu.solven.matmul.Verifier.isExactNonCubic(child);
			System.out.printf("%nmaster ⟨6,18,26⟩ claim=1726 -> replayed r=%d  verifies=%s%n",
					child.r, ok);
		} catch (RuntimeException e) {
			System.out.printf("%nmaster ⟨6,18,26⟩ replay FAILED (phantom): %s%n", e);
		}
	}

	private static void probe(String label, LineageReplayer rep, File f) {
		try {
			NonCubicBilinearAlgorithm parent = rep.replayFromFile(f);
			System.out.printf("%n%s -> replayed r=%d%n", label, parent.r);
			int[] margins = ProjectionSearch.axisMargins(parent);
			System.out.printf("  axis margins (n,m,p) = [%d,%d,%d]%n",
					margins[0], margins[1], margins[2]);
			Optional<ProjectionSearch.Hit> hit =
					ProjectionSearch.bestFor(6, 18, 26, List.of(parent), 1730, 1);
			if (hit.isPresent()) {
				System.out.printf("  best ⟨6,18,26⟩ projection = %d  (master=1726, branch concat=1728)%n",
						hit.get().scheme().r);
			} else {
				System.out.println("  no projection of this parent reaches ⟨6,18,26⟩ under 1730");
			}
		} catch (RuntimeException e) {
			System.out.printf("%n%s -> replay/probe FAILED: %s%n", label, e);
		}
	}
}
