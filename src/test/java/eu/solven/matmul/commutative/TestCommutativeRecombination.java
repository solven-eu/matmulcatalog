package eu.solven.matmul.commutative;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.util.Optional;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.recombination.BlockSplitSearch;
import eu.solven.matmul.commutative.CommutativeBounds;
import eu.solven.matmul.recombination.Recombination;
import eu.solven.matmul.papers.rosowski2019.RosowskiBound;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Reproduces DIS11 Table 4 (commutative cubic ranks) via the same recombination
 * machinery used for non-commutative work, just with a commutative sota
 * ({@link CommutativeBounds}).
 *
 * <p>Setup: NC Strassen ⟨2,2,2⟩=7 outer × commutative sub-ranks. This is
 * SUB-OPTIMAL vs a fully-commutative recombination (which would use
 * Hopcroft/Winograd 1971 ⟨2,2,2⟩=6 as outer — not yet implemented), but
 * already produces interesting improvements on top of Rosowski/DIS09 for
 * many sizes.</p>
 */
public class TestCommutativeRecombination {

	@Test
	public void scan_commutative_cubic_via_strassen_recombine() throws Exception {
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		CommutativeBounds cmt = new CommutativeBounds();
		Recombination.SotaResolver sota = cmt.asSotaResolver();

		System.out.println();
		System.out.printf("%8s | %10s | %10s | %14s | best vs Rosowski/DIS09%n",
				"⟨n,n,n⟩", "Rosowski", "DIS09 T4", "Strassen-rec");
		System.out.println("-".repeat(75));
		int beats = 0, ties = 0, loses = 0;
		for (int n = 4; n <= 20; n++) {
			Optional<Long> ros = RosowskiBound.commutativeBound(n, n, n);
			Optional<Long> direct = cmt.bestRank(n, n, n);
			int directInt = direct.isPresent() ? direct.get().intValue() : Integer.MAX_VALUE / 100;
			// Forbid self-lookup in the sota to force a non-trivial split.
			int tn = n;
			Recombination.SotaResolver pureSota = (a, b, c) -> {
				if (a == tn && b == tn && c == tn) return Integer.MAX_VALUE / 100;
				return sota.getRank(a, b, c);
			};
			Optional<BlockSplitSearch.NonCubicSplitCandidate> best =
					BlockSplitSearch.findBestSplitNonCubic(n, n, n, strassen, pureSota);
			if (best.isEmpty()) continue;
			long recomb = best.get().rank();
			String marker;
			if (recomb < directInt) { marker = "✓ (Δ=" + (directInt - recomb) + ")"; beats++; }
			else if (recomb == directInt) { marker = "= (matches)"; ties++; }
			else { marker = "✗ (+" + (recomb - directInt) + ")"; loses++; }
			System.out.printf("%8s | %10d | %10d | %14d | %s%n",
					"⟨" + n + "," + n + "," + n + "⟩",
					ros.orElse(-1L),
					directInt < Integer.MAX_VALUE / 100 ? directInt : -1,
					recomb,
					marker);
		}
		System.out.println();
		System.out.println("Summary: " + beats + " beat direct, " + ties + " match, " + loses + " worse than direct.");
		assertThat(beats + ties + loses).isGreaterThan(0);
	}
}
