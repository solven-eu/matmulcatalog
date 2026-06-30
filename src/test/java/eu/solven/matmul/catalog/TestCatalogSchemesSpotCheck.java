package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.Optional;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Tag;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.FieldAwareLookup;

/**
 * Spot-check every inner scheme the FieldAwareLookup picks at sizes
 * involved in the failing recombinations. If any catalog scheme itself
 * fails the spot-check, that's the root cause (we've been recombining
 * a broken inner scheme).
 */
@Tag("catalog-iterating")
public class TestCatalogSchemesSpotCheck {

	@Test
	public void check_inner_shapes_for_failing_recombines() {
		FieldAwareLookup lookup = new FieldAwareLookup("R");
		// Shapes that recombine touches for failing cases.
		int[][] shapes = {
				{8, 8, 8}, {3, 8, 8}, {4, 8, 8}, {5, 8, 8}, {6, 8, 8}, {7, 8, 8},
				{9, 9, 9}, {10, 10, 10}, {11, 11, 11},
				{9, 10, 10}, {9, 11, 11},
				// Cover all shapes touched by [9,10]³ and [9,11]³ Strassen recombines.
				{9, 9, 10}, {9, 10, 9}, {10, 9, 9}, {10, 10, 9}, {10, 9, 10}, {9, 10, 10},
				{9, 9, 11}, {9, 11, 9}, {11, 9, 9}, {11, 11, 9}, {11, 9, 11}, {9, 11, 11},
		};
		int failures = 0;
		for (int[] s : shapes) {
			Optional<NonCubicBilinearAlgorithm> opt = lookup.find(s[0], s[1], s[2]);
			if (opt.isEmpty()) {
				System.out.printf("⟨%d,%d,%d⟩  MISSING%n", s[0], s[1], s[2]);
				continue;
			}
			NonCubicBilinearAlgorithm alg = opt.get();
			boolean ok = Verifier.passesRandomMatmulSpotCheck(alg);
			System.out.printf("⟨%d,%d,%d⟩ found rank=%d   %s%n",
					s[0], s[1], s[2], alg.r, ok ? "PASS" : "FAIL");
			if (!ok) failures++;
		}
		assertThat(failures).as("catalog inner schemes that fail spot-check").isZero();
	}
}
