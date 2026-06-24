package eu.solven.matmul.research;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.CatalogLimits;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.recombination.BlockSplitSearch;
import eu.solven.matmul.search.RecursiveClosureSota;
import eu.solven.matmul.search.RecursiveMaterialiser;

/**
 * Multi-field recursive materialisation sweep with field-inclusion-aware
 * pruning. Processes one shape at a time across all fields {C, R, Q, Z, F2},
 * skipping narrower fields when the broadest field's materialised rank
 * already saturates the narrower catalog.
 *
 * <p><strong>Pruning rule.</strong> For a shape, let {@code m_C} be the
 * rank materialised under field C (broadest characteristic-0 field —
 * accepts every Z, Q, R, C scheme as a leaf). For any narrower field
 * {@code F ∈ {R, Q, Z}}: the materialiser at {@code F} can use only a
 * subset of leaves, so {@code m_F ≥ m_C}. If the catalog already has an
 * F-valid scheme at rank {@code c_F ≤ m_C}, then {@code m_F ≥ m_C ≥ c_F}
 * — no improvement possible — and we skip F entirely.</p>
 *
 * <p>F2 is independent (characteristic 2) and runs as its own pass.</p>
 *
 * <p>Compared to running per-field sweeps sequentially, this approach
 * (a) walks the shape grid once instead of once per field;
 * (b) shares the recursive-closure SOTA cache;
 * (c) typically skips 70-90% of narrow-field work because the catalog's
 * Z entries already cover the materialisable range.</p>
 */
// TODO: migrate ad-hoc progress reporting to eu.solven.matmul.util.ProgressMonitor.
public final class MaterialiseMultiFieldSweep {

	public static void main(String[] args) throws Exception {
		// Order: broadest characteristic-0 → narrowest, then F2 separately.
		Field[] charZeroChain = { Field.C, Field.R, Field.Q, Field.Z };

		// Per-field setup. Each FieldAwareLookup wraps the same disk root but
		// applies its own field-chain filter on candidate files.
		FieldAwareLookup[] lookups = new FieldAwareLookup[charZeroChain.length];
		RecursiveMaterialiser[] mats = new RecursiveMaterialiser[charZeroChain.length];
		Path writeRoot = Path.of("src/main/resources/schemes");
		for (int i = 0; i < charZeroChain.length; i++) {
			lookups[i] = new FieldAwareLookup(charZeroChain[i]);
			List<BlockSplitSearch.NamedBase> pool = BlockSplitSearch.defaultPool();
			RecursiveClosureSota sota = new RecursiveClosureSota(lookups[i], pool, true, true);
			mats[i] = new RecursiveMaterialiser(lookups[i], pool, sota, writeRoot,
					/*writeNew=*/ true, /*balancedOnly=*/ true);
		}
		FieldAwareLookup f2Lookup = new FieldAwareLookup(Field.F2);
		List<BlockSplitSearch.NamedBase> f2Pool = BlockSplitSearch.defaultPool();
		RecursiveClosureSota f2Sota = new RecursiveClosureSota(f2Lookup, f2Pool, true, true);
		RecursiveMaterialiser f2Mat = new RecursiveMaterialiser(
				f2Lookup, f2Pool, f2Sota, writeRoot, true, true);

		int MAX_CUBIC = Math.min(32, CatalogLimits.MAX_DIM);
		int MAX_NC = Math.min(32, CatalogLimits.MAX_DIM);

		long t0 = System.nanoTime();
		long[] tallies = new long[3 * charZeroChain.length + 3];
		// tallies layout per chain-field: [wins, skipped (catalog already saturates), evaluated]
		// then [f2_wins, f2_skipped, f2_evaluated], then total shapes seen.
		int totalIdx = tallies.length - 1;

		List<int[]> shapes = new ArrayList<>();
		for (int n = 4; n <= MAX_CUBIC; n++) shapes.add(new int[] { n, n, n });
		for (int n = 3; n <= MAX_NC; n++) {
			for (int m = n; m <= MAX_NC; m++) {
				for (int p = m; p <= MAX_NC; p++) {
					if (n == m && m == p) continue;
					shapes.add(new int[] { n, m, p });
				}
			}
		}
		System.out.printf("Processing %d shapes across %d char-0 fields + F2%n",
				shapes.size(), charZeroChain.length);

		for (int[] s : shapes) {
			int n = s[0], m = s[1], p = s[2];

			// Char-0 chain: materialise at broadest (C) first, then prune.
			Integer mC = null;  // m_C — materialised rank at C
			for (int i = 0; i < charZeroChain.length; i++) {
				Field f = charZeroChain[i];
				int catalogF = lookups[i].find(n, m, p).map(a -> a.r).orElse(Integer.MAX_VALUE);
				if (mC != null && catalogF <= mC) {
					// Catalog already saturates this narrower field; skip.
					tallies[i * 3 + 1]++;
					continue;
				}
				tallies[i * 3 + 2]++;
				Optional<RecursiveMaterialiser.Result> r = mats[i].materialise(n, m, p);
				if (r.isPresent() && !r.get().fromDisk()) {
					int matRank = r.get().alg().r;
					if (matRank < catalogF) {
						tallies[i * 3]++;
						if (i == 0) mC = matRank;
					} else if (i == 0) {
						mC = matRank;
					}
				} else if (r.isPresent() && i == 0) {
					mC = r.get().alg().r;
				}
			}

			// F2: independent.
			int catalogF2 = f2Lookup.find(n, m, p).map(a -> a.r).orElse(Integer.MAX_VALUE);
			tallies[3 * charZeroChain.length + 2]++;
			Optional<RecursiveMaterialiser.Result> rF2 = f2Mat.materialise(n, m, p);
			if (rF2.isPresent() && !rF2.get().fromDisk()
					&& rF2.get().alg().r < catalogF2) {
				tallies[3 * charZeroChain.length]++;
			}

			tallies[totalIdx]++;
			if (tallies[totalIdx] % 100 == 0) {
				long elapsedMs = (System.nanoTime() - t0) / 1_000_000L;
				StringBuilder sb = new StringBuilder();
				sb.append(String.format("  [progress] %4d/%d  %dms", tallies[totalIdx], shapes.size(), elapsedMs));
				for (int i = 0; i < charZeroChain.length; i++) {
					sb.append(String.format("  %s:[%dw/%dskip/%deval]",
							charZeroChain[i], tallies[i*3], tallies[i*3+1], tallies[i*3+2]));
				}
				sb.append(String.format("  F2:[%dw/%deval]",
						tallies[3 * charZeroChain.length], tallies[3 * charZeroChain.length + 2]));
				System.out.println(sb);
			}
		}

		System.out.println("-".repeat(80));
		System.out.printf("Total: %d shapes%n", tallies[totalIdx]);
		for (int i = 0; i < charZeroChain.length; i++) {
			System.out.printf("  %s: %d wins, %d skipped (catalog saturates), %d evaluated%n",
					charZeroChain[i], tallies[i*3], tallies[i*3+1], tallies[i*3+2]);
		}
		System.out.printf("  F2: %d wins, %d evaluated%n",
				tallies[3 * charZeroChain.length], tallies[3 * charZeroChain.length + 2]);
	}
}
