package eu.solven.matmul.docs.verify;

import eu.solven.matmul.recombination.Recombination;

import java.nio.file.Path;
import java.util.List;
import java.util.Optional;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.recombination.Recombination.SotaResolver;
import eu.solven.matmul.recombination.BlockSplitSearch;
import eu.solven.matmul.recombination.BlockSplitSearch.NamedBase;
import eu.solven.matmul.recombination.BlockSplitSearch.NonCubicStrategy;
import eu.solven.matmul.search.PoolConfig;
import eu.solven.matmul.search.RecursiveClosureSota;
import eu.solven.matmul.search.RecursiveMaterialiser;

/**
 * Verify the 5 "composition beats catalog" candidates surfaced by the
 * 2026-06-04 band 11-16 evaluate sweep. For each shape:
 * <ol>
 *   <li>reproduce the evaluate composition rank via {@code findBestStrategy}
 *       (catalog SOTA, unbalanced) — should match the sweep's number;</li>
 *   <li>materialise the ACTUAL recombined scheme (unbalanced, no disk write) and
 *       run {@link Verifier#isExactNonCubic} — proves it genuinely computes the
 *       matmul;</li>
 *   <li>compare the verified scheme's rank to the on-disk catalog entry.</li>
 * </ol>
 * A candidate is a real win iff: Verifier PASSES and verified-rank &lt; catalog.
 */
public final class VerifyBand11to16Candidates {

	private VerifyBand11to16Candidates() {}

	private static SotaResolver catalogSota(FieldAwareLookup lk) {
		return (p, q, r) -> {
			if (p == 0 || q == 0 || r == 0) return 0;
			if (p == 1) return q * r;
			if (q == 1) return p * r;
			if (r == 1) return p * q;
			int v = lk.findRank(p, q, r);
			return v >= Recombination.SotaResolver.UNKNOWN_RANK ? p * q * r : v;
		};
	}

	public static void main(String[] args) throws Exception {
		int[][] cands = {
				{ 2, 10, 15 }, { 2, 10, 16 }, { 2, 12, 16 }, { 3, 3, 14 }, { 3, 3, 15 } };

		FieldAwareLookup lk = new FieldAwareLookup("R");
		List<NamedBase> pool = BlockSplitSearch.buildPool(PoolConfig.simple());
		SotaResolver sota = catalogSota(lk);
		RecursiveClosureSota recSota = new RecursiveClosureSota(lk, pool, true, true);
		Path tmp = Path.of("target/verify-tmp");
		// writeNewSchemes=false (no catalog mutation), balancedOnly=false (match the
		// unbalanced evaluate path that produced the wins).
		RecursiveMaterialiser mat = new RecursiveMaterialiser(lk, pool, recSota, tmp, false, false);

		System.out.printf("%-12s  %-9s  %-8s  %-9s  %-9s  %-8s  %s%n",
				"shape", "evalComp", "catalog", "matRank", "verified", "fromDisk", "verdict");
		System.out.println("-".repeat(86));
		for (int[] c : cands) {
			int n = c[0], m = c[1], p = c[2];
			String shape = String.format("⟨%d,%d,%d⟩", n, m, p);

			Optional<NonCubicStrategy> strat =
					BlockSplitSearch.findBestStrategy(n, m, p, pool, sota, false);
			long evalComp = strat.map(NonCubicStrategy::rank).orElse(-1L);
			String stratLabel = strat.map(NonCubicStrategy::label).orElse("none");
			System.out.println("  " + shape + " best strategy: " + stratLabel + " = " + evalComp);

			int catalog = lk.findRank(n, m, p);
			boolean haveCat = catalog < Recombination.SotaResolver.UNKNOWN_RANK;

			long matRank = -1; boolean verified = false, fromDisk = false; String note = "";
			try {
				Optional<RecursiveMaterialiser.Result> r = mat.materialise(n, m, p);
				if (r.isPresent()) {
					NonCubicBilinearAlgorithm alg = r.get().alg();
					matRank = alg.r;
					fromDisk = r.get().fromDisk();
					verified = Verifier.isExactNonCubic(alg);
				} else {
					note = "materialise-empty";
				}
			} catch (Exception e) {
				note = "EXC:" + e.getClass().getSimpleName();
			}

			String verdict;
			if (!verified) verdict = "✗ NOT VERIFIED " + note;
			else if (!haveCat) verdict = "verified, no catalog entry";
			else if (matRank < catalog) verdict = "✔ REAL WIN (−" + (catalog - matRank) + ")";
			else if (matRank == catalog) verdict = "tie (eval win not reproduced by materialiser)";
			else verdict = "materialiser worse than catalog";

			System.out.printf("%-12s  %-9d  %-8s  %-9d  %-9s  %-8s  %s%n",
					shape, evalComp, haveCat ? String.valueOf(catalog) : "none",
					matRank, verified, fromDisk, verdict);
		}
	}
}
