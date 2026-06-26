package eu.solven.matmul.docs.verify;

import java.nio.file.Path;
import java.util.List;
import java.util.Optional;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.recombination.BlockSplitSearch;
import eu.solven.matmul.recombination.BlockSplitSearch.NamedBase;
import eu.solven.matmul.recombination.Recombination.SotaResolver;
import eu.solven.matmul.search.RecombinationPoolConfig;
import eu.solven.matmul.search.RecursiveClosureSota;
import eu.solven.matmul.search.RecursiveMaterialiser;

/**
 * Build-verify arbitrary shapes (no disk write). Used to resolve the {@code find()}-blind >16
 * sub-blocks of the ⟨18,30,31⟩=8970 candidate: if each materialises + {@link Verifier} passes at a
 * rank ≤ its {@code findRank}, the catalog sub-rank is a real (buildable) bound, not a phantom stub.
 */
public final class VerifySubBlocks {
	private VerifySubBlocks() {}

	private static SotaResolver catalogSota(FieldAwareLookup lk) {
		return (p, q, r) -> {
			if (p == 0 || q == 0 || r == 0) return 0;
			if (p == 1) return q * r;
			if (q == 1) return p * r;
			if (r == 1) return p * q;
			int v = lk.findRank(p, q, r);
			return v >= SotaResolver.UNKNOWN_RANK ? p * q * r : v;
		};
	}

	public static void main(String[] args) throws Exception {
		FieldAwareLookup lk = new FieldAwareLookup("R");
		List<NamedBase> pool = BlockSplitSearch.buildPool(RecombinationPoolConfig.thorough());
		SotaResolver sota = catalogSota(lk);
		RecursiveClosureSota recSota = new RecursiveClosureSota(lk, pool, false, true);
		RecursiveMaterialiser mat = new RecursiveMaterialiser(lk, pool, recSota, Path.of("target/verify-tmp"), false, false);

		int[][] shapes = args.length >= 3
				? new int[][] { { Integer.parseInt(args[0]), Integer.parseInt(args[1]), Integer.parseInt(args[2]) } }
				: new int[][] { { 9, 18, 18 }, { 9, 18, 13 }, { 9, 12, 18 }, { 9, 12, 13 } };
		System.out.printf("%-12s %9s %9s %9s  %s%n", "shape", "findRank", "matRank", "verified", "verdict");
		System.out.println("-".repeat(70));
		for (int[] s : shapes) {
			int n = s[0], m = s[1], p = s[2];
			int fr = lk.findRank(n, m, p);
			long matRank = -1; boolean verified = false; String note = "";
			try {
				Optional<RecursiveMaterialiser.Result> r = mat.materialise(n, m, p);
				if (r.isPresent()) { matRank = r.get().alg().r; verified = Verifier.isExactNonCubic(r.get().alg()); }
				else note = "empty";
			} catch (Throwable t) { note = "EXC:" + t.getClass().getSimpleName(); }
			String verdict = !verified ? "✗ NOT BUILDABLE " + note
					: matRank <= fr ? "✔ buildable at " + matRank + (matRank < fr ? " (BETTER than findRank!)" : " (= findRank)")
					: "buildable but worse (" + matRank + " > " + fr + " → findRank was a PHANTOM)";
			System.out.printf("⟨%d,%d,%d⟩ %9d %9d %9s  %s%n", n, m, p, fr, matRank, verified, verdict);
		}
	}
}
