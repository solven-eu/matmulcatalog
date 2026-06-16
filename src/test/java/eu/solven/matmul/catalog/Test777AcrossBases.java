package eu.solven.matmul.catalog;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Optional;
import java.util.stream.Stream;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.Recombination;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.catalog.Recombination.SotaResolver;

/**
 * Diagnostic — for the {@code ⟨7,7,7⟩} target, evaluate the formula-derived
 * rank using OUTER bases of various sizes (⟨2,2,2⟩, ⟨3,3,3⟩, ⟨4,4,4⟩,
 * ⟨5,5,5⟩, …). For each base, enumerate all non-degenerate allocations
 * summing to 7 and pick the min. Prints a table for inclusion in
 * TRILINEAR_AGGREGATION.md.
 */
public class Test777AcrossBases {

	private static Recombination.AlgorithmLookup catalogLookup() {
		return (n, m, p) -> {
			int[] sorted = { n, m, p };
			java.util.Arrays.sort(sorted);
			String prefix = sorted[0] + "x" + sorted[1] + "x" + sorted[2];
			Path root = Path.of("src/main/resources/schemes");
			try (Stream<Path> s = Files.walk(root)) {
				return s.filter(p_ -> {
					String name = p_.getFileName().toString();
					return name.endsWith(".json")
							&& name.matches(".*[_-]" + prefix + "_[rm].*")
							&& !name.contains("F2") && !name.contains("Z2");
				}).findFirst().flatMap(p_ -> {
					try {
						NonCubicBilinearAlgorithm alg = SchemeIO.readBilinear(p_.toFile());
						return alg.orientAs(n, m, p);
					} catch (Exception e) {
						return Optional.empty();
					}
				});
			} catch (IOException e) {
				return Optional.empty();
			}
		};
	}

	@Test
	public void enumerate_777_via_various_outer_bases() throws IOException {
		Recombination.AlgorithmLookup lookup = catalogLookup();
		// Forbid direct ⟨7,7,7⟩ catalog lookups so degenerate allocations don't trivially win.
		// Trivial fallbacks for ⟨1,n,m⟩-family (matrix-vector product, rank n·m).
		SotaResolver sota = (a, b, c) -> {
			if (a == 0 || b == 0 || c == 0) return 0;
			if (a == 7 && b == 7 && c == 7) return Integer.MAX_VALUE / 100;
			if (a == 1) return b * c;
			if (b == 1) return a * c;
			if (c == 1) return a * b;
			return lookup.find(a, b, c).map(alg -> alg.r).orElse(Integer.MAX_VALUE / 100);
		};

		int[] baseDims = { 2, 3, 4, 5, 6 };
		System.out.println();
		System.out.println("⟨7,7,7⟩ rank via OUTER base ⟨b,b,b⟩ + best non-degenerate allocation [a1..ab]³:");
		System.out.println();
		System.out.printf("%6s | %5s | %35s | %20s%n",
				"base", "rank", "best allocation [a1..ab]³", "best rank for ⟨7,7,7⟩");
		System.out.println("-".repeat(75));
		// The choice of WHICH rank-7 ⟨2,2,2⟩ algorithm matters. The canonical Strassen
		// file (sign/permutation conventions matching Strassen 1969) yields rank 250 at
		// [4,3]³; the alphatensor-Z variant (different equivalent algorithm) yields 255.
		// Same rank, different U/V/W layouts → different sub-shape distribution under
		// min-reduction. So always pin the specific file when measuring.
		NonCubicBilinearAlgorithm strassenCanonical = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		for (int b : baseDims) {
			Optional<NonCubicBilinearAlgorithm> baseOpt = (b == 2)
					? Optional.of(strassenCanonical)
					: lookup.find(b, b, b);
			if (baseOpt.isEmpty()) {
				System.out.printf("%6d | %5s | %35s | %20s%n", b, "—", "—", "(missing in catalog)");
				continue;
			}
			NonCubicBilinearAlgorithm base = baseOpt.get();
			List<int[]> allocs = Recombination.blockFillings(b, 7);
			int bestRank = Integer.MAX_VALUE;
			int[] bestAlloc = null;
			for (int[] alloc : allocs) {
				// Non-degenerate: every block must be > 0.
				boolean ok = true;
				for (int x : alloc) if (x == 0) { ok = false; break; }
				if (!ok) continue;
				Recombination.Result r = Recombination.recombineWithAllocation(base, sota, alloc, alloc, alloc);
				if (r.totalRank < bestRank) {
					bestRank = (int) r.totalRank;
					bestAlloc = alloc;
				}
			}
			System.out.printf("%6s | %5d | %35s | %20s%n",
					"⟨" + b + "," + b + "," + b + "⟩", base.r,
					bestAlloc == null ? "(none)" : java.util.Arrays.toString(bestAlloc) + "³",
					bestRank == Integer.MAX_VALUE ? "—" : bestRank);
		}
		System.out.println();
		System.out.println("Note: catalog's direct ⟨7,7,7⟩=250 (Perminov/Sedoglavic) — listed as");
		System.out.println("    reference. ⟨2,2,2⟩ Strassen with [4,3]³ matches this exactly.");
	}
}
