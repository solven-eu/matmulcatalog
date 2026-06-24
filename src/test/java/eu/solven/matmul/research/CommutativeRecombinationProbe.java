package eu.solven.matmul.research;

import java.io.File;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.recombination.Recombination;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.papers.waksman1970.WaksmanBound;

/**
 * Probe: for each ⟨n,n,n⟩ target the user cares about, compute the rank
 * a Strassen-recombination would yield if the inner sub-products used
 * the best-known <strong>commutative</strong> sub-ranks, and compare to
 * DIS09 Table 4 (commutative).
 *
 * <p>Outer layer is non-commutative Strassen (7 sub-products per level).
 * The inner sub-blocks are scalars in the outer-block ring, which is a
 * non-commutative matrix ring — so technically commutative inner
 * algorithms only apply when the OUTERMOST level uses scalars. Treat the
 * numbers below as a "fantasy bound": what rank we'd reach if the
 * Strassen recursion could host commutative leaves.</p>
 *
 * <p>This does NOT write a scheme — the result isn't representable in
 * our bilinear scheme format. Output is rank only.</p>
 */
public final class CommutativeRecombinationProbe {

	private CommutativeRecombinationProbe() {}

	/** Best known commutative rank for ⟨a,b,c⟩, used as the SOTA resolver. */
	private static int commutativeRank(int a, int b, int c) {
		if (a == 0 || b == 0 || c == 0) return 0;
		// Sort a ≤ b ≤ c canonical (commutative rank is permutation-invariant under tensor symmetries).
		int[] s = { a, b, c };
		java.util.Arrays.sort(s);
		int x = s[0], y = s[1], z = s[2];

		// Known explicit commutative results (publish year / source).
		// Cubic:
		if (x == y && y == z) {
			switch (x) {
				case 1: return 1;
				case 2: return 6;    // Hopcroft-Kerr 1971 / Winograd 1971
				case 3: return 21;   // Hopcroft-Kerr 1971 / Rosowski 2019 Corollary 1
				case 4: return 46;   // Waksman 1970
				case 5: return 93;   // Waksman 1970
				case 6: return 141;  // Waksman 1970
				case 7: return 235;  // Waksman 1970
				case 8: return 316;  // Waksman 1970
				case 9: return 472;  // DIS09 Table 4 (via Hopcroft 3×3=21 composition)
				case 10: return 595; // Waksman 1970
				case 11: return 825; // DIS09 Table 4
				case 12: return 987; // DIS09 Table 4
				case 13: return 1318; // DIS09
				case 14: return 1525; // DIS09
				case 15: return 1941; // DIS09
				default: return (int) WaksmanBound.forShape(x, x, x);
			}
		}
		// Non-cubic ⟨n,3,3⟩ Rosowski 2019 Algorithm 1: 6n+3
		if (y == 3 && z == 3) return 6 * x + 3;
		// Non-cubic ⟨2,2,n⟩: 3n+1 (cf. Hopcroft-Kerr) — well-known commutative
		if (x == 2 && y == 2) return 3 * z + 1;
		// General: Waksman closed form, best over axis perms.
		return (int) Math.min(
				Math.min(WaksmanBound.forShape(x, y, z), WaksmanBound.forShape(y, x, z)),
				WaksmanBound.forShape(x, z, y));
	}

	private record Probe(int n, int[] alloc, int dis09Comm) {}

	/** Split {@code n} into {@code k} near-balanced parts, e.g. (14,3) → [5,5,4]. */
	private static int[] splitInto(int n, int k) {
		int[] a = new int[k];
		int base = n / k, rem = n % k;
		for (int i = 0; i < k; i++) a[i] = base + (i < rem ? 1 : 0);
		return a;
	}

	public static void main(String[] args) throws Exception {
		// DIS09 Table 4 (commutative) cubic values, where known.
		Probe[] probes = {
				new Probe(14, new int[] { 7, 7 }, 1525),
				new Probe(17, new int[] { 8, 9 }, 2435), // ≈ approx; verify against table
				new Probe(19, new int[] { 9, 10 }, 3463),
				new Probe(21, new int[] { 9, 12 }, 4502),
				new Probe(23, new int[] { 11, 12 }, 5839),
		};

		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		NonCubicBilinearAlgorithm laderman = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section3/laderman_1976-3x3x3_m23_a98.json"));
		Recombination.SotaResolver commSota = CommutativeRecombinationProbe::commutativeRank;

		record OuterBase(String name, NonCubicBilinearAlgorithm alg) {}
		OuterBase[] outers = { new OuterBase("Strassen<2,2,2>=7", strassen),
				new OuterBase("Laderman<3,3,3>=23", laderman) };

		System.out.printf("%-3s  %-22s  %-12s  %-7s  %-8s  %-8s  %-8s%n",
				"n", "outer", "alloc", "rank", "DIS09cmt", "Δ", "verdict");
		for (Probe p : probes) {
			for (OuterBase ob : outers) {
				int dim = ob.alg.n;          // ⟨d,d,d⟩ cubic outer assumed below
				int[] alloc;
				if (dim == 2) alloc = p.alloc;       // [a, b]
				else if (dim == 3) alloc = splitInto(p.n, 3); // [a, b, c]
				else continue;
				Recombination.Result rec = Recombination.recombineWithAllocation(
						ob.alg, commSota, alloc, alloc, alloc);
				long rank = rec.totalRank;
				long delta = rank - p.dis09Comm;
				String verdict = rank < p.dis09Comm ? "WIN" : (rank == p.dis09Comm ? "TIE" : "LOSE");
				System.out.printf("%-3d  %-22s  %-12s  %-7d  %-8d  %+8d  %-8s%n",
						p.n, ob.name, java.util.Arrays.toString(alloc),
						rank, p.dis09Comm, delta, verdict);
			}
		}
		System.out.println();
		System.out.println("Note: ranks above are 'fantasy bounds' — they assume the Strassen");
		System.out.println("recursion can host commutative inner algorithms, which only works");
		System.out.println("at the outermost (scalar) level. Treat as an upper bound on what a");
		System.out.println("dedicated commutative composition framework could yield.");
	}
}
