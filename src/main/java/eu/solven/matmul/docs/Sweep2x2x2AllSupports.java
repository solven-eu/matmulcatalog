package eu.solven.matmul.docs;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.Recombination;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.search.CitedBound;
import eu.solven.matmul.search.LineageReplayer;

/**
 * Sweep using <strong>only ⟨2,2,2⟩ bases — but every distinct support pattern</strong>
 * (Strassen, Winograd, Waksman, AlphaTensor-Z, Perminov, the solven orbit reps, …)
 * in every axis-orientation, recombined at every 2-part-per-axis allocation. Each
 * (support × allocation) realises a <em>different recombination multiset</em>, so
 * this is the constructable way to span the 2×2×2 base's multiset variety (the
 * change-of-basis orbit predicts more, but those aren't buildable without the
 * realising transform).
 *
 * <p>For each cubic target it predicts the rank of every (base, allocation) via
 * {@link Recombination#recombineWithAllocation}, constructs + verifies the best,
 * and compares to the catalog. Reports WIN (&lt; catalog) / TIE (= catalog) /
 * worse. {@code --apply} materialises ties-or-better as {@code derived/} stubs.
 * {@code --min/--max} bound the cubic band (default 2..16).</p>
 */
public final class Sweep2x2x2AllSupports {
	private Sweep2x2x2AllSupports() {}

	private static final String ROOT = "src/main/resources/schemes";

	public static void main(String[] args) throws Exception {
		boolean apply = List.of(args).contains("--apply");
		int minN = intArg(args, "--min", 2), maxN = intArg(args, "--max", 16);

		FieldAwareLookup lookup = new FieldAwareLookup(Field.Q);
		LineageReplayer replayer = LineageReplayer.withDefaultPool(lookup);
		CitedBound sota = new CitedBound(lookup);
		Recombination.AlgorithmLookup recombLookup = (a, b, c) -> {
			try { return Optional.of(replayer.replay(new Lineage.Atom(a + "x" + b + "x" + c))); }
			catch (RuntimeException e) { return Optional.empty(); }
		};

		// Every distinct-support ⟨2,2,2⟩=7 base, each in its 3 cyclic orientations,
		// deduped by content hash. (cyclicShift rotates the axes U→V→W.)
		List<NonCubicBilinearAlgorithm> bases = new ArrayList<>();
		Map<String, Boolean> seen = new LinkedHashMap<>();
		try (var s = Files.walk(Path.of(ROOT))) {
			for (Path p : s.filter(x -> x.getFileName().toString().matches("2x2x2-.*\\.json")).sorted().toList()) {
				NonCubicBilinearAlgorithm a;
				try { a = SchemeIO.read(p.toFile()); } catch (Exception e) { continue; }
				if (a.r != 7) continue;
				NonCubicBilinearAlgorithm cur = a;
				for (int rot = 0; rot < 3; rot++) {
					String h = SchemeIO.contentHash(cur);
					if (seen.putIfAbsent(h, true) == null) bases.add(cur);
					cur = cur.cyclicShift();
				}
			}
		}
		System.out.printf("⟨2,2,2⟩ base variants (distinct support × orientation): %d%n", bases.size());
		System.out.printf("cubic band %d..%d  %s%n%n", minN, maxN, apply ? "[APPLY]" : "[report]");
		System.out.println("shape      best2x2x2   catalog   outcome   (base / alloc)");
		System.out.println("---------------------------------------------------------------");

		int wins = 0, ties = 0;
		for (int N = minN; N <= maxN; N++) {
			long t0 = System.currentTimeMillis();
			long bestPred = Long.MAX_VALUE;
			NonCubicBilinearAlgorithm bestBase = null;
			int[] bA = null, bB = null, bC = null;
			for (NonCubicBilinearAlgorithm base : bases) {
				for (int[] aA : comps(N))
					for (int[] aB : comps(N))
						for (int[] aC : comps(N)) {
							Recombination.Result pr = Recombination.recombineWithAllocation(base, sota, aA, aB, aC);
							if (pr != null && pr.totalRank > 0 && pr.totalRank < bestPred) {
								bestPred = pr.totalRank; bestBase = base; bA = aA; bB = aB; bC = aC;
							}
						}
			}
			long cat = lookup.findRank(N, N, N);
			if (bestBase == null) continue;
			// construct + verify the best candidate
			String outcome; boolean verified = false; int built = -1;
			try {
				NonCubicBilinearAlgorithm alg = Recombination.constructWithAllocation(bestBase, recombLookup, bA, bB, bC);
				built = alg.r;
				verified = alg.n == N && alg.m == N && alg.p == N && alg.r == bestPred
						&& Verifier.passesRandomMatmulSpotCheck(alg);
			} catch (RuntimeException ignored) { /* not realisable */ }
			if (!verified) outcome = "pred-only(" + bestPred + ")";
			else if (built < cat) { outcome = "*** WIN"; wins++; }
			else if (built == cat) { outcome = "= TIE"; ties++; }
			else outcome = "worse";
			long ms = System.currentTimeMillis() - t0;
			System.out.printf("⟨%d,%d,%d⟩  %-9d  %-8d  %-9s base=%s n=%s m=%s p=%s  (%dms)%n",
					N, N, N, bestPred, cat, outcome,
					bestBase.n + "x" + bestBase.m + "x" + bestBase.p,
					java.util.Arrays.toString(bA), java.util.Arrays.toString(bB),
					java.util.Arrays.toString(bC), ms);
		}
		System.out.printf("%n%d wins, %d ties over cubic %d..%d (vs current catalog)%n", wins, ties, minN, maxN);
	}

	/** All ordered 2-part compositions of {@code d} ({@code [a, d-a]}, 1≤a≤d-1). */
	private static List<int[]> comps(int d) {
		List<int[]> out = new ArrayList<>();
		for (int a = 1; a < d; a++) out.add(new int[] { a, d - a });
		return out;
	}

	private static int intArg(String[] args, String key, int def) {
		for (String a : args) if (a.startsWith(key + "=")) return Integer.parseInt(a.substring(key.length() + 1));
		return def;
	}
}
