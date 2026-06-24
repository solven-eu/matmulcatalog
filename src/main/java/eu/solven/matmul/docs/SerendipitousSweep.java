package eu.solven.matmul.docs;

import eu.solven.matmul.recombination.Recombination;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.catalog.SerendipitousBudProduct;

/**
 * Discovery sweep for serendipitous products (#159 closure step). Scans catalog
 * schemes for bud-rich BASES, multiplies each by small second schemes, and
 * reports any target whose serendipitous rank beats the current catalog SOTA —
 * with the result verified by {@link Verifier#isExactNonCubic}.
 *
 * <pre>
 *   mvn -q -o -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.SerendipitousSweep \
 *       -Dexec.args="--baseCap=8 --secondCap=6 --targetCap=16 --field=Q"
 * </pre>
 */
public final class SerendipitousSweep {

	private SerendipitousSweep() {}

	public static void main(String[] args) throws Exception {
		int baseCap = intArg(args, "--baseCap", 8);
		int secondCap = intArg(args, "--secondCap", 6);
		int targetCap = intArg(args, "--targetCap", 16);
		int targetMin = intArg(args, "--targetMin", 2);
		String field = strArg(args, "--field", "Q");
		FieldAwareLookup lk = new FieldAwareLookup(field);

		System.out.printf("[serendipitous] baseCap=%d secondCap=%d targetMin=%d targetCap=%d field=%s%n",
				baseCap, secondCap, targetMin, targetCap, field);

		// 1) Gather bud-rich bases from the catalog.
		List<NonCubicBilinearAlgorithm> bases = new ArrayList<>();
		List<String> baseLabels = new ArrayList<>();
		Path root = Path.of("src/main/resources/schemes");
		List<Path> files;
		try (var s = Files.walk(root)) {
			files = s.filter(p -> p.toString().endsWith(".json")).sorted().toList();
		}
		int scanned = 0, budRich = 0;
		for (Path f : files) {
			NonCubicBilinearAlgorithm a;
			try { a = SchemeIO.read(f.toFile()); } catch (Exception e) { continue; }
			scanned++;
			if (Math.max(a.n, Math.max(a.m, a.p)) > baseCap) continue;
			var dec = SerendipitousBudProduct.findBuds(a);
			if (dec.buds().isEmpty()) continue;
			bases.add(a);
			baseLabels.add(f.getFileName().toString().replace(".json", ""));
			budRich++;
		}
		System.out.printf("[serendipitous] scanned %d files, %d bud-rich bases (maxDim≤%d)%n",
				scanned, budRich, baseCap);

		// 1b) True per-shape SOTA — content-driven via findRank (stub-aware AND
		//     orientation-invariant: it keys on the sorted shape over fileIndex,
		//     which includes maxDim>16 stubs that findWithSource skips). This is
		//     what a serendipitous r_s must beat.
		//     WAS a filename regex (`(\d+)x(\d+)x(\d+)_m(\d+)`) that the 2026-06
		//     `-r{rank}-` rename silently zeroed → empty index → trueSota always
		//     -1 → every candidate skipped → the sweep could NEVER report a win.
		//     Guard: TestSerendipitousSweep.sota_oracle_is_content_driven.

		// 2) Sweep base × second-shape → target; predict, then build+verify wins.
		Map<String, String> wins = new LinkedHashMap<>(); // target → description (best per target)
		Map<String, Long> winRank = new LinkedHashMap<>();
		long t0 = System.currentTimeMillis();
		for (int bi = 0; bi < bases.size(); bi++) {
			NonCubicBilinearAlgorithm base = bases.get(bi);
			var dec = SerendipitousBudProduct.findBuds(base);
			for (int n2 = 1; n2 <= secondCap; n2++)
				for (int m2 = 1; m2 <= secondCap; m2++)
					for (int p2 = 1; p2 <= secondCap; p2++) {
						if (n2 == 1 && m2 == 1 && p2 == 1) continue;
						int tn = base.n * n2, tm = base.m * m2, tp = base.p * p2;
						int tMax = Math.max(tn, Math.max(tm, tp));
						if (tMax > targetCap || tMax < targetMin) continue;
						long pred = predict(dec, lk, n2, m2, p2);
						if (pred < 0) continue;
						// TRUE SOTA: catalog best (incl stubs, any orientation) AND plain
						// Kronecker over all factorisations. Only a genuine improvement
						// passes (fixes the ⟨18,18,18⟩=3200 false positive: 3200 = plain
						// Kron ⟨3,3,6⟩⊗⟨6,6,3⟩, so not a win).
						long sota = trueSota(lk, tn, tm, tp);
						if (sota < 0 || pred >= sota) continue;
						// candidate win → build + verify
						NonCubicBilinearAlgorithm built;
						try {
							built = SerendipitousBudProduct.productViaBuds(base, lk, n2, m2, p2);
						} catch (Exception e) { continue; }
						if (built.r != pred || !Verifier.isExactNonCubic(built)) continue;
						String key = tn + "x" + tm + "x" + tp;
						if (!winRank.containsKey(key) || built.r < winRank.get(key)) {
							winRank.put(key, (long) built.r);
							wins.put(key, String.format("⟨%s⟩ %d < sota %d  via %s ⊗ ⟨%d,%d,%d⟩",
									key, built.r, sota, baseLabels.get(bi), n2, m2, p2));
						}
					}
			if ((bi + 1) % 25 == 0 || bi + 1 == bases.size()) {
				System.out.printf("[progress] %d/%d bases, %d target-wins, %dms%n",
						bi + 1, bases.size(), wins.size(), System.currentTimeMillis() - t0);
			}
		}

		System.out.println("=== serendipitous wins (verified, r_s < catalog SOTA) ===");
		if (wins.isEmpty()) System.out.println("  (none below current SOTA in this range)");
		wins.values().forEach(w -> System.out.println("  " + w));
		System.out.printf("=== done: %d target-wins ===%n", wins.size());
	}

	/** Best known NC rank for {@code ⟨n,m,p⟩} (any orientation, incl stubs) via the
	 *  content-driven {@link FieldAwareLookup#findRank} index; -1 if unknown. */
	static long catalogRank(FieldAwareLookup lk, int n, int m, int p) {
		int r = lk.findRank(n, m, p);
		return r >= Recombination.SotaResolver.UNKNOWN_RANK ? -1L : r;
	}

	/** Best plain-Kronecker rank over ALL proper factorisations; -1 if none buildable. */
	static long kronBest(FieldAwareLookup lk, int n, int m, int p) {
		long best = -1;
		for (int n1 = 1; n1 <= n; n1++) {
			if (n % n1 != 0) continue;
			for (int m1 = 1; m1 <= m; m1++) {
				if (m % m1 != 0) continue;
				for (int p1 = 1; p1 <= p; p1++) {
					if (p % p1 != 0) continue;
					if (n1 * m1 * p1 == 1) continue;
					if ((n / n1) * (m / m1) * (p / p1) == 1) continue;
					long r1 = catalogRank(lk, n1, m1, p1);
					long r2 = catalogRank(lk, n / n1, m / m1, p / p1);
					if (r1 < 0 || r2 < 0) continue;
					long prod = r1 * r2;
					if (best < 0 || prod < best) best = prod;
				}
			}
		}
		return best;
	}

	/** The rank a serendipitous result must strictly beat: min(catalog, plain-Kron). */
	static long trueSota(FieldAwareLookup lk, int n, int m, int p) {
		long cat = catalogRank(lk, n, m, p);
		long kron = kronBest(lk, n, m, p);
		if (cat < 0) return kron;
		if (kron < 0) return cat;
		return Math.min(cat, kron);
	}

	private static long predict(SerendipitousBudProduct.BudDecomposition dec,
			FieldAwareLookup lk, int n2, int m2, int p2) {
		long triv = rank(lk, n2, m2, p2);
		if (triv < 0) return -1;
		long total = (long) dec.trivial().length * triv;
		for (var bud : dec.buds()) {
			int k = bud.terms().length;
			long r = switch (bud.type()) {
				case U -> rank(lk, n2, m2, k * p2);
				case V -> rank(lk, k * n2, m2, p2);
				case W -> rank(lk, n2, k * m2, p2);
			};
			if (r < 0) return -1;
			total += r;
		}
		return total;
	}

	private static long rank(FieldAwareLookup lk, int n, int m, int p) {
		return lk.findWithSource(n, m, p).map(ws -> (long) ws.alg().r).orElse(-1L);
	}

	private static int intArg(String[] a, String key, int def) {
		for (String s : a) if (s.startsWith(key + "=")) return Integer.parseInt(s.substring(key.length() + 1));
		return def;
	}

	private static String strArg(String[] a, String key, String def) {
		for (String s : a) if (s.startsWith(key + "=")) return s.substring(key.length() + 1);
		return def;
	}
}
