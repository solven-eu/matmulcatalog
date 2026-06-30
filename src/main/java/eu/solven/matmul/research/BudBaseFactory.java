package eu.solven.matmul.research;

import eu.solven.matmul.recombination.Recombination;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.Compose;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.catalog.SerendipitousBudProduct;
import lombok.extern.slf4j.Slf4j;

/**
 * Bud-base factory (v1) — explores whether bud-rich schemes give serendipitous
 * rank wins, using the cost-only primitive
 * {@link SerendipitousBudProduct#serendipitousCost}.
 *
 * <p>Three steps:</p>
 * <ol>
 *   <li><b>Bud scan</b>: read every small catalog scheme, score its bud
 *       structure, and report the bud-richest per shape (answers "is there a
 *       ⟨3,3,3⟩ with non-trivial buds?").</li>
 *   <li><b>Serendipitous evaluation</b>: for each bud-rich base and a grid of
 *       inner shapes, compute the predicted serendipitous rank
 *       {@code base ⊗ˢ inner} and compare it to the catalog's best rank at the
 *       target ⟨n₁n₂,m₁m₂,p₁p₂⟩. Report wins (serendipitous &lt; catalog).</li>
 *   <li><b>Factory</b>: manufacture larger bud-rich bases by Kronecker-ing the
 *       top seeds (Kron multiplies bud sizes — exact per LineageBudInference),
 *       and evaluate those as bases too.</li>
 * </ol>
 *
 * <p>The objective is the <em>actual</em> serendipitous cost vs catalog best —
 * no proxy ratio. A win means a bud-rich (possibly non-rank-minimal) base
 * produced a target rank below anything currently in the catalog.</p>
 */
@Slf4j
public final class BudBaseFactory {

	private BudBaseFactory() {}

	private static final Pattern NAME = Pattern.compile(
			"(?<source>.*)[_-](?<n>\\d+)x(?<m>\\d+)x(?<p>\\d+)_(?:r|m)(?<rank>\\d+)[^/]*\\.json");

	/** Σ over independent U/V/W class sizes ≥ 2 — total terms participating in a bud. */
	static int budScore(NonCubicBilinearAlgorithm a) {
		int[][] cs = SerendipitousBudProduct.independentClassSizes(a);
		int s = 0;
		for (int[] classes : cs) {
			for (int x : classes) {
				if (x >= 2) {
					s += x;
				}
			}
		}
		return s;
	}

	private static final int MAX_BASE_DIM = 8;
	private static final int MAX_MANUFACTURED_DIM = 12;
	private static final int[][] INNERS = {
			{ 2, 2, 2 }, { 3, 3, 3 }, { 2, 2, 3 }, { 2, 3, 3 }, { 2, 2, 4 },
			{ 2, 3, 4 }, { 3, 4, 3 }, { 2, 4, 5 }, { 2, 3, 6 }, { 4, 4, 4 }
	};

	public static void main(String[] args) throws Exception {
		FieldAwareLookup lookup = new FieldAwareLookup(Field.Q);

		// ── step 1: bud scan over small catalog schemes ──
		List<Path> files;
		try (var st = Files.walk(new File("src/main/resources/schemes").toPath())) {
			files = st.filter(p -> p.toString().endsWith(".json")).sorted().toList();
		}
		Map<String, NonCubicBilinearAlgorithm> richest = new HashMap<>();
		Map<String, Integer> richestScore = new HashMap<>();
		for (Path pth : files) {
			Matcher m = NAME.matcher(pth.getFileName().toString());
			if (!m.matches()) {
				continue;
			}
			int n = Integer.parseInt(m.group("n"));
			int mm = Integer.parseInt(m.group("m"));
			int p = Integer.parseInt(m.group("p"));
			if (Math.max(n, Math.max(mm, p)) > MAX_BASE_DIM) {
				continue;
			}
			NonCubicBilinearAlgorithm a;
			try {
				a = SchemeIO.read(pth.toFile());
			} catch (Exception ignored) {
				continue;
			}
			int sc = budScore(a);
			if (sc == 0) {
				continue;
			}
			String key = n + "x" + mm + "x" + p;
			if (sc > richestScore.getOrDefault(key, 0)) {
				richestScore.put(key, sc);
				richest.put(key, a);
			}
		}
		log.info("=== bud-rich catalog bases (maxDim ≤ {}): {} shapes ===", MAX_BASE_DIM, richest.size());
		richestScore.entrySet().stream()
				.sorted((x, y) -> y.getValue() - x.getValue())
				.limit(30)
				.forEach(e -> log.info("  {}  bud_score={}  rank={}",
						e.getKey(), e.getValue(), richest.get(e.getKey()).r));

		// Specific probe: every ⟨3,3,3⟩ scheme's bud score.
		log.info("=== ⟨3,3,3⟩ bud probe ===");
		for (Path pth : files) {
			if (!pth.getFileName().toString().contains("3x3x3")) {
				continue;
			}
			try {
				NonCubicBilinearAlgorithm a = SchemeIO.read(pth.toFile());
				if (a.n == 3 && a.m == 3 && a.p == 3) {
					log.info("  r={} bud_score={}  {}", a.r, budScore(a), pth.getFileName());
				}
			} catch (Exception ignored) {
				// skip
			}
		}

		// ── step 3: manufacture Kron candidates from the richest seeds ──
		List<NonCubicBilinearAlgorithm> candidates = new ArrayList<>(richest.values());
		List<NonCubicBilinearAlgorithm> seeds = richestScore.entrySet().stream()
				.sorted((x, y) -> y.getValue() - x.getValue())
				.limit(8)
				.map(e -> richest.get(e.getKey()))
				.toList();
		for (NonCubicBilinearAlgorithm s1 : seeds) {
			for (NonCubicBilinearAlgorithm s2 : seeds) {
				int N = s1.n * s2.n, M = s1.m * s2.m, P = s1.p * s2.p;
				if (Math.max(N, Math.max(M, P)) > MAX_MANUFACTURED_DIM) {
					continue;
				}
				try {
					candidates.add(Compose.kroneckerGeneral(s1, s2));
				} catch (Exception ignored) {
					// skip
				}
			}
		}
		log.info("=== serendipitous evaluation: {} bases × {} inners ===", candidates.size(), INNERS.length);

		// ── step 2: serendipitous cost vs catalog best ──
		record Win(int N, int M, int P, long cost, long best, NonCubicBilinearAlgorithm base,
				int n2, int m2, int p2) {}
		List<Win> winList = new ArrayList<>();
		int evals = 0;
		for (NonCubicBilinearAlgorithm base : candidates) {
			for (int[] inr : INNERS) {
				long cost = SerendipitousBudProduct.serendipitousCost(base, lookup, inr[0], inr[1], inr[2]);
				if (cost >= Long.MAX_VALUE / 4) {
					continue;
				}
				int N = base.n * inr[0], M = base.m * inr[1], P = base.p * inr[2];
				long best = lookup.findRank(N, M, P);
				if (best >= Recombination.SotaResolver.UNKNOWN_RANK) {
					continue;  // no catalog entry to beat → can't claim a win
				}
				evals++;
				if (cost < best) {
					winList.add(new Win(N, M, P, cost, best, base, inr[0], inr[1], inr[2]));
				}
			}
		}
		winList.sort((a, b) -> Long.compare((long) a.N * a.M * a.P, (long) b.N * b.M * b.P));
		log.info("=== {} evaluations, {} PREDICTED serendipitous wins over catalog ===",
				evals, winList.size());
		for (Win w : winList) {
			log.info("  pred ⟨{},{},{}⟩ {} < {}  (base ⟨{},{},{}⟩ r{}, inner ⟨{},{},{}⟩)",
					w.N, w.M, w.P, w.cost, w.best, w.base.n, w.base.m, w.base.p, w.base.r, w.n2, w.m2, w.p2);
		}

		// ── verification: actually BUILD + spot-check the smallest wins, confirming
		// the predicted rank is realized by a valid scheme (optimality discipline:
		// a prediction is a bound until materialised + verified). ──
		log.info("=== verifying smallest wins (build + random matmul spot-check) ===");
		int verified = 0, checked = 0;
		for (Win w : winList) {
			if (checked >= 8) {
				break;
			}
			checked++;
			try {
				NonCubicBilinearAlgorithm built =
						SerendipitousBudProduct.productViaBuds(w.base, lookup, w.n2, w.m2, w.p2);
				boolean ok = eu.solven.matmul.verifiers.Verifier.passesRandomMatmulSpotCheck(built);
				boolean rankOk = built.r == w.cost;
				if (ok && rankOk) {
					verified++;
					log.info("  CONFIRMED ⟨{},{},{}⟩ = {} (was {}) — built+verified",
							w.N, w.M, w.P, built.r, w.best);
				} else {
					log.info("  ⚠ ⟨{},{},{}⟩ built r={} (pred {}) verify={} — NOT confirmed",
							w.N, w.M, w.P, built.r, w.cost, ok);
				}
			} catch (Exception e) {
				log.info("  — ⟨{},{},{}⟩ unverifiable ({}: enlarged inner likely a stub)",
						w.N, w.M, w.P, e.getClass().getSimpleName());
			}
		}
		log.info("=== done: {} predicted wins, {}/{} smallest verified ===",
				winList.size(), verified, checked);
	}
}
