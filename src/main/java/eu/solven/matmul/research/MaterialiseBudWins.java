package eu.solven.matmul.research;

import java.io.File;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.catalog.SerendipitousBudProduct;
import eu.solven.matmul.search.LineageReplayer;
import lombok.extern.slf4j.Slf4j;

/**
 * Materialise the 6 verified bud-base-factory serendipitous wins (2026-06-06)
 * to disk JSON. Each is {@code base ⟨n1,m1,p1⟩ ⊗ˢ inner ⟨n2,m2,p2⟩}. For each:
 * pick the bud-richest base at the base shape, build via
 * {@link SerendipitousBudProduct#productViaBuds}, spot-check, then write a
 * replayable {@code SerendipitousProduct} stub — but only after confirming the
 * stub's lineage replays to the same rank (else fall back to full matrices, so
 * the discovery is never lost to a fragile shape-ref resolution).
 */
@Slf4j
public final class MaterialiseBudWins {
	private MaterialiseBudWins() {}

	private record Win(int n1, int m1, int p1, int n2, int m2, int p2, int expectedRank) {}

	private static final Win[] WINS = {
			new Win(2, 3, 7, 3, 3, 3, 709),   // ⟨6,9,21⟩
			new Win(2, 4, 8, 3, 3, 3, 1029),  // ⟨6,12,24⟩
			new Win(2, 4, 8, 3, 4, 3, 1383),  // ⟨6,16,24⟩
			new Win(2, 5, 5, 3, 4, 3, 1136),  // ⟨6,20,15⟩
			new Win(2, 6, 6, 2, 3, 4, 1119),  // ⟨4,18,24⟩
			new Win(2, 7, 7, 2, 3, 4, 1518),  // ⟨4,21,28⟩
	};

	public static void main(String[] args) throws Exception {
		FieldAwareLookup lookup = new FieldAwareLookup(Field.Q);
		LineageReplayer replayer = LineageReplayer.withDefaultPool(lookup);
		int written = 0;
		for (Win w : WINS) {
			NonCubicBilinearAlgorithm base = richest(lookup, w.n1, w.m1, w.p1);
			if (base == null) {
				log.warn("no base at ⟨{},{},{}⟩ — skip", w.n1, w.m1, w.p1);
				continue;
			}
			NonCubicBilinearAlgorithm built =
					SerendipitousBudProduct.productViaBuds(base, lookup, w.n2, w.m2, w.p2);
			int N = w.n1 * w.n2, M = w.m1 * w.m2, P = w.p1 * w.p2;
			if (built.r != w.expectedRank || !Verifier.passesRandomMatmulSpotCheck(built)) {
				log.warn("⟨{},{},{}⟩ built r={} (exp {}) verify={} — SKIP",
						N, M, P, built.r, w.expectedRank, Verifier.passesRandomMatmulSpotCheck(built));
				continue;
			}
			int adds = Verifier.additionCount(built);
			// Reference the base by HASH-ref "{shape}@{contentHash}" so replay
			// resolves the EXACT bud-rich base, not whatever the shape-ref would
			// pick (a bud-poor sibling) — the replayability fix.
			String baseRef = w.n1 + "x" + w.m1 + "x" + w.p1 + "@" + SchemeIO.contentHash(base);
			Lineage.Node lineage = new Lineage.SerendipitousProduct(
					new Lineage.Atom(baseRef), w.n2, w.m2, w.p2);

			// Does the stub's lineage replay to the same scheme? If so, a stub is safe.
			boolean replayOk;
			try {
				NonCubicBilinearAlgorithm rep = replayer.replay(lineage);
				replayOk = rep.r == built.r && Verifier.passesRandomMatmulSpotCheck(rep);
			} catch (Exception e) {
				replayOk = false;
			}

			File dir = new File("src/main/resources/schemes/derived/section" + Math.max(N, Math.max(M, P)));
			dir.mkdirs();
			File f = new File(dir, SchemeIO.canonicalName(built, "derived_serendipitous"));
			if (replayOk) {
				SchemeIO.writeStub(built, f, lineage);
				log.info("CONFIRMED+stub ⟨{},{},{}⟩ = {} (a{}) → {} (replay verified)",
						N, M, P, built.r, adds, f.getName());
			} else {
				SchemeIO.write(built, f, lineage);
				log.info("CONFIRMED+matrices ⟨{},{},{}⟩ = {} (a{}) → {} (replay NOT reproducible, full matrices)",
						N, M, P, built.r, adds, f.getName());
			}
			written++;
		}
		log.info("=== materialised {} / {} bud-win schemes ===", written, WINS.length);
	}

	private static NonCubicBilinearAlgorithm richest(FieldAwareLookup lookup, int n, int m, int p)
			throws Exception {
		NonCubicBilinearAlgorithm best = null;
		int bestScore = -1;
		for (var path : lookup.findFiles(n, m, p)) {
			try {
				NonCubicBilinearAlgorithm a = SchemeIO.read(path.toFile());
				if (a.n != n || a.m != m || a.p != p) {
					continue;  // exact orientation only
				}
				int sc = BudBaseFactory.budScore(a);
				if (sc > bestScore) {
					bestScore = sc;
					best = a;
				}
			} catch (Exception ignored) {
				// skip
			}
		}
		return best;
	}
}
