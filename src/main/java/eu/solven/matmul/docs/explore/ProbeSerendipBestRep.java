package eu.solven.matmul.docs.explore;

import java.io.File;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.catalog.SerendipitousBudProduct;
import lombok.extern.slf4j.Slf4j;

/**
 * Enumerate every catalog scheme at a base shape/rank and price its
 * serendipitous product against a given inner — to find whether a bud-richer
 * <em>representative</em> already exists (buds are GL-invariant, so only a
 * distinct GL-orbit can carry more σ). Args:
 * {@code bN bM bP baseRank iN iM iP sota}.
 */
@Slf4j
public class ProbeSerendipBestRep {
	public static void main(String[] args) {
		int bN = Integer.parseInt(args[0]), bM = Integer.parseInt(args[1]), bP = Integer.parseInt(args[2]);
		int baseRank = Integer.parseInt(args[3]);
		int n2 = Integer.parseInt(args[4]), m2 = Integer.parseInt(args[5]), p2 = Integer.parseInt(args[6]);
		long sota = Long.parseLong(args[7]);
		FieldAwareLookup q = new FieldAwareLookup(Field.Q);
		List<NonCubicBilinearAlgorithm> reps = new ArrayList<>();
		for (Path pth : q.findFiles(bN, bM, bP)) {
			try {
				NonCubicBilinearAlgorithm a = SchemeIO.read(pth.toFile());
				if (a.r == baseRank) {
					reps.add(a);
				}
			} catch (Exception ignore) {
				// stub / unreadable
			}
		}
		log.info("base ⟨{},{},{}⟩ r={}: {} explicit rep(s); inner ⟨{},{},{}⟩; SOTA(target)={}",
				bN, bM, bP, baseRank, reps.size(), n2, m2, p2, sota);
		long best = Long.MAX_VALUE;
		String bestFile = null;
		for (int i = 0; i < reps.size(); i++) {
			NonCubicBilinearAlgorithm a = reps.get(i);
			long cost = SerendipitousBudProduct.serendipitousCost(a, q, n2, m2, p2);
			long sigma = (long) a.r * q.findRank(n2, m2, p2) - cost;
			String tag = cost < sota ? "  *** BEATS SOTA ***" : (cost == sota ? "  (ties)" : "");
			log.info("  rep {}/{}: predicted={} σ={}{}", i + 1, reps.size(), cost, sigma, tag);
			if (cost < best) {
				best = cost;
			}
		}
		log.info("BEST rep predicts ⟨{},{},{}⟩ ≤ {} (SOTA {}) — {}",
				bN * n2, bM * m2, bP * p2, best, sota,
				best < sota ? "*** IMPROVEMENT AVAILABLE ***" : "no catalog rep beats SOTA");
	}
}
