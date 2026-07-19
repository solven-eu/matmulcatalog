package eu.solven.matmul.docs.explore;

import eu.solven.matmul.recombination.Recombination;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.ProjectionSearch;
import lombok.extern.slf4j.Slf4j;

/**
 * PROJECTION-margin frontier scan (user 2026-06-12: the projection analogue of
 * the serendipity study — "meta-flips can find bases with slightly higher rank
 * but much better projection margin"). For every catalog parent shape P (dims
 * ≤ {@code maxDim}) and every 1-drop target T = P−eᵢ:
 *
 * <pre>  gap(P→T) = projectedRank(best catalog P-scheme → T) − R_cat(T)</pre>
 *
 * gap = 0: projection already matches the catalog (or produced it). gap &gt; 0:
 * a margin-richer P-base at rank r+δ wins T iff it reaches margin μ ≥ r + δ −
 * R_cat(T) + 1 — rank's exchange rate against margin is exactly 1, so EVERY δ
 * is in play (no serendipity-style finite window). "closed" = R_cat(T) equals
 * the flattening lower bound max(nm,mp,np): nothing can project below it.
 *
 * <p>Output: walkable candidates sorted by gap (small gap = cheapest demo for
 * a {@code --objective=project} flip walk), with the needed Δμ at δ=0/1/2.</p>
 *
 * <p>Args: {@code [maxDim=7]}.</p>
 */
@Slf4j
public class ScanProjectionFrontier {

	record Row(int[] parent, int[] target, int r, int proj, int rt, int lb, long gap) {}

	/** Every readable scheme file at the minimal rank among READABLE files —
	 *  the catalog-best rank may be claimed by a lineage stub with no factor
	 *  matrices; the projection study needs concrete representatives. */
	static List<NonCubicBilinearAlgorithm> minRankRepresentatives(FieldAwareLookup q,
			int n, int m, int p) {
		List<NonCubicBilinearAlgorithm> readable = new ArrayList<>();
		int best = Integer.MAX_VALUE;
		for (java.nio.file.Path f : q.findFiles(n, m, p)) {
			try {
				NonCubicBilinearAlgorithm alg = eu.solven.matmul.catalog.SchemeIO.read(f.toFile());
				readable.add(alg);
				best = Math.min(best, alg.r);
			} catch (Exception stubOrUnreadable) {
				// lineage stubs carry no factor matrices — skip
			}
		}
		final int min = best;
		readable.removeIf(a -> a.r != min);
		// Union with the lookup's own resolution: some shapes have no same-shape
		// file at all (stub replay / derivation) — find() is their only concrete
		// representative, and it may also out-rank every file.
		q.find(n, m, p).ifPresent(alg -> {
			if (readable.isEmpty() || alg.r < min) {
				readable.clear();
				readable.add(alg);
			} else if (alg.r == min) {
				readable.add(alg);
			}
		});
		return readable;
	}

	public static void main(String[] args) {
		int maxDim = args.length > 0 ? Integer.parseInt(args[0]) : 7;
		FieldAwareLookup q = new FieldAwareLookup(Field.Q);
		List<Row> rows = new ArrayList<>();
		List<Row> matchedRows = new ArrayList<>();
		int closed = 0;
		int matched = 0;
		for (int n = 2; n <= maxDim; n++) {
			for (int m = 2; m <= maxDim; m++) {
				for (int p = 2; p <= maxDim; p++) {
					// ALL minimal-rank representatives, not just the (arbitrary)
					// catalog-best file: margin is heavily scheme-dependent at
					// fixed rank (rank-23 ⟨3,3,3⟩: kauers_2026 projects to 15,
					// the other six representatives to 16).
					List<NonCubicBilinearAlgorithm> parents = minRankRepresentatives(q, n, m, p);
					if (parents.isEmpty()) {
						continue;
					}
					NonCubicBilinearAlgorithm parent = parents.get(0);
					for (int ax = 0; ax < 3; ax++) {
						int tn = n - (ax == 0 ? 1 : 0);
						int tm = m - (ax == 1 ? 1 : 0);
						int tp = p - (ax == 2 ? 1 : 0);
						if (tn < 2 || tm < 2 || tp < 2) {
							continue;
						}
						int rt = q.findRank(tn, tm, tp);
						if (rt >= Recombination.SotaResolver.UNKNOWN_RANK) {
							continue;
						}
						int lb = Math.max(tn * tm, Math.max(tm * tp, tn * tp));
						long proj = Long.MAX_VALUE;
						for (NonCubicBilinearAlgorithm rep : parents) {
							long pr = ProjectionSearch.projectedRank(rep, tn, tm, tp, 1);
							if (pr >= 0) {
								proj = Math.min(proj, pr);
							}
						}
						if (proj == Long.MAX_VALUE) {
							continue;
						}
						long gap = proj - rt;
						if (gap <= 0) {
							matched++;
							// A projection that ACHIEVES the catalog rank on a
							// target with room above its LB: one extra margin
							// point at the same base rank improves the catalog.
							if (rt > lb) {
								matchedRows.add(new Row(new int[] { n, m, p },
										new int[] { tn, tm, tp }, parent.r, (int) proj,
										rt, lb, gap));
							}
							continue;
						}
						if (rt <= lb) {
							closed++;
							continue;
						}
						rows.add(new Row(new int[] { n, m, p }, new int[] { tn, tm, tp },
								parent.r, (int) proj, rt, lb, gap));
					}
				}
			}
		}
		rows.sort(Comparator.comparingLong(Row::gap));
		for (Row row : rows.subList(0, Math.min(25, rows.size()))) {
			// Needed margin at base rank r+δ to strictly beat the catalog target.
			String needs = String.format("μ≥%d/%d/%d at δ=0/1/2",
					row.r() - row.rt() + 1, row.r() + 1 - row.rt() + 1, row.r() + 2 - row.rt() + 1);
			log.info("gap {} ⟨{},{},{}⟩(r={})→⟨{},{},{}⟩: proj {} vs catalog {} (LB {}, room {}) — {}",
					row.gap(), row.parent()[0], row.parent()[1], row.parent()[2], row.r(),
					row.target()[0], row.target()[1], row.target()[2],
					row.proj(), row.rt(), row.lb(), row.rt() - row.lb(), needs);
		}
		matchedRows.sort(Comparator.comparingInt(r -> -(r.rt() - r.lb())));
		for (Row row : matchedRows.subList(0, Math.min(25, matchedRows.size()))) {
			log.info("MATCHED ⟨{},{},{}⟩(r={})→⟨{},{},{}⟩: proj {} = catalog {} (LB {}, room {})"
					+ " — Δμ=+1 at δ=0 would improve the catalog",
					row.parent()[0], row.parent()[1], row.parent()[2], row.r(),
					row.target()[0], row.target()[1], row.target()[2],
					row.proj(), row.rt(), row.lb(), row.rt() - row.lb());
		}
		log.info("projection frontier (dims ≤ {}): {} matched-or-better ({} with room shown), "
				+ "{} closed-at-LB, {} open gaps ({} shown)", maxDim, matched,
				Math.min(25, matchedRows.size()), closed, rows.size(),
				Math.min(25, rows.size()));
		// Gap histogram: how far does the projected rank sit above the catalog SOTA
		// of the target, across ALL open gaps (not just the 25 displayed)?
		java.util.TreeMap<Long, Integer> hist = new java.util.TreeMap<>();
		for (Row row : rows) {
			hist.merge(row.gap(), 1, Integer::sum);
		}
		StringBuilder sb = new StringBuilder("gap histogram (proj − R_cat(T)) over ")
				.append(rows.size()).append(" open gaps: ");
		hist.forEach((g, c) -> sb.append(g).append("→").append(c).append("  "));
		log.info(sb.toString().trim());
	}
}
