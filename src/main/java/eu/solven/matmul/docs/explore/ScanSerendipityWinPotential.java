package eu.solven.matmul.docs.explore;

import eu.solven.matmul.recombination.Recombination;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import lombok.extern.slf4j.Slf4j;

/**
 * Base-AGNOSTIC serendipity win-potential screen — the "drive a directed search"
 * layer (companion to {@code ScanSerendipityGapMap}, which only prices the
 * sub-additivity; and to {@code docs.SerendipitousSweep}, which only tests buds
 * EXISTING bases already carry).
 *
 * <p>For every catalog target {@code T} and every factorization
 * {@code T = ⟨n₁,m₁,p₁⟩ ⊗ ⟨n₂,m₂,p₂⟩}, asks: <em>what bud structure would a base
 * of shape ⟨n₁,m₁,p₁⟩ need so that the serendipitous product beats SOTA(T)?</em>
 * A base of rank {@code r₁ = R(base)} with {@code b} disjoint size-{@code k} buds
 * on one axis predicts
 *
 * <pre>
 *   r_s = r₁·R(inner) − b·σ(k),    σ(k) = k·R(inner) − R(k-enlarged inner)
 * </pre>
 *
 * so the WIN ({@code r_s < SOTA(T)}) needs at least
 *
 * <pre>
 *   b* = ⌊(r₁·R(inner) − SOTA(T)) / σ(k)⌋ + 1     (σ(k) &gt; 0)
 * </pre>
 *
 * size-{@code k} buds — FEASIBLE only if they fit the base
 * ({@code b*·k ≤ r₁}) and stay within the per-axis capacity bound
 * ({@code b*·k ≤ 2·(r₁ − d_axis)}, {@code d_U=n₁m₁, d_V=m₁p₁, d_W=n₁p₁}; see
 * {@code references/BUD_STRUCTURE_THEORY.md}).
 *
 * <p>Every feasible row is a <strong>pre-certified meta-flip goal</strong>:
 * "find a rank-{@code r₁} ⟨n₁,m₁,p₁⟩ scheme with {@code b*} size-{@code k}
 * {axis}-buds, and the serendipitous product is GUARANTEED to beat ⟨T⟩'s SOTA."
 * The hit is certified by arithmetic; only the realizability (does such a base
 * exist at that rank?) is left to the directed flip / ALS search — the rigidity
 * wall. Small {@code b*} = easiest realizability = best target.</p>
 *
 * <p><strong>Honesty:</strong> SOTA is {@code min(findRank(T), best plain-Kron)}
 * — so if the catalog already books this serendipity (entry ≤ plain-Kron), the
 * required σ is huge and no small-{@code b*} row appears; that absence is the
 * "already-banked" verdict, not a gap. Args: {@code maxDim=32 kMax=4 maxBuds=6}.</p>
 */
@Slf4j
public class ScanSerendipityWinPotential {

	enum BudAxis { U, V, W }

	record Goal(int tn, int tm, int tp, long sota, int bn, int bm, int bp, int r1,
			int in, int im, int ip, int rInner, BudAxis axis, int k, long sigma,
			int bStar, long rs, boolean fits, boolean withinCapacity,
			long rsRealistic, boolean survivesRankCost) {
		long margin() {
			return sota - rs;
		}
	}

	public static void main(String[] args) {
		int maxDim = args.length > 0 ? Integer.parseInt(args[0]) : 32;
		int kMax = args.length > 1 ? Integer.parseInt(args[1]) : 4;
		int maxBuds = args.length > 2 ? Integer.parseInt(args[2]) : 6;
		Field field = Field.Q;
		FieldAwareLookup lk = new FieldAwareLookup(field);
		final int UNKNOWN = Recombination.SotaResolver.UNKNOWN_RANK;

		// Best feasible goal (smallest b*) per target, plus all feasible goals.
		List<Goal> feasible = new ArrayList<>();
		long targetsKnown = 0;
		long targetsWithGoal = 0;

		for (int tn = 2; tn <= maxDim; tn++) {
			for (int tm = tn; tm <= maxDim; tm++) {
				for (int tp = tm; tp <= maxDim; tp++) {
					long sota = trueSota(lk, tn, tm, tp, UNKNOWN);
					if (sota < 0) {
						continue;
					}
					targetsKnown++;
					boolean any = false;
					// Factorize each axis: base dim divides target dim.
					for (int bn : divisors(tn))
						for (int bm : divisors(tm))
							for (int bp : divisors(tp)) {
								int in = tn / bn, im = tm / bm, ip = tp / bp;
								if (bn * bm * bp == 1 || in * im * ip == 1) {
									continue; // trivial base or trivial inner
								}
								int r1 = lk.findRank(bn, bm, bp);
								int rInner = lk.findRank(in, im, ip);
								if (r1 >= UNKNOWN || rInner >= UNKNOWN) {
									continue;
								}
								long plainKron = (long) r1 * rInner;
								if (plainKron <= sota) {
									continue; // even with 0 buds this factorization can't strictly win;
									          // buds only help if plain-Kron is ABOVE sota
								}
								for (int k = 2; k <= kMax; k++) {
									any |= tryGoal(lk, feasible, BudAxis.U, tn, tm, tp, sota,
											bn, bm, bp, r1, in, im, k * ip, in, im, ip, rInner, k,
											maxBuds, UNKNOWN);
									any |= tryGoal(lk, feasible, BudAxis.V, tn, tm, tp, sota,
											bn, bm, bp, r1, k * in, im, ip, in, im, ip, rInner, k,
											maxBuds, UNKNOWN);
									any |= tryGoal(lk, feasible, BudAxis.W, tn, tm, tp, sota,
											bn, bm, bp, r1, in, k * im, ip, in, im, ip, rInner, k,
											maxBuds, UNKNOWN);
								}
							}
					if (any) {
						targetsWithGoal++;
					}
				}
			}
		}

		long survivors = feasible.stream().filter(Goal::survivesRankCost).count();

		log.info("=== Serendipity win-potential ({}), targets≤{}, k≤{}, b*≤{} ===",
				field, maxDim, kMax, maxBuds);
		log.info("known targets: {} | targets with a NAIVE goal (bud-at-optimal): {} | naive goals: {}",
				targetsKnown, targetsWithGoal, feasible.size());
		log.info("goals surviving the 1-rank/bud creation charge (REALISTIC): {}", survivors);

		log.info("");
		log.info("--- NAIVE goals, easiest first — ALL assume buds at OPTIMAL rank ---");
		log.info("    (at optimal rank schemes are bud-sterile; ⟨2,2,2⟩=7 provably so)");
		feasible.stream()
				.sorted(Comparator.<Goal>comparingInt(Goal::bStar)
						.thenComparing(Comparator.comparingLong(Goal::margin).reversed()))
				.limit(12)
				.forEach(g -> log.info(
						"  ⟨{},{},{}⟩ SOTA={} ← r_s={} (−{}) | base ⟨{},{},{}⟩ r1={} ⊗ ⟨{},{},{}⟩ R={}"
								+ " | {}×{}-bud k={} σ={}/bud | rs_real(+rank)={} {}",
						g.tn(), g.tm(), g.tp(), g.sota(), g.rs(), g.margin(),
						g.bn(), g.bm(), g.bp(), g.r1(), g.in(), g.im(), g.ip(), g.rInner(),
						g.bStar(), g.axis(), g.k(), g.sigma(), g.rsRealistic(),
						g.survivesRankCost() ? "← SURVIVES" : "✗ dies (≥SOTA)"));

		log.info("");
		if (survivors == 0) {
			log.info("--- VERDICT: 0 goals survive bud rank-cost ---");
			log.info("    Every naive goal needs a bud at OPTIMAL rank. Charging even 1 rank/bud");
			log.info("    (σ < R(inner) by the 28% sub-additivity ceiling) pushes r_s ≥ plain-Kron ≥ SOTA.");
			log.info("    Serendipity self-improvement is CLOSED over Q≤32 unless a base is bud-rich");
			log.info("    AT its optimal rank — and SerendipitousSweep found 0 such (de Groote-style");
			log.info("    rigidity at optimum). This is the provable form of that empirical 0.");
		} else {
			log.info("--- {} REALISTIC goals (survive rank-cost) — directed meta-flip targets ---", survivors);
			feasible.stream().filter(Goal::survivesRankCost)
					.sorted(Comparator.comparingLong(g -> g.rsRealistic() - g.sota()))
					.limit(30)
					.forEach(g -> log.info(
							"  ⟨{},{},{}⟩ SOTA={} ← rs_real={} | base ⟨{},{},{}⟩ r1={} ⊗ ⟨{},{},{}⟩"
									+ " | {}×{}-bud k={} σ={}",
							g.tn(), g.tm(), g.tp(), g.sota(), g.rsRealistic(),
							g.bn(), g.bm(), g.bp(), g.r1(), g.in(), g.im(), g.ip(),
							g.bStar(), g.axis(), g.k(), g.sigma()));
		}
	}

	private static boolean tryGoal(FieldAwareLookup lk, List<Goal> out, BudAxis axis,
			int tn, int tm, int tp, long sota, int bn, int bm, int bp, int r1,
			int en, int em, int ep, int in, int im, int ip, int rInner, int k,
			int maxBuds, int unknown) {
		int enlarged = lk.findRank(en, em, ep);
		if (enlarged >= unknown) {
			return false;
		}
		long sigma = (long) k * rInner - enlarged;
		if (sigma <= 0) {
			return false; // this bud saves nothing against this inner
		}
		long plainKron = (long) r1 * rInner;
		// b* = minimal #buds to push r_s strictly below sota.
		long bStar = (plainKron - sota) / sigma + 1;
		if (bStar > maxBuds) {
			return false;
		}
		long usedTerms = bStar * k;
		boolean fits = usedTerms <= r1;
		if (!fits) {
			return false; // can't even place that many size-k buds in r1 terms
		}
		int dAxis = switch (axis) {
			case U -> bn * bm;
			case V -> bm * bp;
			case W -> bn * bp;
		};
		boolean withinCapacity = usedTerms <= 2L * (r1 - dAxis);
		long rs = plainKron - bStar * sigma;
		// REALISTIC cost: the win above assumes the bud structure exists at the
		// base's OPTIMAL rank r1 — but at optimal rank schemes are bud-sterile
		// (capacity → 0; ⟨2,2,2⟩=7 PROVABLY so by de Groote uniqueness). Creating
		// b buds costs rank: charge the minimal 1 rank/bud (an optimistic lower
		// bound — a size-k class realistically costs more). Then
		//   rs_real = (r1 + bStar)·R(inner) − bStar·σ = plainKron + bStar·(R(inner) − σ).
		// Since the 28% sub-additivity ceiling forces σ ≤ ~0.28·k·R(inner) < R(inner)
		// for k≤3 (and ≈R(inner) at the k=4 ceiling), R(inner) − σ ≥ 0, so charging
		// even ONE rank pushes rs_real ≥ plainKron ≥ sota → the win evaporates.
		long rsRealistic = plainKron + bStar * (rInner - sigma);
		boolean survivesRankCost = rsRealistic < sota;
		out.add(new Goal(tn, tm, tp, sota, bn, bm, bp, r1, in, im, ip, rInner, axis, k,
				sigma, (int) bStar, rs, true, withinCapacity, rsRealistic, survivesRankCost));
		return true;
	}

	private static List<Integer> divisors(int n) {
		List<Integer> d = new ArrayList<>();
		for (int i = 1; i <= n; i++) {
			if (n % i == 0) {
				d.add(i);
			}
		}
		return d;
	}

	/** min(catalog rank, best plain-Kron over factorizations); -1 if unknown. */
	static long trueSota(FieldAwareLookup lk, int n, int m, int p, int unknown) {
		int cat = lk.findRank(n, m, p);
		long best = cat >= unknown ? -1 : cat;
		for (int n1 : divisors(n))
			for (int m1 : divisors(m))
				for (int p1 : divisors(p)) {
					if (n1 * m1 * p1 == 1 || (n / n1) * (m / m1) * (p / p1) == 1) {
						continue;
					}
					int r1 = lk.findRank(n1, m1, p1);
					int r2 = lk.findRank(n / n1, m / m1, p / p1);
					if (r1 >= unknown || r2 >= unknown) {
						continue;
					}
					long prod = (long) r1 * r2;
					if (best < 0 || prod < best) {
						best = prod;
					}
				}
		return best;
	}
}
