package eu.solven.matmul.docs.verify;

import java.util.ArrayList;
import java.util.List;

import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.recombination.Recombination;

/**
 * Find the SMALLEST matmul shape where the current catalog rank is NOT midpoint-convex along an
 * axis — i.e. {@code 2·R(center) > R(left)+R(right)} for {@code center=(a,b,c)},
 * {@code left=(a,b,c-1)}, {@code right=(a,b,c+1)} (the middle "bumps up" above the chord, so a
 * small+big split is cheaper than two mediums — the convexity heuristic fails there). Only triples
 * whose three ranks are all KNOWN (not the cubic {@code a·b·c} fallback) count. Iterating every
 * {@code (a,b,c)} and testing the c-axis covers all axes by relabeling. Smallest = least
 * {@code (a+b+c, max, lexicographic)}.
 *
 * <p>Run: {@code mvn -q -ntp exec:java
 * -Dexec.mainClass=eu.solven.matmul.docs.verify.FindRankNonConvexity}</p>
 */
public final class FindRankNonConvexity {
	private FindRankNonConvexity() {}

	public static void main(String[] args) {
		FieldAwareLookup lookup = new FieldAwareLookup(Field.R);
		int MAX = args.length > 0 ? Integer.parseInt(args[0]) : 12;

		record Violation(int a, int b, int c, int rl, int rc, int rr, int slack) {}
		List<Violation> hits = new ArrayList<>();
		for (int a = 1; a <= MAX; a++)
			for (int b = 1; b <= MAX; b++)
				for (int c = 2; c <= MAX; c++) {
					int rl = lookup.findRank(a, b, c - 1), rc = lookup.findRank(a, b, c), rr = lookup.findRank(a, b, c + 1);
					int U = Recombination.SotaResolver.UNKNOWN_RANK;
					if (rl == U || rc == U || rr == U) continue;
					int slack = 2 * rc - (rl + rr); // > 0 ⇒ non-convex (bump)
					if (slack > 0) hits.add(new Violation(a, b, c, rl, rc, rr, slack));
				}

		hits.sort((x, y) -> {
			int sx = x.a + x.b + x.c, sy = y.a + y.b + y.c;
			if (sx != sy) return Integer.compare(sx, sy);
			int mx = Math.max(x.a, Math.max(x.b, x.c)), my = Math.max(y.a, Math.max(y.b, y.c));
			if (mx != my) return Integer.compare(mx, my);
			if (x.a != y.a) return Integer.compare(x.a, y.a);
			if (x.b != y.b) return Integer.compare(x.b, y.b);
			return Integer.compare(x.c, y.c);
		});

		System.out.printf("Rank non-convexity scan over R⟨a,b,c⟩, dims ≤ %d, known ranks only.%n", MAX);
		System.out.printf("A violation at center c means: R⟨a,b,c-1⟩ + R⟨a,b,c+1⟩ < 2·R⟨a,b,c⟩"
				+ " (split into small+big beats two mediums).%n%n");
		if (hits.isEmpty()) { System.out.println("No non-convexity found in range (all known triples convex)."); return; }
		System.out.printf("SMALLEST: ⟨%d,%d,%d⟩  R⟨%d,%d,%d⟩=%d + R⟨%d,%d,%d⟩=%d = %d  <  2·R⟨%d,%d,%d⟩=%d  (slack %d)%n%n",
				hits.get(0).a, hits.get(0).b, hits.get(0).c,
				hits.get(0).a, hits.get(0).b, hits.get(0).c - 1, hits.get(0).rl,
				hits.get(0).a, hits.get(0).b, hits.get(0).c + 1, hits.get(0).rr, hits.get(0).rl + hits.get(0).rr,
				hits.get(0).a, hits.get(0).b, hits.get(0).c, hits.get(0).rc, hits.get(0).slack);
		System.out.println("first 15 violations (smallest-first):");
		int shown = 0;
		for (var v : hits) {
			if (shown++ >= 15) break;
			System.out.printf("  ⟨%d,%d,%d⟩: %d + %d = %d < 2·%d = %d  (slack +%d)%n",
					v.a, v.b, v.c, v.rl, v.rr, v.rl + v.rr, v.rc, 2 * v.rc, v.slack);
		}
		System.out.printf("%ntotal violations found: %d%n", hits.size());
	}
}
