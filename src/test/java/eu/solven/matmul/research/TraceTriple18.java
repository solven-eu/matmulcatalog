package eu.solven.matmul.research;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Optional;
import java.util.Set;

import eu.solven.matmul.recombination.BlockSplitSearch;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.recombination.Recombination;
import eu.solven.matmul.search.RecursiveClosureSota;

/**
 * Trace the recursive-closure path to {@code R(⟨18,18,18⟩) ≤ 3200}
 * reported by {@link StrategySurvey}: it comes from
 * {@code concat-p[9,9]} which requires {@code R(⟨18,18,9⟩) ≤ 1600},
 * which is itself a recursive prediction. This walker recursively
 * expands the chosen strategy at each shape until it bottoms out at
 * a direct catalog hit, exposing every link in the chain.
 *
 * <p>Goal: determine whether {@code R(⟨18,18,9⟩) = 1600} is backed by
 * an actual scheme on disk or is a phantom — a prediction whose
 * sub-shapes don't exist as schemes.</p>
 */
public final class TraceTriple18 {

	public static void main(String[] args) throws Exception {
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		List<BlockSplitSearch.NamedBase> pool = BlockSplitSearch.defaultPool();
		RecursiveClosureSota sota = new RecursiveClosureSota(lookup, pool, true, true);

		int[][] toTrace = {
				{18, 18, 18},
				{18, 18, 9},
				{18, 9, 18},
				{9, 18, 18},
		};

		Set<String> visited = new HashSet<>();
		for (int[] s : toTrace) {
			trace(s[0], s[1], s[2], 0, lookup, pool, sota, visited);
			System.out.println();
		}
	}

	private static void trace(int a, int b, int c, int depth,
			FieldAwareLookup lookup, List<BlockSplitSearch.NamedBase> pool,
			Recombination.SotaResolver sota, Set<String> visited) {
		String key = canon(a, b, c);
		String indent = "  ".repeat(depth);
		int direct = lookup.find(a, b, c).map(x -> x.r).orElse(-1);
		int closure = sota.getRank(a, b, c);
		System.out.printf("%s⟨%d,%d,%d⟩  closure=%d  catalog=%s%n",
				indent, a, b, c, closure,
				direct < 0 ? "MISSING" : Integer.toString(direct));

		if (depth > 4) {
			System.out.printf("%s  …depth cap%n", indent);
			return;
		}
		if (!visited.add(key)) {
			System.out.printf("%s  …seen, stop%n", indent);
			return;
		}
		if (direct >= 0 && direct == closure) {
			System.out.printf("%s  ✓ catalog scheme materialises this rank%n", indent);
			return;
		}

		Optional<BlockSplitSearch.NonCubicStrategy> picked =
				BlockSplitSearch.findBestStrategy(a, b, c, pool, sota, true);
		if (picked.isEmpty()) {
			System.out.printf("%s  ⚠ no strategy%n", indent);
			return;
		}
		BlockSplitSearch.NonCubicStrategy ps = picked.get();
		System.out.printf("%s  → %s  rank=%d%n", indent, ps.label(), ps.rank());

		// Expand sub-shapes by strategy type.
		List<int[]> children = subShapes(ps, a, b, c);
		for (int[] ch : children) {
			trace(ch[0], ch[1], ch[2], depth + 1, lookup, pool, sota, visited);
		}
	}

	private static List<int[]> subShapes(BlockSplitSearch.NonCubicStrategy ps, int a, int b, int c) {
		List<int[]> out = new ArrayList<>();
		if (ps.kronecker() != null) {
			var k = ps.kronecker();
			out.add(new int[] { k.n1(), k.m1(), k.p1() });
			out.add(new int[] { k.n2(), k.m2(), k.p2() });
		} else if (ps.concat() != null) {
			var cc = ps.concat();
			int leftSize = cc.leftSize();
			int rightSize = cc.rightSize();
			int axis = cc.axis();  // 0 = n-axis, 2 = p-axis
			if (axis == 0) {
				out.add(new int[] { leftSize, b, c });
				out.add(new int[] { rightSize, b, c });
			} else {
				out.add(new int[] { a, b, leftSize });
				out.add(new int[] { a, b, rightSize });
			}
		} else if (ps.recombination() != null) {
			var rec = ps.recombination();
			// allocA/B/C: each is array of block sizes summing to a/b/c.
			// The sub-shape consulted is the max-block product (Strassen base × inner shape).
			int[] aa = rec.allocA(), bb = rec.allocB(), cc = rec.allocC();
			System.out.print("    allocA=" + java.util.Arrays.toString(aa)
					+ " allocB=" + java.util.Arrays.toString(bb)
					+ " allocC=" + java.util.Arrays.toString(cc) + " ");
			System.out.println("base=" + rec.baseLabel());
			// emit max block as primary child
			int maxA = max(aa), maxB = max(bb), maxC = max(cc);
			out.add(new int[] { maxA, maxB, maxC });
		}
		return out;
	}

	private static int max(int[] xs) {
		int m = Integer.MIN_VALUE;
		for (int x : xs) if (x > m) m = x;
		return m;
	}

	private static String canon(int a, int b, int c) {
		int[] s = { a, b, c };
		java.util.Arrays.sort(s);
		return s[0] + "x" + s[1] + "x" + s[2];
	}
}
