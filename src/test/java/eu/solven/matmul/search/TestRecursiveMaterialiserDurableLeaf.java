package eu.solven.matmul.search;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.stream.Stream;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Regression guard for the durable building-block leaf ref
 * ({@code RecursiveMaterialiser.durableLeafRef}).
 *
 * <p>The silent bug: a building block sourced from a lineage-only STUB was pinned
 * by {@code contentHash(orientedInMemoryAlg)} — a hash stamped on NO file — so
 * {@code resolveLeaf} could never resolve it, silently fell back to shape-best, the
 * replayed rank drifted, and the write-guard discarded the win (~46% of the non-5
 * sweep's candidates were refused, all stub-sourced). The fix pins the SOURCE file's
 * native shape + stamped hash, wrapped in {@link Lineage.OrientAs} for reorientation.
 *
 * <p>This test asserts the foundation the fix relies on: such a ref resolves
 * bit-exactly in BOTH the native and a permuted orientation. Were the resolution to
 * regress to the shape-best fallback, the permuted replay would drift off the pinned
 * rank and the assertion would fail.
 */
public class TestRecursiveMaterialiserDurableLeaf {

	/** First non-cubic stub (maxDim&gt;16, lineage-only, distinct dims) on disk. */
	private static Path findNonCubicStub() throws IOException {
		Path root = Path.of("src/main/resources/schemes");
		try (Stream<Path> files = Files.walk(root)) {
			return files.filter(p -> p.toString().endsWith(".json"))
					.filter(TestRecursiveMaterialiserDurableLeaf::isDistinctDimStub)
					.findFirst()
					.orElseThrow(() -> new AssertionError("no non-cubic stub found under " + root));
		}
	}

	private static boolean isDistinctDimStub(Path p) {
		try {
			var root = SchemeIO.parseJson(p.toFile());
			if (root.has("u") || root.has("u_sparse")) return false; // must be lineage-only
			if (!root.has("lineage") || SchemeIO.readHash(root) == null) return false;
			var n = root.get("n");
			if (n == null || n.size() != 3) return false;
			int a = n.get(0).asInt(), b = n.get(1).asInt(), c = n.get(2).asInt();
			return Math.max(a, Math.max(b, c)) > 16 && !(a == b && b == c);
		} catch (Exception e) {
			return false;
		}
	}

	@Test
	public void stubLeaf_pinnedBySourceHash_replaysBitExactInAnyOrientation() throws Exception {
		Path stub = findNonCubicStub();
		var root = SchemeIO.parseJson(stub.toFile());
		String hash = SchemeIO.readHash(root);
		int sn = root.get("n").get(0).asInt();
		int sm = root.get("n").get(1).asInt();
		int sp = root.get("n").get(2).asInt();
		int rank = root.get("m").asInt();
		String shapeRef = sn + "x" + sm + "x" + sp;

		LineageReplayer replayer =
				LineageReplayer.withDefaultPool(new FieldAwareLookup("R"));

		// (1) native orientation: bare source-hash atom must resolve to the exact stub.
		NonCubicBilinearAlgorithm nativeAlg =
				replayer.replay(new Lineage.Atom(shapeRef + "@" + hash));
		assertThat(new int[] { nativeAlg.n, nativeAlg.m, nativeAlg.p })
				.as("native replay shape").containsExactly(sn, sm, sp);
		assertThat(nativeAlg.r).as("native replay rank").isEqualTo(rank);

		// (2) permuted orientation via OrientAs: must reconstruct the SAME scheme at the
		// permuted shape and the SAME rank — i.e. resolve the pin, never shape-best drift.
		Lineage.Node reoriented =
				new Lineage.OrientAs(new Lineage.Atom(shapeRef + "@" + hash), sm, sn, sp);
		NonCubicBilinearAlgorithm permAlg = replayer.replay(reoriented);
		assertThat(new int[] { permAlg.n, permAlg.m, permAlg.p })
				.as("permuted replay shape").containsExactly(sm, sn, sp);
		assertThat(permAlg.r).as("permuted replay rank (no shape-best drift)").isEqualTo(rank);
	}
}
