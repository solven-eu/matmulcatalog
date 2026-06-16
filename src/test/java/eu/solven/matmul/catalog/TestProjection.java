package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;

/**
 * Projection operator (#159 / ROADMAP): restrict a larger scheme to kept
 * indices + DCE. Pins the drop-position/axis dependence on Laderman
 * ⟨3,3,3⟩=23 and checks every projection computes the smaller matmul exactly.
 */
public class TestProjection {

	private static NonCubicBilinearAlgorithm laderman() throws Exception {
		return SchemeIO.read(eu.solven.matmul.catalog.SchemeResolver.byHint(
				"src/main/resources/schemes/known/section3/laderman_1976-3x3x3_m23_a98_b0.json"));
	}

	@Test
	public void n_axis_projection_is_position_dependent() throws Exception {
		NonCubicBilinearAlgorithm a = laderman();
		int[] all = { 0, 1, 2 };
		// Drop n-index 0 → 18; drop 1 or 2 → 15 (= optimum R(⟨2,3,3⟩)).
		NonCubicBilinearAlgorithm d0 = Compose.project(a, Compose.keepExcept(3, 0), all, all);
		NonCubicBilinearAlgorithm d1 = Compose.project(a, Compose.keepExcept(3, 1), all, all);
		NonCubicBilinearAlgorithm d2 = Compose.project(a, Compose.keepExcept(3, 2), all, all);
		System.out.printf("n-drop ranks: 0→%d 1→%d 2→%d%n", d0.r, d1.r, d2.r);
		// Position-dependent: dropping index 0 leaves 18, dropping 1 or 2 leaves 16.
		assertThat(d0.r).isEqualTo(18);
		assertThat(d1.r).isEqualTo(16);
		assertThat(d2.r).isEqualTo(16);
		for (NonCubicBilinearAlgorithm s : new NonCubicBilinearAlgorithm[] { d0, d1, d2 }) {
			assertThat(s.n).isEqualTo(2);
			assertThat(s.m).isEqualTo(3);
			assertThat(s.p).isEqualTo(3);
			assertThat(Verifier.isExactNonCubic(s)).as("projection must be exact ⟨2,3,3⟩").isTrue();
		}
	}

	@Test
	public void inner_axis_projection_differs_from_outer() throws Exception {
		NonCubicBilinearAlgorithm a = laderman();
		int[] all = { 0, 1, 2 };
		// Contraction axis m: drop 0→16, 1→18, 2→16 — output ⟨3,2,3⟩, not optimum 15.
		NonCubicBilinearAlgorithm m0 = Compose.project(a, all, Compose.keepExcept(3, 0), all);
		NonCubicBilinearAlgorithm m1 = Compose.project(a, all, Compose.keepExcept(3, 1), all);
		System.out.printf("m-drop ranks: 0→%d 1→%d%n", m0.r, m1.r);
		assertThat(m0.r).isEqualTo(16);
		assertThat(m1.r).isEqualTo(18);
		assertThat(m0.n).isEqualTo(3);
		assertThat(m0.m).isEqualTo(2);
		assertThat(m0.p).isEqualTo(3);
		assertThat(Verifier.isExactNonCubic(m0)).as("inner projection must be exact ⟨3,2,3⟩").isTrue();
	}

	@Test
	public void search_finds_best_projection() throws Exception {
		// Exhaustively project Laderman ⟨3,3,3⟩ down to ⟨2,3,3⟩; best is 16.
		var hit = ProjectionSearch.bestFor(2, 3, 3,
				java.util.List.of(laderman()), /*upperBound*/ 1000, /*maxDelta*/ 1);
		assertThat(hit).isPresent();
		System.out.printf("best projection ⟨2,3,3⟩ rank=%d, keepN=%s%n",
				hit.get().rank(), java.util.Arrays.toString(hit.get().keepN()));
		assertThat(hit.get().rank()).isEqualTo(16L);
		assertThat(Verifier.isExactNonCubic(hit.get().scheme())).isTrue();
	}

	@Test
	public void lineage_node_roundtrips_and_replays() throws Exception {
		// Build a Project lineage node over ⟨3,3,3⟩ PINNED to Laderman 23 (hash ref) —
		// a bare "3x3x3" shape-ref resolves to catalog-best, which is arbitrary among the
		// several rank-23 schemes (the 2026 rename re-ordered them). Pinning keeps the
		// projection deterministic. Render → JSON, parse back, replay. Drop n-index 1 → 16.
		Lineage.Node node = new Lineage.Project(
				new Lineage.Atom("3x3x3@b173cf2"),
				Compose.keepExcept(3, 1), new int[] { 0, 1, 2 }, new int[] { 0, 1, 2 });
		String json = Lineage.toJson(node);
		assertThat(json).contains("\"op\":\"Project\"").contains("\"keepN\":[0,2]");

		File tmp = File.createTempFile("project-lineage", ".json");
		java.nio.file.Files.writeString(tmp.toPath(), "{\"lineage\":" + json + "}");
		var parsed = SchemeIO.readLineage(tmp).orElseThrow();
		assertThat(parsed).isInstanceOf(Lineage.Project.class);

		FieldAwareLookup lk = new FieldAwareLookup("Q");
		var replayer = eu.solven.matmul.search.LineageReplayer.withDefaultPool(lk);
		NonCubicBilinearAlgorithm out = replayer.replay(parsed);
		System.out.printf("replayed projection ⟨%d,%d,%d⟩ r=%d%n", out.n, out.m, out.p, out.r);
		assertThat(out.n).isEqualTo(2);
		assertThat(out.m).isEqualTo(3);
		assertThat(out.p).isEqualTo(3);
		assertThat(out.r).isEqualTo(16);
		assertThat(Verifier.isExactNonCubic(out)).isTrue();
	}

	@Test
	public void dis09_lemma4_stub_is_replayable() throws Exception {
		// The dis09 even cubes (projection's preferred parents) are stored as
		// stubs whose lineage is just Atom("DIS09Lemma4(n=N)"). Replaying that ref
		// must reconstruct the ⟨N,N,N⟩ scheme via PanTrilinearAggregation.build.
		FieldAwareLookup lk = new FieldAwareLookup("Q");
		var replayer = eu.solven.matmul.search.LineageReplayer.withDefaultPool(lk);

		NonCubicBilinearAlgorithm direct =
				replayer.replay(new eu.solven.matmul.catalog.Lineage.Atom("DIS09Lemma4(n=20)"));
		assertThat(direct.n).isEqualTo(20);
		assertThat(direct.m).isEqualTo(20);
		assertThat(direct.p).isEqualTo(20);
		assertThat(direct.r).isEqualTo(4340); // matches dis09_Q-20x20x20_m4340 on disk
		assertThat(Verifier.passesRandomMatmulSpotCheck(direct)).isTrue();

		// And resolving the shape (as projection's resolveParent does) must pick the
		// best-rank dis09 stub and replay it — not a worse derived_recursive file.
		NonCubicBilinearAlgorithm viaShape =
				replayer.replay(new eu.solven.matmul.catalog.Lineage.Atom("20x20x20"));
		assertThat(viaShape.r).isEqualTo(4340);
	}

	@Test
	public void empty_projection_is_identity() throws Exception {
		NonCubicBilinearAlgorithm a = laderman();
		NonCubicBilinearAlgorithm same =
				Compose.project(a, new int[] { 0, 1, 2 }, new int[] { 0, 1, 2 }, new int[] { 0, 1, 2 });
		assertThat(same.r).isEqualTo(a.r); // nothing dropped → same rank
		assertThat(Verifier.isExactNonCubic(same)).isTrue();
	}
}
