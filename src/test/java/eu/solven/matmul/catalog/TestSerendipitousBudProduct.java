package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;

/**
 * Serendipitous product via bud decomposition (#159). Rediscovers ⟨4,12,8⟩=272
 * (≅ ⟨4,8,12⟩) from a committed ⟨2,3,4⟩=20 base — exact, combinatorial, no
 * import of the result, no border rank / ALS / SAT.
 */
public class TestSerendipitousBudProduct {

	/** Committed AT-Z ⟨2,3,4⟩=20 base — has 4 U-buds (same structure as HK 1971). */
	private static NonCubicBilinearAlgorithm base234() throws Exception {
		return SchemeIO.read(eu.solven.matmul.catalog.SchemeResolver.byHint(
				"src/main/resources/schemes/known/section4/alphatensor_Z-2x3x4_m20_a88.json"));
	}

	@Test
	public void base_234_has_four_u_buds() throws Exception {
		NonCubicBilinearAlgorithm a = base234();
		assertThat(Verifier.isExactNonCubic(a)).isTrue();
		var dec = SerendipitousBudProduct.findUBuds(a);
		System.out.println("U-buds: " + dec.uBuds().size() + ", trivial: " + dec.trivial().length);
		assertThat(dec.uBuds()).hasSize(4);
		assertThat(dec.uBuds()).allMatch(b -> b.length == 2);
		assertThat(dec.trivial()).hasSize(12);
	}

	@Test
	public void bud_summary_is_readable() throws Exception {
		var bs = SerendipitousBudProduct.summarise(base234());
		System.out.println("bud summary: " + bs.summary());
		assertThat(bs.hasBuds()).isTrue();
		assertThat(bs.uBuds()).isEqualTo(4);
		assertThat(bs.summary()).isEqualTo("4×U⟨1,1,2⟩ + 12×⟨1,1,1⟩");
	}

	@Test
	public void rediscover_4x12x8_272() throws Exception {
		NonCubicBilinearAlgorithm a = base234();
		FieldAwareLookup lk = new FieldAwareLookup("Q");
		NonCubicBilinearAlgorithm out = SerendipitousBudProduct.productViaBuds(a, lk, 2, 4, 2);
		System.out.printf("serendipitous ⟨%d,%d,%d⟩ rank=%d (naive Kron=%d)%n",
				out.n, out.m, out.p, out.r, a.r * lk.findWithSource(2, 4, 2).orElseThrow().alg().r);
		assertThat(out.n).isEqualTo(4);
		assertThat(out.m).isEqualTo(12);
		assertThat(out.p).isEqualTo(8);
		assertThat(out.r).isEqualTo(272);
		assertThat(Verifier.isExactNonCubic(out)).isTrue();
	}

	@Test
	public void w_bud_path() throws Exception {
		// cyclicShift rotates U→W, so cyclicShift(base) has 4 W-buds.
		NonCubicBilinearAlgorithm shifted = base234().cyclicShift();
		var dec = SerendipitousBudProduct.findBuds(shifted);
		long w = dec.buds().stream().filter(b -> b.type() == SerendipitousBudProduct.BudType.W).count();
		assertThat(w).isEqualTo(4);
		FieldAwareLookup lk = new FieldAwareLookup("Q");
		NonCubicBilinearAlgorithm out = SerendipitousBudProduct.productViaBuds(shifted, lk, 2, 4, 2);
		System.out.printf("W-bud product ⟨%d,%d,%d⟩ r=%d%n", out.n, out.m, out.p, out.r);
		assertThat(Verifier.isExactNonCubic(out)).isTrue();
	}

	@Test
	public void v_bud_path() throws Exception {
		// cyclicShift² rotates U→V, so cyclicShift²(base) has 4 V-buds.
		NonCubicBilinearAlgorithm shifted = base234().cyclicShift().cyclicShift();
		var dec = SerendipitousBudProduct.findBuds(shifted);
		long v = dec.buds().stream().filter(b -> b.type() == SerendipitousBudProduct.BudType.V).count();
		assertThat(v).isEqualTo(4);
		FieldAwareLookup lk = new FieldAwareLookup("Q");
		NonCubicBilinearAlgorithm out = SerendipitousBudProduct.productViaBuds(shifted, lk, 2, 4, 2);
		System.out.printf("V-bud product ⟨%d,%d,%d⟩ r=%d%n", out.n, out.m, out.p, out.r);
		assertThat(Verifier.isExactNonCubic(out)).isTrue();
	}

	@Test
	public void lineage_node_roundtrips_and_replays() throws Exception {
		// Build the SerendipitousProduct lineage node, render → JSON, parse it back
		// via SchemeIO, and replay via LineageReplayer (base shape-ref → catalog
		// best ⟨2,3,4⟩ → re-bud → product). Validates the full node plumbing.
		eu.solven.matmul.catalog.Lineage.Node node =
				new eu.solven.matmul.catalog.Lineage.SerendipitousProduct(
						new eu.solven.matmul.catalog.Lineage.Atom("2x3x4"), 2, 4, 2);
		String json = eu.solven.matmul.catalog.Lineage.toJson(node);
		assertThat(json).contains("SerendipitousProduct").contains("\"n2\":2");

		File tmp = File.createTempFile("serendip-lineage", ".json");
		java.nio.file.Files.writeString(tmp.toPath(), "{\"lineage\":" + json + "}");
		var parsed = SchemeIO.readLineage(tmp).orElseThrow();
		assertThat(parsed).isInstanceOf(
				eu.solven.matmul.catalog.Lineage.SerendipitousProduct.class);

		FieldAwareLookup lk = new FieldAwareLookup("Q");
		var replayer = eu.solven.matmul.search.LineageReplayer.withDefaultPool(lk);
		NonCubicBilinearAlgorithm out = replayer.replay(parsed);
		System.out.printf("replayed serendipitous ⟨%d,%d,%d⟩ r=%d%n", out.n, out.m, out.p, out.r);
		assertThat(out.n).isEqualTo(4);
		assertThat(out.m).isEqualTo(12);
		assertThat(out.p).isEqualTo(8);
		assertThat(Verifier.isExactNonCubic(out)).isTrue();
		assertThat(out.r).isLessThanOrEqualTo(280); // ≤ naive Kronecker; 272 if base is bud-rich
	}

	@Test
	public void bud_ordering_recovers_8x9x9_430() throws Exception {
		// ⟨8,9,9⟩=430 = (⟨4,3,3⟩=29 − 3)⊗⟨2,3,3⟩ + ⟨6,3,3⟩=40 is, in the bud model,
		// a size-3 V-bud of ⟨4,3,3⟩ (V-bud size k → inner ⟨k·n2,m2,p2⟩ = ⟨6,3,3⟩=40
		// vs 3·⟨2,3,3⟩=45, saving 5 over the 435 plain Kronecker). The default
		// U→V→W greedy masks it into a size-2 bud (→434); only trying all type
		// orderings recovers the size-3 V-bud (→430). Regression guard for that fix.
		FieldAwareLookup lk = new FieldAwareLookup("Q");
		java.util.List<NonCubicBilinearAlgorithm> bases = new java.util.ArrayList<>();
		boolean v3 = false;
		for (java.nio.file.Path path : lk.findFiles(4, 3, 3)) {
			NonCubicBilinearAlgorithm b;
			try { b = SchemeIO.read(path.toFile()); } catch (Exception e) { continue; }
			var or = b.orientAs(4, 3, 3);
			if (or.isEmpty() || or.get().r != 29) continue;
			bases.add(or.get());
			// a size-3 V-bud is reachable under a V-first ordering but not U-first
			var vFirst = SerendipitousBudProduct.findBuds(or.get(), new SerendipitousBudProduct.BudType[] {
					SerendipitousBudProduct.BudType.V, SerendipitousBudProduct.BudType.U,
					SerendipitousBudProduct.BudType.W });
			v3 |= vFirst.buds().stream().anyMatch(
					x -> x.type() == SerendipitousBudProduct.BudType.V && x.terms().length == 3);
		}
		// PRESENCE-FIRST: the 430 is only reachable if a bud-RICH ⟨4,3,3⟩=29 base is on disk.
		// Check that explicitly BEFORE the search, so a catalog that lost the bud-rich base
		// (e.g. a churn that replaced it with an equal-rank, bud-poorer variant) fails HERE
		// with a clear "base missing" message — not later with a cryptic "rank 432 ≠ 430".
		assertThat(bases).as("at least one ⟨4,3,3⟩=29 base must be present in the catalog").isNotEmpty();
		assertThat(v3).as("a ⟨4,3,3⟩=29 base with HIGH bud structure (a size-3 V-bud under V-first "
				+ "ordering) must be present — without it the serendipitous ⟨8,9,9⟩=430 is unreachable")
				.isTrue();

		var hit = SerendipitousSearch.bestFor(8, 9, 9, bases, lk, 999);
		assertThat(hit).isPresent();
		assertThat(hit.get().scheme().r).isEqualTo(430);
		assertThat(Verifier.isExactNonCubic(hit.get().scheme())).isTrue();
	}

	@Test
	public void search_discovers_272() throws Exception {
		NonCubicBilinearAlgorithm a = base234();
		FieldAwareLookup lk = new FieldAwareLookup("Q");
		var hit = SerendipitousSearch.bestFor(4, 12, 8, java.util.List.of(a), lk, 10_000);
		assertThat(hit).isPresent();
		System.out.printf("search hit: ⟨4,12,8⟩ rank=%d via ⟨%d,%d,%d⟩ ⊗ ⟨%d,%d,%d⟩%n",
				hit.get().rank(), hit.get().base().n, hit.get().base().m, hit.get().base().p,
				hit.get().n2(), hit.get().m2(), hit.get().p2());
		assertThat(hit.get().rank()).isEqualTo(272L);
		assertThat(Verifier.isExactNonCubic(hit.get().scheme())).isTrue();
	}
}
