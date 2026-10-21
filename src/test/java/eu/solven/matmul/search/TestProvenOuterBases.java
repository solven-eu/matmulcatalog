package eu.solven.matmul.search;

import static org.assertj.core.api.Assertions.assertThat;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.recombination.BlockSplitSearch;

/**
 * Pool = cap + proven-useful list (2026-10-21). The default recombination pool caps
 * non-cubic bases at dim 5, yet ~500 committed stubs stood on dim-6/7/8 outer bases
 * (⟨3,4,6⟩, ⟨2,5,6⟩, ⟨3,4,7⟩, ⟨2,6,6⟩, ⟨3,5,8⟩, …) that only one-off {@code --base=} runs
 * or older, wider sweeps had offered. The closure sweep could not re-derive what the
 * catalog held, so its "no improvement" fixpoint over those shapes said nothing — a
 * silent failure. Now every file the catalog's lineages name as a recombination base is
 * admitted whatever its dimension, and this guard fails the moment a committed lineage
 * stands on a base the default pool does not offer.
 */
public class TestProvenOuterBases {

	private static final String DPS_347 = "known/section7/3x4x7-r63-dumas_pernet_sedoglavic_2025-ac0e1ad.json";
	private static final String DPS_347_HASH7 = "ac0e1ad";

	/** The hashes every pool entry pins (through its origin lineage), full length. */
	private static Set<String> poolHashes(List<BlockSplitSearch.NamedBase> pool) {
		Set<String> out = new HashSet<>();
		for (BlockSplitSearch.NamedBase nb : pool) {
			if (nb.originLineage() == null) continue;
			Lineage.Node leaf = ProvenOuterBases.unwrap(nb.originLineage());
			if (leaf instanceof Lineage.Atom a && a.ref().indexOf('@') > 0) {
				out.add(a.ref().substring(a.ref().indexOf('@') + 1));
			}
		}
		return out;
	}

	private static boolean representedByTwin(ProvenOuterBases.Proven p, List<List<BlockSplitSearch.NamedBase>> pools,
			List<FieldAwareLookup> lookups) {
		Set<String> keys = new HashSet<>();
		for (List<BlockSplitSearch.NamedBase> pool : pools) {
			for (BlockSplitSearch.NamedBase nb : pool) keys.add(BlockSplitSearch.poolContentKey(nb.base()));
		}
		return lookups.stream().anyMatch(lk -> lk
				.findByHash(p.shape()[0], p.shape()[1], p.shape()[2], p.hashPrefix())
				.map(ws -> keys.contains(BlockSplitSearch.poolContentKey(ws.alg())))
				.orElse(false));
	}

	private static boolean contains(Set<String> hashes, String prefix) {
		for (String h : hashes) {
			if (h.startsWith(prefix) || prefix.startsWith(h)) return true;
		}
		return false;
	}

	private static void stub(Path root, String rel, String baseRef) throws Exception {
		Path f = root.resolve(rel);
		Files.createDirectories(f.getParent());
		Files.writeString(f, "{\n"
				+ "  \"lineage_compact\": \"R[3x4x7; 4,4,4 | 4,4,4,4 | 4,4,4,4,4,4,3]\",\n"
				+ "  \"lineage\": {\n"
				+ "    \"op\": \"Recombination\",\n"
				+ "    \"base\": {\"op\": \"Atom\", \"ref\": \"" + baseRef + "\"},\n"
				+ "    \"allocA\": [4, 4, 4],\n"
				+ "    \"allocB\": [4, 4, 4, 4],\n"
				+ "    \"allocC\": [4, 4, 4, 4, 4, 4, 3],\n"
				+ "    \"leaves\": [{\"op\": \"Atom\", \"ref\": \"4x4x4\"}]\n"
				+ "  },\n"
				+ "  \"n\": [12, 16, 27],\n"
				+ "  \"m\": 2964,\n"
				+ "  \"fields\": [\"Q\", \"R\", \"C\"],\n"
				+ "  \"scheme_type\": \"stub\"\n"
				+ "}\n", StandardCharsets.UTF_8);
	}

	private static void leaf(Path root) throws Exception {
		Path dst = root.resolve(DPS_347);
		Files.createDirectories(dst.getParent());
		Files.copy(BlockSplitSearch.DEFAULT_SCHEMES_ROOT.resolve(DPS_347), dst);
	}

	/** A dim-7 leaf is outside the cap-5 default pool … */
	@Test
	public void a_dim7_leaf_alone_stays_outside_the_capped_pool(@TempDir Path root) throws Exception {
		leaf(root);
		List<BlockSplitSearch.NamedBase> pool =
				BlockSplitSearch.buildPool(RecombinationPoolConfig.defaultLarge(), "Q", root);
		assertThat(contains(poolHashes(pool), DPS_347_HASH7)).as("no lineage stands on it: capped out").isFalse();
		assertThat(pool).noneMatch(BlockSplitSearch::isProven);
	}

	/** … until a committed lineage pins it as an outer base: then it is admitted above the cap. */
	@Test
	public void a_stub_pinned_to_a_dim7_base_admits_that_file(@TempDir Path root) throws Exception {
		leaf(root);
		stub(root, "derived/section27/12x16x27-r2964-derived-fd854e5.json", "3x4x7@" + DPS_347_HASH7);
		ProvenOuterBases proven = ProvenOuterBases.of(root);
		assertThat(proven.all()).hasSize(1);
		assertThat(proven.all().get(0).shapeKey()).isEqualTo("3x4x7");
		assertThat(proven.unresolved()).isEmpty();

		List<BlockSplitSearch.NamedBase> pool =
				BlockSplitSearch.buildPool(RecombinationPoolConfig.defaultLarge(), "Q", root);
		assertThat(contains(poolHashes(pool), DPS_347_HASH7)).as("the proven dim-7 base is in the pool").isTrue();
		List<BlockSplitSearch.NamedBase> provenEntries = pool.stream().filter(BlockSplitSearch::isProven).toList();
		assertThat(provenEntries).isNotEmpty();
		assertThat(provenEntries).allMatch(nb -> nb.label().contains("3x4x7-r63-dumas_pernet_sedoglavic_2025"));
		// The cap still binds everything else: no other dim-7 base, no ⟨3,4,7⟩ twin.
		assertThat(pool.stream().filter(nb -> !BlockSplitSearch.isProven(nb))
				.noneMatch(nb -> Math.max(nb.base().n, Math.max(nb.base().m, nb.base().p)) > 5
						&& !(nb.base().n == nb.base().m && nb.base().m == nb.base().p))).isTrue();
		// An explicit narrowing is respected: cubicOnly drops the rectangular proven base.
		RecombinationPoolConfig cubicOnly = new RecombinationPoolConfig(true, false,
				eu.solven.matmul.SymmetryTransforms.InternalOrbitMode.AXIS_FLIP, 5, 8, false,
				RecombinationPoolConfig.UNBOUNDED_IMBALANCE, RecombinationPoolConfig.UNBOUNDED_COMBINATIONS,
				RecombinationPoolConfig.NO_PADDING, 8);
		assertThat(BlockSplitSearch.buildPool(cubicOnly, "Q", root)).noneMatch(BlockSplitSearch::isProven);
	}

	/**
	 * A pre-pinning stub names its base by the pool label of the day,
	 * {@code ext[<stem>.json] / <shape>=rank :: CANONICAL}; the stem's hash7 names the file.
	 */
	@Test
	public void a_legacy_ext_label_names_its_stem_file(@TempDir Path root) throws Exception {
		leaf(root);
		stub(root, "derived/section27/12x16x27-r2964-derived-0000000.json",
				"ext[3x4x7-r63-dumas_pernet_sedoglavic_2025-ac0e1ad.json] / <3,4,7>=63 :: CANONICAL");
		ProvenOuterBases proven = ProvenOuterBases.of(root);
		assertThat(proven.admits("ac0e1ad89643fa7e0bacda21c38e173a36c63334ab067ca1cce61a5b1b794914")).isTrue();
		assertThat(proven.unresolved()).isEmpty();
	}

	/**
	 * THE guard: every recombination base the committed catalog names is offered by the
	 * default pool, so the closure can always re-derive what it holds. Char-0 lineages are
	 * checked against the C pool (its chain covers Z/Q/R), F₂/F₃ ones against their own.
	 */
	@Test
	public void every_outer_base_of_a_committed_lineage_is_in_the_default_pool() {
		Path root = BlockSplitSearch.DEFAULT_SCHEMES_ROOT;
		ProvenOuterBases proven = ProvenOuterBases.of(root);
		assertThat(proven.lineageFiles()).as("the whole catalog was scanned").isGreaterThan(4000);
		assertThat(proven.all().size()).isGreaterThan(20);
		// The historical reason for this guard: dim-6/7/8 bases that one-off runs introduced.
		assertThat(proven.aboveCap(5, 8)).extracting(ProvenOuterBases.Proven::shapeKey)
				.contains("3x4x7", "3x4x6", "2x5x6");

		RecombinationPoolConfig cfg = RecombinationPoolConfig.defaultLarge();
		List<List<BlockSplitSearch.NamedBase>> pools = new ArrayList<>();
		List<FieldAwareLookup> lookups = new ArrayList<>();
		pools.add(BlockSplitSearch.buildPool(cfg, "C", root));
		lookups.add(new FieldAwareLookup(Field.C, root));
		Set<String> hashes = poolHashes(pools.get(0));
		List<ProvenOuterBases.Proven> missing = new ArrayList<>();
		for (ProvenOuterBases.Proven p : proven.all()) {
			if (!contains(hashes, p.hashPrefix())) missing.add(p);
		}
		// The pool dedups bases by product-support multiset (two such files tile every target
		// identically), so a proven file may be represented by a support twin — e.g.
		// AlphaTensor's ⟨3,5,5⟩=58 by the derived ⟨5,3,5⟩ twin that also carries a pin.
		if (!missing.isEmpty()) missing.removeIf(p -> representedByTwin(p, pools, lookups));
		if (!missing.isEmpty()) {
			// F₂ / F₃ universes: a base valid only there is absent from the C pool by design.
			for (Field f : new Field[] { Field.F2, Field.F3 }) {
				List<BlockSplitSearch.NamedBase> pool = BlockSplitSearch.buildPool(cfg, f.tag(), root);
				pools.add(pool);
				lookups.add(new FieldAwareLookup(f, root));
				Set<String> h = poolHashes(pool);
				missing.removeIf(p -> contains(h, p.hashPrefix()) || representedByTwin(p, pools, lookups));
			}
		}
		List<String> report = new ArrayList<>();
		for (ProvenOuterBases.Proven p : missing) {
			report.add("⟨" + p.shapeKey() + "⟩@" + p.hashPrefix().substring(0, 7) + " (" + p.uses()
					+ " lineage(s), e.g. " + p.example().getFileName() + ")");
		}
		assertThat(report).as("recombination bases the catalog stands on but the default pool does not offer").isEmpty();
	}

	/** A base ref that resolves to nothing is a dangling lineage, which the audit owns; none may hide here. */
	@Test
	public void no_committed_outer_base_is_unresolvable() {
		ProvenOuterBases proven = ProvenOuterBases.of(BlockSplitSearch.DEFAULT_SCHEMES_ROOT);
		List<String> report = new ArrayList<>();
		for (ProvenOuterBases.Unresolved u : proven.unresolved()) {
			report.add(u.ref() + "   [in " + u.file().getFileName() + "]");
		}
		assertThat(report).isEmpty();
	}
}
