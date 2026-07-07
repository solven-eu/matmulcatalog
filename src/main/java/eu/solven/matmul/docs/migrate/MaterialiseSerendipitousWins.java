package eu.solven.matmul.docs.migrate;

import java.nio.file.Path;
import java.util.List;
import java.util.Optional;
import java.util.Set;

import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.recombination.BlockSplitSearch;
import eu.solven.matmul.search.RecursiveClosureSota;
import eu.solven.matmul.search.RecursiveMaterialiser;
import lombok.extern.slf4j.Slf4j;

/**
 * One-shot: persist serendipitous-bud-product wins that the normal sweep cannot
 * reach because the target carries a dense ≤16 import (the sweep skips it as
 * "direct", never composing). Runs an IMPROVE-mode, serendipitous-only
 * materialiser with writes enabled, so {@code compose()} → {@code trySerendipitous}
 * actually derives the win and persists it when it strictly beats the catalog.
 *
 * <p>Canonical case ⟨8,9,9⟩=430 = (⟨4,3,3⟩=29−3)⊗⟨2,3,3⟩+⟨6,3,3⟩=40: the disk best
 * was a perminov-432 import, and the σ-base-selection fix (2026-06-23,
 * {@code budRichestAt}→{@code budBasesAt}) makes {@code trySerendipitous} feed the
 * size-3 V-bud base to {@code SerendipitousSearch.bestFor} → 430. This driver lands
 * that on disk (and cascades: ⟨17,17,17⟩'s recombination uses an ⟨8,9,9⟩ leaf).</p>
 *
 * <pre>mvn -q -o exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.MaterialiseSerendipitousWins</pre>
 */
@Slf4j
public final class MaterialiseSerendipitousWins {
	private MaterialiseSerendipitousWins() {}

	/** Serendipitous targets to (re-)derive and persist if strictly better than the catalog. */
	private static final int[][] TARGETS = {
			{ 8, 9, 9 },
	};

	/**
	 * Parses {@code --shapes=NxMxP,NxMxP,…} into a target list; without the flag,
	 * falls back to the hardcoded {@link #TARGETS} (the historical ⟨8,9,9⟩ case).
	 * Lets a SerendipitousSweep run feed its win list straight into persistence.
	 */
	private static int[][] targetsFrom(String[] args) {
		for (String a : args) {
			if (a.startsWith("--shapes=")) {
				String[] shapes = a.substring("--shapes=".length()).split(",");
				int[][] out = new int[shapes.length][];
				for (int i = 0; i < shapes.length; i++) {
					String[] d = shapes[i].trim().toLowerCase().split("x");
					if (d.length != 3) throw new IllegalArgumentException("bad shape: " + shapes[i]);
					out[i] = new int[] { Integer.parseInt(d[0]), Integer.parseInt(d[1]), Integer.parseInt(d[2]) };
				}
				return out;
			}
		}
		return TARGETS;
	}

	public static void main(String[] args) {
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		List<BlockSplitSearch.NamedBase> pool = BlockSplitSearch.defaultPool();
		RecursiveClosureSota sota = new RecursiveClosureSota(lookup, pool, true, true);
		Path root = Path.of("src/main/resources/schemes");
		// writeRoot=root, writeNewSchemes=true, balancedOnly=false, improveExisting=true, deriveBest=false.
		RecursiveMaterialiser mat =
				new RecursiveMaterialiser(lookup, pool, sota, root, true, false, true, false);
		// Serendipitous strategy only: it's the mechanism these targets need, and skipping the
		// recombination B&B keeps the run fast and the lineage clean.
		mat.setStrategies(Set.of(RecursiveMaterialiser.STRAT_SERENDIPITOUS));

		for (int[] t : targetsFrom(args)) {
			int n = t[0], m = t[1], p = t[2];
			int before = lookup.findRank(n, m, p);
			Optional<RecursiveMaterialiser.Result> r = mat.materialise(n, m, p);
			if (r.isEmpty()) {
				log.warn("⟨{},{},{}⟩ — no serendipitous result (catalog best stays r={})", n, m, p, before);
				continue;
			}
			int got = r.get().alg().r;
			log.info("⟨{},{},{}⟩ — catalog was r={}, serendipitous derived r={} ({})",
					n, m, p, before, got, got < before ? "WIN, persisted" : "no improvement, not written");
		}
		log.info("Done. Re-run GenerateCatalogManifest to refresh docs/catalog.json.");
	}
}
