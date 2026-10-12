package eu.solven.matmul.commutative;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;
import java.util.Optional;

import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.papers.rosowski2019.RosowskiBound;
import eu.solven.matmul.recombination.Recombination;
import lombok.extern.slf4j.Slf4j;
import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;

/**
 * Best-known COMMUTATIVE matmul rank bounds, combining:
 *
 * <ul>
 *   <li><strong>Rosowski 2019/2020</strong> closed-form formulas
 *       ({@link RosowskiBound#commutativeBound}).</li>
 *   <li><strong>DIS09 Table 4</strong> commutative entries
 *       (loaded from {@code references/dis09-cubic-tables.json}).</li>
 *   <li><strong>Trivial</strong> {@code ⟨1,n,m⟩ = nm} (matrix-vector, holds
 *       commutatively).</li>
 *   <li><strong>The on-disk catalog</strong> via
 *       {@link FieldAwareLookup#findRankAllowCommutative} — see below.</li>
 * </ul>
 *
 * <p>For each requested {@code ⟨a,b,c⟩}, returns the min over all axis
 * permutations of all sources.</p>
 *
 * <p><strong>Why the catalog floor</strong>: every non-commutative scheme is
 * also valid over a commutative ring, so {@code R_c ≤ R}. Consulting only the
 * commutative formulas therefore <em>under-reports</em> wherever the catalog
 * beats Rosowski/DIS09 — it did so on <strong>209 of the 5456 shapes</strong>
 * we hold over {@code Q} (smallest: ⟨9,18,18⟩, formulas 1624 vs 1600 on disk;
 * largest gap: ⟨32,32,32⟩, 17392 vs 15079). Since this class is the SOTA
 * resolver for commutative recombination, that inflated baseline made
 * "improvements" look larger than they were. Pass {@code false} to
 * {@link #CommutativeBounds(boolean)} for formula-only behaviour (used to
 * exhibit the gap in tests).</p>
 *
 * <p>The floor consults schemes we actually hold on disk. Cited-only bounds
 * ({@code docs/cited-bounds.json}, no factor matrices) are NOT folded in —
 * a further tightening, deliberately left open.</p>
 *
 * <p>NOT a constructive lookup — only ranks (no factor matrices). Strassen's
 * rank-7 {@code ⟨2,2,2⟩} construction is already optimal even when scalars
 * commute, by Winograd's 1971 lower bound.</p>
 */
@Slf4j
public final class CommutativeBounds {

	private final Map<String, Integer> dis09Cubic;

	/** Whether {@link #bestRank} floors at the on-disk catalog (see class javadoc). */
	private final boolean useCatalogFloor;

	/** Lazily built — the index behind it is statically cached by FieldAwareLookup. */
	private FieldAwareLookup lookup;

	public CommutativeBounds() {
		this(true);
	}

	public CommutativeBounds(boolean useCatalogFloor) {
		this.dis09Cubic = loadDis09CubicCommutative();
		this.useCatalogFloor = useCatalogFloor;
	}

	/** Best-known commutative rank for {@code ⟨a, b, c⟩}, or empty. */
	public Optional<Long> bestRank(int a, int b, int c) {
		if (a < 1 || b < 1 || c < 1) return Optional.empty();
		if (a == 1) return Optional.of((long) b * c);
		if (b == 1) return Optional.of((long) a * c);
		if (c == 1) return Optional.of((long) a * b);

		long best = Long.MAX_VALUE;
		// Try Rosowski over all 6 axis permutations.
		int[][] perms = {
				{ a, b, c }, { a, c, b }, { b, a, c },
				{ b, c, a }, { c, a, b }, { c, b, a }
		};
		for (int[] p : perms) {
			Optional<Long> r = RosowskiBound.commutativeBound(p[0], p[1], p[2]);
			if (r.isPresent() && r.get() < best) best = r.get();
		}
		// DIS09 Table 4 (cubic only).
		if (a == b && b == c) {
			Integer d = dis09Cubic.get(String.valueOf(a));
			if (d != null && d < best) best = d;
		}
		// Catalog floor: R_c ≤ R, so any catalogued scheme — commutative-only OR
		// non-commutative — caps the commutative bound. Field Q per CLAUDE.md's
		// "sweeps default to --field=Q" (its fallback chain admits Z + Q ingredients).
		long fromCatalog = catalogFloor(a, b, c);
		if (fromCatalog < best) best = fromCatalog;
		return best == Long.MAX_VALUE ? Optional.empty() : Optional.of(best);
	}

	/** Best on-disk rank usable commutatively, or {@link Long#MAX_VALUE} when unknown. */
	private long catalogFloor(int a, int b, int c) {
		if (!useCatalogFloor) return Long.MAX_VALUE;
		if (lookup == null) lookup = new FieldAwareLookup(Field.Q);
		int r = lookup.findRankAllowCommutative(a, b, c);
		return r >= Recombination.SotaResolver.UNKNOWN_RANK ? Long.MAX_VALUE : r;
	}

	/** SotaResolver adapter — used by {@link Recombination#recombineWithAllocation}. */
	public Recombination.SotaResolver asSotaResolver() {
		return (a, b, c) -> {
			Optional<Long> r = bestRank(a, b, c);
			if (r.isEmpty()) return Recombination.SotaResolver.UNKNOWN_RANK;
			long v = r.get();
			return v > Recombination.SotaResolver.UNKNOWN_RANK ? Recombination.SotaResolver.UNKNOWN_RANK : (int) v;
		};
	}

	private static Map<String, Integer> loadDis09CubicCommutative() {
		Map<String, Integer> out = new HashMap<>();
		Path src = Path.of("references/dis09-cubic-tables.json");
		try {
			JsonNode root = JsonMapper.builder().build().readTree(Files.readString(src));
			JsonNode arr = root.get("commutative");
			if (arr != null && arr.isArray()) {
				for (JsonNode row : arr) {
					int n = row.get("format").get(0).asInt();
					int rank = row.get("rank").asInt();
					out.put(String.valueOf(n), rank);
				}
			}
		} catch (IOException e) {
			throw new RuntimeException("loading dis09 commutative table", e);
		}
		return out;
	}

	private Integer dis09For(int n) {
		return dis09Cubic.get(String.valueOf(n));
	}

	/** CLI: print cubic commutative bounds n ∈ [2, 30], showing Rosowski / DIS09 / best. */
	public static void main(String[] args) {
		CommutativeBounds cb = new CommutativeBounds();
		log.info(String.format("%6s | %10s | %10s | %10s | source%n", "n", "Rosowski", "DIS09 T4", "best"));
		log.info("-".repeat(60));
		for (int n = 2; n <= 30; n++) {
			Optional<Long> ros = RosowskiBound.commutativeBound(n, n, n);
			Integer d = cb.dis09For(n);
			Optional<Long> best = cb.bestRank(n, n, n);
			String src = "";
			if (best.isPresent() && ros.isPresent() && best.get() == ros.get().longValue()) src = "Rosowski";
			else if (best.isPresent() && d != null && best.get() == d.longValue()) src = "DIS09";
			else if (best.isPresent()) src = "catalog (NC scheme, valid commutatively)";
			log.info(String.format("%6d | %10s | %10s | %10s | %s%n",
					n,
					ros.isPresent() ? ros.get() : "—",
					d == null ? "—" : d,
					best.isPresent() ? best.get() : "—",
					src));
		}
	}
}
