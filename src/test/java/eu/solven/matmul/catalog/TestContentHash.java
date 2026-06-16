package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;

/**
 * Content hash = a precise, representation-invariant scheme reference. The
 * load-bearing case: two ⟨3,3,3⟩=23 schemes of equal rank/additions but
 * different content (Laderman bud=0 vs Perminov bud-rich) must hash differently,
 * and {@code findByHash} must resolve to the exact one — fixing the shape-ref
 * ambiguity that forced ⟨6,20,15⟩ to be stored as full matrices.
 */
public class TestContentHash {

	private static final String LADERMAN =
			"src/main/resources/schemes/known/section3/laderman_1976-3x3x3_m23_a98_b0.json";
	private static final String PERMINOV =
			"src/main/resources/schemes/known/section3/perminov_c88_ZT-3x3x3_m23_a88_b6.json";

	private static NonCubicBilinearAlgorithm read(String path) throws Exception {
		return SchemeIO.read(eu.solven.matmul.catalog.SchemeResolver.byHint(path));
	}

	@Test
	public void hash_is_deterministic() throws Exception {
		assertThat(SchemeIO.contentHash(read(LADERMAN)))
				.isEqualTo(SchemeIO.contentHash(read(LADERMAN)));
	}

	@Test
	public void same_shape_rank_adds_but_different_content_hash_differently() throws Exception {
		// Both are ⟨3,3,3⟩=23; the canonical (shape,rank) key would collide.
		NonCubicBilinearAlgorithm lad = read(LADERMAN);
		NonCubicBilinearAlgorithm per = read(PERMINOV);
		assertThat(lad.r).isEqualTo(per.r).isEqualTo(23);
		assertThat(SchemeIO.contentHash(lad)).isNotEqualTo(SchemeIO.contentHash(per));
	}

	@Test
	public void findByHash_resolves_the_exact_scheme() throws Exception {
		FieldAwareLookup lookup = new FieldAwareLookup(Field.Q);
		String perHash = SchemeIO.contentHash(read(PERMINOV));

		var hit = lookup.findByHash(3, 3, 3, perHash);
		assertThat(hit).isPresent();
		// The resolved scheme must hash to the requested hash (i.e. the Perminov
		// one), NOT some other ⟨3,3,3⟩=23 the shape-ref would have picked.
		assertThat(SchemeIO.contentHash(hit.get().alg())).isEqualTo(perHash);
	}

	@Test
	public void short_prefix_also_resolves() throws Exception {
		FieldAwareLookup lookup = new FieldAwareLookup(Field.Q);
		String shortH = SchemeIO.shortHash(read(PERMINOV));  // 12-hex prefix
		var hit = lookup.findByHash(3, 3, 3, shortH);
		assertThat(hit).isPresent();
		assertThat(SchemeIO.contentHash(hit.get().alg())).startsWith(shortH);
	}
}
