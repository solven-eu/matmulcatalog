package eu.solven.matmul.docs.explore;

import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import lombok.extern.slf4j.Slf4j;

/**
 * One-shot σ-table for the ⟨6,8,9⟩ demonstration: what each bud type/size is
 * worth for base ⟨2,4,3⟩ ⊗ inner ⟨3,2,3⟩ — σ(type,k) = k·R(inner) − R(enlarged).
 */
@Slf4j
public class ProbeSerendipitySigmaTable {

	public static void main(String[] args) {
		FieldAwareLookup q = new FieldAwareLookup(Field.Q);
		int inner = q.findRank(3, 2, 3);
		log.info("R⟨3,2,3⟩ = {}  (plain Kron ⟨6,8,9⟩ = 20·{} = {})", inner, inner, 20 * inner);
		for (int k = 2; k <= 4; k++) {
			log.info("k={}: σ_U = {}·{} − R⟨3,2,{}⟩={} → {} | σ_V = −R⟨{},2,3⟩={} → {}"
					+ " | σ_W = −R⟨3,{},3⟩={} → {}",
					k, k, inner,
					3 * k, q.findRank(3, 2, 3 * k), k * inner - q.findRank(3, 2, 3 * k),
					3 * k, q.findRank(3 * k, 2, 3), k * inner - q.findRank(3 * k, 2, 3),
					2 * k, q.findRank(3, 2 * k, 3), k * inner - q.findRank(3, 2 * k, 3));
		}
	}
}
