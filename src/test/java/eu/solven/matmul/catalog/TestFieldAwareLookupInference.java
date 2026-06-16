package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.List;
import java.util.Optional;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.algebra.Field;

/**
 * Unit tests for {@link FieldAwareLookup#inferFieldFromLineage} and the
 * lattice-combine helper. Synthetic lineage trees built from {@link Lineage}
 * primitives; leaves resolve against the real on-disk catalog where the
 * filename suffix discipline is well-tested
 * ({@code _Z_}, {@code _Q_}, {@code _R_}, etc.).
 */
public class TestFieldAwareLookupInference {

	private final FieldAwareLookup lookup = new FieldAwareLookup(Field.C);

	@Test
	public void combine_max_in_char0_lattice() {
		// Z < Q < R < C — combining picks the MAX (smallest field containing both).
		assertThat(FieldAwareLookup.combine(Optional.of(Field.Z), Optional.of(Field.Z)))
				.contains(Field.Z);
		assertThat(FieldAwareLookup.combine(Optional.of(Field.Z), Optional.of(Field.Q)))
				.contains(Field.Q);
		assertThat(FieldAwareLookup.combine(Optional.of(Field.Q), Optional.of(Field.R)))
				.contains(Field.R);
		assertThat(FieldAwareLookup.combine(Optional.of(Field.R), Optional.of(Field.C)))
				.contains(Field.C);
		assertThat(FieldAwareLookup.combine(Optional.of(Field.Z), Optional.of(Field.C)))
				.contains(Field.C);
	}

	@Test
	public void combine_f2_dominates_char0() {
		assertThat(FieldAwareLookup.combine(Optional.of(Field.F2), Optional.of(Field.Z)))
				.contains(Field.F2);
		assertThat(FieldAwareLookup.combine(Optional.of(Field.Q), Optional.of(Field.F2)))
				.contains(Field.F2);
	}

	@Test
	public void combine_f2_vs_f3_is_unknown() {
		assertThat(FieldAwareLookup.combine(Optional.of(Field.F2), Optional.of(Field.F3)))
				.isEmpty();
	}

	@Test
	public void combine_unknown_propagates() {
		assertThat(FieldAwareLookup.combine(Optional.empty(), Optional.of(Field.Z)))
				.isEmpty();
		assertThat(FieldAwareLookup.combine(Optional.of(Field.R), Optional.empty()))
				.isEmpty();
	}

	@Test
	public void leaf_named_base_is_Z() {
		// Strassen<2,2,2>=7 is the canonical hand-coded integer base.
		Lineage.Node n = new Lineage.Atom("Strassen<2,2,2>=7");
		FieldAwareLookup.InferredField inf = lookup.inferFieldFromLineage(n);
		assertThat(inf.field()).contains(Field.Z);
		assertThat(inf.unknownLeaves()).isEmpty();
	}

	@Test
	public void leaf_unknown_ref_marks_unknown() {
		Lineage.Node n = new Lineage.Atom("does-not-exist-anywhere_99x99x99");
		FieldAwareLookup.InferredField inf = lookup.inferFieldFromLineage(n);
		assertThat(inf.field()).isEmpty();
		assertThat(inf.unknownLeaves()).contains("does-not-exist-anywhere_99x99x99");
	}

	@Test
	public void kronproduct_of_two_named_bases_is_Z() {
		Lineage.Node n = new Lineage.KronProduct(
				new Lineage.Atom("Strassen<2,2,2>=7"),
				new Lineage.Atom("Strassen<2,2,2>=7"));
		FieldAwareLookup.InferredField inf = lookup.inferFieldFromLineage(n);
		assertThat(inf.field()).contains(Field.Z);
	}

	@Test
	public void kronchain_of_named_bases_is_Z() {
		Lineage.Node n = new Lineage.KronChain(List.of(
				new Lineage.Atom("Strassen<2,2,2>=7"),
				new Lineage.Atom("Laderman<3,3,3>=23"),
				new Lineage.Atom("Strassen<2,2,2>=7")));
		FieldAwareLookup.InferredField inf = lookup.inferFieldFromLineage(n);
		assertThat(inf.field()).contains(Field.Z);
	}

	@Test
	public void transpose_preserves_field() {
		Lineage.Node n = new Lineage.Transpose(
				new Lineage.Atom("Strassen<2,2,2>=7"), "NMP->MPN");
		FieldAwareLookup.InferredField inf = lookup.inferFieldFromLineage(n);
		assertThat(inf.field()).contains(Field.Z);
	}

	@Test
	public void dce_preserves_field() {
		Lineage.Node n = new Lineage.Dce(new Lineage.Atom("Strassen<2,2,2>=7"));
		FieldAwareLookup.InferredField inf = lookup.inferFieldFromLineage(n);
		assertThat(inf.field()).contains(Field.Z);
	}

	@Test
	public void concat_right_combines_leaves() {
		Lineage.Node n = new Lineage.ConcatCols(
				new Lineage.Atom("Strassen<2,2,2>=7"),
				new Lineage.Atom("Strassen<2,2,2>=7"));
		FieldAwareLookup.InferredField inf = lookup.inferFieldFromLineage(n);
		assertThat(inf.field()).contains(Field.Z);
	}

	@Test
	public void recombination_combines_base_and_leaves() {
		Lineage.Node n = new Lineage.RecombinationN(
				new Lineage.Atom("Strassen<2,2,2>=7"),
				new int[] { 1, 1 }, new int[] { 1, 1 }, new int[] { 1, 1 },
				List.of(new Lineage.Atom("Strassen<2,2,2>=7"),
						new Lineage.Atom("Strassen<2,2,2>=7")));
		FieldAwareLookup.InferredField inf = lookup.inferFieldFromLineage(n);
		assertThat(inf.field()).contains(Field.Z);
	}

	@Test
	public void audit_real_catalog_finds_mismatches_or_at_least_processes() throws Exception {
		// Walk a sample of lineage-bearing derived schemes and run lineage
		// inference; assert we processed at least a few files and surface any
		// mismatch. Derived schemes live under derived/; filenames are pure
		// labels (post-2026 rename), so we select by folder + the JSON's own
		// lineage, never a filename prefix.
		java.nio.file.Path root = java.nio.file.Path.of("src/main/resources/schemes/derived");
		if (!java.nio.file.Files.isDirectory(root)) return;  // packaged JAR test → skip
		int processed = 0, mismatches = 0;
		java.util.Map<String, Integer> byPair = new java.util.TreeMap<>();
		java.util.List<String> sample = new java.util.ArrayList<>();
		try (java.util.stream.Stream<java.nio.file.Path> walk = java.nio.file.Files.walk(root)) {
			java.util.Iterator<java.nio.file.Path> it = walk.iterator();
			while (it.hasNext() && processed < 500) {
				java.nio.file.Path p = it.next();
				String name = p.getFileName().toString();
				if (!name.endsWith(".json")) continue;
				try {
					tools.jackson.databind.JsonNode jroot = SchemeIO.parseJson(p.toFile());
					Optional<Lineage.Node> ln = SchemeIO.readLineage(jroot);
					if (ln.isEmpty()) continue;
					processed++;
					FieldAwareLookup.InferredField inf = lookup.inferFieldFromLineage(ln.get());
					// Content-driven baseline (scheme's own fields[]), never the filename.
					Field byContent = FieldAwareLookup.fieldFromContent(p);
					if (inf.field().isPresent() && inf.field().get() != byContent) {
						mismatches++;
						byPair.merge(byContent.tag() + "->" + inf.field().get().tag(), 1, Integer::sum);
						if (sample.size() < 5) sample.add(name + " : " + byContent + "->" + inf.field().get());
					}
				} catch (Exception ignored) { /* skip parse failures */ }
			}
		}
		System.out.println("[audit] processed=" + processed + " mismatches=" + mismatches + " byPair=" + byPair);
		for (String s : sample) System.out.println("[audit] " + s);
		// We just want to confirm the inference ran on real catalog files;
		// finding zero mismatches is fine (means the catalog is consistent).
		assertThat(processed).isGreaterThan(0);
	}

}
