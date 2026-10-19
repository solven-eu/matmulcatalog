package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

/**
 * A lineage DAG with a shared subtree is written with one definition
 * ({@code "id": "L0"}) and {@code {"op":"@ref","id":"L0"}} back-references. Reading it
 * back must resolve every back-reference to the same node instance, whatever the order
 * in which the parser visits a node's fields. Until 2026-10-19 the Recombination case
 * parsed {@code leaves} before {@code base} while the writer renders base first, so a
 * nested recombination's back-reference to the shared base became a literal
 * {@code Atom("@ref?:L0")} — a lossy round trip that a bulk stub rewrite persisted in 40
 * catalog files.
 */
public class TestLineageJsonSharedNodes {

	private static List<String> atomRefs(Lineage.Node n, List<String> out) {
		switch (n) {
			case Lineage.Atom a -> out.add(a.ref());
			case Lineage.RecombinationN r -> {
				atomRefs(r.base(), out);
				for (Lineage.Node l : r.leaves()) atomRefs(l, out);
			}
			default -> throw new IllegalStateException("unexpected node " + n);
		}
		return out;
	}

	@Test
	void a_back_reference_from_a_leaf_to_the_shared_base_survives_the_round_trip(@TempDir Path dir)
			throws Exception {
		Lineage.Node base = new Lineage.Atom("2x2x2@db11bcc");
		Lineage.Node inner = new Lineage.RecombinationN(base, new int[] { 1, 1 }, new int[] { 1, 1 },
				new int[] { 1, 1 }, List.of(new Lineage.Atom("naive-1x1x1")));
		Lineage.Node outer = new Lineage.RecombinationN(base, new int[] { 2, 2 }, new int[] { 2, 2 },
				new int[] { 2, 2 }, List.of(new Lineage.Atom("3x3x3@00460b6"), inner));

		String json = Lineage.toJson(outer);
		assertThat(json).as("the shared base is written once and back-referenced").contains("\"op\":\"@ref\"");

		Path f = dir.resolve("stub.json");
		Files.writeString(f, "{\n  \"n\": [4, 4, 4],\n  \"m\": 49,\n  \"lineage\": " + json + "\n}");
		Lineage.Node parsed = SchemeIO.readLineage(f.toFile()).orElseThrow();

		assertThat(atomRefs(parsed, new ArrayList<>())).noneMatch(r -> r.startsWith("@ref?:"));
		Lineage.RecombinationN po = (Lineage.RecombinationN) parsed;
		Lineage.RecombinationN pi = (Lineage.RecombinationN) po.leaves().get(1);
		assertThat(pi.base()).as("the nested base is the very same node as the outer base").isSameAs(po.base());
		assertThat(((Lineage.Atom) pi.base()).ref()).isEqualTo("2x2x2@db11bcc");
		assertThat(Lineage.toJson(parsed)).as("writing the parsed DAG again is stable").isEqualTo(json);
	}
}
