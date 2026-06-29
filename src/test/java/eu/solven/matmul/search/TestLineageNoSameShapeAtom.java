package eu.solven.matmul.search;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

import org.junit.jupiter.api.Test;

import tools.jackson.databind.JsonNode;

import eu.solven.matmul.catalog.SchemeIO;

/**
 * Catalog-integrity guard: a DERIVED scheme's lineage must never contain an {@code Atom}
 * whose ORDERED shape equals the root's shape — that is a degenerate self-derivation
 * ("build ⟨2,2,3⟩ from a ⟨2,2,3⟩ scheme"), e.g. the corrupt
 * {@code ConcatCols(naive-⟨2,2,1⟩, Project(⟨2,2,3⟩→⟨2,2,2⟩))} entry.
 *
 * <p>The existing cycle guards ({@link LineageReplayer} replay-descent, the
 * RecursiveMaterialiser write-time acyclicity check) catch only GRAPH cycles (infinite
 * recursion → StackOverflow) and self-HASH references. A same-shape atom that points at a
 * DIFFERENT concrete scheme terminates fine, so it slips past those guards — yet it is
 * still logically circular (a shape must be built from STRICTLY different shapes). This
 * test encodes that stronger invariant over the committed catalog so the corruption can
 * never silently return. An orientation/transpose atom (same multiset, different order —
 * e.g. ⟨2,3,2⟩ from ⟨2,2,3⟩) is legitimate and intentionally NOT flagged.</p>
 */
public class TestLineageNoSameShapeAtom {

	private static final Path DERIVED = Path.of("src/main/resources/schemes/derived");
	private static final Pattern SHAPE = Pattern.compile("(\\d+)x(\\d+)x(\\d+)");

	@Test
	public void no_derived_lineage_references_an_atom_of_its_own_shape() throws IOException {
		List<String> violations = new ArrayList<>();
		try (Stream<Path> files = Files.walk(DERIVED)) {
			files.filter(p -> p.toString().endsWith(".json")).forEach(p -> {
				Matcher fm = SHAPE.matcher(p.getFileName().toString());
				if (!fm.lookingAt()) return;
				int[] root = { Integer.parseInt(fm.group(1)), Integer.parseInt(fm.group(2)),
						Integer.parseInt(fm.group(3)) };
				JsonNode lineage;
				try {
					lineage = SchemeIO.parseJson(p.toFile()).get("lineage");
				} catch (IOException e) {
					throw new RuntimeException("unreadable " + p, e);
				}
				if (lineage == null || lineage.isNull()) return;
				if (hasSameShapeAtom(lineage, root)) {
					violations.add(p.toString().replace(DERIVED + "/", ""));
				}
			});
		}
		assertThat(violations)
				.as("derived lineages that build a shape from an Atom of the SAME ordered shape "
						+ "(degenerate self-derivation) — purge + re-derive these")
				.isEmpty();
	}

	/** True iff any descendant {@code op:"Atom"} ref decodes to the same ordered shape as {@code root}. */
	private static boolean hasSameShapeAtom(JsonNode node, int[] root) {
		if (node == null) return false;
		if (node.isObject()) {
			if ("Atom".equals(node.path("op").asText())) {
				Matcher m = SHAPE.matcher(node.path("ref").asText());
				if (m.find() && Integer.parseInt(m.group(1)) == root[0]
						&& Integer.parseInt(m.group(2)) == root[1]
						&& Integer.parseInt(m.group(3)) == root[2]) {
					return true;
				}
			}
			for (JsonNode child : node) {
				if (hasSameShapeAtom(child, root)) return true;
			}
		} else if (node.isArray()) {
			for (JsonNode child : node) {
				if (hasSameShapeAtom(child, root)) return true;
			}
		}
		return false;
	}
}
