package eu.solven.matmul.docs.verify;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.recombination.BlockSplitSearch;
import eu.solven.matmul.search.LineageReplayer;

import tools.jackson.databind.JsonNode;

/**
 * The lineage-ref audit must agree with the replayer, ref by ref: a ref it calls
 * resolvable replays, a ref it calls dangling / unresolvable throws. Until 2026-10-18
 * the audit kept its own regexes and reported 17 859 "UNRESOLVABLE" refs that all
 * replayed (canonical filename stems, {@code TA_lita(n=N)} — an underscore its
 * parametric pattern rejected — and {@code Strassen<2,2,2>=7 :: CANONICAL} labels),
 * plus four "dangling" pins that were simply F₂ files invisible to its C lookup.
 */
public class TestAuditLineageRefs {

	private static final Pattern CANONICAL = Pattern.compile("^\\d+x\\d+x\\d+-r\\d+-.+-([0-9a-f]{7})\\.json$");

	private static FieldAwareLookup q;
	private static LineageReplayer replayer;

	@BeforeAll
	static void setup() {
		q = new FieldAwareLookup("Q");
		replayer = LineageReplayer.withDefaultPool(q);
	}

	/** A canonical {@code {shape}-r{rank}-{note}-{hash7}} file at the shape (explicit, char-0). */
	private static Path canonicalFileAt(FieldAwareLookup lookup, int n, int m, int p) {
		return lookup.findFiles(n, m, p).stream().filter(f -> {
			if (!CANONICAL.matcher(f.getFileName().toString()).matches()) return false;
			try {
				JsonNode r = SchemeIO.parseJson(f.toFile());
				return !SchemeIO.isStub(r) && !SchemeIO.isComplex(r) && !SchemeIO.isNonBilinear(r);
			} catch (Exception e) {
				return false;
			}
		}).findFirst().orElseThrow();
	}

	private static String hash7Of(Path f) {
		Matcher m = CANONICAL.matcher(f.getFileName().toString());
		m.matches();
		return m.group(1);
	}

	@Test
	void the_audit_agrees_with_the_replayer_on_every_kind_of_ref() {
		Path strassen = canonicalFileAt(q, 2, 2, 2);
		String stem337 = canonicalFileAt(q, 3, 3, 7).getFileName().toString().replaceAll("\\.json$", "");
		// A pool label whose shape has no unit axis (an AxisSplit<2,1,1>=2 label would be NAIVE).
		String named = BlockSplitSearch.defaultPool().stream().map(BlockSplitSearch.NamedBase::label)
				.filter(l -> l.contains("<2,2,2>")).findFirst().orElseThrow();

		Map<String, String> expected = new LinkedHashMap<>();
		expected.put(stem337, "LEGACY_STEM"); // a pre-2026-10-18 canonical stem (419 distinct ones in the catalog)
		expected.put("naive-1x1x1", "NAIVE");
		expected.put("1x4x4", "NAIVE"); // a width-1 axis is built in code, no file needed
		expected.put("2x3x3", "BARE");
		expected.put("TA_lita_npz(n=13,r=1379)", "PARAMETRIC"); // underscores in the constructor name
		expected.put(named + " :: CANONICAL", "LEGACY_STEM"); // 2 347 files carry such pool labels
		expected.put(named, "NAMED");
		expected.put("2x2x2@" + hash7Of(strassen), "PINNED_OK");
		expected.put("2x2x2@0000000000000000", "PINNED_DANGLING");
		expected.put("97x98x99", "UNRESOLVABLE");
		expected.put("no-shape-in-this-ref", "UNRESOLVABLE");

		for (Map.Entry<String, String> e : expected.entrySet()) {
			String ref = e.getKey();
			String cls = AuditLineageRefs.classify(ref, replayer);
			assertThat(cls).as(ref).isEqualTo(e.getValue());
			boolean resolvable = !cls.equals("PINNED_DANGLING") && !cls.equals("UNRESOLVABLE");
			if (resolvable) {
				assertThat(replayer.replay(new Lineage.Atom(ref))).as("the replayer resolves " + ref).isNotNull();
			} else {
				assertThatThrownBy(() -> replayer.replay(new Lineage.Atom(ref))).as("the replayer rejects " + ref)
						.isInstanceOf(IllegalStateException.class);
			}
		}
	}

	@Test
	void a_pin_to_an_f2_only_file_is_dangling_for_a_char0_lookup_but_fine_under_the_files_own_field() {
		FieldAwareLookup f2 = new FieldAwareLookup("F2");
		// An F₂-only explicit ⟨4,4,4⟩ (AlphaTensor's 47): no char-0 tag at all.
		Path f2only = f2.findFiles(4, 4, 4).stream().filter(f -> {
			if (!CANONICAL.matcher(f.getFileName().toString()).matches()) return false;
			try {
				JsonNode r = SchemeIO.parseJson(f.toFile());
				List<String> tags = SchemeIO.fieldTags(r);
				return !SchemeIO.isStub(r) && tags.contains("F2") && !tags.contains("Z") && !tags.contains("Q")
						&& !tags.contains("R") && !tags.contains("C");
			} catch (Exception e) {
				return false;
			}
		}).findFirst().orElseThrow();
		String ref = "4x4x4@" + hash7Of(f2only);

		assertThat(AuditLineageRefs.classify(ref, replayer)).as("under Q").isEqualTo("PINNED_DANGLING");
		assertThat(AuditLineageRefs.classify(ref, AuditLineageRefs.replayerFor(Field.F2))).as("under F2")
				.isEqualTo("PINNED_OK");
		// The audit tries the file's own fields: both lookups together resolve it.
		assertThat(AuditLineageRefs.classify(ref, List.of(replayer, AuditLineageRefs.replayerFor(Field.F2))))
				.isEqualTo("PINNED_OK");

		assertThat(AuditLineageRefs.fieldsToTry(List.of("F2"))).containsExactly(Field.F2);
		assertThat(AuditLineageRefs.fieldsToTry(List.of("F2", "F3", "Z", "Q", "R", "C")))
				.containsExactly(Field.C, Field.F2, Field.F3);
		assertThat(AuditLineageRefs.fieldsToTry(List.of("Q", "R", "C"))).containsExactly(Field.C);
		assertThat(AuditLineageRefs.fieldsToTry(List.of())).as("no fields[] = data bug: try everything")
				.containsExactly(Field.C, Field.F2, Field.F3);
	}
}
