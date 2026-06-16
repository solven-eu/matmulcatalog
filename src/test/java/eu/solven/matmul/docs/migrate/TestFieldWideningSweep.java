package eu.solven.matmul.docs.migrate;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import eu.solven.matmul.docs.migrate.FieldWideningSweep.CoefficientProfile;
import eu.solven.matmul.docs.migrate.FieldWideningSweep.FieldTag;
import eu.solven.matmul.docs.migrate.FieldWideningSweep.Proposal;

/**
 * Round-trip tests for {@link FieldWideningSweep}. We never iterate the
 * full catalog here — that's the driver's job. Instead we hand-pick a
 * couple of real scheme files and a few synthetic fixtures so the test
 * runs in well under a second.
 */
public class TestFieldWideningSweep {

	/** Hand-picked alphatensor-Z scheme — coefficients are all in {-1,0,+1} ∪ {2} (Z). */
	private static final Path ALPHATENSOR_Z_333 =
			Path.of("src/main/resources/schemes/known/section3/alphatensor_Z-3x3x3_m23_a110.json");

	/** Strassen — coefficients are all in {-1,0,+1}, F2-compatible after sign-strip. */
	private static final Path STRASSEN_222 =
			Path.of("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json");

	@Test
	public void smallest_denominator_recognises_integers() {
		assertThat(FieldWideningSweep.smallestDenominator(2.0)).isEqualTo(1);
		assertThat(FieldWideningSweep.smallestDenominator(-3.0)).isEqualTo(1);
		assertThat(FieldWideningSweep.smallestDenominator(0.5)).isEqualTo(2);
		assertThat(FieldWideningSweep.smallestDenominator(0.125)).isEqualTo(8);
		// "Irrational" / unbounded should give -1 within the cap.
		assertThat(FieldWideningSweep.smallestDenominator(Math.PI)).isEqualTo(-1);
	}

	@Test
	public void profile_pure_pm1_classifies_as_F2_pm1_or_pm1() {
		CoefficientProfile pm1 = new CoefficientProfile();
		for (double x : new double[] { -1, 0, 1, 0, -1, 1 }) {
			pm1.accumulateReal(x);
		}
		assertThat(pm1.allInPM1).isTrue();
		assertThat(pm1.f2Compatible).isFalse(); // has -1 entries
		assertThat(pm1.toFieldTag()).isEqualTo(FieldTag.PM1);

		CoefficientProfile f2 = new CoefficientProfile();
		for (double x : new double[] { 0, 1, 0, 1, 1, 0 }) {
			f2.accumulateReal(x);
		}
		assertThat(f2.f2Compatible).isTrue();
		assertThat(f2.toFieldTag()).isEqualTo(FieldTag.F2_PM1);
	}

	@Test
	public void profile_integer_outside_pm1_classifies_as_Z() {
		CoefficientProfile pr = new CoefficientProfile();
		for (double x : new double[] { 2, -1, 0, 1, 0 }) {
			pr.accumulateReal(x);
		}
		assertThat(pr.allRealInZ).isTrue();
		assertThat(pr.allInPM1).isFalse();
		assertThat(pr.toFieldTag()).isEqualTo(FieldTag.Z);
	}

	@Test
	public void profile_rational_classifies_as_Q_with_finite_denom() {
		CoefficientProfile pr = new CoefficientProfile();
		for (double x : new double[] { 0.5, -0.25, 1, 0, 0.125 }) {
			pr.accumulateReal(x);
		}
		assertThat(pr.allRealInZ).isFalse();
		assertThat(pr.maxDenom).isEqualTo(8);
		assertThat(pr.toFieldTag()).isEqualTo(FieldTag.Q);
	}

	@Test
	public void profile_zero_imag_complex_is_real_or_stricter() {
		// A "complex" matrix with imag = 0 everywhere and real ∈ {-1, 0, 1}
		// should classify down to PM1 (not stay at C).
		CoefficientProfile pr = new CoefficientProfile();
		double[] reals = { 1, -1, 0, 0, 1 };
		double[] imags = { 0, 0, 0, 0, 0 };
		for (int i = 0; i < reals.length; i++) {
			pr.accumulateComplex(reals[i], imags[i]);
		}
		assertThat(pr.hasComplex).isFalse();
		assertThat(pr.toFieldTag()).isEqualTo(FieldTag.PM1);
	}

	@Test
	public void profile_genuinely_complex_stays_C() {
		CoefficientProfile pr = new CoefficientProfile();
		pr.accumulateComplex(0.5, 0.5);
		pr.accumulateComplex(0, 1);
		assertThat(pr.hasComplex).isTrue();
		assertThat(pr.toFieldTag()).isEqualTo(FieldTag.C);
	}

	@Test
	public void strassen_222_detects_as_PM1() throws Exception {
		Proposal p = FieldWideningSweep.analyze(STRASSEN_222);
		assertThat(p).isNotNull();
		// Strassen entries are exactly {-1, 0, +1} → PM1 (within Z).
		assertThat(p.detected).isEqualTo(FieldTag.PM1);
	}

	@Test
	public void alphatensor_Z_3x3x3_detects_as_Z_not_pm1() throws Exception {
		// This scheme has at least one entry equal to 2 (and -2) — so it's
		// integer but NOT all-in-{-1,0,+1}. Detector must say Z, not PM1.
		Proposal p = FieldWideningSweep.analyze(ALPHATENSOR_Z_333);
		assertThat(p).isNotNull();
		assertThat(p.detected).isEqualTo(FieldTag.Z);
		// Current tag is already Z — so this should NOT be a widening.
		assertThat(p.isWidening()).isFalse();
	}

	@Test
	public void synthetic_Q_tagged_with_integer_entries_proposes_Z(@TempDir Path tmp) throws Exception {
		// Hand-crafted "Q-tagged" 2x2x2 file whose entries are actually all
		// integers. The driver should propose narrowing Q -> Z.
		String json = "{\n"
				+ "  \"field\": \"Q\",\n"
				+ "  \"u\": [[1,0,0,1],[0,0,1,1],[1,0,0,0],[0,0,0,1],[1,1,0,0],[-1,0,1,0],[0,1,0,-1]],\n"
				+ "  \"v\": [[1,0,0,1],[1,0,0,0],[0,1,0,-1],[-1,0,1,0],[0,0,0,1],[1,1,0,0],[0,0,1,1]],\n"
				+ "  \"w\": [[1,0,0,1],[0,1,0,-1],[0,0,1,1],[1,1,0,0],[-1,0,1,0],[0,0,0,1],[1,0,0,0]],\n"
				+ "  \"m\": 7,\n"
				+ "  \"n\": [2,2,2]\n"
				+ "}\n";
		Path file = tmp.resolve("synthetic_Q-2x2x2_m7_a18.json");
		Files.writeString(file, json);

		Proposal p = FieldWideningSweep.analyze(file);
		assertThat(p).isNotNull();
		assertThat(p.current).isEqualTo(FieldTag.Q);
		// All entries are in {-1, 0, +1} so detector goes PM1, which is
		// strictly inside Z — widening proposal triggers.
		assertThat(p.detected).isEqualTo(FieldTag.PM1);
		assertThat(p.isWidening()).isTrue();

		// Sanity: re-verification still passes (we didn't break the scheme).
		assertThat(FieldWideningSweep.reverifies(file)).isTrue();
	}

	@Test
	public void synthetic_complex_with_zero_imag_proposes_R_or_stricter(@TempDir Path tmp) throws Exception {
		// A "complex"-tagged 2x2x2 scheme whose imag parts are all zero —
		// the detector should propose R (or stricter, if integers).
		// Encode Strassen as complex with imag=0; the imag layer is purely
		// padding. Use the dronperminov column-major layout for w.
		String json = "{\n"
				+ "  \"complex\": true,\n"
				+ "  \"field\": \"C\",\n"
				+ "  \"u\": [\n"
				+ "    [[1,0],[0,0],[0,0],[1,0]],\n"
				+ "    [[0,0],[0,0],[1,0],[1,0]],\n"
				+ "    [[1,0],[0,0],[0,0],[0,0]],\n"
				+ "    [[0,0],[0,0],[0,0],[1,0]],\n"
				+ "    [[1,0],[1,0],[0,0],[0,0]],\n"
				+ "    [[-1,0],[0,0],[1,0],[0,0]],\n"
				+ "    [[0,0],[1,0],[0,0],[-1,0]]\n"
				+ "  ],\n"
				+ "  \"v\": [\n"
				+ "    [[1,0],[0,0],[0,0],[1,0]],\n"
				+ "    [[1,0],[0,0],[0,0],[0,0]],\n"
				+ "    [[0,0],[1,0],[0,0],[-1,0]],\n"
				+ "    [[-1,0],[0,0],[1,0],[0,0]],\n"
				+ "    [[0,0],[0,0],[0,0],[1,0]],\n"
				+ "    [[1,0],[1,0],[0,0],[0,0]],\n"
				+ "    [[0,0],[0,0],[1,0],[1,0]]\n"
				+ "  ],\n"
				+ "  \"w\": [\n"
				+ "    [[1,0],[0,0],[0,0],[1,0]],\n"
				+ "    [[0,0],[0,0],[1,0],[-1,0]],\n"
				+ "    [[0,0],[1,0],[0,0],[1,0]],\n"
				+ "    [[1,0],[0,0],[1,0],[0,0]],\n"
				+ "    [[-1,0],[1,0],[0,0],[0,0]],\n"
				+ "    [[0,0],[0,0],[0,0],[1,0]],\n"
				+ "    [[1,0],[0,0],[0,0],[0,0]]\n"
				+ "  ],\n"
				+ "  \"m\": 7,\n"
				+ "  \"n\": [2,2,2]\n"
				+ "}\n";
		Path file = tmp.resolve("synthetic_padded_C-2x2x2_m7.json");
		Files.writeString(file, json);

		Proposal p = FieldWideningSweep.analyze(file);
		assertThat(p).isNotNull();
		assertThat(p.current).isEqualTo(FieldTag.C);
		// All real entries are {-1, 0, +1}, all imag zero → PM1.
		assertThat(p.detected).isEqualTo(FieldTag.PM1);
		assertThat(p.isWidening()).isTrue();
	}

	@Test
	public void field_tag_lattice_is_ordered() {
		// PM1 ⊂ Z ⊂ Q ⊂ R ⊂ C — permissiveness should be monotone.
		assertThat(FieldTag.PM1.permissiveness)
				.isLessThan(FieldTag.Z.permissiveness);
		assertThat(FieldTag.Z.permissiveness)
				.isLessThan(FieldTag.Q.permissiveness);
		assertThat(FieldTag.Q.permissiveness)
				.isLessThan(FieldTag.R.permissiveness);
		assertThat(FieldTag.R.permissiveness)
				.isLessThan(FieldTag.C.permissiveness);
	}
}
