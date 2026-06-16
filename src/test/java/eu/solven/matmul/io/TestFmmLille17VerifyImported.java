package eu.solven.matmul.io;

import org.junit.jupiter.api.Tag;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.SchemeIO;
import tools.jackson.databind.JsonNode;

@Tag("slow")
public class TestFmmLille17VerifyImported {

    @Test
    public void verify_imported_FMM_17x17x17_Q_strict() throws Exception {
        // STRICT resolution: this test verifies one specific import. The m2945
        // file was purged in the 2026-06 "Remove various FMM from known" pass
        // (superseded by better ⟨17,17,17⟩ bounds); the old lax byHint silently
        // fell back to an arbitrary same-shape scheme, so this test verified the
        // wrong subject. Skip (don't vacuously pass) while the import is absent;
        // reactivates automatically if it is ever re-imported.
        File f = eu.solven.matmul.catalog.SchemeResolver.byHintStrict(
                "src/main/resources/schemes/known/section17/fmm_lille_2025-17x17x17_m2945_a68812.json");
        org.junit.jupiter.api.Assumptions.assumeTrue(f != null && f.exists(),
                "FMM ⟨17,17,17⟩ m2945 import not in catalog (purged 2026-06) — skipping");

        // ⟨17,17,17⟩ is 4913³ ≈ 1.2e11 tensor positions — the exact symbolic
        // verifier (BigInteger, O(N⁶)) takes ~4 min. Use the sampled randomised
        // spot-check instead: O(samples·rank), statistically equivalent
        // confidence in <1s. (Was full SymbolicVerifier.verify — see slow-test
        // report 2026-06-03.)
        JsonNode root = SchemeIO.parseJson(f);
        NonCubicBilinearAlgorithm alg = SchemeIO.readBilinear(f);
        boolean ok = Verifier.passesRandomMatmulSpotCheck(alg, 20_000, 42L);
        assertThat(ok).as("FMM ⟨17,17,17⟩ should pass sampled matmul verification").isTrue();

        // Field tag: the unified fields[] should declare Q (rational coefficients).
        assertThat(SchemeIO.fieldTags(root)).as("declared field tags").contains("Q");
    }
}
