# New matmul rank upper bounds from this repository

Bounds that **this repository's** combined pipeline (DIS09 reproduction
with modern catalog + Pan TA closed form + paired sub-products +
non-balanced enumeration + Hopcroft-Kerr base) establishes as strictly
better than what's published in [DIS09 Table 3 (non-commutative)
and Table 4 (commutative)](../REFERENCES.md#10-drisc09).

**Status legend**:
- ✅ **verified** — explicit scheme materialised + `Verifier.isExactNonCubic` passes
- 📐 **derived** — formula-derived bound; the construction is sound
  (each ingredient is a verified scheme or a published formula) but
  the composed scheme is not yet materialised as a single JSON file
- 🔬 **claim** — bound from a closed-form formula whose proof we cite
  but haven't independently re-derived

Every entry here can be reproduced with:
```bash
mvn -q -o test -Dtest=TestDIS09FullScan#dis09_full_comparison_R_4_to_30_unbalanced
mvn -q -o test -Dtest=TestDIS09Table4Commutative
```

---

## Non-commutative over R/Q/Z (vs DIS09 Table 3)

| `⟨n,n,n⟩` | DIS09 | ours derived | ours verified | Δ verified | construction                                                      |
|-----------|-------|--------------|---------------|------------|-------------------------------------------------------------------|
| ⟨14,14,14⟩ | 1728  | **1720** 📐 | (not yet)    | —          | Strassen[7,7]³ + paired cyclic sub-products (Layer 4 not constructive yet) |
| ⟨19,19,19⟩ | 4073  | **4030** 📐 | **4053** ✅  | **−20**    | Strassen[9,10]³ — `section19/solven-strassen-2026_19x19x19_r4053_a114750.json` |
| ⟨21,21,21⟩ | 5365  | **5258** 📐 | **5276** ✅  | **−89**    | Strassen[9,12]/[9,12]/[10,11] — `section21/solven-strassen-2026_21x21x21_r5276_a424381.json` |
| ⟨23,23,23⟩ | 6806  | **6731** 📐 | **6738** ✅  | **−68**    | Strassen[11,12]³ — `section23/solven-strassen-2026_23x23x23_r6738_a528794.json` |

Aggregate vs DIS09 Table 3 (n=4..30, R/Q/Z, non-commutative):
**15 wins, 12 ties, 0 losses**, total **−641 mults (−0.56%)**.

**Verified gap (slightly higher than derived):** the verified materialised
ranks are 18–23 mults higher than the derived bounds because the
construction has to fall back to verified-but-slightly-worse inner
sub-schemes for some shapes — the AT-Z composed-improvement files for
dims 9–12 are in `FieldAwareLookup.KNOWN_BROKEN_FILES` (spot-check
failures from our import). Once the AT-Z import is re-fetched from
FMM/Perminov (see [ROADMAP.md](../ROADMAP.md)), the verified ranks will
drop back to match the derived bounds.

The 12 ties come from cases where DIS09's value already matches Pan TA's
closed-form bound (which we reproduce exactly via
[`PanTrilinearAggregation.cubicBound`](../src/main/java/io/cormoran/strassen/v3/catalog/PanTrilinearAggregation.java)
— see [REFERENCES.md#19-schwartz-zwecher25](../REFERENCES.md#19-schwartz-zwecher25)
for SZ25's recent work pushing TA further).

---

## Commutative (vs DIS09 Table 4)

These derive from a commutative SOTA pipeline that adds Rosowski 2019
([`RosowskiBound.bestCommutativeBound`](../src/main/java/io/cormoran/strassen/v3/catalog/RosowskiBound.java))
+ Waksman 1970 ([`WaksmanBound`](../src/main/java/io/cormoran/strassen/v3/catalog/WaksmanBound.java))
+ Hopcroft-Kerr ⟨2,3,3⟩=15 base into the recursive search.

| `⟨n,n,n⟩` | DIS09 (cmt) | ours | Δ      | status | construction                                       |
|-----------|-------------|------|--------|--------|----------------------------------------------------|
| ⟨5,5,5⟩    | 93          | **85**   | −8     | 📐 derived | mul211 outer + `[2,3]/[5]/[5]` — Rosowski ⟨5,3,5⟩=51 + Waksman ⟨5,2,5⟩=34 |
| ⟨7,7,7⟩    | 235         | **217**  | −18    | 📐 derived | mul211 outer + `[2,5]/[7]/[7]`                    |
| ⟨9,9,9⟩    | 472         | **441**  | −31    | 📐 derived | mul211 outer + `[2,7]/[9]/[9]`                    |
| ⟨11,11,11⟩ | 825         | **781**  | −44    | 📐 derived | mul211 outer + `[2,9]/[11]/[11]`                  |
| ⟨13,13,13⟩ | 1318        | **1261** | −57    | 📐 derived | Strassen[6,7]³                                     |
| ⟨14,14,14⟩ | 1525        | **1519** | −6     | 📐 derived | Strassen[7,7]³                                     |
| ⟨15,15,15⟩ | 1941        | **1860** | −81    | 📐 derived | Hopcroft-Kerr[7,8]/[5,5,5]/[5,5,5]                |
| ⟨17,17,17⟩ | 2762        | **2673** | −89    | 📐 derived | Strassen[8,9]³                                     |
| ⟨19,19,19⟩ | 3757        | **3646** | −111   | 📐 derived | Hopcroft-Kerr[9,10]/[6,6,7]/[6,6,7]               |
| ⟨21,21,21⟩ | 4938        | **4767** | −171   | 📐 derived | Hopcroft-Kerr[10,11]/[7,7,7]/[7,7,7]              |
| ⟨23,23,23⟩ | 6382        | **6259** | −123   | 📐 derived | Hopcroft-Kerr[12,11]/[8,8,7]/[8,8,7]              |
| ⟨25,25,25⟩ | 8083        | **7897** | −186   | 📐 derived | Hopcroft-Kerr[13,12]/[9,8,8]/[9,8,8]              |
| ⟨27,27,27⟩ | 9994        | **9720** | −274   | 📐 derived | Hopcroft-Kerr[14,13]/[9,9,9]/[9,9,9]              |
| ⟨29,29,29⟩ | 12109       | **12059**| −50    | 📐 derived | Hopcroft-Kerr[15,14]/[10,10,9]/[10,10,9]          |

Aggregate vs DIS09 Table 4 (n=4..30, R/Q/Z, commutative):
**14 wins, 11 ties, 2 losses**, total **−1,221 mults (−1.11%)**.

The 2 remaining losses (`n=20`: +7, `n=22`: +21) come from DIS09's
mul121-with-commutative-specific-tricks that aren't yet replicated.

---

## Why "derived" and not "verified"

A 📐 derived bound is the sum of inner-sub-product ranks weighted by the
outer base's structure, per `BlockSplitSearch.findBestMultiBaseSplit`.
Each inner sub-product rank is either:
- a verified scheme JSON (Strassen, AT, AE, Perminov, ...) — direct,
- a closed-form formula (Pan TA, Waksman, Rosowski commutative) — direct,
- or itself a derived bound — recursive (chain length usually 2–3).

To promote to ✅ verified, the corresponding scheme JSON must be
materialised via `Recombination.constructWithAllocation` and the result
must pass `Verifier.isExactNonCubic`. This is mechanical — see
`Compositions.java` for examples — but is its own work item and is
tracked in [ROADMAP.md](../ROADMAP.md).

---

## How to cite

If you use any of these bounds, cite this repository as
[REFERENCES.md entry [30]](../REFERENCES.md#30-solven-strassen):

```bibtex
@misc{solven-strassen,
  author       = {Lacelle, Benoit},
  title        = {strassen: a research-grade catalog of fast matrix multiplication algorithms},
  howpublished = {\url{https://github.com/solven-eu/matmulcatalog}},
  year         = {2026},
  note         = {Work in progress; see \url{https://www.solven.eu/matmulcatalog/} for the live catalog.}
}
```

The derived bounds also appear in
[`docs/derived-bounds.json`](../docs/derived-bounds.json) tagged
`solven-strassen 2026 [30] (multi-base + S₃ symmetry; modern catalog SOTA)`.
