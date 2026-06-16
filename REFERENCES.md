# References

Central bibliography for the matmul-rank work in this repo.

**Entries are numbered `[1]..[N]` chronologically by publication year**, so
they can be cited compactly from other docs:

- Link the number: `[\[8\]](REFERENCES.md#8-cw90)` (e.g. for ω timelines).
- Link the key: `[strassen69](REFERENCES.md#1-strassen69)` (preferred in
  prose where the name carries meaning).

Each entry has:
1. A **BibTeX block** (paste-ready for LaTeX).
2. A short annotation about why it matters here.
3. A **Cited in** line listing the MD files that use the entry.

> **Confidence on BibTeX detail**: entries with DOIs / arXiv IDs are
> verified against the canonical source. Entries marked `[?]` need
> volume/page cross-check — pull requests welcome.

## Scheme provenance — sources we pull from

This catalog ingests schemes from a small number of upstream sources.
Each scheme file is named with its source prefix
(`alphatensor-*`, `perminov-*`, `fmm-lille_*`, `solven-strassen-*`,
`smirnov-*`, etc.) so the lineage is grep-able.

| Source | What it provides | Notes |
| --- | --- | --- |
| **Perminov** ([github.com/dronperminov/FastMatrixMultiplication](https://github.com/dronperminov/FastMatrixMultiplication)) | Per-field rank summaries via `schemes/status.json`; ZT-arithmetic scheme JSONs; mirrors of FMM Maple files | Authoritative for field classification (Z / Q / ZT). His own `schemes/results/ZT/` is the largest single contributor to our catalog. |
| **FMM (Université de Lille)** ([fmm.univ-lille.fr](https://fmm.univ-lille.fr/)) | Maple tensor files for Q-arithmetic schemes from various papers (Smirnov, Pan, Sedoglavic, …) | Perminov references these as upstream for Q-rank rows; we fetch via Perminov's plain-text mirror by default, fall back to FMM bz2. See `tools/import_fmm_maple.py:FORCE_FMM_UPSTREAM` for shapes where Perminov's mirror is known stale. |
| **AlphaTensor (DeepMind)** ([github.com/google-deepmind/alphatensor](https://github.com/google-deepmind/alphatensor)) | F₂ and R `factorizations_{f2,r}.npz` | Imported via `tools/import_alphatensor_*.py`. |
| **AlphaEvolve 2025** ([arXiv:2506.13131](https://arxiv.org/abs/2506.13131); [HTML](https://arxiv.org/html/2506.13131)) | ⟨4,4,4⟩=48 over C (and a handful of others) | Imported via `tools/import_alphaevolve.py`. |
| **Sedoglavic FMM digest** ([github.com/sedoglavic/fmm_digest](https://github.com/sedoglavic/fmm_digest)) | A manually-curated `fmm_sota.json` covering 5,456 shapes with rank + tensor-exponent. Complements FMM Université de Lille — includes shapes neither Perminov nor FMM tabulates (e.g. cubic ⟨18,18,18⟩=3200). | Cached at `references/sedoglavic-fmm-sota.json`. Periodic refresh via `eu.solven.matmul.docs.explore.ScanUpstreamSources` (`mvn exec:java`) + GitHub Action `.github/workflows/scan-upstream-sources.yml`. |
| **arbenson/fast-matmul** ([github.com/arbenson/fast-matmul](https://github.com/arbenson/fast-matmul)) | Reference implementations (Python + C++) of fast-matmul algorithms with explicit `U/V/W` factor matrices: Makarov, Hopcroft-Kerr, Pan TA, Strassen-Winograd variants. Especially valuable for *cross-validation* of reconstructions where we struggled to nail down formula details (e.g. HK page-10 same-method, Pan TA constructive recipe). | Not yet automated; consult per shape via the repo's `algorithms/` directory. Use as ground-truth when our derivations disagree with sympy. |
| **mkauers/matrix-multiplication** ([github.com/mkauers/matrix-multiplication](https://github.com/mkauers/matrix-multiplication)) | Manuel Kauers' (JKU Linz) GitHub data repo of explicit fast-matmul **schemes** — the machine-readable backing for the group's SAT-search [[84]](#84-heule-kauers-seidl-2019), flip-graph / meta-flip-graph, and structure-exploiting [[83]](#83-kauers-2026-structure) work (Kauers–Moosbauer–Wood). Many `⟨n,m,p⟩` decompositions across fields, including the bud-rich / structured ones our serendipitous engine (`#159`) feeds on. Partly already in our catalog as `kauers_2026-*`. | **Not yet systematically imported** — register as a primary source to mine for bud-rich bases and orbit diversity; the concrete data location for the JKU page [[85]](#85-linz-mm-catalog). Cross-check against Perminov/FMM before importing (dedup by content hash). |
| **DPS 2025** ([arXiv:2506.13242](https://arxiv.org/abs/2506.13242); [HTML](https://arxiv.org/html/2506.13242)) | Non-complex ⟨4,4,4⟩=48 (improves AlphaEvolve to R-arithmetic) | Pulled via FMM Maple. |
| **`solven-strassen-2026`** | Our own materialisations via `Recombination` over imported leaves | Tagged with year-and-author for repo attribution. |
| **plinopt — Dumas/Pernet/Sedoglavic** ([github.com/jgdumas/plinopt](https://github.com/jgdumas/plinopt)) | C++ toolkit for *addition-count optimization* of bilinear algorithms: given U/V/W factor matrices, searches for an equivalent factorization with fewer additions (independent of rank, which is fixed). References Dumas–Pernet–Sedoglavic papers on common-subexpression elimination + Kaporin-style optimizations. Read as: "we have R(⟨3,3,3⟩)=23, but what's the smallest `a` such that 23 mults + `a` additions suffice?" Adjacent in concept to the recent `_a59` / `_a58` / `_a60` ⟨3,3,3⟩ entries in our catalog. | Not integrated into our pipeline (additions-only optimization, no rank improvement). Worth running on our `solven-strassen-2026_*` recombination outputs to tighten their addition counts. |
| **GroupNames — Tim Dokchitser** ([people.maths.bris.ac.uk/~matyd/GroupNames/](https://people.maths.bris.ac.uk/~matyd/GroupNames/)) | A catalog of small finite groups (orders 1–500) with structure / character / subgroup tables. Tangentially relevant: some fast-matmul schemes are obtained by exploiting the group structure of the symmetry orbit on the matmul tensor (Strassen's invariance group, AlphaEvolve's continuous symmetries, etc.) — GroupNames is a quick reference for identifying which abstract group is acting. | Not currently linked to any catalog scheme metadata. Useful when classifying a newly-found scheme's symmetry group (e.g. "is this Strassen orbit S₃ × S₃ × S₃?"). |

### Hopcroft-Kerr family `⟨a,2,c⟩=(3ac+max(a,c))/2`: caveat

Hopcroft-Kerr 1971 (ref [2] below) introduces a **meta-algorithm**,
not a fixed scheme. The paper proves that `⌈(3pn+max(p,n))/2⌉`
multiplications *suffice* for `⟨p,2,n⟩` matmul without commutativity,
but the construction is parameterised by:

1. **Linear functionals** (paper Lemma 1) — multiple valid choices,
   none canonical.
2. **Method sequence** for diagonal elements (paper Lemma 3) — picks
   one of three Strassen-style decompositions per row, subject to a
   colouring constraint on consecutive rows.

As a consequence, **no upstream — neither Perminov nor FMM — publishes
a Hopcroft-Kerr scheme as explicit factor matrices**. For shapes like
`⟨2,10,15⟩=233`, FMM lists the *rank claim* (matching the HK formula
bound), but its Maple file contains a different (sub-optimal)
algorithm; Perminov tabulates the bound in `status.json` but his
`schemes.Q[]` array is empty for those entries. The Maple files he
*does* host for `⟨2,k,k⟩` shapes contain HK-equivalent schemes (we
imported them as `fmm-lille_2x*_r*` files).

For the genuinely-missing HK schemes, we have:

- A **closed-form rank bound** registered as a family entry in
  `docs/cited-bounds.json` via `GenerateCitedBounds.familyJson(…)`
  and rendered in the SPA via `expandFamilyEntry(c)` — so the catalog
  *table* shows the bound for any `⟨2,b,c⟩` up to `MAX_DIM=32`.
- **The constructor itself — DONE (2026-06-11/12, tasks #7/#10/#11)**:
  `GenerateHk2npConstructed` emits the full family as explicit integer
  factor matrices — 465 schemes (`3 ≤ p ≤ 32, p ≤ n ≤ 32`), 459 at the
  exact formula, registered under `schemes/constructed/` with
  `hk71`-tagged filenames and dual-certified. We are, to our knowledge,
  the first to publish HK-attaining explicit factor matrices for the
  shapes above (e.g. `⟨2,10,15⟩=233/ℤ`). Recipe + the one open residue
  (six `g ≥ 6` shapes at +1..+3, task #9):
  `research/hopcroft-kerr-2np/CONSTRUCTIVE_METHOD.md`.

The PDF is archived at
[`references/papers/hopcroft_kerr_1971_2bc_2n2.pdf`](references/papers/hopcroft_kerr_1971_2bc_2n2.pdf).

**Independent confirmation of the formula**: Hopcroft & Musinski 1973
([\[81\]](#81-hopcroft-musinski-1973), draft p. 76) restate the HK bound
`⌈(3pn+max(n,p))/2⌉` for `p×2 by 2×n` (citing HK as their ref [3]) and prove the
**S₃-duality** (Thm 6 / Cor 7) that carries it to all dual formats — the
cleanest secondary source. For small cases, Alekseev & Smirnov 2013
([\[80\]](#80-alekseev-smirnov-2013), in Russian) compute the exact and border
ranks of `⟨4,2,2⟩` and `⟨2,2,2⟩`.

> ⚠️ **Not BCS97.** Contrary to an earlier note here, Bürgisser–Clausen–Shokrollahi
> 1997 ([\[65\]](#65-bcs97)) does **not** reproduce the HK `⟨a,2,c⟩` *bilinear*
> formula. Its 2×n material (Ex. 14.20, Notes §14.9) is the **commutative**
> Waksman/Winograd bound `L(⟨e,2,e⟩) ≤ e²+2e−1` plus Feig's `L(⟨n,2,2⟩) ≤ 3n+1`;
> it attributes the duality to Hopcroft–*Musinski* and cites Hopcroft–*Kerr*
> [252] only for the **substitution method** (lower bounds). Use [81]/[80], not
> [65], when sourcing the HK rank formula.

---

## Local PDF archive

Open-access PDFs of foundational and recently-cited papers are saved
under [`references/papers/`](references/papers/) so this repo remains
self-contained for offline / archival use. Each entry below that has a
local file is marked with a **Local PDF** line giving the relative
path.

Currently archived:
- `winograd_1971_2x2x2_lb7.pdf` ([3]) — proves `R(⟨2,2,2⟩) ≥ 7` (and `R_c ≥ 7` too).
- `laderman_1976_3x3x3_r23.pdf` ([4])
- `sedoglavic_2017_5x5x5_r99_hal01562131.pdf` (Sedoglavic 2017, `⟨5,5,5⟩=99` — NC, follows the same structured-padding recipe as his `⟨7,7,7⟩=250` paper)
- `drevet_islam_schost_2009_DrIsSc09.pdf` ([10])
- `drevet_schost_2010_MITACS_poster.pdf` ([20])
- `sedoglavic_2017_7x7x7_r250.pdf` ([21])
- `alphaevolve_2025_novikov_arxiv2506.13131.pdf` ([14])
- `wang_2026_F2_3x3x3_lb20_arxiv2603.07280.pdf` ([15])
- `stapleton_2025_3x3x3_r23_a60_arxiv2508.03857.pdf` (=[add60])
- `martensson_stankovski_stapleton_2025_3x3x3_r23_a59_arxiv2601.05272.pdf` ([27])
- `perminov_2025_3x3x3_r23_a58_arxiv2512.21980.pdf` ([79] / `[add58]`, by Andrew I. Perminov — the **same** researcher who maintains the [18] `FastMatrixMultiplication` repo; "Andrey" = transliteration of "Andrew")
- `rosowski_2019_commutative_matmul_arxiv1904.07683.pdf` ([28], commutative `R_c(⟨3,3,3⟩) ≤ 21` + generic formula for `⟨l,n,m⟩`)
- `alekseev_smirnov_2013_RU_4x2_2x2_bilinear_complexity_spm47.pdf` ([80], 🇷🇺 **in Russian** — exact + border bilinear ranks of ⟨4,2,2⟩ and ⟨2,2,2⟩; independent confirmation of the Hopcroft-Kerr `⟨a,2,c⟩` family. Math-Net.Ru, freely distributed)
- `schwartz_zwecher_2025_feasible_matmul_arxiv2508.01748.pdf` ([19], feasible ω ≤ 2.773203 via trilinear aggregation + de Groote)
- `alphatensor_2022_fawzi_nature.pdf` ([12], full Nature paper with Extended Data Table 1 used for AT provenance audit)
- `mezzarobba_2007_mthesis_dfinite.pdf` ([26], MPRI master's thesis — appendix tabulates commutative small-matmul ranks up to ⟨28,28,28⟩, the table DIS09 Table 4 cites as its commutative baseline; Fig. 3 gives the explicit Waksman ⟨3,3,3⟩=23 formula)
- `islam_2009_msc_optim_matmul.pdf` ([68], 96-page MSc thesis (Schost supervised) — Ch 6 generalizes Waksman to non-square ⟨m,n,p⟩; Chs 4-5 contain the addition-aware decomposition material that became DIS09 §3)
- `pan_2014_trilinear_apa_arxiv1412.1145.pdf` ([66], Pan's survey of trilinear decompositions + APA + summation; explicit TA recipes for small disjoint MM)
- `pan_2014b_history_arxiv1411.1972.pdf` ([67], "Better Late Than Never" — 6-page historical narrative + concrete disjoint-MM recipes)
- `dumas_pernet_sedoglavic_2025_4x4x4_r48_arxiv2506.13242.pdf` (Dumas-Pernet-Sedoglavic 2025, non-complex `R(⟨4,4,4⟩) ≤ 48` — improves AlphaEvolve 2025 by removing the complex-coefficient requirement)
- Dumas–Pernet–Sedoglavic 2025 (HAL [hal-05121550](https://hal.science/hal-05121550), June 2025) — "A non-commutative algorithm for multiplying a 3×4 matrix by a 4×7 matrix using 63 non-complex multiplications." Q/R/Z counterpart to AlphaEvolve's `⟨3,4,7⟩=63` over C. **Not yet imported** (see task #81). The catalog currently elects `perminov-ZT_3x4x7_r64_a454.json` (64) for Q targets and falls back to that for our search prediction — once the DPS scheme is imported, the last non-cubic FMM-Lille gap closes.

**TODO** (need manual download — paywalled or anti-bot protected):
- `strassen69` (Numerische Mathematik vol 13, Springer): no known open
  PDF — paywalled. The ALGORITHM ITSELF is reproduced in countless
  open sources (CLRS, Wikipedia, classroom notes), and we have it
  hardcoded in [`Strassen7.java`](src/main/java/io/cormoran/strassen/v3/Strassen7.java)
  + as scheme file `strassen_2x2x2_r7_a18.json`. The 1969 PDF would
  add provenance/history value only — the math is fully captured.
  Save to `references/papers/strassen_1969_gaussian_elimination_not_optimal.pdf`
  if accessible.
- `hk71` (Hopcroft-Kerr 1971): no known open PDF. SIAM J. Applied Math
  vol 20 — same era as Winograd 1971, gives independent LB proof.
- `makarov87` ([23], Makarov 1987, `USSR Comput. Math. Math. Phys.` vol 27;
  Russian original `Zh. Vychisl. Mat. i Mat. Fiz.` 27): try [mathnet.ru](https://www.mathnet.ru/eng/zvmmf)
  (manual browse — automated search blocked by JS) or
  ScienceDirect's translated edition
  ([USSR Comp. Math. Math. Phys. vol 10/5](https://www.sciencedirect.com/journal/ussr-computational-mathematics-and-mathematical-physics/vol/10/issue/5)).
- `makarov86` ([29], same journal, vol 26 no 2 pp 293-294, 1986): try
  [mathnet.ru](https://www.mathnet.ru/eng/zvmmf) or
  [USSR Comp. Math. Math. Phys. vol 26/2](https://www.sciencedirect.com/journal/ussr-computational-mathematics-and-mathematical-physics/vol/26/issue/2).
  Likely paperid in the 4000-range based on era. Russian-original 2-page
  algorithm announcement.
- `smirnov2017` ([70], "Several Bilinear Algorithms for Matrix
  Multiplication Problems ⟨3, P, Q⟩"): hosted on ResearchGate
  ([publication 313064941](https://www.researchgate.net/publication/313064941_Several_bilinear_algorithms_for_matrix_multiplication))
  and Academia.edu ([36266758](https://www.academia.edu/36266758/Several_Bilinear_Algorithms_for_Matrix_Multiplication));
  both return 403 to curl/WebFetch. Save to
  `references/papers/smirnov_2017_3pq_several_bilinear.pdf` after a
  manual browser download. Contains explicit factor matrices for
  ⟨3,4,6⟩=56, ⟨3,5,5⟩=58, ⟨3,4,7⟩=66, ⟨3,4,8⟩=75, ⟨3,5,7⟩=82.

When adding a new reference, save a copy locally if open-access and
add the **Local PDF** line linking the file path.

---

## Index

| # | key | year | type | topic |
|---|---|---|---|---|
| [1](#1-strassen69) | strassen69 | 1969 | algorithm | `⟨2,2,2⟩ = 7`; field-agnostic recursion |
| [2](#2-hk71) | hk71 | 1971 | LB + algorithm | `R(⟨2,2,2⟩) = 7` tight; `R(⟨2,2,3⟩) ≥ 11`, `R(⟨2,3,3⟩) ≥ 15` |
| [3](#3-winograd71) | winograd71 | 1971 | LB proof | `R(⟨2,2,2⟩) ≥ 7` (matches Strassen UB; also `R_c ≥ 7` — no commutative advantage at size 2); complex mult in 3 |
| [4](#4-lad76) | lad76 | 1976 | algorithm | `R(⟨3,3,3⟩) ≤ 23` — still best |
| [5](#5-pan78) | pan78 | 1978 | algorithm + ω | Trilinear aggregation, ω ≤ 2.78 |
| [6](#6-bini79) | bini79 | 1979 | border rank | `R̃(⟨2,2,2⟩) ≤ 5`; ω ≤ 2.7799 |
| [7](#7-sch81) | sch81 | 1981 | ω | τ-theorem; ω ≤ 2.522 |
| [8](#8-cw90) | cw90 | 1990 | ω | Coppersmith–Winograd framework; ω ≤ 2.376 |
| [9](#9-blaser03) | blaser03 | 2003 | LB | `R(⟨3,3,3⟩) ≥ 19` (substitution method) |
| [10](#10-drisc09) | drisc09 | 2009/2011 | survey | Cat up to ⟨30,30,30⟩; commutative vs non-commutative |
| [11](#11-smirnov) | smirnov | 2013 | catalog | Small-format `{-1,0,+1}` algorithms |
| [12](#12-alphatensor) | alphatensor | 2022 | algorithm (RL) | `R_{F₂}(⟨4,4,4⟩) ≤ 47`; many F₂ improvements |
| [13](#13-williams2024) | williams2024 | 2024 | ω | ω ≤ 2.371552 (current best) |
| [14](#14-alphaevolve) | alphaevolve | 2025 | algorithm (LLM-search) | `R_C(⟨4,4,4⟩) ≤ 48`; complex-specific |
| [15](#15-wang26) | wang26 | 2026 | LB | `R_{F₂}(⟨3,3,3⟩) ≥ 20` |
| [16](#16-phialsbasement) | phialsbasement | 2025 | verification | AlphaEvolve `⟨4,4,4⟩=48` factor matrices |
| [17](#17-fmm-lille) | fmm-lille | (active) | catalog | Université de Lille hosted catalog, 5,426 entries, includes 2025 AlphaEvolve/Kauers–Wood/Perminov |
| [18](#18-perminov) | perminov | (active) | catalog | Community fast-matmul scheme database |
| [19](#19-schwartz-zwecher25) | schwartz-zwecher25 | 2025 | ω (feasible) | `O(n^{2.773203})` — improves Pan 1982 for practical base cases ≥ 28 |
| [20](#20-drevet-schost-poster) | drevet-schost-poster | 2010 | survey / poster | MITACS poster — early presentation of the work that became Drevet–Islam–Schost 2011 |
| [22](#22-waksman70) | waksman70 | 1970 | algorithm (commutative) | Commutative matmul family used in DIS09 Table 4 |
| [23](#23-makarov87) | makarov87 | 1987 [?] | algorithm | NC `⟨5,5,5⟩=100` |
| [29](#29-makarov86) | makarov86 | 1986 | algorithm (commutative) | Cmt `⟨3,3,3⟩=22` (different paper from [23]; both by O. M. Makarov) |
| [24](#24-probert-fischer-80) | probert-fischer-80 | 1980 | survey | Computer-aided baseline cited in DIS09 Table 3 |
| [25](#25-smith-2002) | smith-2002 | 2002 | survey | Computer-aided NC baseline cited in DIS09 Table 3 |
| [26](#26-mezzarobba-2007) | mezzarobba-2007 | 2007 | survey | Commutative baseline cited in DIS09 Table 4 |
| [27](#27-mss25) | mss25 | 2025 | algorithm (additions) | `⟨3,3,3⟩ r=23` with **59 additions** (no basis change) — current best |
| [28](#28-rosowski19) | rosowski19 | 2019/2020 | algorithm (commutative) | **`R_c(⟨3,3,3⟩) ≤ 21`** beats Makarov's 22; general formula for `⟨l,n,m⟩` |
| [21](#21-sedoglavic17) | sedoglavic17 | 2017 | algorithm + methodology | `R(⟨7,7,7⟩) ≤ 250` via padded Kronecker; foundation of the fmm-lille catalog |
| [30](#30-solven-strassen) | solven-strassen | 2026 (WIP) | catalog + reproduction | This repo: DIS09 reproduction with modern catalog; multi-base + S₃ symmetry; improves 12 of 27 DIS09 ⟨n,n,n⟩ NC bounds for n=4..30 |
| [31](#31-schachtel78) | schachtel78 | 1978 | algorithm | NC `⟨5,5,5⟩=103`; predates Smirnov / DIS09 by decades |
| [65](#65-bcs97) | bcs97 | 1997 | textbook | Comprehensive treatment of algebraic-complexity theory; §14–16 cover bilinear-form complexity, fast matmul, trilinear aggregation, partial matrix mult |
| [66](#66-pan2014) | pan2014 | 2014 | survey | Pan's own survey of trilinear decompositions, APA algorithms, and summation tricks — explicit construction patterns for TA used by DIS09 |
| [67](#67-pan2014b) | pan2014b | 2014 | history | "Better Late Than Never" — Pan's 6-page narrative of the void he's filling, with concrete TA recipes for small disjoint MM cases |
| [68](#68-islam2009) | islam2009 | 2009 | thesis | MSc thesis: generalizes Waksman's commutative algorithm to non-square ⟨m,n,p⟩; superseded for ⟨n,3,3⟩ by Rosowski 2019 but historically the first non-square commutative bound. Also chapter-summarises DIS09's pre-publication material |
| [69](#69-dumas-pernet-sedoglavic-2025) | dumas-pernet-sedoglavic-2025 | 2025 | algorithm | `Q⟨3,4,7⟩:m=63` — non-complex rational realisation of AlphaEvolve's complex `⟨3,4,7⟩=63` via a Klein-four isotropy. Closes the last non-cubic FMM-Lille gap. |
| [70](#70-smirnov2017) | smirnov2017 | 2017 | catalog (NC, ⟨3,P,Q⟩) | Follow-up to [11] Smirnov 2013 — explicit factor matrices for ⟨3,4,6⟩=56, ⟨3,5,5⟩=58, ⟨3,4,7⟩=66, ⟨3,4,8⟩=75, ⟨3,5,7⟩=82 over ℤ |
| [71](#71-sedoglavic-smirnov-2021) | sedoglavic-smirnov-2021 | 2021 | border-rank | A. Sedoglavic & Alexey Vladimirovich Smirnov, "The tensor rank of 5x5 matrices multiplication is bounded by 98 and its border rank by 89" — tightens ⟨5,5,5⟩ border-rank UB to 89 (still SOTA); exact-rank 98 now dominated by AlphaEvolve 93 |
| [72](#72-sedoglavic-2017-fmm-methodology) | sedoglavic-2017-fmm-methodology | 2017 | **METHODOLOGY** | A. Sedoglavic, "A non-commutative algorithm for multiplying (7×7) matrices using 250 multiplications" (hal-01572046v2). **The article behind the FMM-Lille catalog for cubic decompositions** — Sedoglavic Prop 1: `⟨u+v,u+v,u+v⟩ ≤ ⟨u,u,u⟩ + 3⟨u,u,v⟩ + 3⟨v,v,u⟩` for u>v. Yields ⟨7⟩³=250, ⟨11⟩³=873, ⟨19⟩³=4044 and more via single closed form. |
| [73](#73-burichenko-2014) | burichenko-2014 | 2014 | symmetry | V. P. Burichenko, "On symmetries of the Strassen algorithm" (arXiv:1408.6273). Computes the discrete stabilizer of Strassen's ⟨2,2,2⟩=7 (order 36 = S₃ × Z₂³). Cited by our rank-7 ⟨2,2,2⟩ orbit enumeration. |
| [74](#74-landsberg-2008) | landsberg-2008 | 2008 | LB survey | J. M. Landsberg, "Geometry and the complexity of matrix multiplication" (Bull. AMS 45(2): 247–284). Geometric survey of matmul complexity; Theorem 3.8.4 cites **Bläser's LB R(⟨m,m,m⟩) ≥ (5/2)m² − 3m** plus the sharper R(⟨3,3,3⟩) ≥ 19. |
| [75](#75-perminov-2026-arxiv) | perminov-2026-arxiv | 2026 | catalog | A. I. Perminov, *Fast Matrix Multiplication in Small Formats: Discovering New Schemes with an Open-Source Flip Graph Framework*, arXiv:2603.02398 (2026). The umbrella paper behind the `FastMatrixMultiplication` repo — describes the open-source flip-graph framework that produced most `perminov-*` schemes we import. Cite this when our catalog imports a `perminov-*` scheme.  [Abs](https://arxiv.org/abs/2603.02398) |
| [76](#76-smirnov-2022) | smirnov-2022 | 2022 | algorithm | Alexey Vladimirovich Smirnov, "Bilinear algorithm for matrix multiplication ⟨4,4,9⟩=104. An irreducibly irrational solution of the Brent system?" (Oct 2022). Source for ⟨4,4,9⟩=104 imported via FMM-Lille Maple mirror — re-attributed from fmm-lille → Smirnov 2022 per 2026-06-03 audit. |
| [77](#77-perminov-2025-metaflip) | perminov-2025-metaflip | 2025 | method | A. I. Perminov, *Fast Matrix Multiplication via Ternary Meta Flip Graphs*, arXiv:2511.20317 (2025). The ternary meta-flip-graph search method underlying many `perminov-*` rank results. [Abs](https://arxiv.org/abs/2511.20317) |
| [78](#78-perminov-2025-additive) | perminov-2025-additive | 2025 | method | A. I. Perminov, *Parallel Heuristic Exploration for Additive Complexity Reduction in Fast Matrix Multiplication*, arXiv:2512.13365 (2025). Addition-count (`_a…`) reduction heuristic — source of the low-addition variants of `perminov-*` schemes. [Abs](https://arxiv.org/abs/2512.13365) |
| [79](#79-perminov-2025-3x3-a58) | perminov-2025-3x3-a58 | 2025 | algorithm | A. I. Perminov, *A 58-Addition, Rank-23 Scheme for General 3×3 Matrix Multiplication*, arXiv:2512.21980 (Dec 2025). The `R(⟨3,3,3⟩)=23` scheme with only 58 additions and no basis change (also cited as `[add58]` in [RANK_KNOWLEDGE.md](paper/theory/RANK_KNOWLEDGE.md)). [Abs](https://arxiv.org/abs/2512.21980) |
| [80](#80-alekseev-smirnov-2013) | alekseev-smirnov-2013 | 2013 | algorithm (🇷🇺 Russian) | V. B. Alekseev & Alexey Vladimirovich Smirnov, "On the exact and approximate bilinear complexities of multiplication of 4×2 and 2×2 matrices" (**in Russian**), *Sovr. Probl. Mat.* 17:135–152. Independent confirmation of the exact + border ranks for the inner-dimension-2 Hopcroft-Kerr family (`⟨4,2,2⟩`, `⟨2,2,2⟩`). DOI [10.4213/spm47](https://doi.org/10.4213/spm47); [ResearchGate (RU)](https://www.researchgate.net/publication/284529885_On_the_exact_and_approximate_bilinear_complexities_of_multiplication_of_42_and_22_matrices). |
| [81](#81-hopcroft-musinski-1973) | hopcroft-musinski-1973 | 1973 | duality / algorithm | J. E. Hopcroft & J. Musinski, "Duality applied to the complexity of matrix multiplication and other bilinear forms," *SIAM J. Comput.* 2(3):159–173. S₃-invariance of matmul rank (Thm 6 / Cor 7); states the Hopcroft-Kerr `⌈(3pn+max(n,p))/2⌉` bound for `p×2 by 2×n` (draft p. 76) and propagates it by duality. Local draft PDF + [cr.yp.to mirror](https://cr.yp.to/bib/1973/hopcroft-duality-draft.pdf). |
| [82](#82-deza-2023-cp) | deza-2023-cp | 2023 | method (exact / infeasibility) | A. Deza, C. Liu, P. Vaezipoor & E. B. Khalil, *Fast Matrix Multiplication Without Tears: A Constraint Programming Approach*, CP 2023 (LIPIcs vol. 280, art. 26), arXiv:2306.01097. CSP/CP formulation of the Brent equations over `{−1,0,1}` with **symmetry-breaking constraints + valid inequalities** enabling *exact* search and **infeasibility proofs** up to ⟨3,3,3⟩=23 — the complement to our heuristic upper-bound construction. Directly comparable to our axis-flip / allocation canonicalisation. Code: [khalil-research/Matrix-Mult-CP](https://github.com/khalil-research/Matrix-Mult-CP). [Abs](https://arxiv.org/abs/2306.01097) |
| [83](#83-kauers-2026-structure) | kauers-2026-structure | 2026 | method (**serendipitous / buds**) | M. Kauers, J. Moosbauer & I. Wood, *Exploiting the Structure in Tensor Decompositions for Matrix Multiplication*, arXiv:2602.11041 (May 2026; submitted to Elsevier). **The published formalisation of exactly our bud / serendipitous-product idea**: recursive calls that *share an output* (or whose output feeds multiple positions) are merged into a *single larger* matrix multiplication, giving an effective exponent **lower than the tensor rank suggests** — improving 6×6 from ω 2.8075→**2.8019**. Direct prior/concurrent art for our `#159` serendipitous engine, `SerendipitousBudProduct`, and `BudBaseFactory`; cite in the paper. [Abs](https://arxiv.org/abs/2602.11041) |
| [84](#84-heule-kauers-seidl-2019) | heule-kauers-seidl-2019 | 2019 | method (SAT search) / algorithm | M. J. H. Heule, M. Kauers & M. Seidl, *New Ways to Multiply 3×3-Matrices*, arXiv:1905.10192 (2019); J. Symbolic Computation **104** (2021) 899–916. SAT/heuristic search that produced **many distinct rank-23 ⟨3,3,3⟩ schemes** over ℤ and small fields, mapping out the solution variety (not a rank improvement — ⟨3,3,3⟩=23 is Laderman 1976 — but a *diversity* result, directly relevant to our orbit/uniqueness work and to bud-rich base selection). PDF: [cs.cmu.edu/~mheule/publications/CCA19.pdf](https://www.cs.cmu.edu/~mheule/publications/CCA19.pdf); [arXiv](https://arxiv.org/abs/1905.10192). |
| [85](#85-linz-mm-catalog) | linz-mm-catalog | (active) | catalog (⟨3,3,3⟩, value TBD) | Kauers group (JKU Linz), *Matrix Multiplication* research data page — online catalog of ⟨3,3,3⟩=23 schemes from the SAT search [[84]](#84-heule-kauers-seidl-2019) and later flip-graph work. **Not yet imported / evaluated** — register as a source to mine for bud-rich ⟨3,3,3⟩ bases and orbit diversity; its incremental value over our existing ⟨3,3,3⟩ entries is **undetermined** (see ROADMAP). Web page: [algebra.uni-linz.ac.at/research/matrix-multiplication](http://www.algebra.uni-linz.ac.at/research/matrix-multiplication/index.html); machine-readable schemes: [github.com/mkauers/matrix-multiplication](https://github.com/mkauers/matrix-multiplication) (see the provenance table at the top). |
| [86](#86-perminov-2026-serendipitous) | perminov-2026-serendipitous | 2026 | method (**serendipitous / buds**) | A. I. Perminov, *Meta Flip Graph meets Serendipitous Product: new Fast Matrix Multiplication results*, arXiv:2606.02480 (June 2026). Combines the meta-flip-graph search [[77]](#77-perminov-2025-metaflip) with the **serendipitous product** (bud fusion; cf. Smith 2002 eq. (69) and Kauers–Moosbauer–Wood [[83]](#83-kauers-2026-structure)) — improving **207 rectangular formats** (≤16×16×16) and the **17–32 band** our catalog cites. **The source for the `Perminov 2026 (serendipitous)` rank claims** in `docs/cited-bounds.json` and `references/catalogs/perminov-serendipitous-catalog.json`. A distinct paper from the umbrella framework [[75]](#75-perminov-2026-arxiv) (arXiv:2603.02398). [Abs](https://arxiv.org/abs/2606.02480) |
| [87](#87-kaporin-2024-brent) | kaporin-2024-brent | 2024 | algorithm (**C⟨4,4,4⟩=48**) | I. E. Kaporin, *Semi-analytical solution of Brent equations*, Doklady Mathematics **518**(1):29–34 (2024), DOI [10.31857/S2686954324040056](https://doi.org/10.31857/S2686954324040056). A parametrisation of the Brent equations (cyclic symmetry, several-fold fewer unknowns) solved numerically; yields explicit **complex** designs **(4,4,4;48)** and **(2,4,5;32)**. The (4,4,4;48) scheme — verified to floating-point tolerance in the author's companion `test444r48.for` — is an **independent C-coefficient 48** that **predates AlphaEvolve 2025** [[14]](#14-alphaevolve) (existence of `r<49` over C was conjectured by Li–Zhang–Ke 2023). Loaded as `4x4x4-r48-kaporin_2024-*.json`. [Article](https://journals.rcsi.science/2686-9543/article/view/269374) |

---

## [1] strassen69

```bibtex
@article{strassen69,
  author  = {Strassen, Volker},
  title   = {Gaussian elimination is not optimal},
  journal = {Numerische Mathematik},
  volume  = {13},
  number  = {4},
  pages   = {354--356},
  year    = {1969},
  doi     = {10.1007/BF02165411}
}
```

The original `⟨2,2,2⟩` rank-7 algorithm. Field-agnostic (`±1` coefficients);
recurses to matrix entries. Implies ω ≤ log₂7 ≈ 2.807.
Cited in: [RANK_KNOWLEDGE.md](paper/theory/RANK_KNOWLEDGE.md) §1,
[SMALL_MATMUL_CATALOG.md](paper/theory/SMALL_MATMUL_CATALOG.md) §4.1,
[OMEGA_HISTORY.md](paper/theory/OMEGA_HISTORY.md).

## [2] hk71

The published SIAM 1971 paper was preceded by Cornell Tech Report
69-44 (September 1969) — same month as Strassen's announcement of
⟨2,2,2⟩=7. HK targeted Strassen's algorithm specifically, which is
why their lower bound `R(⟨2,2,2⟩) ≥ 7` so cleanly closes Strassen's
upper bound. The 1969 tech report is the version archived locally
(34 pages; the SIAM article condenses to 7 pages).

```bibtex
@article{hk71,
  author  = {Hopcroft, John E. and Kerr, Leslie R.},
  title   = {On minimizing the number of multiplications necessary for matrix multiplication},
  journal = {SIAM Journal on Applied Mathematics},
  volume  = {20},
  number  = {1},
  pages   = {30--36},
  year    = {1971},
  doi     = {10.1137/0120004}
}

@techreport{hk1969-cornell-tr69-44,
  author      = {Hopcroft, John E. and Kerr, Leslie R.},
  title       = {On minimizing the number of multiplications necessary for matrix multiplication},
  institution = {Cornell University, Department of Computer Science},
  number      = {TR 69-44},
  month       = sep,
  year        = {1969}
}
```

Proves `R(⟨2,2,2⟩) = 7` (tight, char-0 fields) and `R(⟨2,2,3⟩) ≥ 11`,
`R(⟨2,3,3⟩) ≥ 15`. Cited in:
[RANK_KNOWLEDGE.md](paper/theory/RANK_KNOWLEDGE.md) §1, §5,
[SMALL_MATMUL_CATALOG.md](paper/theory/SMALL_MATMUL_CATALOG.md) §4–§5.

## [3] winograd71

```bibtex
@article{winograd71,
  author  = {Winograd, Shmuel},
  title   = {On multiplication of $2 \times 2$ matrices},
  journal = {Linear Algebra and its Applications},
  volume  = {4},
  number  = {4},
  pages   = {381--388},
  year    = {1971},
  doi     = {10.1016/0024-3795(71)90009-7}
}
```

Proves (i) **`R(⟨2,2,2⟩) ≥ 7`** — the lower bound matching Strassen 1969's
upper bound, establishing `R(⟨2,2,2⟩) = 7` exactly. Critically, Winograd's
proof goes through **even when commutativity is assumed**, so
`R_c(⟨2,2,2⟩) ≥ 7` too — commutativity gives no advantage at size 2.
Also proves (ii) `R_C(⟨1,2,2⟩) = 3` (complex multiplication via 3 real
multiplications), the classic `(a+bi)(c+di)` algorithm. Together with
Hopcroft-Kerr 1971 ([2]), the foundational tight LB for matmul.
Together with Strassen 1969 ([1]), establishes the historical baseline
`R(⟨2,2,2⟩) = R_c(⟨2,2,2⟩) = 7` that every subsequent matmul-rank result
builds on. **Local PDF**: [`references/papers/winograd_1971_2x2x2_lb7.pdf`](references/papers/winograd_1971_2x2x2_lb7.pdf).
Cited in: [RANK_KNOWLEDGE.md](paper/theory/RANK_KNOWLEDGE.md) §0, §1.2,
[SMALL_MATMUL_CATALOG.md](paper/theory/SMALL_MATMUL_CATALOG.md) §4.1.

## [4] lad76

```bibtex
@article{lad76,
  author  = {Laderman, Julian D.},
  title   = {A noncommutative algorithm for multiplying $(3 \times 3)$ matrices using 23 multiplications},
  journal = {Bulletin of the American Mathematical Society},
  volume  = {82},
  number  = {1},
  pages   = {126--128},
  year    = {1976},
  doi     = {10.1090/S0002-9904-1976-13988-2}
}
```

`R(⟨3,3,3⟩) ≤ 23` — still the best known after 49 years. Coefficients
`±1`, valid over any field.
**Local PDF**: [`references/papers/laderman_1976_3x3x3_r23.pdf`](references/papers/laderman_1976_3x3x3_r23.pdf)
(via AMS mirror — Project Euclid blocks curl with Incapsula, but
[ams.org direct PDF](https://www.ams.org/journals/bull/1976-82-01/S0002-9904-1976-13988-2/S0002-9904-1976-13988-2.pdf)
serves cleanly). Cited in:
[RANK_KNOWLEDGE.md](paper/theory/RANK_KNOWLEDGE.md) §1,
[SMALL_MATMUL_CATALOG.md](paper/theory/SMALL_MATMUL_CATALOG.md) §4.2.

## [5] pan78

```bibtex
@inproceedings{pan78,
  author    = {Pan, Victor Ya.},
  title     = {Strassen's algorithm is not optimal: Trilinear technique of aggregating, uniting and canceling for constructing fast algorithms for matrix operations},
  booktitle = {19th Annual Symposium on Foundations of Computer Science (FOCS 1978)},
  pages     = {166--176},
  year      = {1978},
  doi       = {10.1109/SFCS.1978.34}
}
```

Trilinear-aggregation framework underpinning sub-`7^k` asymptotic
matmul (ω ≤ 2.78). Cited in: [RANK_KNOWLEDGE.md](paper/theory/RANK_KNOWLEDGE.md) §1,
[SMALL_MATMUL_CATALOG.md](paper/theory/SMALL_MATMUL_CATALOG.md) §4–§5,
[OMEGA_HISTORY.md](paper/theory/OMEGA_HISTORY.md).

## [6] bini79

```bibtex
@article{bini79,
  author  = {Bini, Dario and Capovani, Milvio and Romani, Francesco and Lotti, Grazia},
  title   = {$O(n^{2.7799})$ complexity for $n \times n$ approximate matrix multiplication},
  journal = {Information Processing Letters},
  volume  = {8},
  number  = {5},
  pages   = {234--235},
  year    = {1979},
  doi     = {10.1016/0020-0190(79)90113-3}
}
```

First result that **border rank** can beat exact rank for matmul:
`R̃(⟨2,2,2⟩) ≤ 5`, giving ω ≤ 2.7799 via approximation.
Cited in: [RANK_KNOWLEDGE.md](paper/theory/RANK_KNOWLEDGE.md) §1, §5,
[OMEGA_HISTORY.md](paper/theory/OMEGA_HISTORY.md).

## [7] sch81

```bibtex
@article{sch81,
  author  = {Sch\"{o}nhage, Arnold},
  title   = {Partial and total matrix multiplication},
  journal = {SIAM Journal on Computing},
  volume  = {10},
  number  = {3},
  pages   = {434--455},
  year    = {1981},
  doi     = {10.1137/0210032}
}
```

τ-theorem / asymptotic-sum machinery; ω ≤ 2.522.
Cited in: [RANK_KNOWLEDGE.md](paper/theory/RANK_KNOWLEDGE.md) §1, §4,
[OMEGA_HISTORY.md](paper/theory/OMEGA_HISTORY.md).

## [8] cw90

```bibtex
@article{cw90,
  author  = {Coppersmith, Don and Winograd, Shmuel},
  title   = {Matrix multiplication via arithmetic progressions},
  journal = {Journal of Symbolic Computation},
  volume  = {9},
  number  = {3},
  pages   = {251--280},
  year    = {1990},
  doi     = {10.1016/S0747-7171(08)80013-2}
}
```

The Coppersmith–Winograd framework that underpinned all ω improvements
from 1987 to ~2010. ω ≤ 2.376.
Cited in: [RANK_KNOWLEDGE.md](paper/theory/RANK_KNOWLEDGE.md) §4,
[OMEGA_HISTORY.md](paper/theory/OMEGA_HISTORY.md).

## [9] blaser03

```bibtex
@article{blaser03,
  author  = {Bl\"{a}ser, Markus},
  title   = {On the complexity of the multiplication of matrices of small formats},
  journal = {Journal of Complexity},
  volume  = {19},
  number  = {1},
  pages   = {43--60},
  year    = {2003},
  doi     = {10.1016/S0885-064X(02)00007-9}
}
```

`R(⟨3,3,3⟩) ≥ 19` via the substitution method; standing LB over
arbitrary fields until 2026. Cited in:
[RANK_KNOWLEDGE.md](paper/theory/RANK_KNOWLEDGE.md) §1, §5,
[SMALL_MATMUL_CATALOG.md](paper/theory/SMALL_MATMUL_CATALOG.md) §4.2.

## [10] drisc09

**Technique class**: *sparseness-aware recursive padding* — pad the target to a
multiple of the base pattern, then prune the recursive sub-products that
acquire zero rows/columns or feed only a few needed outputs, computing each at
its *effective* smaller size (§2; the ⟨3,3,3⟩→pad-⟨4,4,4⟩ worked example drops
49→25 mults). This is NOT the Schönhage "serendipitous / partial-product"
construction (that is borrow-and-correct; see `KnownTauIdentities.BorrowAndCorrect`
and [65] BCS97 §14–16). Our `Recombination` peel/padding (#87) + output-side
zero masks (#86, "Islam Ch. 4 γ5 reduction") already implement this technique.

*Local PDF + extracted tables.*

**Local PDF**: [`references/papers/drevet_islam_schost_2009_DrIsSc09.pdf`](references/papers/drevet_islam_schost_2009_DrIsSc09.pdf)
(full paper, 32 pages).
**Extracted tables**: [`references/dis09-cubic-tables.json`](references/dis09-cubic-tables.json)
(Tables 3 & 4: per-cubic-format ranks for n ∈ [2, 30] in both
non-commutative and commutative cases, with provenance vs
Probert-Fischer 1980 / Smith 2002 / Mezzarobba 2007 baselines).


```bibtex
@article{drisc09,
  author  = {Drevet, Charles-{\'E}ric and Islam, Md. Nazrul and Schost, {\'E}ric},
  title   = {Optimization techniques for small matrix multiplication},
  journal = {Theoretical Computer Science},
  volume  = {412},
  number  = {22},
  pages   = {2219--2236},
  year    = {2011},
  doi     = {10.1016/j.tcs.2010.10.012},
  note    = {Preprint dated 2009; published 2011.}
}
```

Historical synthesis going up to `⟨30,30,30⟩`, with separate tables for
commutative vs non-commutative ranks. Pre-dates AlphaTensor and
AlphaEvolve. Direct PDF:
[cs.uwaterloo.ca/~eschost/publications/DrIsSc09.pdf](https://cs.uwaterloo.ca/~eschost/publications/DrIsSc09.pdf).

**⚠ Errata**: Appendix A.1 Lemma 4 (the explicit `t(Ã, B̃, C̃)` formula
that encodes the trilinear-aggregation construction) has a
copy-paste typo in the 7th `u`-correction term — slot 3 reads
`−C̃^{1,2}` where it should be `−C̃^{2,1}` per Islam's 2009 MSc thesis
page 50 (same lemma). See [`references/typos.md`](references/typos.md#drevet-islam-schost-2009--lemma-4-7th-u-correction).
Islam's Magma validation code, originally hosted at
`http://www.csd.uwo.ca/~mislam63/`, is no longer reachable but is
recovered via the Wayback Machine
([snapshot](https://web.archive.org/web/20120223044300/http://www.csd.uwo.ca:80/~mislam63/TA.mgm))
and archived locally at
[`references/islam2009/magma/TA.mgm`](references/islam2009/magma/TA.mgm).
The Magma source backs the thesis version of the formula.
Cited in: [RANK_KNOWLEDGE.md](paper/theory/RANK_KNOWLEDGE.md) §1.2,
[SMALL_MATMUL_CATALOG.md](paper/theory/SMALL_MATMUL_CATALOG.md) §2.

## [11] smirnov

```bibtex
@article{smirnov,
  author  = {Smirnov, Alexey Vladimirovich},
  title   = {Bilinear complexity and practical algorithms for matrix multiplication},
  journal = {Computational Mathematics and Mathematical Physics},
  volume  = {53},
  number  = {12},
  pages   = {1781--1795},
  year    = {2013},
  doi     = {10.1134/S0965542513120129}
}
```

Comprehensive catalog of integer / `{-1, 0, +1}` algorithms for many
small formats. Pre-dates AlphaTensor 2022, so its Z/2 entries are
superseded for some formats. **Attribution-relevant small-shape
results** that Sedoglavic 2017 (`⟨7,7,7⟩=250`) credits to Smirnov:

| shape | rank | additions in Smirnov 2013 |
| --- | --: | --: |
| `⟨3,3,4⟩` | 29 | 36 |
| `⟨3,3,5⟩` | 36 | 45 |
| `⟨3,3,6⟩` | 40 | 54 |
| `⟨3,4,4⟩` | 38 | 48 |
| `⟨3,4,5⟩` | 48 | 60 |

When these shapes appear as leaves of a composed scheme, the
attribution chain should resolve down to Smirnov 2013 (not to a later
re-importer like Perminov or AlphaTensor that merely re-encoded them).

**Access**:
- English translation **paywalled** at Springer:
  <https://link.springer.com/article/10.1134/S0965542513120129>
- Russian original (free, archived locally):
  [`references/papers/smirnov_2013_RU_zvmmf9955.pdf`](references/papers/smirnov_2013_RU_zvmmf9955.pdf)
  — the explicit factor matrices in the Russian PDF are language-
  neutral and usable for cross-checking imported schemes.

Cited in: [RANK_KNOWLEDGE.md](paper/theory/RANK_KNOWLEDGE.md) §1,
[SMALL_MATMUL_CATALOG.md](paper/theory/SMALL_MATMUL_CATALOG.md) §2, §4.4.

## [12] alphatensor

```bibtex
@article{alphatensor,
  author  = {Fawzi, Alhussein and Balog, Matej and Huang, Aja and Hubert, Thomas and Romera-Paredes, Bernardino and Barekatain, Mohammadamin and Novikov, Alexander and Ruiz, Francisco J. R. and Schrittwieser, Julian and Swirszcz, Grzegorz and Silver, David and Hassabis, Demis and Kohli, Pushmeet},
  title   = {Discovering faster matrix multiplication algorithms with reinforcement learning},
  journal = {Nature},
  volume  = {610},
  number  = {7930},
  pages   = {47--53},
  year    = {2022},
  doi     = {10.1038/s41586-022-05172-4}
}
```

Headline result: `R_{F₂}(⟨4,4,4⟩) ≤ 47` (improves Strassen²=49 over
GF(2)), plus `~50` non-cubic improvements (e.g. `R_{F₂}(⟨4,5,5⟩) ≤ 76`).
**`F₂`-specific** — does *not* lift to `R` (still 49 there). Open-access
PDF (free): [Nature open](https://www.nature.com/articles/s41586-022-05172-4.pdf).
**Local PDF**: [`references/papers/alphatensor_2022_fawzi_nature.pdf`](references/papers/alphatensor_2022_fawzi_nature.pdf)
(Nature open-access; includes Extended Data Table 1 used for AT
provenance audit — see CLAUDE.md "Distinguish discoveries from
re-discoveries").
Source code + data: [github.com/deepmind/alphatensor](https://github.com/deepmind/alphatensor).
Cited in: [RANK_KNOWLEDGE.md](paper/theory/RANK_KNOWLEDGE.md) §1, §3, §4,
[SMALL_MATMUL_CATALOG.md](paper/theory/SMALL_MATMUL_CATALOG.md) §4.3, §4.4, §4.5, §5,
[OMEGA_HISTORY.md](paper/theory/OMEGA_HISTORY.md).

## [13] williams2024

```bibtex
@inproceedings{williams2024,
  author    = {Williams, Virginia Vassilevska and Xu, Yinzhan and Xu, Zixuan and Zhou, Renfei},
  title     = {New Bounds for Matrix Multiplication: from Alpha to Omega},
  booktitle = {Proceedings of the 2024 ACM-SIAM Symposium on Discrete Algorithms (SODA)},
  pages     = {3792--3835},
  year      = {2024},
  doi       = {10.1137/1.9781611977912.134}
}
```

Current best `ω ≤ 2.371552`. Cited in:
[RANK_KNOWLEDGE.md](paper/theory/RANK_KNOWLEDGE.md) §4,
[OMEGA_HISTORY.md](paper/theory/OMEGA_HISTORY.md).

## [14] alphaevolve

```bibtex
@article{alphaevolve,
  title   = {Alpha{E}volve: A coding agent for scientific and algorithmic discovery},
  author  = {Novikov, Alexander and V\~{u}, Ng\^{a}n and Eisenberger, Marvin and Dupont, Emilien and Huang, Po-Sen and Wagner, Adam Zsolt and Shirobokov, Sergey and Kozlovskii, Borislav and Ruiz, Francisco J. R. and Mehrabian, Abbas and Kumar, M. Pawan and See, Abigail and Chaudhuri, Swarat and Holland, George and Davies, Alex and Nowozin, Sebastian and Kohli, Pushmeet and Balog, Matej},
  year    = {2025},
  journal = {arXiv preprint arXiv:2506.13131}
}
```

`R_C(⟨4,4,4⟩) ≤ 48` over complex matrices — *"first improvement, after
56 years, over Strassen's algorithm in this setting"*. **`C`-specific**:
does NOT improve `R(⟨4,4,4⟩)` over `R` (still 49).
**Local PDF**: [`references/papers/alphaevolve_2025_novikov_arxiv2506.13131.pdf`](references/papers/alphaevolve_2025_novikov_arxiv2506.13131.pdf).
Direct PDF (online):
[arxiv.org/pdf/2506.13131](https://arxiv.org/pdf/2506.13131) · [HTML](https://arxiv.org/html/2506.13131). Companion data:
[alphaevolve_results/mathematical_results.ipynb](https://github.com/google-deepmind/alphaevolve_results/blob/main/mathematical_results.ipynb)
(imported by `tools/import_alphaevolve.py`). Cited in:
[RANK_KNOWLEDGE.md](paper/theory/RANK_KNOWLEDGE.md) §1,
[SMALL_MATMUL_CATALOG.md](paper/theory/SMALL_MATMUL_CATALOG.md) §4.3, §4.5,
[OMEGA_HISTORY.md](paper/theory/OMEGA_HISTORY.md).

## [15] wang26

```bibtex
@article{wang26,
  author  = {Wang, Chengu},
  title   = {Automated Lower Bounds for Small Matrix Multiplication Complexity over Finite Fields},
  journal = {arXiv preprint arXiv:2603.07280},
  year    = {2026}
}
```

`R_{F₂}(⟨3,3,3⟩) ≥ 20` — tightens [9] over GF(2). Uses orbit
classification + dynamic programming + verifiable proof certificates;
**Local PDF**: [`references/papers/wang_2026_F2_3x3x3_lb20_arxiv2603.07280.pdf`](references/papers/wang_2026_F2_3x3x3_lb20_arxiv2603.07280.pdf).

methodology directly aligns with this codebase's SAT pipeline. Direct
PDF: [arxiv.org/pdf/2603.07280](https://arxiv.org/pdf/2603.07280) · [HTML](https://arxiv.org/html/2603.07280).
Cited in: [RANK_KNOWLEDGE.md](paper/theory/RANK_KNOWLEDGE.md) §1, §1.2bis, §5,
[SMALL_MATMUL_CATALOG.md](paper/theory/SMALL_MATMUL_CATALOG.md) §4.2, §7.1.

## [16] phialsbasement

```bibtex
@misc{phialsbasement,
  author       = {{PhialsBasement}},
  title        = {{AlphaEvolve-MatrixMul-Verification}},
  howpublished = {\url{https://github.com/PhialsBasement/AlphaEvolve-MatrixMul-Verification}},
  year         = {2025},
  note         = {Accessed: 2026-05-26}
}
```

Independent verification of [14]'s `⟨4,4,4⟩=48` over `C`. README has
the factor matrices in a near-Perminov format — useful as a
cross-check source. Cited in:
[SMALL_MATMUL_CATALOG.md](paper/theory/SMALL_MATMUL_CATALOG.md).

## [21] sedoglavic17

```bibtex
@techreport{sedoglavic17,
  author       = {Sedoglavic, Alexandre},
  title        = {A non-commutative algorithm for multiplying $(7 \times 7)$ matrices using 250 multiplications},
  number       = {hal-01572046},
  institution  = {HAL},
  year         = {2017},
  month        = aug,
  url          = {https://hal.science/hal-01572046v2}
}
```

`R(⟨7,7,7⟩) ≤ 250` via **padded Kronecker** (= trilinear aggregation):
embed `⟨7,7,7⟩` into `⟨8,8,8⟩` Strassen³ and drop multiplications whose
input factors are entirely zero-padded. The construction technique is
the foundation of the entire [\[17\]](REFERENCES.md#17-fmm-lille)
fmm-lille catalog. Direct PDF (HAL):
[hal.science/hal-01572046v2](https://hal.science/hal-01572046v2/file/RFC1708.pdf).
**Local copy**: [`references/papers/sedoglavic_2017_7x7x7_r250.pdf`](references/papers/sedoglavic_2017_7x7x7_r250.pdf).
Cited in: [TRILINEAR_AGGREGATION.md](paper/theory/TRILINEAR_AGGREGATION.md),
[REFERENCES.md](REFERENCES.md) §17.

## [17] fmm-lille

```bibtex
@misc{fmm-lille,
  title        = {Collection of fast matrix multiplication algorithms},
  howpublished = {\url{https://fmm.univ-lille.fr/}},
  organization = {Universit\'e de Lille},
  note         = {Accessed: 2026-05-27}
}
```

Active web-hosted catalog (5,426 entries) of fast matmul algorithms from
`⟨2,2,2⟩` up to `⟨32,32,32⟩` — substantially broader than Perminov's
GitHub catalog. Includes recent 2025 results: AlphaEvolve, Kauers & Wood,
Perminov. No structured download (browse HTML); catalog construction
methodology in Sedoglavic (2017, hal-01572046). Cited in:
[RANK_KNOWLEDGE.md](paper/theory/RANK_KNOWLEDGE.md),
[SMALL_MATMUL_CATALOG.md](paper/theory/SMALL_MATMUL_CATALOG.md) §2.

## [19] schwartz-zwecher25

```bibtex
@article{schwartz-zwecher25,
  author  = {Schwartz, Oded and Zwecher, Eyal},
  title   = {Towards Faster Feasible Matrix Multiplication by Trilinear Aggregation},
  journal = {arXiv preprint arXiv:2508.01748},
  year    = {2025}
}
```

`O(n^{2.773203})` for **practically feasible** matrix multiplication
(small base cases), improving Pan 1982's `O(n^{2.773372})`. Claims
"the fastest matrix multiplication algorithm with base case smaller
than 1000"; best asymptotic complexity for many small base cases
starting at `n_0 = 28`. Technique: trilinear aggregation + de Groote
equivalence + sparse decomposition. Direct PDF:
[arxiv.org/pdf/2508.01748](https://arxiv.org/pdf/2508.01748) · [HTML](https://arxiv.org/html/2508.01748).
**Local PDF**: [`references/papers/schwartz_zwecher_2025_feasible_matmul_arxiv2508.01748.pdf`](references/papers/schwartz_zwecher_2025_feasible_matmul_arxiv2508.01748.pdf).
**Note**:
this is the *feasible* ω line (small-n algorithms that actually run),
distinct from the asymptotic ω race ([\[13\]](REFERENCES.md#13-williams2024))
which requires astronomically large `n`. Cited in:
[RANK_KNOWLEDGE.md](paper/theory/RANK_KNOWLEDGE.md) §4,
[OMEGA_HISTORY.md](paper/theory/OMEGA_HISTORY.md).

## [20] drevet-schost-poster

*Local PDF.*

**Local file:** [`references/papers/drevet_schost_2010_MITACS_poster.pdf`](references/papers/drevet_schost_2010_MITACS_poster.pdf)


```bibtex
@misc{drevet-schost-poster,
  author       = {Drevet, Charles-{\'E}ric and Schost, {\'E}ric},
  title        = {Optimization techniques for small matrix multiplication},
  howpublished = {\url{https://www.cecm.sfu.ca/~pborwein/MITACS/posters/EricEric10.pdf}},
  year         = {2010},
  note         = {MITACS poster; precursor to the journal version published as Drevet--Islam--Schost 2011 ([10]).}
}
```

MITACS poster — the early-stage presentation of the work that became the
[\[10\]](REFERENCES.md#10-drisc09) paper (Drevet–Islam–Schost 2011, journal
version). Useful as a graphical at-a-glance summary of the same per-format
upper-bound table. Cited in:
[SMALL_MATMUL_CATALOG.md](paper/theory/SMALL_MATMUL_CATALOG.md) §2,
[REFERENCES.md](REFERENCES.md) (cross-link with [10]).

## [18] perminov

*Andrew I. Perminov's `FastMatrixMultiplication` community catalog (the
software artefact). The methods and results are written up in the
companion papers [\[75\]](#75-perminov-2026-arxiv),
[\[77\]](#77-perminov-2025-metaflip), [\[78\]](#78-perminov-2025-additive),
[\[79\]](#79-perminov-2025-3x3-a58),
[\[86\]](#86-perminov-2026-serendipitous) — cite the relevant paper, not just the
repo, when a `perminov-*` scheme's rank or addition-count is attributable
to a specific method.*

```bibtex
@misc{perminov,
  author       = {Perminov, Andrew I.},
  title        = {{FastMatrixMultiplication}: A catalog of fast matrix multiplication algorithms},
  howpublished = {\url{https://github.com/dronperminov/FastMatrixMultiplication}},
  note         = {Accessed: 2026-05-26}
}
```

Cite the artefact as **Perminov, Andrew I.** — `dronperminov` is the
GitHub username hosting the work, not a separate author. (Earlier notes
in this repo treated "Andrew I. Perminov" and "Andrey Perminov" as
different people; they are the **same** researcher — "Andrey" is just the
transliteration of the given name.)

Active community catalog of fast matmul schemes up to `⟨16,16,16⟩`,
structured as per-algorithm JSON files. **This is the format `SchemeIO`
reads and writes** — spec at the
[example section](https://github.com/dronperminov/FastMatrixMultiplication#example).
Internal file-name prefixes (e.g. `dronperminov-cr60_…`) and source
labels in `catalog.json` retain the `dronperminov` token to preserve
the provenance trail back to the GitHub artefact; human-facing
references should always say "Perminov".

Cited in: [RANK_KNOWLEDGE.md](paper/theory/RANK_KNOWLEDGE.md),
[SMALL_MATMUL_CATALOG.md](paper/theory/SMALL_MATMUL_CATALOG.md) §2, §4–§5.

---

## [75] <a name="75-perminov-2026-arxiv"></a>perminov-2026-arxiv

```bibtex
@article{perminov2026fast,
  title  = {Fast Matrix Multiplication in Small Formats: Discovering New Schemes with an Open-Source Flip Graph Framework},
  author = {Perminov, Andrew I.},
  journal = {arXiv preprint arXiv:2603.02398},
  url    = {https://arxiv.org/abs/2603.02398},
  year   = {2026}
}
```

The umbrella write-up of the `FastMatrixMultiplication` repo [18]: the
open-source flip-graph framework and the schemes it discovered across
small formats. **This is the paper to cite when our catalog imports a
generic `perminov-*` scheme** whose discovery isn't attributable to one
of the more specific method papers below.

> ⚠️ Earlier revisions of this repo cited `arXiv:2606.02480` for **this**
> umbrella work — that was wrong: the umbrella *Small Formats* paper is
> **arXiv:2603.02398**. `arXiv:2606.02480` is a *separate*, real paper —
> Perminov's June-2026 serendipitous-product work, catalogued as
> [[86]](#86-perminov-2026-serendipitous). Don't conflate the two.

Cited in: [docs/index.html](docs/index.html) (ω &lt; Strassen filter
tooltip), [docs/catalog.js](docs/catalog.js).

---

## [77] <a name="77-perminov-2025-metaflip"></a>perminov-2025-metaflip

```bibtex
@article{perminov2025fast,
  title  = {Fast Matrix Multiplication via Ternary Meta Flip Graphs},
  author = {Perminov, Andrew I.},
  journal = {arXiv preprint arXiv:2511.20317},
  url    = {https://arxiv.org/abs/2511.20317},
  year   = {2025}
}
```

The **ternary meta flip-graph** search method behind many of the rank
results in the catalog. Cite this for `perminov-*` schemes whose *rank*
came from the meta-flip-graph search (as opposed to a low-addition
re-optimisation, which is [78]).

---

## [78] <a name="78-perminov-2025-additive"></a>perminov-2025-additive

```bibtex
@article{perminov2025parallel,
  title  = {Parallel Heuristic Exploration for Additive Complexity Reduction in Fast Matrix Multiplication},
  author = {Perminov, Andrew I.},
  journal = {arXiv preprint arXiv:2512.13365},
  url    = {https://arxiv.org/abs/2512.13365},
  year   = {2025}
}
```

The **additive-complexity reduction** heuristic — the source of the
low-addition (`…_a{N}`) variants of `perminov-*` schemes in our catalog.
Cite this when a `perminov-*` scheme is elected for its *addition count*
rather than (or in addition to) its rank.

---

## [79] <a name="79-perminov-2025-3x3-a58"></a>perminov-2025-3x3-a58

```bibtex
@article{perminov202558,
  title  = {A 58-Addition, Rank-23 Scheme for General 3x3 Matrix Multiplication},
  author = {Perminov, Andrew I.},
  journal = {arXiv preprint arXiv:2512.21980},
  url    = {https://arxiv.org/abs/2512.21980},
  year   = {2025}
}
```

The concrete `R(⟨3,3,3⟩)=23` scheme with only **58 additions** and **no
basis change** — the current addition-count record for general 3×3 over
a commutative ring, beating Stapleton 2025 (60) and
Mårtensson–Stankovski–Stapleton 2025 (59). Also referenced as `[add58]`
in [RANK_KNOWLEDGE.md](paper/theory/RANK_KNOWLEDGE.md). Local PDF:
[`references/papers/perminov_2025_3x3x3_r23_a58_arxiv2512.21980.pdf`](references/papers/perminov_2025_3x3x3_r23_a58_arxiv2512.21980.pdf).

---

## [82] <a name="82-deza-2023-cp"></a>deza-2023-cp

```bibtex
@inproceedings{deza2023fmmcp,
  title     = {Fast Matrix Multiplication Without Tears: A Constraint Programming Approach},
  author    = {Deza, Arnaud and Liu, Chang and Vaezipoor, Pashootan and Khalil, Elias B.},
  booktitle = {29th International Conference on Principles and Practice of Constraint Programming (CP 2023)},
  series    = {LIPIcs},
  volume    = {280},
  pages     = {26:1--26:15},
  year      = {2023},
  doi       = {10.4230/LIPIcs.CP.2023.26},
  url       = {https://arxiv.org/abs/2306.01097}
}
```

A **Constraint Programming** formulation of fast matrix multiplication: the
Brent equations are posed as a CSP over `U,V,W ∈ {−1,0,1}` and solved with CP
Optimizer. Two contributions are directly relevant to this repo:

- **Symmetry-breaking constraints + valid inequalities** — they break the same
  symmetry group we exploit (sign/permutation/orbit equivalence of low-rank
  decompositions, after Strassen's uniqueness [de Groote]). Useful cross-reference
  for our axis-flip / allocation-canonicalisation work and the content-hash
  dedup (`SchemeIO.contentHash`).
- **Infeasibility proofs** — CP can *prove* no rank-R scheme exists for a shape,
  the complement to our heuristic **upper-bound** construction. They reach exact
  results up to ⟨3,3,3⟩=23 and frame the open `19 ≤ R(⟨3,3,3⟩) ≤ 22` question
  (they cite a known `R ≥ 19` lower bound). This is the methodology to cite when
  we discuss *proven-optimal* vs *bound* (optimality discipline).

Field: `{−1,0,1}` (ternary), non-commutative. No new record scheme (3×3 R=23 is
known) — a methods/tooling paper, not a scheme to import. Code:
[khalil-research/Matrix-Mult-CP](https://github.com/khalil-research/Matrix-Mult-CP).
Local PDF:
[`references/papers/deza_2023_matmul_constraint_programming_arxiv2306.01097.pdf`](references/papers/deza_2023_matmul_constraint_programming_arxiv2306.01097.pdf).

**Cited in**: `REFERENCES.md` (this entry). Candidate cross-refs:
`SOLVING_STRATEGIES.md` (CP/SAT exact-search route), the symmetry/orbit notes.

---

## [83] <a name="83-kauers-2026-structure"></a>kauers-2026-structure

```bibtex
@article{kauers2026structure,
  title   = {Exploiting the Structure in Tensor Decompositions for Matrix Multiplication},
  author  = {Kauers, Manuel and Moosbauer, Jakob and Wood, Isaac},
  journal = {arXiv preprint arXiv:2602.11041},
  note    = {Submitted to Elsevier},
  url     = {https://arxiv.org/abs/2602.11041},
  year    = {2026}
}
```

**The published formalisation of the technique we independently built as
"buds" + the "serendipitous product".** A tensor decomposition whose recursive
calls have *special structure* — some recursive calls **share an output**, or an
output is **used in multiple positions** — lets those calls be treated as a
**single matrix multiplication of larger size**. Applied repeatedly, this lowers
the *effective* exponent **below what the tensor rank alone would suggest**,
without reducing the number of multiplications in the base case. Concretely they
improve 6×6 from the Moosbauer–Poole exponent ω 2.8075 to **2.8019** (vs
Strassen's 2.8073), with a reasonable leading coefficient. They connect it to
Schönhage's 1981 ASI special-property, Romani's ASI generalisation, and
Schwartz–Zwecher 2025.

Relationship to this repo (load-bearing — this is direct prior/concurrent art):
- "Recursive calls sharing an output / output in multiple positions → merge into
  a larger MM" **is** our **bud-fusion**: rank-one terms sharing a `u`/`v`/`w`
  vector fuse the inner block into `⟨n, m, k·p⟩` (U-bud), etc. — see
  `SerendipitousBudProduct`, `LineageBudInference`, `BudBaseFactory`, and the
  `#159` serendipitous engine. We also cite Smith 2002 eq. (69) and Perminov's
  draft for the same identity; this Kauers–Moosbauer–Wood paper is the rigorous,
  exponent-level treatment and **must be cited in `paper/` (`#129`)** as the
  state of the art our serendipitous work relates to.
- Their "effective exponent < rank suggests" is exactly our "(rank, buds) Pareto"
  intuition that a bud-richer higher-rank scheme can be the better building block.

Field: general (non-commutative). Methods/exponent paper — not a single scheme to
import, though their improved 6×6 construction is a candidate to replicate. Local
PDF:
[`references/papers/kauers_moosbauer_wood_2026_tensor_structure_arxiv2602.11041.pdf`](references/papers/kauers_moosbauer_wood_2026_tensor_structure_arxiv2602.11041.pdf).

**Cited in**: `REFERENCES.md` (this entry). Candidate cross-refs: the
serendipitous/bud notes, `paper/` related-work, `ROADMAP.md` (serendipitous).

---

## [86] <a name="86-perminov-2026-serendipitous"></a>perminov-2026-serendipitous

```bibtex
@article{perminov2026serendipitous,
  title   = {Meta Flip Graph meets Serendipitous Product: new Fast Matrix Multiplication results},
  author  = {Perminov, Andrew I.},
  journal = {arXiv preprint arXiv:2606.02480},
  url     = {https://arxiv.org/abs/2606.02480},
  year    = {2026}
}
```

Perminov's **June-2026** paper combining the ternary meta-flip-graph search
[[77]](#77-perminov-2025-metaflip) with the **serendipitous product** — the same
bud-fusion identity we build in `SerendipitousBudProduct` / the `#159` engine
(cf. Smith 2002 eq. (69) and the exponent-level Kauers–Moosbauer–Wood treatment
[[83]](#83-kauers-2026-structure)). Reports rank improvements over **207
rectangular formats** up to 16×16×16, 84 newly-ternary formats, and 23 new
schemes with ω &lt; log₂7. **This is the source for the serendipitous 17–32
band** our catalog cites — every `"source": "Perminov 2026 (serendipitous)"`
row in [`docs/cited-bounds.json`](docs/cited-bounds.json) and
[`references/catalogs/perminov-serendipitous-catalog.json`](references/catalogs/perminov-serendipitous-catalog.json)
attributes here (with a `source_scheme_url` pointing at the base scheme's file in
the [`FastMatrixMultiplication`](https://github.com/dronperminov/FastMatrixMultiplication)
repo [18]). The paper is also cited in `paper/` (serendipitous section, after Smith
2002 — `paper/sections/strategies.tex`, key `perminov2026serendipitous`).

> ⚠️ **Not** the same as [[75]](#75-perminov-2026-arxiv) (arXiv:2603.02398, the
> umbrella *Small Formats* framework paper). Earlier repo revisions wrongly
> treated `arXiv:2606.02480` as a mistaken ID for 2603.02398; it is in fact this
> distinct serendipitous-product paper. The two are different works — cite [75]
> for generic `perminov-*` imports, [86] for the serendipitous 17–32 bounds.

**Cited in**: `REFERENCES.md` (this entry, index row [86]), `docs/cited-bounds.json`,
`references/catalogs/perminov-serendipitous-catalog.json`, `paper/refs.bib`
(`perminov2026serendipitous`), `paper/sections/strategies.tex`.

---

## [87] <a name="87-kaporin-2024-brent"></a>kaporin-2024-brent

```bibtex
@article{kaporin2024brent,
  title   = {Semi-analytical solution of Brent equations},
  author  = {Kaporin, I. E.},
  journal = {Doklady Mathematics},
  volume  = {518},
  number  = {1},
  pages   = {29--34},
  year    = {2024},
  doi     = {10.31857/S2686954324040056}
}
```

Kaporin parametrises the Brent equations so that the cyclic (A→B→C) symmetry of
the matmul tensor collapses the unknowns several-fold, then solves the reduced
trilinear/cubic systems numerically. Among the resulting fast-matmul designs are
explicit **complex** schemes **(4,4,4;48)** and **(2,4,5;32)** — "many known
values of rank are reproduced and even improved".

Relevance to this repo:
- **C⟨4,4,4⟩=48** — an explicit, **independently-constructed** complex 48, **a year
  before AlphaEvolve 2025** [[14]](#14-alphaevolve). (Existence of `(4,4,4;r<49)`
  over C was conjectured by Li–Zhang–Ke 2023, arXiv:2310.11686.) This sharpens
  the catalog's canonical field-discipline example: 48/C is **not** original to
  AlphaEvolve. The scheme is **loaded** as
  `src/main/resources/schemes/known/section4/4x4x4-r48-kaporin_2024-8a17320.json`
  (`source: "Kaporin 2024"`, `fields: ["C"]`, `complex: true`), extracted from the
  author's companion verification program and re-verified exactly
  (`Verifier.isExactComplex`, residual ~4.7e-15). It is a **numerical** scheme
  (complex floats, not exact rationals) and **differs from AlphaEvolve's C=48
  original**. Whether it rationalises (à la Dumas–Pernet–Sedoglavic 2025) is open.
- **C⟨2,4,5⟩=32** — a below-Hopcroft–Kerr rank (cf. [[project_below_hk_uses_half_symmetrization]]);
  not yet loaded.

Local copies (downloaded with the user's authorisation):
[`references/papers/kaporin_2024_brent_equations_doklady518.pdf`](references/papers/kaporin_2024_brent_equations_doklady518.pdf)
(Russian full text) and the scheme source
[`references/papers/kaporin_2024_test444r48.for`](references/papers/kaporin_2024_test444r48.for)
(the author's `(4,4,4;48)` verification program, from the paper's ref [8]:
`https://cloud.mail.ru/public/Yfij/ErDxopqBh`). Note: `fmm-lille-biblio.json` also
lists a *different* Kaporin 2024 paper (*Finding complex-valued solutions of Brent
equations using nonlinear least squares*, Comp. Math. & Math. Phys. **64**(9)) —
do not conflate the two.

**Cited in**: `REFERENCES.md` (this entry, index row [87]), the loaded scheme
JSON, `paper/refs.bib` (`kaporin2024brent`), `paper/sections/intro.tex`.

---

## Historical contributors cited in DIS09 (no separate dedicated entry yet)

The entries below are referenced by [10] (DIS09) and surface as
attribution in `docs/cited-bounds.json` rows. They aren't yet
expanded into full numbered entries above — pull requests welcome to
add them as canonical `[N]` entries with full BibTeX.

### [22] waksman70

```bibtex
@article{waksman70,
  author  = {Waksman, A.},
  title   = {On Winograd's Algorithm for Inner Products},
  journal = {IEEE Transactions on Computers},
  volume  = {C-19},
  pages   = {360--361},
  year    = {1970},
  doi     = {10.1109/T-C.1970.222926}
}
```

**Commutative** algorithm family for matmul: rank
`(b · (ac + a + c − 1)) / 2` for `⟨a, b, c⟩` (per DIS09 Table 2). Used
in DIS09's Table 4 commutative-case rankings for `⟨n,n,n⟩`,
`n ∈ {4..10}`. Not useful for recursive matmul (commutative-only),
but a baseline for scalar-only matmul. No factor matrices imported.

### [23] makarov87

*Makarov 1987 — non-commutative `⟨5,5,5⟩=100`.*

```bibtex
@article{makarov87,
  author  = {Makarov, O. M.},
  title   = {A non-commutative algorithm for multiplying ($5\times5$) matrices using one hundred multiplications},
  journal = {USSR Computational Mathematics and Mathematical Physics},
  volume  = {27},
  year    = {1987},
  note    = {Russian original: Zh. Vychisl. Mat. i Mat. Fiz. 27 (1987). Volume/issue/pages [?] — needs cross-check against mathnet.ru.}
}
```

Source for the **non-commutative `⟨5,5,5⟩ = 100`** algorithm cited in
DIS09 Table 3 (labelled "Makarov"). No factor matrices imported.

> **Correction (2026-06)**: earlier revisions of this entry dated the
> scheme `Makarov 1970` with a linear-equations title — that was wrong.
> The non-commutative `⟨5,5,5⟩=100` result is Makarov's **1987** paper;
> there is no relevant 1970 Makarov matmul article. (The entry keeps the
> number [23] for link stability even though it is now out of strict
> chronological order.)
>
> **Disambiguation**: Two distinct Makarov papers appear in this
> bibliography:
> - **Makarov 1987** ([this entry], ⟨5,5,5⟩=100 non-commutative)
> - **Makarov 1986** ([\[29\]](#29-makarov86), ⟨3,3,3⟩=22 commutative)
>
> Both are by O. M. Makarov but address different formats and modes.
> DIS09's labels `Makarov` and `Makarov333` refer to these two papers
> respectively.

### [29] makarov86

*Makarov 1986 — `⟨3,3,3⟩=22`, commutative.*

```bibtex
@article{makarov86,
  author  = {Makarov, O. M.},
  title   = {An algorithm for multiplication of $3 \times 3$ matrices},
  journal = {Zh. Vychisl. Mat. i Mat. Fiz.},
  volume  = {26},
  number  = {2},
  pages   = {293--294, 320},
  year    = {1986}
}
```

Source for the **commutative `⟨3,3,3⟩ = 22`** algorithm cited in DIS09
Table 4 (labelled "Makarov333"). This was the best published
commutative rank for ⟨3,3,3⟩ until [Rosowski 2019/2020](#28-rosowski19)
improved it to 21. The 22-mult result is notable as the first published
rank STRICTLY BELOW Laderman's non-commutative 23 — but only valid for
scalar matmul, not recursive.

🔒 **External PayWall (English translation):** [sciencedirect.com/science/article/abs/pii/004155538690203X](https://www.sciencedirect.com/science/article/abs/pii/004155538690203X)
(USSR Computational Mathematics and Mathematical Physics, 26(2):293–294, 1986).
DOI: [10.1016/0041-5553(86)90203-X](https://doi.org/10.1016/0041-5553(86)90203-X).
🔗 **External Portal (Russian, mathnet.ru):** [mathnet.ru/php/archive.phtml?wshow=paper&jrnid=zvmmf&paperid=4056&option_lang=eng](https://www.mathnet.ru/php/archive.phtml?wshow=paper&jrnid=zvmmf&paperid=4056&option_lang=eng)
(English-language portal entry for the Russian paper, with bibliographic
metadata + access links).
🔗 **External PDF (Russian original):** [mathnet.ru/links/00395e3617837f6ba51655ed9fbc8f2f/zvmmf4056.pdf](https://www.mathnet.ru/links/00395e3617837f6ba51655ed9fbc8f2f/zvmmf4056.pdf)
(*Zh. Vychisl. Mat. i Mat. Fiz.* 26(2):293–294 + 320).
📄 **Local PDF:** [`references/papers/makarov_1986_RU_zvmmf4056.pdf`](references/papers/makarov_1986_RU_zvmmf4056.pdf).
**Factor matrices materialised:**
[`section3/makarov-1986_3x3x3_r22_a105_commutative.json`](src/main/resources/schemes/section3/makarov-1986_3x3x3_r22_a105_commutative.json)
(rank 22, 105 additions). Constructor:
`eu.solven.matmul.papers.makarov1986.Makarov22` /
`eu.solven.matmul.research.MaterializeMakarov1986`.

**Errata in secondary transcriptions:** Islam 2009 MSc thesis §3.3.1
transcribes Makarov's algorithm with a **single-index typo in γ_18**:
the second factor reads `b_{3,2}` where the Russian original (M₁₈ in
the image of `zvmmf4056`) gives `k_6 = b_{2,3}`. This typo causes 4
of 9 output combinations to fail by residual ≈ 4.9 if blindly
transcribed. Cross-checking against the Russian original (linked
above) was required to locate the error. The fix is documented in
`Makarov22.java`'s comment for γ_18.

### [24] probert-fischer-80

```bibtex
@techreport{probert-fischer-80,
  author      = {Probert, R. L. and Fischer, P. C.},
  title       = {Some uses of computers in the development of efficient algorithms for matrix multiplication},
  institution = {University of Waterloo (CS-80-43)},
  year        = {1980}
}
```

Computer-aided synthesis of best-known matmul algorithms for square
sizes up to `⟨40,40,40⟩`. DIS09 Table 3 uses Probert-Fischer 1980 as
its non-commutative baseline. Many of their bounds were beaten by
DIS09 and subsequently by AlphaTensor / AlphaEvolve.
**Local PDF**: [`references/papers/probert_fischer_1980_cs8043.pdf`](references/papers/probert_fischer_1980_cs8043.pdf)
— 43-page scanned tech report from UWaterloo (no OCR; viewable but
text not searchable).

### [25] smith-2002

Computer-aided exploration of small matmul algorithms (cubic and
non-cubic, up to `⟨28,28,28⟩`). DIS09 Table 3 cites Smith 2002 as a
second non-commutative baseline. Title: *"Fast matrix multiplication
formulae — report of the prospectors"*, Warren D. Smith, 2002 (tech
report). **PDF held by the user, NOT committed (copyright).**

**§9 = our composition operators.** Smith §9.1 "How to combine
constructions" enumerates exactly the operators this repo implements:
*hexality* `Rk⟨a,b,c⟩=Rk⟨b,c,a⟩` (our S₃ symmetry / shape
canonicalisation), *additive combination*
`Rk⟨a,b,c+k⟩ ≤ Rk⟨a,b,c⟩+Rk⟨a,b,k⟩` (axis concatenation), *zero padding*
`Rk⟨a,b,c⟩ ≤ Rk⟨a,b,c+1⟩` (peeling), and *multiplicative/recursive*
`Rk⟨aA,bB,cC⟩ ≤ Rk⟨a,b,c⟩·Rk⟨A,B,C⟩` (Kronecker). §9.3 adds
*"serendipitous equalities"* — eq. (69)
`Rk⟨3a,3b,3c⟩ ≤ 19·Rk⟨a,b,c⟩ + 2·Rk⟨2a,b,c⟩` from the
Johnson–McLoughlin `⟨3,3,3⟩=23` scheme (→ `⟨9,9,9⟩≤527`) — our **buds**
(`SerendipitousBudProduct`, Perminov Def 2.12 generalises it).

### [28] rosowski19

*Rosowski 2019/2020 — commutative bounds.*

```bibtex
@misc{rosowski19,
  author       = {Rosowski, Andreas},
  title        = {Fast Commutative Matrix Algorithm},
  howpublished = {arXiv:1904.07683 (v2: 2020)},
  year         = {2020}
}
```

**Commutative** `R_c(⟨3,3,3⟩) ≤ 21` — strictly improves Makarov's 22
and beats every other commutative cubic result published before. Also
gives closed-form generic formulas for `⟨l, n, m⟩` over a commutative
ring:

- **n odd, m odd**: `R_c(⟨l,n,m⟩) ≤ n(lm + l + m − 1) / 2`
- **n odd, m even**: `R_c(⟨l,n,m⟩) ≤ (n(lm + l + m − 1) + l − 1) / 2`
- For even n: Rosowski gives a divisions-free variant.

These formulas reduce Waksman/Islam (used in DIS09 Table 4) by `lm` in
both odd and even cases — material savings for `n ≥ 3`. The paper also
introduces "approximate non-bilinear algorithms" — non-bilinear
recursive schemes that beat bilinear for small `n < 20`. Notable
example: `⟨5,5,5⟩` approximate non-bilinear with **89 multiplications**.

These results are **commutative-only** (don't lift to recursive matmul
over non-commutative entries). They're reflected in
`docs/cited-bounds.json` and surfaced in Pages with the `cmt` tag.

**Local PDF**: [`references/papers/rosowski_2019_commutative_matmul_arxiv1904.07683.pdf`](references/papers/rosowski_2019_commutative_matmul_arxiv1904.07683.pdf).
**Algorithms extracted** to
[`references/rosowski-algorithms.md`](references/rosowski-algorithms.md):
Algorithm 1 (⟨n,3,3⟩ in 6n+3) and Corollary 1 (⟨3,3,3⟩ in 21
multiplications), with full product list. These are NON-BILINEAR
(products mix A and B entries in both factors) and so don't round-trip
through standard `SchemeIO` — see the markdown for why and what a
future non-bilinear scheme format would look like.
Cited in: [RANK_KNOWLEDGE.md](paper/theory/RANK_KNOWLEDGE.md) §1.2.

### [27] mss25

*Mårtensson–Stankovski Wagner–Stapleton 2025.*

```bibtex
@misc{mss25,
  author       = {M{\aa}rtensson, Erik and {Stankovski Wagner}, Paul and Stapleton, Joshua},
  title        = {A Rank 23 Algorithm for Multiplying 3 × 3 Matrices with an Arithmetic Complexity of 59},
  howpublished = {arXiv preprint arXiv:2601.05272},
  year         = {2025},
  note         = {Submitted Dec 18, 2025; technical companion to be published.}
}
```

`⟨3,3,3⟩` with **rank 23 + 59 additions**, **without basis change** —
current best addition-count for rank-23 `⟨3,3,3⟩` algorithms in the
no-basis-change setting. Combines methods of Mårtensson-Stankovski
Wagner (62 adds) and Stapleton (60 adds). Note that Karstadt-Schwartz
achieve 61 adds with basis change, and earlier Schwart et al. went
lower with basis change too — the no-basis-change line is the practical
constraint here. Joins the lineage of recent addition-count records
([add58], [add60]) tracked in
[RANK_KNOWLEDGE.md](paper/theory/RANK_KNOWLEDGE.md) §6.
**Local PDF**: [`references/papers/martensson_stankovski_stapleton_2025_3x3x3_r23_a59_arxiv2601.05272.pdf`](references/papers/martensson_stankovski_stapleton_2025_3x3x3_r23_a59_arxiv2601.05272.pdf).
Cited in: [RANK_KNOWLEDGE.md](paper/theory/RANK_KNOWLEDGE.md) §6.

### [26] mezzarobba-2007

```bibtex
@mastersthesis{mezzarobba07,
  author = {Mezzarobba, Marc},
  title  = {G\'en\'eration automatique de proc\'edures num\'eriques
            pour les fonctions D-finies},
  school = {Master parisien de recherche en informatique (MPRI)},
  year   = {2007},
  note   = {Available at \url{https://marc.mezzarobba.net/m2/Mezzarobba_MScThesisMPRI2007-1.2.pdf}}
}
```

MPRI master's thesis. Although the topic is automatic generation
of numerical procedures for D-finite (holonomic) functions, an
appendix tabulates **commutative** small-matmul ranks for square
sizes up to `⟨28,28,28⟩`, which is the table DIS09 §4 (Table 4)
cites as its commutative baseline.
**Fig. 3** (Waksman ⟨3,3,3⟩) gives the explicit `r=23`
commutative formula via the `(a+b)(c+d) - (a-b)(c-d) = 2(ac+bd)`
identity (requires `1/2 ∈ K`, so the field must be of
characteristic ≠ 2).
**Caveat**: the 2007 table apparently did NOT include Makarov
1986's commutative `⟨3,3,3⟩=22` (it lists 23 instead), so DIS09
Table 4 marks the `22` cell as "new" — but the actual discovery
is Makarov 1986, not DIS09. (See attribution logic in
`GenerateCitedBounds.attributionLabel`.)
**Local PDF**: [`references/papers/mezzarobba_2007_mthesis_dfinite.pdf`](references/papers/mezzarobba_2007_mthesis_dfinite.pdf).
Cited in: `docs/cited-bounds.json` (Waksman ⟨3,3,3⟩=23 cmt),
[SMALL_MATMUL_CATALOG.md](paper/theory/SMALL_MATMUL_CATALOG.md).

---

## [30] solven-strassen

*This repository (WIP).*

```bibtex
@misc{solven-strassen,
  author       = {Lacelle, Benoit},
  title        = {strassen: a research-grade catalog of fast matrix multiplication algorithms},
  howpublished = {\url{https://github.com/solven-eu/matmulcatalog}},
  year         = {2026},
  note         = {Work in progress; see \url{https://www.solven.eu/matmulcatalog/} for the live catalog.}
}
```

This repository — a self-citation tracked as a numbered reference so
that derived/cited bounds produced by our own pipeline can be
attributed back without ambiguity (vs older sources like DIS09 [10] or
fmm-lille [17]).

**What's in scope** (as of 2026-05-28):
- Full catalog ingestion from [11] Smirnov, [12] AlphaTensor, [14] AlphaEvolve,
  [17] fmm-lille, [18] Perminov / FastMatrixMultiplication, with provenance
  preserved.
- DIS09 reproduction pipeline (`BlockSplitSearch.findBestMultiBaseSplit`)
  layered as:
  - **Layer 1** — multi-base recursion: Strassen, Laderman, axis-split
    bases (mul211 / mul121 / mul112).
  - **Layer 3** — S₃ symmetry transforms of each base
    (`SymmetryTransforms.s3Orbit`).
  - SOTA resolver pulls inner sub-shape ranks from the **modern catalog**
    (Sedoglavic [21], AlphaTensor, AlphaEvolve), so the reproduction is
    "DIS09 algorithm × 2026 catalog".
- Layers 4–5 from DIS09 §3 (sub-product pairing, full Pan TA
  parameterisation) are **not yet implemented** — this is why our
  reproduction loses to DIS09 for `n ≥ 18` where DIS09 attributes the
  best rank to "TA".

**Concrete results** (NC `⟨n,n,n⟩` over R/Q/Z, vs DIS09 Table 3):
- 12 strict improvements at `n ∈ {5, 7, 9, 10, 11, 12, 13, 15, 16, 17, 19, 23}`.
- 3 ties at `n ∈ {4, 6, 8}`.
- 12 losses at `n ∈ {14, 18, 20, 21, 22, 24, 25, 26, 27, 28, 29, 30}` —
  all but `n = 14` and `n = 21` correspond to rows where DIS09 used Pan
  TA (Layers 4–5).

See `src/test/java/io/cormoran/strassen/v3/catalog/TestDIS09FullScan.java`
for the runnable comparison. **No publication yet** — entries that
strictly improve DIS09 are tagged with this `[30]` self-reference in
`docs/cited-bounds.json` so a future write-up has clean attribution.

Cited in: `docs/cited-bounds.json`, `docs/derived-bounds.json`,
[NEW_BOUNDS.md](research/NEW_BOUNDS.md),
[DIS09_REPRODUCTION_PLAN.md](research/trilinear-aggregation/DIS09_REPRODUCTION_PLAN.md).

---

## [31] schachtel78

*Schachtel 1978 — non-commutative ⟨5,5,5⟩=103.*

```bibtex
@article{schachtel78,
  author  = {Schachtel, G.},
  title   = {A non-commutative algorithm for multiplying $5 \times 5$ matrices using 103 multiplications},
  journal = {Information Processing Letters},
  volume  = {7},
  number  = {4},
  pages   = {180--182},
  year    = {1978},
  doi     = {10.1016/0020-0190(78)90069-X}
}
```

The original 1978 non-commutative {@code ⟨5,5,5⟩=103} result.
Notable as the first non-trivial improvement over the naive 125 for
5×5 matmul over a general ring, predating Smirnov 2013 (also 103) and
Makarov 1987 (100 — also NC but Russian original is hard to access).
No factor matrices imported (paywalled IPL article); registered as a
cited bound in `docs/cited-bounds.json`. Cited in:
[SMALL_MATMUL_CATALOG.md](paper/theory/SMALL_MATMUL_CATALOG.md) §5.

---

## [65] bcs97

The canonical textbook treatment of algebraic complexity, including
the most thorough exposition of fast matrix multiplication outside
the primary papers themselves.

```bibtex
@book{bcs97,
  author    = {B{\"u}rgisser, Peter and Clausen, Michael and Shokrollahi, Mohammad Amin},
  title     = {Algebraic Complexity Theory},
  series    = {Grundlehren der mathematischen Wissenschaften},
  volume    = {315},
  publisher = {Springer},
  year      = {1997},
  doi       = {10.1007/978-3-662-03338-8}
}
```

Relevant for this catalog: §14–16 reproduce Strassen, Pan
trilinear-aggregation, and partial-matrix-multiplication (Schönhage)
constructions with explicit factor decompositions — often the
clearest source when the original papers are paywalled or terse.

> **Correction (2026-06)**: an earlier version of this note also claimed BCS97
> "reproduces Hopcroft-Kerr". It does **not** carry the HK `⟨a,2,c⟩` *bilinear*
> rank formula. Its 2×n results (Ex. 14.20, Notes §14.9, p. 372–374) are the
> **commutative** Waksman/Winograd bound `L(⟨e,2,e⟩) ≤ e²+2e−1` and Feig's
> `L(⟨n,2,2⟩) ≤ 3n+1`; it attributes matmul-rank S₃-duality to
> Hopcroft–Musinski [253], and cites Hopcroft–Kerr [252] only for the
> substitution method (lower bounds, ≈ pp. 143–160). For the HK rank formula
> itself use [\[81\]](#81-hopcroft-musinski-1973) (Hopcroft–Musinski) or
> [\[80\]](#80-alekseev-smirnov-2013) (Alekseev–Smirnov).

**Local PDF**: held privately outside this repo at
`Dropbox/Solven/JEI/Dossier R&D/Matrix Multiplication/1997 Algebraic
Complexity Theory.pdf` (Springer copyright — must NOT be added to
the repository).

Cited in: REFERENCES.md (as fallback source for HK §3 construction,
Pan TA closed forms, and Pan pair-product recipe).

---

## [66] pan2014

Victor Pan's own survey of his trilinear aggregation (TA) techniques
and APA (any-precision-approximate) algorithms — likely the clearest
free source for Pan's TA constructions used implicitly by DIS09.

```bibtex
@article{pan2014,
  author  = {Pan, Victor Y.},
  title   = {Matrix Multiplication, Trilinear Decompositions, APA Algorithms, and Summation},
  journal = {arXiv preprint},
  number  = {1412.1145},
  year    = {2014},
  eprint  = {1412.1145},
  archivePrefix = {arXiv},
  primaryClass = {cs.CC}
}
```

**Local PDF**: [`references/papers/pan_2014_trilinear_apa_arxiv1412.1145.pdf`](references/papers/pan_2014_trilinear_apa_arxiv1412.1145.pdf).

Relevance: Pan's TA cubic closed form `(n³+12n²+11n)/3` (and the
`(n³+15n²+14n−6)/3` odd-n branch) is implemented as
{@code PanTrilinearAggregation.cubicBound}, but our constructor for
the matching scheme files is still TODO (task #42 / #44). This paper
should contain the explicit recipe for both the cubic aggregation and
the pair-product trick.

---

## [67] pan2014b

Pan's companion piece to [66] — a 6-page historical narrative
explaining what he considers the "void" in the literature about
trilinear aggregation, with concrete recipes for small disjoint MM
cases (which is exactly what the pair-product technique needs).

```bibtex
@misc{pan2014b,
  author        = {Pan, Victor Y.},
  title         = {Better Late Than Never: Filling a Void in the History of Fast Matrix Multiplication and Tensor Decompositions},
  year          = {2014},
  eprint        = {1411.1972},
  archivePrefix = {arXiv},
  primaryClass  = {cs.CC}
}
```

**Local PDF**: [`references/papers/pan_2014b_history_arxiv1411.1972.pdf`](references/papers/pan_2014b_history_arxiv1411.1972.pdf).

---

## [68] islam2009

Master's thesis of Md. Nazrul Islam (Western Ontario, supervisor: Éric
Schost). Chapter 6 generalises Waksman's commutative algorithm to
arbitrary {@code ⟨m,n,p⟩}; Chapters 4–5 contain the addition-aware
peeling-vs-padding analysis and the trilinear-aggregation reformulation
that became the published DIS09 §3 material.

```bibtex
@mastersthesis{islam2009,
  author = {Islam, Md. Nazrul},
  title  = {Optimization Techniques For Matrix Multiplication},
  school = {The University of Western Ontario},
  year   = {2009},
  type   = {{Master of Science Thesis}},
  note   = {Supervisor: Éric Schost. 96 pages.}
}
```

**Local PDF**: [`references/papers/islam_2009_msc_optim_matmul.pdf`](references/papers/islam_2009_msc_optim_matmul.pdf).

**Local PDF**: see line above. Note unlike BCS 1997 [65], Islam's
thesis is openly distributable (UWO MSc theses are public via
Scholarship@Western); we archive it locally for self-contained
reproducibility.

**Proposition 3** (Chapter 6, generalized Waksman):
- ⟨m, n, p⟩ commutative, {@code n} even, {@code n' = n/2}:
  rank ≤ `mn'p + mn' + n'p − n'`
- {@code n} odd, {@code n' = (n−1)/2}:
  rank ≤ `mn'p + mn' + n'p + pm − n'`

Both match Waksman for square ⟨n,n,n⟩, but Islam's presentation extends
to arbitrary rectangular shapes. Superseded for {@code ⟨n,3,3⟩} by
Rosowski 2019 (`6n+3`), but worth registering as the historical first
non-square commutative bound.


---

## [69] dumas-pernet-sedoglavic-2025

Rational-coefficient (non-complex) algorithm for `Q⟨3,4,7⟩:m=63`.
Derived from AlphaEvolve's complex-coefficient `C⟨3,4,7⟩:m=63`
([14]) by applying an isotropy from the Klein-four stabiliser of the
original tensor decomposition; the resulting orbit point lives over Q
(coefficients in {0, ±1, ±1/2}), removing the dependency on a square
root of -1 while preserving tensor rank. Closes the last non-cubic gap
versus FMM-Lille at ⟨3,4,7⟩.

```bibtex
@unpublished{dumas-pernet-sedoglavic-2025,
  title  = {A non-commutative algorithm for multiplying a $3\times 4$ matrix by a $4\times 7$ matrix using 63 non-complex multiplications},
  author = {Dumas, Jean-Guillaume and Pernet, Cl{\'e}ment and Sedoglavic, Alexandre},
  year   = {2025},
  month  = {July},
  note   = {HAL preprint hal-05121550v2, submitted 15 Jul 2025},
  url    = {https://hal.science/hal-05121550v2}
}
```

**Local PDF**: [`references/papers/dumas_pernet_sedoglavic_2025_3x4x7_r63_hal05121550.pdf`](references/papers/dumas_pernet_sedoglavic_2025_3x4x7_r63_hal05121550.pdf).
**Scheme JSON**: [`src/main/resources/schemes/section7/dumas-pernet-sedoglavic-2025_3x4x7_r63_a588.json`](src/main/resources/schemes/section7/dumas-pernet-sedoglavic-2025_3x4x7_r63_a588.json) — 63 rank-one triads, 588 additions, field Q (1/2 required). Imported from FMM-Lille's Maple mirror at [`https://fmm.univ-lille.fr/3x4x7_tensor.mpl.bz2`](https://fmm.univ-lille.fr/3x4x7.html) via `tools/import_fmm_maple.py` + manual rename / attribution. Verified by `TestDPS2025_3x4x7`.

**Attribution**: the rank `C⟨3,4,7⟩=63` is due to AlphaEvolve
([14]); DPS-2025 contributes the Q-coefficient realisation, not a
new rank record. Catalog JSON sets `"discovery": true` for the
rational decomposition specifically, with
`"attribution_for_rank"` pointing back to AlphaEvolve.


---

## [70] <a name="70-smirnov2017"></a>smirnov2017

```bibtex
@techreport{smirnov2017,
  author      = {Smirnov, Alexey Vladimirovich},
  title       = {Several Bilinear Algorithms for Matrix Multiplication Problems $\langle 3, P, Q\rangle$},
  institution = {Self-published (ResearchGate / Academia.edu preprint)},
  year        = {2017},
  month       = jan,
  note        = {ResearchGate publication ID 313064941; also mirrored at academia.edu/36266758. Russian Federal Center of Forensic Examination (Moscow).},
  url         = {https://www.researchgate.net/publication/313064941_Several_bilinear_algorithms_for_matrix_multiplication}
}
```

Follow-up to [[11] Smirnov 2013](#11-smirnov), focused on rectangular
shapes ⟨3, P, Q⟩. Lists **explicit non-commutative bilinear algorithms
over ℤ** for the formats:

| shape | rank claimed |
| --- | --: |
| `R⟨3,4,6⟩` | 56 |
| `R⟨3,5,5⟩` | 58 |
| `R⟨3,4,7⟩` | 66 |
| `R⟨3,4,8⟩` | 75 |
| `R⟨3,5,7⟩` | 82 |

The paper applies the same computer-aided search methodology as
Smirnov 2013 (numerical optimisation followed by rational rounding into
the `{-1, 0, +1}` family). The bounds are **non-commutative** and lift
to recursive matmul over any ring.

**Attribution-relevant context**: these are the historical first
appearance of the listed ranks in print (as best the catalog can
verify). Several of them have been improved since:

| shape | Smirnov 2017 | current best NC | improver |
| --- | --: | --: | --- |
| `R⟨3,4,6⟩` | 56 | 54 (AlphaEvolve 2025); 56 also re-found by Moosbauer flip-graph | [14], jakobmoosbauer_flips |
| `R⟨3,5,5⟩` | 58 | 58 (re-proved by Sedoglavic–Smirnov 2021 + AlphaTensor 2022) | [2101.12568](https://arxiv.org/abs/2101.12568) · [HTML](https://arxiv.org/html/2101.12568), [12] |
| `R⟨3,4,7⟩` | 66 | 63 (DPS 2025 over Q; AlphaEvolve 2025 over C); 64 (meta-flip-graph / Perminov) | [69], [14] |
| `R⟨3,4,8⟩` | 75 | 73 (Kauers–Wood meta-flip-graph 2025); 74 (AlphaEvolve 2025) | [arXiv:2510.19787](https://arxiv.org/abs/2510.19787) · [HTML](https://arxiv.org/html/2510.19787), [14] |
| `R⟨3,5,7⟩` | 82 | 79 (Kauers–Wood meta-flip-graph 2025); 80 (AlphaEvolve 2025) | [arXiv:2510.19787](https://arxiv.org/abs/2510.19787) · [HTML](https://arxiv.org/html/2510.19787), [14] |

Per the project policy ("Distinguish discoveries from re-discoveries"),
on-disk schemes that hit these rank values but were imported under a
later source (AlphaTensor, Perminov, jakobmoosbauer_flips, …) should
have `attribution_for_rank` resolved down to **Smirnov 2017** for the
formats that Smirnov 2017 first reached, NOT to the importer. The
audit table is tracked in
[`references/smirnov_3pq_attribution_audit.md`](references/smirnov_3pq_attribution_audit.md).

**Access**:
- ResearchGate (free, may require login):
  <https://www.researchgate.net/publication/313064941_Several_bilinear_algorithms_for_matrix_multiplication>
- Academia.edu mirror (free, requires login):
  <https://www.academia.edu/36266758/Several_Bilinear_Algorithms_for_Matrix_Multiplication>
- Local PDF: **not archived** — both hosts return HTTP 403 to
  unauthenticated curl/WebFetch. Mark as TODO in the "Local PDF
  archive" §; manual browser download required.

Cited in: [SMALL_MATMUL_CATALOG.md](paper/theory/SMALL_MATMUL_CATALOG.md) (when the
⟨3,P,Q⟩ rows resolve attribution to this paper).


## [71] <a name="71-sedoglavic-smirnov-2021"></a>sedoglavic-smirnov-2021

```bibtex
@unpublished{sedoglavic-smirnov-2021,
  author = {Sedoglavic, Alexandre and Smirnov, Alexey Vladimirovich},
  title  = {The tensor rank of $5\times 5$ matrices multiplication is bounded by 98 and its border rank by 89},
  year   = {2021},
  note   = {Cited as [6] in Smirnov 2017 references; preprint / technical report}
}
```

Two distinct upper bounds at ⟨5,5,5⟩:

| quantity | bound | status |
| --- | --: | --- |
| `R(⟨5,5,5⟩)` (tensor rank, NC) | ≤ 98 | dominated by AlphaTensor 96 (2022) and AlphaEvolve 93 (2025) |
| `R̃(⟨5,5,5⟩)` (border rank, NC) | ≤ 89 | **still SOTA** as of 2026 |

Border rank `R̃` is the limit of rank as a parametric scheme's
auxiliary indeterminate `x → 0`; valid algorithms admit
Laurent-polynomial coefficients that cancel in the limit. Only
relates to exact rank `R` via the Schönhage τ-theorem at the
asymptotic level — so border-rank improvements are tracked but
do not directly close exact-rank catalog rows.

**Access**: pointed to by Smirnov 2017 references; haven't located a
public-host link as of 2026-06. Self-published preprint, similar
hosting profile to Smirnov 2017.

Cited in: [docs/cited-bounds.json](docs/cited-bounds.json) under the
⟨5,5,5⟩ shape; also surfaces a `kind: "border"` companion row tagged
`border_rank: 89`.


## [72] <a name="72-sedoglavic-2017-fmm-methodology"></a>sedoglavic-2017-fmm-methodology

```bibtex
@unpublished{sedoglavic-2017-fmm-methodology,
  author = {Sedoglavic, Alexandre},
  title  = {A non-commutative algorithm for multiplying $(7\times 7)$ matrices using 250 multiplications},
  year   = {2017},
  month  = dec,
  note   = {HAL preprint hal-01572046v2, last revised 2019-01-18. THE methodology paper behind the FMM-Lille catalog for cubic decompositions.},
  url    = {https://hal.science/hal-01572046v2}
}
```

**Why this entry matters**: Sedoglavic's FMM-Lille website
(<https://fmm.univ-lille.fr/>) answers "How this table was made?" with
a pointer to this paper. It is THE methodology behind the cubic-shape
decompositions in the FMM catalog.

**Proposition 1** (the entire machinery in one line):

> `⟨u+v, u+v, u+v⟩ ≤ ⟨u,u,u⟩ + 3·⟨u,u,v⟩ + 3·⟨v,v,u⟩` when `u > v`.

Plug in (u,v) and read off bounds — the closed form generates many
of FMM's cubic decompositions in one shot:

| (u, v) | target | implied rank | match? |
| --- | --- | --: | --- |
| (4, 3) | ⟨7, 7, 7⟩ | 250 with Strassen² ⟨4,4,4⟩=49 (paper title bound) | ✓ |
| (4, 3) | ⟨7, 7, 7⟩ | 249 with AlphaEvolve ⟨4,4,4⟩=48 (current Q) | ✓ |
| (5, 4) | ⟨9, 9, 9⟩ | 520 with Smirnov 2013 ⟨5,5,4⟩=70 (cited in paper) | ✓ |
| (6, 5) | ⟨11, 11, 11⟩ | **873** = ⟨6,6,6⟩=153 + 3·⟨6,6,5⟩=130 + 3·⟨5,5,6⟩=110 | ✓ FMM-Lille bound |
| (9, 8) | ⟨17, 17, 17⟩ | 2940 with current Q atoms (matches our solven-closure-2026) | ✓ |
| (10, 9) | ⟨19, 19, 19⟩ | 4044 (matches catalog) | ✓ |

**Doubling extension** (paper extends to u = v = k via Pan TA pairing):

> `⟨2k, 2k, 2k⟩ ≤ ⟨k,k,k⟩ + 3·pair_cost(k,k,k)`, where `pair_cost(k,k,k) = k³ + 3k²`.

| k | target | implied rank | match? |
| --- | --- | --: | --- |
| 7 | ⟨14,14,14⟩ | 249 + 3·490 = **1719** | ✓ matches catalog (solven-strassen-2026) |
| 11 | ⟨22,22,22⟩ | 873 + 3·1694 = 5955 | ✓ |

**Implementation status**: encoded as identities in
[`KnownTauIdentities.java`](src/main/java/eu/solven/matmul/search/KnownTauIdentities.java)
(verified against the live catalog). The generic `SedoglavicProp1`
constructive method that enumerates all (u, v) pairs at a target
shape and emits the minimum-predicted rank with explicit lineage is
queued as part of #161.

**Access**: HAL open archive, public domain CC0 — local PDF at
`references/papers/sedoglavic_2017_fmm_methodology_hal01572046.pdf`.

Cited in: [STRATEGIES.md](STRATEGIES.md), `KnownTauIdentities` Javadoc.


## [73] <a name="73-burichenko-2014"></a>burichenko-2014

```bibtex
@misc{burichenko-2014,
  author = {Burichenko, Vladimir P.},
  title  = {On symmetries of the Strassen algorithm},
  year   = {2014},
  eprint = {1408.6273},
  archivePrefix = {arXiv},
  primaryClass = {math.GR},
  url    = {https://arxiv.org/abs/1408.6273}
}
```

Computes the **discrete stabilizer** of Strassen's ⟨2,2,2⟩=7 algorithm
under the natural symmetry group on bilinear algorithms. The stabilizer
has order 36 — isomorphic to S₃ × Z₂³ (i.e. permutation of the three
axes × independent axis-flips). Useful for:

1. Characterising rank-7 ⟨2,2,2⟩ orbits: every Strassen-equivalent
   scheme decomposes into orbits of the order-36 group acting on the
   tensor. This bounds how many "distinct" rank-7 schemes can exist
   up to symmetry.
2. The complement question of whether Strassen's orbit is the ONLY
   rank-7 orbit (the "uniqueness of Strassen" question) is touched
   on but not fully resolved here — De Groote 1978a/b is the source
   for that, and Heun 1994 catalogues normal forms.

**Access**: arXiv:1408.6273 — local PDF at
`references/papers/burichenko_2014_strassen_symmetries_arxiv1408.6273.pdf`.

Cited in: [docs/notes/enumerating-rank-7-2x2x2-schemes.md](docs/notes/enumerating-rank-7-2x2x2-schemes.md),
[docs/notes/heun-1994.md](docs/notes/heun-1994.md).


## [74] <a name="74-landsberg-2008"></a>landsberg-2008

```bibtex
@article{landsberg-2008,
  author = {Landsberg, J. M.},
  title  = {Geometry and the complexity of matrix multiplication},
  journal = {Bulletin of the American Mathematical Society},
  volume = {45},
  number = {2},
  pages  = {247--284},
  year   = {2008},
  doi    = {10.1090/S0273-0979-08-01176-2},
  url    = {https://pubs.ams.org/journals/bull/2008-45-02/S0273-0979-08-01176-2/}
}
```

Survey article: geometric perspective on the complexity of matrix
multiplication. Key concrete bounds (Theorem 3.8.4, attributed to
**Markus Bläser**):

> **Bläser**: `R(M_{m,m,m}) ≥ (5/2)·m² − 3m`.
>
> A sharper result of Bläser: `R(M_{3,3,3}) ≥ 19`.

Numerical values from the generic Bläser bound:

| m | (5/2)m² − 3m | ceil | known UB | gap |
| --- | --: | --: | --: | --: |
| 2 |  4 |  4 | 7 | 3 |
| 3 | 13.5 | **14** (sharper: 19 per Bläser) | 23 (Laderman) | 4 |
| 4 | 28 | 28 | 47 (AT-F2) / 48 (AE) / 49 (Strassen²) | 19 |
| 5 | 47.5 | 48 | 93 (AE) | 45 |
| 7 | 101.5 | 102 | 249 | 147 |
| 11 | 269.5 | 270 | 873 | 603 |

The gap widens with m, reflecting the famously hard "ω lower-bound"
problem. None of these bounds is tight for m ≥ 3 — the conjectured
ω = 2 corresponds to `R(⟨m,m,m⟩) ~ m²⁺ᵒ⁽¹⁾`, but the proven lower
bounds are `Θ(m²)` with only the constant improving over decades
(Strassen's 1983 `3m² − o(m²)` → Landsberg–Michałek's `3m² − √m³`
in 2018 etc.).

**Access**: Bull. AMS, open access — local PDF at
`references/papers/landsberg_2008_geometry_complexity_matmul_bull_ams.pdf`.

Cited in: [docs/lower-bounds.json](docs/lower-bounds.json) (Bläser's
generic LB for ⟨m,m,m⟩, m = 2..32).


---

## [76] <a name="76-smirnov-2022"></a>smirnov-2022

*Alexey Vladimirovich Smirnov 2022 — non-commutative `⟨4,4,9⟩=104`.*

```bibtex
@techreport{smirnov-2022,
  author      = {Smirnov, Alexey Vladimirovich},
  title       = {Bilinear algorithm {$\langle 4,4,9\rangle = 104$} for matrix
                 multiplication: an irreducibly irrational solution of the Brent system?},
  institution = {ResearchGate (technical report)},
  month       = oct,
  year        = {2022}
}
```

Source for the catalog's `⟨4,4,9⟩ = 104` scheme, which we ingested via the
FMM-Lille Maple mirror (filename `fmm_lille-4x4x9_m104_a1425.json`). Per the
FMM-Lille per-format page <https://fmm.univ-lille.fr/4x4x9.html>, the scheme is
originally **Smirnov 2022**, not an FMM discovery — the on-disk `source` stays
`fmm-lille` (the importer) and the origin is recorded in the scheme's `notes`
(we no longer carry a separate `attribution_for_rank` field).

> **Disambiguation**: there are two distinct Smirnov Oct-2022 ResearchGate
> reports — this `⟨4,4,9⟩=104` one, and a separate `⟨4,4,6⟩=73` one
> (`Smirnov:2022aa` in `references/fmm-lille-biblio.json`). Don't conflate them.
>
> Exact ResearchGate URL/DOI [?] — to cross-check; the FMM-Lille page above is
> the verified pointer.

Cited in: the `⟨4,4,9⟩=104` scheme's `notes` (origin pointer).

---

## [80] <a name="80-alekseev-smirnov-2013"></a>alekseev-smirnov-2013

*Alekseev & Smirnov 2013 (**in Russian**) — exact and border ("approximate")
bilinear complexity of `⟨4,2,2⟩` and `⟨2,2,2⟩`.*

```bibtex
@article{alekseev-smirnov-2013,
  author   = {Alekseev, V. B. and Smirnov, Alexey Vladimirovich},
  title    = {On the exact and approximate bilinear complexities of
              multiplication of {$4\times2$} and {$2\times2$} matrices},
  journal  = {Sovremennye Problemy Matematiki (Contemporary Problems of Mathematics)},
  number   = {17},
  pages    = {135--152},
  year     = {2013},
  language = {russian},
  doi      = {10.4213/spm47},
  note     = {In Russian. Math-Net.Ru: mathnet.ru/spm47. English translation
              in Proceedings of the Steklov Institute of Mathematics
              [volume/pages [?] — verify].}
}
```

Independent confirmation of the **exact** (and the border / "approximate")
bilinear ranks for the inner-dimension-2 formats `⟨4,2,2⟩` and `⟨2,2,2⟩` —
i.e. specific small cases of the Hopcroft-Kerr `⟨a,2,c⟩` family (see the
[Hopcroft-Kerr family](#hopcroft-kerr-family-a2c3acmaxac2-caveat) caveat
above). Useful as a second source whenever we cite the HK closed-form
`⌈(3ac+max(a,c))/2⌉` for an inner-2 shape.

> 🇷🇺 **Language: the paper is in Russian.** The DOI / Math-Net / ResearchGate
> links below all point to the **Russian original**; an English translation
> appears in *Proceedings of the Steklov Institute of Mathematics*.

📄 **Local PDF (Russian)**:
[`references/papers/alekseev_smirnov_2013_RU_4x2_2x2_bilinear_complexity_spm47.pdf`](references/papers/alekseev_smirnov_2013_RU_4x2_2x2_bilinear_complexity_spm47.pdf)
— archived from Math-Net.Ru (freely distributed, like
[\[29\]](#29-makarov86) Makarov 1986).

**Access** (Russian original):
- DOI: [10.4213/spm47](https://doi.org/10.4213/spm47)
- Math-Net.Ru: <https://www.mathnet.ru/eng/spm47>
- ResearchGate (Russian): <https://www.researchgate.net/publication/284529885_On_the_exact_and_approximate_bilinear_complexities_of_multiplication_of_42_and_22_matrices>

Cited in: REFERENCES.md (Hopcroft-Kerr `⟨a,2,c⟩` family caveat).

---

## [81] <a name="81-hopcroft-musinski-1973"></a>hopcroft-musinski-1973

*Hopcroft & Musinski 1973 — S₃-duality of bilinear/matmul rank; states and
propagates the Hopcroft-Kerr `⟨p,2,n⟩` formula.*

```bibtex
@article{hopcroft-musinski-1973,
  author  = {Hopcroft, John E. and Musinski, Jean},
  title   = {Duality applied to the complexity of matrix multiplication and other bilinear forms},
  journal = {SIAM Journal on Computing},
  volume  = {2},
  number  = {3},
  pages   = {159--173},
  year    = {1973},
  doi     = {10.1137/0202013}
}
```

Establishes that the (non-commutative) bilinear rank of matrix multiplication is
**invariant under the action of S₃** on the three dimensions (Theorem 6 +
Corollary 7) — i.e. `R(⟨m,n,p⟩)` is the same for all six permutations of
`(m,n,p)`. This is the duality the catalog already leans on (e.g. the
`⟨2,3,3⟩=15` sub-products) and that BCS97 [\[65\]](#65-bcs97) attributes to
"Pan and independently Hopcroft and Musinski".

**Where the `⟨2,n,p⟩` formulas are** — **draft page 76** (physical PDF page 4 of
the local copy), in the "Matrix Multiplication" section right after Corollary 7:

> "… in [3] it is shown that `⌈(3pn + max(n,p))/2⌉` multiplications is sufficient
> for p×2 by 2×n matrix multiplication. It follows that `⌈(3pn + max(n,p))/2⌉`
> multiplications is sufficient for 2×p by p×n matrix multiplication. Since
> `⌈7n/2⌉` multiplications are necessary and sufficient for 2×2 by 2×n matrix
> multiplication [3] it follows that `⌈7n/2⌉` multiplications are necessary and
> sufficient for 2×n by n×2 matrix multiplication. Similarly since 15
> multiplications are necessary and sufficient for 3×2 by 2×3 matrix
> multiplication, 15 multiplications are necessary and sufficient for 3×3 by 3×2
> matrix multiplication."

Here HM73's reference `[3]` **is Hopcroft-Kerr 1971** (our `hk71`). So this paper
is an independent primary statement of the HK `⌈(3pn+max(n,p))/2⌉` bound, plus the
S₃-duality that carries it across all inner-dimension-2 (and dual) formats.
Explicit worked dual theorems (`2×2·2×n`, `2×n·n×2`, …) follow on draft pp. 81–82
(PDF pp. 9–10).

📄 **Local PDF (draft)**:
[`references/papers/hopcroft_musinski_1973_duality_bilinear_forms_draft.pdf`](references/papers/hopcroft_musinski_1973_duality_bilinear_forms_draft.pdf)
— the **draft** version (internal pages 73–87; the published SIAM article is
pp. 159–173, so draft p. 76 ≈ published p. 162). Source mirror:
<https://cr.yp.to/bib/1973/hopcroft-duality-draft.pdf>

Cited in: REFERENCES.md (Hopcroft-Kerr `⟨a,2,c⟩` family caveat), FUTURE_WORK.md.

---

## Scheme provenance policy (composition chains)

When a scheme on disk is a simple **Kronecker product** (or other
trivial composition) of previously-known schemes, its `source` field
should resolve to the **underlying contributing papers** — the
"leaves" of the composition tree — not the composer.

When a scheme requires **non-trivial analysis** (paired sub-products,
recursive composition with addition processing, asymmetric padding),
the analyzer is also recorded.

Example: `⟨6,6,12⟩=280` is `⟨2,2,2⟩=7 (Strassen) ⊗ ⟨3,3,6⟩=40 (Smirnov)`.
Source attribution: `Strassen 1969 + Smirnov 2013` (no analyzer
attribution because the Kronecker product is trivial). Distinguish
from `⟨14,14,14⟩=1720` which requires Pan pair-product on 3 cyclic
pairs — that's `Strassen 1969 + Perminov ⟨7,7,7⟩=250 + Pan
pair-product 1978; composition by solven-strassen 2026 (Pan pair
realisation in `PanPairProduct`)`.

For complex Sedoglavic/DIS09-style multi-step compositions (e.g.
⟨18,18,18⟩=3200 reported by Sedoglavic's `fmm_digest` and reproducible
via Strassen[9,9]³ + Pan pair on cyclic ⟨9,9,9⟩ triples), the
analyzer is attributed (Sedoglavic / DIS 2009 / solven-strassen) in
addition to the leaf bases.

Concretely the catalog JSON should grow two fields:
- `leaf_sources`: list of contributing leaf-scheme attributions
- `composition_analyzer`: optional paper/tool that did the
  non-trivial assembly

(Schema migration tracked separately; see task #44 "construction
lineage" which subsumes this provenance work.)

---

## How to add a new reference

1. Append to the index table above with the **next available number**.
2. Add a `## [N] <a name="N-key"></a>key` section below in chronological
   order.
3. Paste a `bibtex` code block — required for every entry.
4. Annotate: what the work proves / contains, and **what field/setting
   it applies to** (every matmul-rank claim must name its field — see
   `memory: feedback-state-field-explicitly`).
5. Add a "Cited in" line listing the consuming MDs.

If the BibTeX detail is uncertain (you have only an arXiv ID, the
volume isn't published yet), mark the affected fields with `[?]` —
better an honest placeholder than fabricated data.
