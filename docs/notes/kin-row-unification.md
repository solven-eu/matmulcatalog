# Kin-row unification

A note on what "kin-row unification" is, where it comes from, and how it
spares multiplications. Written 2026-06-08 while investigating why
FMM-Lille's `⟨27,28,28⟩=10413` sits below our projection's `10442`.

## TL;DR

**Kin-row unification** is a *construction-time* reduction in bilinear
matrix-multiplication algorithms: two products are **kin** when they have
**identical lines in `U` and `V`** (the same left A-form *and* the same
right B-form), i.e. they compute the *same* bilinear product `L·R`.
Uniting kin products computes that product **once** instead of twice,
sparing one multiplication per united pair. The cleverness is not the
uniting (trivial once the kinship exists) but **deliberately creating the
kinship** via a de Groote equivalence transform on a sub-block, so that a
sub-algorithm's products line up with products the aggregation step
already produces.

## Where the term comes from (references)

- **Primary.** Oded Schwartz & Eyal Zwecher, *"Towards Faster Feasible
  Matrix Multiplication by Trilinear Aggregation"*, **arXiv:2508.01748**
  (Aug 2025). §3 (the TA-New25 family) and §4 (TA-New25b) define and use
  kin-row unification; Theorems 2.22 and 3.4 are the load-bearing
  statements. The paper's headline is an `O(n^2.773203)` *feasible*
  algorithm (base case from `n₀=28`), beating Pan 1982's `O(n^2.773372)`
  for practical input sizes.
- **The equivalence used to *create* kin terms.** de Groote's equivalence
  of bilinear algorithms (de Groote 1978, *"On varieties of optimal
  algorithms for the computation of bilinear mappings"*). A sub-part of
  the algorithm equivalent to a small matmul `⟨2,2,2;7⟩` is replaced by an
  *equivalent* `⟨2,2,2;7⟩` algorithm — parameterised by invertible
  `K_U, K_V ∈ F^{2×2}` (Thm 2.22) — chosen so its products become kin to
  the aggregation's products.
- **The aggregation it sits on — "implicit canceling".** The base scheme
  is built by *trilinear aggregation with implicit canceling*, originally
  Victor Pan, *"Trilinear aggregating with implicit canceling for a new
  acceleration of matrix multiplication"* (Comp. & Math. with Appl., 1982),
  and re-formalised constructively (via *Generating Tables* + linear
  transformations) by **Tor Hadas & Oded Schwartz**, *"Towards Practical
  Fast Matrix Multiplication based on Trilinear Aggregation"*, **ISSAC
  2023**. See the aside below for what "implicit canceling" means. Our
  `PanTrilinearAggregation` emits this bound and the base scheme that
  kin-row unification then reduces.

### Aside: what is "implicit canceling"?

Trilinear aggregation multiplies *pre-summed* inputs. A single aggregated
product like `(a_i + a_j)·(b_k + b_l)` expands to
`a_i b_k + a_i b_l + a_j b_k + a_j b_l` — it yields the two **wanted**
cross-terms (`a_i b_l`, `a_j b_k`) *plus* two **unwanted "garbage" terms**
(`a_i b_k`, `a_j b_l`). The whole method only pays off if the garbage is
removed cheaply:

- **Explicit canceling** subtracts dedicated *correction products* to kill
  the garbage — but each correction costs a multiplication, eating the
  saving.
- **Implicit canceling** instead *arranges the aggregation* (signs,
  pairing, the Generating Table layout) so the garbage terms produced by
  **different** aggregated products **cancel one another** in the final
  sum — for **free**, with no correction multiplications. That free
  cancellation is the source of trilinear aggregation's speed.

Kin-row unification (above) then squeezes out a further constant factor on
top of an implicitly-canceled aggregation scheme.
- **Where we meet it in practice.** The **FMM-Lille** catalog
  (Sedoglavic) ships `_raw.mpl` schemes whose product count is the *pre*-
  unification number, while the headline rank is *post*-unification
  (e.g. `⟨17,17,17⟩` 2934→2931; `⟨27,28,28⟩` 10442→10413). See
  `references/fmm-lille/*/`.

## How it spares multiplications (mechanism)

1. **Aggregation.** Pan's trilinear aggregation builds `⟨n,n,n⟩` by summing
   several outer products through a shared bilinear core — already fewer
   multiplications than naïve recursion. The result is organised as
   *disjoint-sum components* (blocks of the recursive construction).
2. **Find an equivalent sub-block.** Identify a part of the algorithm that
   is itself equivalent to a small matrix multiplication (the `⟨2,2,2;7⟩`
   Strassen-type core).
3. **Re-coordinate it (de Groote).** Replace that sub-block by an
   *equivalent* algorithm under de Groote's equivalence — a basis change
   `K_U, K_V` on its inputs. Equivalence preserves correctness but changes
   the explicit `U`/`V` lines.
4. **Engineer kinship.** Choose `K_U, K_V` so that the transformed
   sub-block's products have **identical `U`- and `V`-lines** to products
   the aggregation step already computes — they become *kin*.
5. **Unite.** Two kin products are the *same* multiplication `L·R`; keep
   one and redirect every output use of the other to it (the `W`/output
   side absorbs the merge). Each united kin pair removes **one
   multiplication** → the rank drops.

Net effect: a scheme of `R` products with `k` kin pairs realisable
collapses to `R − k`. For `⟨27,28,28⟩`, `k=29` (10442 → 10413); for
`⟨17,17,17⟩`, `k=3` (2934 → 2931).

### Why it is NOT a trivial post-hoc dedup

Crucial subtlety, confirmed empirically on FMM's `27x28x28_raw.mpl`
(`AnalyzeFmmKinRows.java`): in the **raw, un-transformed** scheme there
are **zero** products sharing both their `U`- and `V`-lines, **zero**
products sharing any two of the three factor directions, and the
shared-left-factor groups (the 27 size-3 groups — note `27` = the reduced
dimension) all have **full-rank** `Σ v_k⊗w_k` slices. So you cannot
recover the 29 by inspecting/merging the raw scheme alone. The kinship is
**latent**: it only appears *after* applying the right de Groote transform
to the equivalent sub-block. This is why kin-row unification is a
**construction strategy**, not a generic scheme-minimiser, and why our
downward projection reaches FMM's *published* `10442` but not the headline
`10413` (the latter needs the SZ construction, not a tensor post-process).

## Status in this repo

- We **import** SZ2025 schemes (n=20..32) with
  `"derivation_task": "TBD-SZ2025-kin-row-constructor"` — we carry the
  matrices but do **not** yet synthesise them.
- The constructor is a tracked follow-up (ROADMAP "Schwartz-Zwecher 2025 —
  kin-row unification constructor"): implement arXiv:2508.01748 §3/§4 from
  Pan's aggregation tables + a chosen `⟨2,2,2;7⟩` algorithm + the
  Hadas–Schwartz implicit-canceling transform.
- Related fusion code paths that already spare muls by cross-product
  sharing: `PairFusedRecombination`, `PanTrilinearAggregation` (see
  `docs/notes/materialisation-and-overlap.md`).

## Relation to flip-graph search (load-bearing)

Kin-row unification and **flip-graph reductions** are the *same final
move* — unite two products once they coincide — reached two ways:

- **Flip graph** (Kauers–Moosbauer 2022; Perminov's *ternary meta flip
  graphs*, arXiv:2511.20317; the open-source framework arXiv:2603.02398):
  a sequence of rank-**preserving** *flips* shuffles the decomposition. A
  flip takes two rank-1 products sharing one factor and rewrites
  `a⊗b⊗c + a⊗b'⊗c'  ⟶  a⊗b⊗(c−c') + a⊗(b+b')⊗c'` (same sum, still two
  products). Repeated flips can drive a factor to zero or make **two
  products identical**, at which point a **reduction** drops/merges them →
  rank −1. That coincidence is a *kin pair*; the reduction is the
  unification. **So flips GENERATE kin; reductions UNITE them.**
- **SZ kin-row unification** (arXiv:2508.01748): a *deterministic* de
  Groote transform `K_U,K_V` on an equivalent `⟨2,2,2;7⟩` sub-block
  engineers the identity in one shot, then unites.

i.e. **search (flips) vs construction (de Groote)** — same endpoint.

**Practical consequence.** Our `AnalyzeFmmKinRows` finding — FMM's raw
`⟨27,28,28⟩=10442` has *no manifest kin* — means the 29 reductions sit
several rank-preserving flips away. That is exactly what a flip-graph
search explores. So an alternative to implementing the SZ de Groote
constructor (`TBD-SZ2025-kin-row-constructor`) is to **run flip-graph /
meta-flip search on the materialised 10442 schemes** and let it chase the
reductions to 10413. Caveat: flip-graph is **incomplete/stochastic** — not
guaranteed to find what SZ proves constructively; SZ's transform is a known
shortcut to the reduced configuration.

## de Groote vs flip graph: construction vs search

The two are different *kinds* of move, and the distinction is the whole
story of "deterministic warp" vs "random walk".

| | de Groote equivalence | flip graph |
| --- | --- | --- |
| What it is | action of the matmul tensor's **isotropy (symmetry) group** — sandwiching `A↦XAY⁻¹, B↦YBZ⁻¹, C↦ZCX⁻¹`, plus product permute/scale | a **discrete, local** rank-preserving rewrite `a⊗b⊗c + a⊗b'⊗c' ⟶ a⊗b⊗(c−c') + a⊗(b+b')⊗c'` |
| Nature | continuous (Lie group), global, **algebraic** | discrete, local, **combinatorial** |
| How it reaches kin | **solve** for the group element that creates the coincidence | **search** (random walk) until a coincidence happens |
| Cost to reach a known kin | one algebraic solve (closed form for the SZ family) | exponential blind walk; may stall |
| Needs structure? | **yes** — must know *where* (an equivalent sub-block) and set up the kin equations | **no** — walks any scheme |

**Why de Groote is the faster generator (when it applies).** SZ don't
search; they *solve*. They locate a sub-block equivalent to `⟨2,2,2;7⟩`,
write the kinship condition "transformed product line = an aggregation
product line" as equations in the isotropy parameters `K_U,K_V ∈ F^{2×2}`,
and solve them in closed form (rational `K`, denominators dividing
`n/2+1`). The kin configuration is reached **directly**, not stumbled
upon. This is exactly the "deterministic warp" — de Groote knows the
destination; flips wander toward it.

**The catch — not a universal fast kin-generator.** The shortcut needs the
*structure* (an identifiable equivalent sub-block + an aggregation to
write kin equations against). For an arbitrary scheme (AlphaTensor output,
a random rank-`r` decomposition), finding *which* isotropy element creates
kin is itself a hard algebraic problem. de Groote gives the transform
cheaply only once you know *where and what to solve for*; flip-graph needs
none of that, which is why it is the workhorse for unstructured schemes.

**Different spaces — neither subsumes the other.** A generic de Groote
element is *continuous*; it generally **cannot** be realised by any finite
sequence of (discrete, rational) flips — so de Groote can reach kin
configurations flips cannot express. Conversely, flips are **not confined
to one isotropy orbit**: distinct rank-`r` decompositions of the same
matmul tensor can lie in different de Groote orbits, and flips can cross
between them — so flips reach configurations a fixed-orbit de Groote
cannot. Complementary generating sets, not nested.

**Data point — even for buds, the state of the art only *searches*.**
Kauers–Moosbauer–Wood 2026 (arXiv:2602.11041), the bud / serendipitous
("divide less, conquer more") paper, obtain bud-rich decompositions by
(i) a **flip-graph search** maximising the count of shareable
`⟨1,1,k⟩`/`⟨1,ℓ,1⟩`/`⟨m,1,1⟩` sub-blocks under a non-overlap constraint,
then (ii) **random** de Groote orbit elements for support reduction. Their
6×6 `ω≤2.8019` came from *analysing* a Moosbauer–Poole rank-153 scheme and
re-searching for one richer in buds. They report **no general
construction** of the structure. So the de-Groote-as-deterministic-kin /
bud constructor exists today **only** for SZ's narrow kin-row family;
for buds it is open — and applying de Groote *randomly* (as KMW do) is
"flip-like" in spirit.

**Open question / experiment.** Could a deterministic de Groote *solve*
(à la SZ) replace the random search for **buds** too? A concrete test on a
case we already hold: take one of our PanTA-derived `10442` `⟨27,28,28⟩`
schemes and (a) try to *construct* the `K_U,K_V` that surfaces SZ's kin
directly, vs (b) run a flip search on the same scheme — and compare
whether/how fast each reaches `10413`. That pins down "warp vs walk"
empirically.

## Don't conflate with

- **Projection / DCE** (drop an axis index, delete dead products) — a
  *downward* operator; gives FMM's *raw* ranks, not the kin-united ones.
- **Pairwise rank-1 merge** (two products that are scalar multiples) — a
  special, trivial case; the SZ kinship is *engineered*, not found.
- **Bud / serendipitous merging** (Kauers–Moosbauer–Wood 2026,
  arXiv:2602.11041) — merges *shared-output recursive calls* into a larger
  matmul; a different "sharing" axis than identical-product kin merging.
- **τ-theorem disjoint-sum identities** — a different large-shape gap (see
  the `⟨17,17,17⟩` memory).
