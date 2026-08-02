# HK71 for ⟨2,p,n⟩ — the logic in one pass, and where reachability goes fuzzy

A synthetic companion to [OVERVIEW.md](OVERVIEW.md) (status) and
[CONSTRUCTIVE_METHOD.md](CONSTRUCTIVE_METHOD.md) (full recipe + proofs).
This note answers two questions only:

1. **Why does the Hopcroft–Kerr construction give
   `R(⟨2,p,n⟩) ≤ ⌈(3pn + max(p,n))/2⌉`?** — the causal chain, not the algebra.
2. **Why is it *unclear* whether some larger shapes can reach that formula
   at all?** — what exactly is unknown, and at which layer.

Everything here is a re-reading of results already established in this repo;
no new claim. All ranks are ℤ-field claims (`Z ⇒ Q,R,C`, and F₂/F₃ by
reduction) and all bounds are **upper-bound tier** — lower bounds are parked
(HK proved optimality only for `⟨2,2,n⟩` and `⟨2,3,3⟩`).

---

## Part 1 — the logic

### The one structural fact: with a `2` in the shape, cells are cheap in *pairs*

`Y = Ā·X` with `Ā ∈ K^{n×2}`, `X ∈ K^{2×n}`. Every output cell `y_ij` is a
**2-term inner product** `a_{i1}x_{1j} + a_{i2}x_{2j}` — naively 2 products.
HK's whole construction is an accounting exercise over three cell prices:

| Object | Naive | HK price | Mechanism |
|---|---:|---:|---|
| diagonal cell `y_ii` | 2 | **2** | but both products stay **reusable** (methods 1/2/3) |
| off-diagonal pair `(y_ij, y_ji)` | 4 | **3** | Lemma 2, *provided* `c(i) ≠ c(j)` |
| a lone off-diagonal cell | 2 | 2 | no saving |

A "method" is just which 2 products you spend on the diagonal:
`1:{a₂(x₁+x₂), (a₁−a₂)x₁}`, `2:{(a₁−a₂)x₂, a₁(x₁+x₂)}`, `3:{a₂x₂, a₁x₁}`.
All three compute `y_ii`; they differ in *what residue they leave behind*.
Lemma 2 says: if rows `i` and `j` were paid with **different** methods, their
four leftovers plus **3** new products reconstruct both `y_ij` and `y_ji`.
That single saved multiplication, replicated over the whole tiling, *is* the
theorem.

### Lemma 1 — why you only need a *band*, and where `max(p,n)` comes from

The real problem is `p×2 · 2×n`, not `n×2 · 2×n`. HK augments: pick
`M ∈ ℤ^{n×p}` (first `p` rows = identity) such that **every cyclic `p`-window
of rows is nonsingular**, and work with `Ā = M·A`. Then any `p` *contiguous*
cells of an output column determine the whole column by one linear solve
(`A·x_j = M_j⁻¹ · (ĀX)_{band(j),j}`) — no multiplications, only the
back-substitution.

Two consequences that drive everything downstream:

- **You only pay for a band of `p` cells per column**, not `n`. The `n`-vs-`p`
  asymmetry in the formula (`max(p,n)`) is exactly this: the internal problem
  is sized `n`, the band width is `p`.
- **The scheme's field is decided here.** Dense ±1 augmentation rows can never
  be unimodular for `n−p ≥ 2` (every 2×2 ±1-minor is even) — that is where the
  historical denominators 2/4/8/16 came from. The Euclidean comb construction
  (`LemmaOneAugmentation.buildUnimodular`) makes every window determinant ±1,
  so every window inverse is integral and all 465 emitted schemes are over ℤ.
  This axis is *orthogonal* to the rank question below.

### The budget is exact — that's the whole difficulty

**Odd `p = 2k+1`.** Band per column = the symmetric window `[j−k, j+k]`, so
each column needs its diagonal plus `k` cells on each side — i.e. all cyclic
pairs at distances `1..k`. Cost:

```
2n  (diagonals)  +  3 · k · n  (n pairs per distance, 3 each)  =  n(3k+2)
                                                              =  (3pn + n)/2      ✓ formula
```

**Even `p = 2k+2`.** Distances `1..k` give `2k+1 = p−1` band cells per column;
each column needs **one more**, at distance `k+1`. The formula allows exactly
`3n/2` further products — i.e. `n/2` pairs at 3 products each, and **nothing
else**.

> This is the load-bearing observation: **the HK formula has zero slack.**
> It is not "a bound with room"; it is the exact cost of paying every band
> cell at the 3-products-per-pair rate. The only slack in the whole
> construction is the single unit the ceiling `⌈·⌉` grants when `n` is odd.
> Consequently *any* cell that has to be computed naively costs `+1` over
> formula. Attaining the bound = never being forced into a naive cell.

### So the entire problem reduces to a colouring/scheduling problem

Lemma 2 only fires on **different-method** pairs. With the alternating
colouring `1,2,1,2,…` on the cyclic row order, pairs at odd distance are
different-method (free), pairs at **even** distance are same-method (blocked).
Half the required pairs are blocked. Two devices unblock them:

1. **Bridges** (HK, for same-method pairs). Route `(i,j)` through a third
   position `b` carrying a different method, reusing the `(i,b)` and `(b,j)`
   material; still 3 products.
2. **Bridge selection** (ours — the key unblocking step). Processing pairs in
   order of increasing distance, `(i,j)` at distance `d` may use **any
   arc-interior** `b = i+e (mod n)`, `1 ≤ e < d`, because both `(i,b)` and
   `(b,j)` are already computed. With the alternating colouring an
   opposite-method interior position always exists (`e=1`, or `e=2` when
   `i+1` is the lone method-3 slot). **Hence the historically stuck
   `(2,2,bridge-3)` configuration never arises in the schedule.**

That closes odd `p` at the formula, everywhere.

### Even `p`: the leftover-column mechanism (this is where the trouble lives)

The distance-`(k+1)` pair `(i, i+k+1)` supplies the missing cell of column `i`
*and* of column `i+k+1` — 3 products for two columns, exactly on budget. So
the question is: **can the columns be perfectly matched** under
`i ↔ i+(k+1) (mod n)`? That circulant splits into

```
g := gcd(n, k+1)   cycles, each of length n/g
```

Even-length cycles match perfectly. Each **odd**-length cycle leaves one
column unmatched — `g` leftover columns in the bad case. Budget check: the
formula funds `n/2` pairs (plus one naive cell when `n` is odd, absorbed by the
ceiling). Every unrescued leftover column must be finished naively at 2
products instead of its 1.5-product share.

**The `Z`-trick** (HK's Step 3, repaired by us) rescues leftover columns **two
at a time**: with `α = a_{i₂}+a_{i₃}` and `β = a_{i₁}+a_{i₄}` treated as
*virtual* diagonal rows, whose "diagonal products" already exist as the band
pairs' cross-products, the pair `(α,β)` is an ordinary Lemma-2 pair — 3
products for 2 columns, back on budget. Two leftover columns rescued = 3
products; the same two done naively = 4. **Each dropped `Z`-pair is exactly
`+1`.**

And here is the catch that propagates into Part 2: **each `Z`-pair consumes one
method-3 row**, and method-3 rows must be **pairwise more than `k` apart** —
otherwise the band itself contains a `(3,3)` same-method pair, which (theorem,
below) has no 3-product completion. Attaining the formula therefore requires
placing `⌊g/2⌋` method-3 rows on a cycle of `n` positions with pairwise gaps
`≥ k+2`.

---

## Part 2 — why reachability becomes unclear at higher shapes

### The geometric wall (arc-sum)

Placing `m` rows on `C_n` pairwise `≥ k+2` apart needs `m(k+2) ≤ n`, i.e.

```
m_max = ⌊ n / (k+2) ⌋        method-3 rows available
m_need = ⌊ g / 2 ⌋           method-3 rows required   (g = gcd(n, k+1))
```

For every even-`p` band shape in the swept range `3 ≤ p ≤ 32, p ≤ n ≤ 32`,
`m_max = 2`. So the formula is attainable iff `⌊g/2⌋ ≤ 2`, i.e. **`g ≤ 5`** —
and the excess is the number of `Z`-pairs that could not be placed:

```
excess  =  max(0, ⌊g/2⌋ − ⌊n/(k+2)⌋)
```

This reproduces every recorded deficit exactly (`p = 2k+2`):

| Shape | `k+1` | `g = gcd(n,k+1)` | `⌊g/2⌋` needed | `⌊n/(k+2)⌋` available | excess | emitted |
|---|---:|---:|---:|---:|---:|---:|
| ⟨2,10,15⟩ | 5 | 5 | 2 | 2 | 0 | 233 = formula |
| ⟨2,12,18⟩ | 6 | 6 | 3 | 2 | **+1** | 334 |
| ⟨2,14,21⟩ | 7 | 7 | 3 | 2 | **+1** | 453 |
| ⟨2,24,30⟩ | 12 | 6 | 3 | 2 | **+1** | 1096 |
| ⟨2,16,24⟩ | 8 | 8 | 4 | 2 | **+2** | 590 |
| ⟨2,18,27⟩ | 9 | 9 | 4 | 2 | **+2** | 745 |
| ⟨2,20,30⟩ | 10 | 10 | 5 | 2 | **+3** | 918 |

*(The two ingredients — the arc-sum argument and "+1 per dropped `Z`-pair" —
are both established in CONSTRUCTIVE_METHOD.md; the closed form above is just
their composition, offered because it reproduces all six observed excesses.)*

**Why this bites harder at higher shapes.** `m_need` grows with `g`, which
grows along the divisibility families (`n = 3(k+1)`, `n = 2(k+1)`, …), while
`m_max` stays pinned near `n/(k+2) ≈ 2` because `n ≤ 2p−1 ≈ 4k`. Two curves
that diverge: the demand for method-3 rows scales with `gcd(n, k+1)`, the
supply does not scale at all. So the deficit is not a fixed defect — it
**widens** (`+1 → +3` already inside `p ≤ 32`), and there is no reason to
expect it to stop widening beyond the swept range.

### Now: what exactly is unclear — three distinct layers

The word "unclear" hides three different kinds of not-knowing. Keeping them
apart is the point of this note.

**Layer 1 — inside the framework: settled, negative.**
For these shapes the `+1..+3` is *forced* by two stacked results: the arc-sum
argument above, and the `(3,3,bridge-1/2)` theorem
(`sympy/derive_33bridge_general.py`): **no 3-product completion by ANY rank-1
atoms local to `{i,b,j}` exists over the emitter's 9-product reusable set** —
the required spill would have to be a rank-1 element of `(U₀⊗V₀) ∩ span(shared)`,
which is 1-dimensional and generated by the **rank-2** virtual diagonal
`(a_i+a_b)⊗(x_i+x_b)`. Not a search failure; a proof.

**Layer 2 — the scope of that theorem: genuinely open, and known to be fragile.**
The theorem is scoped to *rank-1 atoms, local to `{i,b,j}`, over a 9-product
reusable set*. That scope has **already been broken once, in this repo**, in
the analogous case:

> `(2,2,bridge-3)` was proved impossible over the narrow **6-product** set
> `S = {Ci,Di,Cj,Dj,Ep,Fp}` (`references/hopcroftkerr1971/README.md`, a
> correct char-0 Gröbner theorem) — and is **solvable** over the emission's
> true **12-product** set, with an explicit identity now shipped in
> `emitSameMethodPair_22_bridge3`.

The same enlargement for `(3,3,·)` (add the `(b,j)`-pair's three products with
free weights) is *partially* analysed: the block forced-structure argument
survives under projection, but the spill-collapse weakens from proportionality
to *projected* proportionality, leaving a finite quadratic system that has not
been eliminated. Unexplored beyond that: non-local atoms (rows/cols outside
`{i,b,j}`), and 4-product `(3,3)` completions traded against a saving
elsewhere (budget bookkeeping says one product saved anywhere pays for one
`+1`). **So "unreachable" here means "unreachable by the atoms we have proved
about", and the precedent says the proof boundary is not the truth boundary.**

**Layer 3 — the bound itself: never actually proved at these shapes.**
This is the deepest source of the fuzziness, and it is HK's, not ours. HK's
Theorem 1 relies on a method sequence (their Lemma 3) that **forces `(3,3)`
same-method pairs** in the band, and the paper waves the case off ("the other
cases follow by symmetry" — demonstrably false: the natural involutions
`a₁↔a₂`, `x₁↔x₂` swap methods 1↔2 but *fix* method 3, so `(1,1,bridge-2)`
does not map onto the method-3 cases under any single involution). So at
exactly the shapes where the `Z`-trick forces clustered method-3 rows,
`⌈(3pn+max(p,n))/2⌉` is a **claim with no published construction behind it**,
not a theorem whose construction we merely failed to reproduce.

Two independent witnesses that this is real rather than an artefact of our
implementation:

- FMM-Lille publishes **every** `⟨2,n,n⟩` for `n = 9..16` *strictly above* the
  HK formula (`+1..+5`), all credited "hopcroft:1971". Perminov shows the same
  pattern. Nobody has ever published a scheme attaining the formula where the
  gap bites.
- The excesses match `+δ`-per-unclosed-case bookkeeping exactly.

And on the other side: our own construction **beat** those catalogs at 219 of
the swept shapes (e.g. `⟨2,10,15⟩ = 233` vs 234 everywhere, `⟨2,24,27⟩ = 986`
vs 1015) — including at `g = 5`, where the catalogs had also stalled. Which is
the honest reading of the situation: the gap where the catalogs sat was partly
a *derivation* gap that yielded to better scheduling, and what is left at
`g ≥ 6` may be the same kind of gap or may be a genuine obstruction. **We
cannot currently tell those apart**, because:

- we have no lower bound (parked by explicit decision — everything emitted is
  upper-bound tier), so we cannot certify `+1` as optimal; and
- we have no impossibility proof at full scope (Layer 2), so we cannot certify
  `+1` as forced beyond local rank-1 atoms.

The six shapes sit in that bracket. They are, note, **already below every
published catalog** — the `+1..+3` is a gap against our own formula target,
not against the field.

### One-line summary

The HK bound is the exact cost of paying every band cell at 3-products-per-pair;
attaining it is a *scheduling* problem in which same-method pairs must be
bridged and leftover columns must be `Z`-paired; each `Z`-pair burns a
method-3 row; method-3 rows must be `> k` apart; and for `g = gcd(n,k+1) ≥ 6`
the ring is too small to hold enough of them. Whether that wall is an artefact
of the atom family we can prove about (it was, once, for `(2,2,bridge-3)`) or a
real obstruction to the formula itself (which HK never proved at these shapes)
is the open question — task #9.

### If you want to resolve it, in order of promise

1. **Eliminate the 12-product `(3,3,bridge)` quadratic system** (Gröbner after
   the structural substitution). A solution closes all six shapes; a second
   impossibility makes `+1` framework-optimal in the strongest local sense.
2. **Non-local atoms** — rows/cols beyond `{i,b,j}`; no analysis yet.
3. **4-product `(3,3)` completions traded against a saving elsewhere.**
4. **A different tiling / colouring altogether** — the residual quotient `Q`
   that all these impossibility proofs are about is a function of the method
   assignment; change the assignment and none of the theorems transfer.
5. **Lower bounds** (currently parked) — the only route to turning "we can't
   reach it" into "it can't be reached".
