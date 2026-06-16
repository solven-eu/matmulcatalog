# Unimodular Lemma-1 augmentation (task #11) — discovery artifacts

Goal: an n×p Lemma-1 matrix `[I_p; B]` with EVERY cyclic p-window
determinant ±1, so the HK ⟨2,p,n⟩ back-substitution stays integer and the
constructed schemes are over ℤ. **Solved** — the production construction is
`LemmaOneAugmentation.buildUnimodular` (Euclidean recursion); these scripts
are the discovery/verification trail, in chronological order:

1. `search_intervals.py` + `intervals.out` — consecutive-ones (totally
   unimodular interval) rows: TU guarantees minors ∈ {0,±1}, but placement
   fails for many (p, m) → dead end, kept as the negative result.
2. `comb_seam_coldfs.py` + `comb_seam.out` — the key reduction (window dets
   of `[I;B]` ⇔ three minor families of B: leading / trailing / sliding) +
   period-m comb body (proven ±1 when m | p) + column-DFS for the seam when
   m ∤ p. Found seams for all p ≤ 11; printing the solutions exposed the
   closed structure.
3. `euclidean.py` + `euclidean.out` — the final construction: comb body +
   tail = `B(m, p mod m)ᵀ` (the seam must solve the SAME three-family
   problem transposed at (m, r) — gcd descent, terminates at the pure
   comb). **Zero failures over the entire band range 3 ≤ p ≤ 32.**
4. `closed_form.py` + `closed_form.out` — a wrong intermediate guess
   (`[I|ones]`-over-`[I]` tail), kept as the falsification record.
5. `comb_seam.py` — earlier row-wise DFS (superseded by the column-wise 2).
6. `hk2np_full_sweep_2026-06-11.log` — log of the pre-unimodular (ℚ)
   certified emission, for provenance.

Key negative worth remembering: **{−1,1}-dense rows can never be
unimodular for m ≥ 2** (every 2×2 ±1-minor is even) — the origin of the
denominators 2/4/8/16 in the first (ℚ) emission.

Full writeup: `research/hopcroft-kerr-2np/CONSTRUCTIVE_METHOD.md`
(ingredient 6); paper: `paper/sections/hk71.tex` block (v).
