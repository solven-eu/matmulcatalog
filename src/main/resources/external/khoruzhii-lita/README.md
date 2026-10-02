# Khoruzhii–Serafin–Gelß–Pokutta 2026 — LITA cube archives (`TA_lita_npz`)

Upstream: <https://github.com/khoruzhii/lita>, directory `schemes/`, repository
state **2026-09-24** (last commit `21b3c73` "Improve LITA7 odd construction.").
Downloaded 2026-09-30 from
`https://raw.githubusercontent.com/khoruzhii/lita/master/schemes/{N}x{N}x{N}_r{R}.npz`.
Reference: REFERENCES.md [82] (`khoruzhii2026lita`).

| what | value |
|---|---|
| files | 20 archives `{N}x{N}x{N}_r{R}.npz`, 13 ≤ N ≤ 32, 12 MB in total (each ≤ 2.1 MB) |
| format | `fmmp.qcsr.v1` — a zip of NumPy `.npy` arrays: per axis `u`/`v`/`w` the CSR arrays `{axis}_indptr` (int64), `{axis}_indices` (int32), `{axis}_numerators` / `{axis}_denominators` (int64); `metadata_json` (tensor dims, rank, coefficient field). Rows = products, columns = ROW-MAJOR flattened matrix entries. |
| ranks | even N: `N³/3 + 3N² + 37N/6 + 5` (14³=1594 … 32³=14197); odd N: `N³/3 + 7N²/2 + 14N/3 − 11/2` (13³=1379 … 31³=13433) |
| field | ℚ (rational; every scheme has even denominators → not reducible mod 2) |

## How the catalog uses them

These files are **not** catalog schemes (the catalog index only reads `*.json`).
Each cube is registered as a hash-stamped stub
`schemes/known/section{N}/{N}x{N}x{N}-r{R}-khoruzhii_2026-{hash7}.json` whose lineage
is the atom `TA_lita_npz(n=N)`; `LineageReplayer` resolves that atom through
`papers.khoruzhii2026.LitaNpzCubes.build(N)` → `io.NpzSchemeReader` (pure Java, no
NumPy). The explicit factors are rebuilt on demand, exactly like the
`DIS09Lemma4(n)` / `TA_lita(n)` formula cubes — they are far too large to store as
catalog JSON (32³ is 117 MB canonical; all 20 total 491 MB).

**Do not replace or re-encode these archives**: projection stubs across the 13–32
band pin their cube parent by `NxNxN@contentHash`, and that hash is a function of
the exact matrices in these files. A newer upstream scheme for some N is a NEW
file (new rank → new stub), never an in-place edit.

Re-import / verify: `docs.migrate.ImportKhoruzhiiLita <dir-of-npz> [--execute]`
(exact BigInteger proof for 13³/15³/17³, spot-check tier above; see the class
Javadoc). Guard: `TestSweepSpotsSota.retains_kgp_lita_cubes`,
`TestLitaNpzCubes`.

## Licence / redistribution — OPEN QUESTION

The upstream repository carries **no licence file** (checked 2026-09-30), and an
"associated manuscript is in preparation". The schemes are mathematical facts
that FMM-Lille already redistributes (as Maple files, credited to the authors),
and this catalog cites them with full attribution — but redistributing the
authors' own archive files is a decision for the repository owner, ideally
settled by asking the authors. The clean long-term exit is the declared
follow-up: port `scripts/lita.py` / `scripts/lita_odd.py` so `TA_lita_npz`
becomes a true generator and these files can be dropped (it must reproduce the
same matrices bit-for-bit, or the dependent projection stubs must be re-derived).
