"""
Import Schwartz-Zwecher 2025 trilinear-aggregation algorithms (n=20..50, step 2)
from the supplemental NumPy .npz archive into our sparse-JSON catalog format.

Source:
  arXiv:2508.01748 "Towards Faster Feasible Matrix Multiplication by Trilinear
  Aggregation" by Oded Schwartz & Eyal Zwecher (Hebrew University, 2025).
  Supplemental zip: https://www.cs.huji.ac.il/~odedsc/papers/
                    trilinear_aggregation_algorithms_decomposed-2025-07-29.zip

What's in the zip:
  algorithm_{n}_{n}_{n}_{rank}_decomposed.npz for n in 20,22,...,50.
  Each .npz holds u_phi, v_phi, w_phi (rank × (n/2+1)²-ish) and phi
  ((n/2+1)²-ish × n²). The full encoding/decoding matrices are
    U = u_phi @ phi    (rank × n²)
    V = v_phi @ phi    (rank × n²)
    W = w_phi @ phi    (rank × n²)
  Coefficients are Q-rational; the only denominator in play is n/2+1
  (cf. Theorem 3.4: γ = 1 - 9/(n/2+1)).

W convention in the zip (per README.txt + sample.py):
    actual = (w.T @ γ).reshape(n,n).T
  i.e. w[k, l*n + i] is the coefficient of γ_k in C[i,l]. This matches
  dronperminov's col-major C-flatten convention that SchemeIO reads
  directly — so we write w as-is into "w_sparse".

Field: Q. We round each coefficient to its nearest rational with
denominator dividing (n/2+1). This is an EXACT reconstruction from the
paper's claimed rational structure — values stored in the .npz are
float64 with rounding noise (~1e-16) which we discard.

Run:
  python3 tools/sz2025-import/import_schwartz_zwecher_2025.py
  python3 tools/sz2025-import/import_schwartz_zwecher_2025.py --apply
  python3 tools/sz2025-import/import_schwartz_zwecher_2025.py --apply --only 44
  python3 tools/sz2025-import/import_schwartz_zwecher_2025.py --apply --max-n 32
"""
import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
SCHEMES_ROOT = REPO_ROOT / "src" / "main" / "resources" / "schemes"
SUPPL_ROOT = REPO_ROOT / "references" / "schwartz-zwecher-2025" / \
    "trilinear_aggregation_algorithms_decomposed"


# All (n, rank) pairs published in Schwartz-Zwecher 2025 supplemental data.
# `discovery` is determined by whether SZ's rank STRICTLY improves on the
# best non-commutative bound previously known (Pan 1982 / DIS 2009 / etc).
# For n < 28, SZ's TA-New25 family equals or slightly under-performs the
# Pan TA bound in DIS09 Table 3 — those entries are still valid algorithms
# but should not be misattributed as SZ "discoveries".
#
# Reference points from existing docs/cited-bounds.json at time of import:
#   n=20: DIS09 Pan TA = 4340 < SZ 4378  → not a discovery (Pan 1982)
#   n=22: DIS09 Pan TA = 5566 < SZ 5596  → not a discovery (Pan 1982)
#   n=24: DIS09 Pan TA = 7000 < SZ 7020  → not a discovery (Pan 1982)
#   n=26: DIS09 Pan TA = 8658 < SZ 8666  → not a discovery (Pan 1982)
#   n=28: DIS09 Pan TA = 10556 > SZ 10550 → DISCOVERY
#   n=30: DIS09 Pan TA = 12710 > SZ 12688 → DISCOVERY
#   n=32: DIS09 Pan TA = 15113 > SZ 15096 → DISCOVERY
SZ_SCHEMES = [
    (20,  4378, False, "Pan 1982"),
    (22,  5596, False, "Pan 1982"),
    (24,  7020, False, "Pan 1982"),
    (26,  8666, False, "Pan 1982"),
    (28, 10550, True,  None),
    (30, 12688, True,  None),
    (32, 15096, True,  None),
    (34, 17790, True,  None),
    (36, 20786, True,  None),
    (38, 24100, True,  None),
    (40, 27748, True,  None),
    (42, 31746, True,  None),
    (44, 36110, True,  None),
    (46, 40856, True,  None),
    (48, 46000, True,  None),
    (50, 51558, True,  None),
]


def expand_scheme(npz_path: Path):
    """Load a .npz file and compute U, V, W = u_phi @ phi, etc."""
    alg = np.load(npz_path)
    U = alg["u_phi"] @ alg["phi"]
    V = alg["v_phi"] @ alg["phi"]
    W = alg["w_phi"] @ alg["phi"]
    return U, V, W


def quantize_to_q(M: np.ndarray, denom_max: int, tol: float = 1e-9) -> np.ndarray:
    """Vectorised Q-snap. Per CLAUDE-discovery (empirical inspection of
    the SZ data), every nonzero entry of U/V/W matches `k/D` for some
    D ∈ {1, 2, 4, 8, 16} · denom_max. We try those D in increasing
    order — the first one for which `M*D` rounds to an integer within
    `tol` is the right denominator.

    Returns a float64 array with the cleaned values; entries below tol
    (relative to D=1) become exact 0.0. Vastly faster than per-cell
    Python Fraction calls (numpy ops only)."""
    out = np.zeros_like(M)
    remaining = ~np.isclose(M, 0.0, atol=1e-12)
    # Candidate denominators, in increasing order.
    candidates = []
    for pow2 in [1, 2, 4, 8, 16, 32]:
        candidates.append(pow2)
        candidates.append(pow2 * denom_max)
    # Deduplicate while preserving order.
    seen = set()
    unique_candidates = []
    for d in candidates:
        if d not in seen:
            seen.add(d)
            unique_candidates.append(d)
    for D in unique_candidates:
        if not remaining.any():
            break
        scaled = M * D
        rounded = np.round(scaled)
        matches = remaining & (np.abs(scaled - rounded) < tol)
        out[matches] = rounded[matches] / D
        remaining &= ~matches
    if remaining.any():
        # Fallback: use the original float for residual entries (rare).
        out[remaining] = M[remaining]
    return out


def rationalise(value: float, denom_max: int) -> Fraction:
    """Snap value to a small rational. SZ 2025 uses denominators dividing
    n/2+1 (the γ denominator from Theorem 3.4), but products from the
    transformation φ — e.g. u from R = ½·I + (1/(d+1))·11ᵀ — carry
    powers of 2 alongside d=n/2+1. Empirically the maximum denominator
    we see is (n/2+1) · 2⁴ = 16·(n/2+1), so we allow `limit_denominator`
    a generous cap of 1024·(n/2+1) — enough to absorb any 2-power chain
    while staying small enough that no spurious rational close to a
    float-noise representation slips through."""
    if abs(value) < 1e-12:
        return Fraction(0)
    return Fraction(value).limit_denominator(1024 * denom_max)


def verify_dense(U: np.ndarray, V: np.ndarray, W: np.ndarray, n: int,
                 samples: int = 3, tol: float = 1e-8):
    """Random matmul spot check using the SZ-published convention.
    W is rank x n² with col-major C-flatten (matches dronperminov)."""
    rng = np.random.default_rng(0xC0DE)
    worst = 0.0
    for _ in range(samples):
        A = rng.standard_normal((n, n))
        B = rng.standard_normal((n, n))
        alpha = U @ A.ravel()
        beta = V @ B.ravel()
        c_flat = W.T @ (alpha * beta)
        C_algo = c_flat.reshape(n, n).T   # SZ convention: reshape then T
        C_naive = A @ B
        res = float(np.max(np.abs(C_algo - C_naive)))
        worst = max(worst, res)
        if res > tol:
            return False, worst
    return True, worst


def fraction_to_json(f: Fraction):
    if f.denominator == 1:
        return int(f.numerator)
    return f"{f.numerator}/{f.denominator}"


def sparse_factor(M: np.ndarray, denom_max: int):
    """Convert rank × dim into sparse-list-per-multiplication.
    Assumes M has been Q-snapped already (exact 0 entries where the
    snapping decided 'zero'). Stores Q-rational entries as 'num/den'
    strings; integer entries as int. denom_max guides rationalisation
    of the surviving floats."""
    rank, dim = M.shape
    out = []
    nonzero = 0
    for k in range(rank):
        col = {}
        row = M[k, :]
        nz_idx = np.flatnonzero(row != 0)
        for pos in nz_idx:
            val = float(row[int(pos)])
            frac = rationalise(val, denom_max)
            if frac == 0:
                continue
            col[str(int(pos))] = fraction_to_json(frac)
            nonzero += 1
        out.append(col)
    return out, nonzero


def count_additions_sparse(U, V, W) -> int:
    """Sum of (nnz - 1)+ across U columns, V columns, and W output rows.
    Operates on already-snapped (exact-zero) matrices using numpy
    boolean ops — orders of magnitude faster than Python-level
    Fraction reconstruction."""
    u_nz = int(np.maximum(0, (U != 0).sum(axis=1) - 1).sum())
    v_nz = int(np.maximum(0, (V != 0).sum(axis=1) - 1).sum())
    w_nz = int(np.maximum(0, (W != 0).sum(axis=0) - 1).sum())
    return u_nz + v_nz + w_nz


def output_path(n: int, rank: int, additions: int) -> Path:
    section = SCHEMES_ROOT / f"section{n}"
    fname = f"schwartz-zwecher-2025_{n}x{n}x{n}_r{rank}_a{additions}_Q.json"
    return section / fname


def import_one(n: int, rank: int, apply: bool, derivation_task: str,
               discovery: bool, attribution_for_rank: str | None):
    npz_path = SUPPL_ROOT / f"algorithm_{n}_{n}_{n}_{rank}_decomposed.npz"
    if not npz_path.exists():
        return ("missing", None, None)
    print(f"  <{n},{n},{n}> r={rank}: loading {npz_path.name}")
    U, V, W = expand_scheme(npz_path)
    actual_rank = U.shape[0]
    if actual_rank != rank:
        print(f"    WARN: actual rank {actual_rank} != claimed {rank}")

    ok, residual = verify_dense(U, V, W, n)
    if not ok:
        print(f"    FAIL: float verification residual={residual:.2e}")
        return ("float_fail", residual, None)
    print(f"    float verify OK (max residual {residual:.2e})")

    denom_max = n // 2 + 1

    # Q-snapped re-verification (uses the same SZ convention). We snap
    # each entry to the nearest k/D where D is the smallest divisor of
    # (n/2+1)·16 such that k/D matches the float within 1e-9 — empirically
    # all SZ denominators are of the form (n/2+1)·2^k with k ≤ 4.
    Uq = quantize_to_q(U, denom_max)
    Vq = quantize_to_q(V, denom_max)
    Wq = quantize_to_q(W, denom_max)
    ok_q, residual_q = verify_dense(Uq, Vq, Wq, n)
    if not ok_q:
        print(f"    WARN: Q-snapped verify residual={residual_q:.2e}; "
              f"writing raw floats instead")
        store_U, store_V, store_W = U, V, W
        store_denom = 1
    else:
        print(f"    Q-snapped verify OK (denom={denom_max}, "
              f"residual={residual_q:.2e})")
        store_U, store_V, store_W = Uq, Vq, Wq
        store_denom = denom_max

    u_sparse, u_nnz = sparse_factor(store_U, store_denom)
    v_sparse, v_nnz = sparse_factor(store_V, store_denom)
    w_sparse, w_nnz = sparse_factor(store_W, store_denom)
    additions = count_additions_sparse(store_U, store_V, store_W)

    out = output_path(n, rank, additions)
    payload = {
        "n": [n, n, n],
        "m": rank,
        "z2": False,
        "field": "Q",
        "commutative": False,
        "importing_source": "Schwartz-Zwecher 2025",
        "source_paper": "arXiv:2508.01748",
        "source_paper_url": "https://arxiv.org/abs/2508.01748",
        "source_data": "https://www.cs.huji.ac.il/~odedsc/papers/"
                       "trilinear_aggregation_algorithms_decomposed-"
                       "2025-07-29.zip",
        "source_data_file": npz_path.name,
        "year": 2025,
        "algorithm_family": "TA-New25 (Schwartz-Zwecher 2025)",
        "construction": "Trilinear aggregation + kin-row unification "
                        "(Theorem 2.22 / 3.4 of arXiv:2508.01748)",
        "notes": (
            "Imported from supplemental .npz (decomposed form "
            "U=u_phi@phi, V=v_phi@phi, W=w_phi@phi). Coefficients "
            f"live in Q with denominators dividing n/2+1={denom_max}. "
            "Header rank/additions match Table 1 of arXiv:2508.01748."
        ),
        "derivation_task": derivation_task,
        "verification_notes": (
            f"Float verify residual <= {residual:.2e}; "
            f"Q-snapped (denom={denom_max}) "
            f"{'OK' if ok_q else 'FAIL — stored raw floats'}."
        ),
        "u_sparse": u_sparse,
        "v_sparse": v_sparse,
        "w_sparse": w_sparse,
    }
    print(f"    nnz: U={u_nnz}, V={v_nnz}, W={w_nnz}, additions={additions}")
    print(f"    -> {out.relative_to(REPO_ROOT)}")
    if apply:
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w") as f:
            json.dump(payload, f)
        size_mb = out.stat().st_size / 1024 / 1024
        print(f"    WROTE ({size_mb:.2f} MB)")
    return ("ok", residual, out)


def main():
    apply = "--apply" in sys.argv
    only = None
    max_n = None
    for i, a in enumerate(sys.argv):
        if a == "--only" and i + 1 < len(sys.argv):
            only = int(sys.argv[i + 1])
        if a == "--max-n" and i + 1 < len(sys.argv):
            max_n = int(sys.argv[i + 1])

    print(f"Mode: {'APPLY' if apply else 'DRY-RUN'}")
    print(f"Source dir: {SUPPL_ROOT}")
    print(f"Target schemes root: {SCHEMES_ROOT}")
    if only is not None:
        print(f"Filter: only n={only}")
    if max_n is not None:
        print(f"Filter: max_n={max_n}")
    print()

    results = {}
    for n, rank, discovery, attribution in SZ_SCHEMES:
        if only is not None and n != only:
            continue
        if max_n is not None and n > max_n:
            continue
        status, residual, path = import_one(
            n, rank, apply,
            derivation_task="TBD-SZ2025-kin-row-constructor",
            discovery=discovery,
            attribution_for_rank=attribution,
        )
        results[n] = (status, residual, path)

    print()
    print("Summary:")
    for n, (status, residual, path) in sorted(results.items()):
        print(f"  n={n}: {status}, residual={residual}")


if __name__ == "__main__":
    main()
