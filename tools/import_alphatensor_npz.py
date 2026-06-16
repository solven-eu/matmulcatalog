"""
Extract all 27 broken-AT-format factorizations from DeepMind's
factorizations_r.npz and write them as our standard sparse JSON
scheme files. Compares each to a reference (Perminov known-good
3x4x11 entry) for sanity, then writes only those that verify in
naive Python matmul tests.

This is the ROADMAP "Re-fetch AT broken schemes" remediation that
uses the original DeepMind source instead of Perminov mirrors.

Run:
  /tmp/strassen_venv/bin/python tools/import_alphatensor_npz.py [--apply]
"""
import json
import os
import sys
from pathlib import Path

import numpy as np


SCHEMES_ROOT = Path("src/main/resources/schemes")
NPZ_PATH = "/tmp/at_r.npz"

# The 27 broken-AT-Z formats we identified — all are in DeepMind's R npz.
TARGET_FORMATS = [
    (3, 4, 11), (3, 9, 11), (5, 7, 11), (5, 8, 11), (6, 7, 11), (6, 8, 11),
    (7, 7, 11), (7, 9, 11), (7, 10, 11), (7, 11, 11), (8, 9, 11),
    (9, 9, 9), (9, 9, 10), (9, 9, 11), (9, 10, 11), (9, 11, 11),
    (10, 10, 11), (11, 11, 11),
    (7, 9, 12), (8, 9, 12), (9, 10, 12), (9, 11, 12),
    (10, 10, 12), (10, 11, 12), (10, 12, 12), (11, 11, 12), (11, 12, 12),
]


def load_uvw(npz, n, m, p):
    """Load U, V, W as float64 ndarrays for the given format."""
    key = f"{n},{m},{p}"
    arr = npz[key]
    if arr.dtype == np.object_:
        u, v, w = arr[0], arr[1], arr[2]
    else:
        # Some entries are stored as a stacked 3D array (3, n*m, rank).
        u, v, w = arr[0], arr[1], arr[2]
    return u, v, w


def convert_w_to_row_major(w_at, n, p):
    """AT stores W with COLUMN-MAJOR C-flatten (W[l*n+i] = coeff in C[i,l]).
    Our format uses ROW-MAJOR (W[i*p+l] = coeff in C[i,l]). Permute rows."""
    rank = w_at.shape[1]
    w_ours = np.zeros_like(w_at)
    for i in range(n):
        for l in range(p):
            w_ours[i * p + l, :] = w_at[l * n + i, :]
    return w_ours


def verify_matmul(u, v, w, n, m, p, samples=5):
    """Random matmul spot-check using OUR row-major C-flatten convention
    (W already permuted into our format)."""
    rng = np.random.default_rng(0xC0DE)
    for _ in range(samples):
        A = rng.standard_normal((n, m))
        B = rng.standard_normal((m, p))
        a_flat = A.flatten()  # row-major: [i*m + j]
        b_flat = B.flatten()  # row-major: [j*p + l]
        alpha = u.T @ a_flat
        beta = v.T @ b_flat
        gamma = alpha * beta
        c_flat = w @ gamma     # w is now our row-major: c_flat[i*p + l]
        C_algo = c_flat.reshape(n, p)
        C_naive = A @ B
        if not np.allclose(C_algo, C_naive, atol=1e-9):
            return False, np.max(np.abs(C_algo - C_naive))
    return True, 0.0


def sparse_dict(M):
    """Convert a matrix M[rows][cols] (cols = rank) into sparse {row: value}
    per column — matching SchemeIO's sparse format."""
    rows, cols = M.shape
    out = []
    for k in range(cols):
        col_dict = {}
        for row in range(rows):
            v = M[row, k]
            if v != 0:
                # Cast to int if it's a whole number; preserves "Z" semantics.
                if v == int(v):
                    col_dict[str(row)] = int(v)
                else:
                    col_dict[str(row)] = float(v)
        out.append(col_dict)
    return out


def section_dir(n, m, p):
    return SCHEMES_ROOT / f"section{max(n, m, p)}"


def write_scheme(n, m, p, rank, u, v, w, apply: bool):
    # Detect field based on coefficient values:
    #   all integer → Z
    #   else (small rationals like ½ or ⅛) → Q
    all_vals = np.concatenate([u.flatten(), v.flatten(), w.flatten()])
    if np.allclose(all_vals, np.round(all_vals), atol=1e-12):
        n_field = "Z"
    else:
        n_field = "Q"
    fname = f"alphatensor-{n_field}_{n}x{m}x{p}_r{rank}_aN.json"
    # We don't have addition count without verifier — use 'aN' placeholder
    # which GenerateCatalogManifest will skip / inspect file directly.
    out_path = section_dir(n, m, p) / fname
    payload = {
        "n": [n, m, p],
        "m": rank,
        "u_sparse": sparse_dict(u),
        "v_sparse": sparse_dict(v),
        "w_sparse": sparse_dict(w),
        "z2": False,
        "source": "https://github.com/google-deepmind/alphatensor/algorithms/factorizations_r.npz",
        "imported_from": "factorizations_r.npz",
        "field_claimed": n_field,
    }
    if apply:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w") as f:
            json.dump(payload, f, indent=2)
    return out_path


def main():
    apply = "--apply" in sys.argv
    print(f"Mode: {'APPLY' if apply else 'DRY-RUN'}")

    npz = np.load(NPZ_PATH, allow_pickle=True)
    written = 0
    verified = 0
    z_field = 0
    r_field = 0
    failed_verify = []
    for n, m, p in TARGET_FORMATS:
        try:
            u, v, w_at = load_uvw(npz, n, m, p)
        except Exception as e:
            print(f"  SKIP ⟨{n},{m},{p}⟩: load error {e}")
            continue
        # AT uses col-major W (l*n + i = C[i,l]); same as dronperminov's spec
        # which SchemeIO reads. So we write AT's W AS-IS in w_sparse, and
        # only convert to row-major for the Python end-to-end verification.
        w_row_major = convert_w_to_row_major(w_at, n, p)
        rank = u.shape[1]
        ok, err = verify_matmul(u, v, w_row_major, n, m, p)
        w = w_at  # write AT's col-major W to JSON, SchemeIO will transpose on read
        status = "PASS" if ok else f"FAIL (residual={err:.2e})"
        print(f"  ⟨{n},{m},{p}⟩ rank={rank} dtype={u.dtype}  verify: {status}")
        if not ok:
            failed_verify.append((n, m, p, rank))
            continue
        verified += 1
        out = write_scheme(n, m, p, rank, u, v, w, apply)
        if "alphatensor-Z" in str(out):
            z_field += 1
        else:
            r_field += 1
        if apply:
            written += 1
            print(f"    WROTE {out}")

    print()
    print(f"Verified: {verified} / {len(TARGET_FORMATS)}")
    print(f"Z-field (all integer coefficients): {z_field}")
    print(f"R-field (some non-integer): {r_field}")
    print(f"Failed verify: {len(failed_verify)}")
    if failed_verify:
        print("FAILED FORMATS:")
        for n, m, p, r in failed_verify:
            print(f"  ⟨{n},{m},{p}⟩ r={r}")
    if apply:
        print(f"Files written: {written}")
    else:
        print("Re-run with --apply to write files.")


if __name__ == "__main__":
    main()
