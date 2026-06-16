"""
F2 sibling of import_alphatensor_npz.py: extracts every factorization
from DeepMind's factorizations_f2.npz, verifies each over GF(2),
writes properly-tagged scheme JSON files.

Same fix as the R-arithmetic import: AT stores W in column-major
C-flatten (W[l·n + i] = coef of C[i,l]); our SchemeIO expects that
and transposes on read. So we write W AS-IS from the npz.

Run:
  /tmp/strassen_venv/bin/python tools/import_alphatensor_f2_npz.py [--apply]
"""
import json
import sys
from pathlib import Path

import numpy as np


SCHEMES_ROOT = Path("src/main/resources/schemes")
NPZ_PATH = "/tmp/at_f2.npz"


def load_uvw(npz, n, m, p):
    key = f"{n},{m},{p}"
    arr = npz[key]
    if arr.dtype == np.object_:
        u, v, w = arr[0], arr[1], arr[2]
    else:
        u, v, w = arr[0], arr[1], arr[2]
    return u.astype(np.int64), v.astype(np.int64), w.astype(np.int64)


def convert_w_to_row_major(w_at, n, p):
    """AT col-major (l·n + i) → our row-major (i·p + l)."""
    w_ours = np.zeros_like(w_at)
    for i in range(n):
        for l in range(p):
            w_ours[i * p + l, :] = w_at[l * n + i, :]
    return w_ours


def verify_f2(u, v, w, n, m, p, samples=5):
    """Verify mod-2 matmul on random {0,1} inputs."""
    rng = np.random.default_rng(0xF2C0DE)
    for _ in range(samples):
        A = rng.integers(0, 2, size=(n, m))
        B = rng.integers(0, 2, size=(m, p))
        a_flat = A.flatten()
        b_flat = B.flatten()
        alpha = (u.T @ a_flat) % 2
        beta = (v.T @ b_flat) % 2
        gamma = (alpha * beta) % 2
        c_flat = (w @ gamma) % 2
        C_algo = c_flat.reshape(n, p)
        C_naive = (A @ B) % 2
        if not np.array_equal(C_algo, C_naive):
            return False
    return True


def sparse_dict(M):
    rows, cols = M.shape
    out = []
    for k in range(cols):
        col_dict = {}
        for row in range(rows):
            v = int(M[row, k])
            if v != 0:
                col_dict[str(row)] = v
        out.append(col_dict)
    return out


def section_dir(n, m, p):
    return SCHEMES_ROOT / f"section{max(n, m, p)}"


def write_scheme(n, m, p, rank, u, v, w, apply: bool):
    fname = f"alphatensor-F2_{n}x{m}x{p}_r{rank}_aN.json"
    out_path = section_dir(n, m, p) / fname
    payload = {
        "n": [n, m, p],
        "m": rank,
        "u_sparse": sparse_dict(u),
        "v_sparse": sparse_dict(v),
        "w_sparse": sparse_dict(w),
        "z2": True,
        "source": "https://github.com/google-deepmind/alphatensor/algorithms/factorizations_f2.npz",
        "imported_from": "factorizations_f2.npz",
        "field_claimed": "F2",
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
    written = verified = failed = 0
    failed_fmts = []
    for key in sorted(npz.files):
        n, m, p = map(int, key.split(","))
        try:
            u, v, w_at = load_uvw(npz, n, m, p)
        except Exception as e:
            print(f"  SKIP ⟨{n},{m},{p}⟩: {e}")
            continue
        # Verify with row-major W in Python, then write col-major W to JSON
        # (which SchemeIO transposes back to row-major on read).
        w_row = convert_w_to_row_major(w_at, n, p)
        ok = verify_f2(u, v, w_row, n, m, p)
        rank = u.shape[1]
        print(f"  ⟨{n},{m},{p}⟩ rank={rank}  verify: {'PASS' if ok else 'FAIL'}")
        if not ok:
            failed += 1
            failed_fmts.append((n, m, p, rank))
            continue
        verified += 1
        out = write_scheme(n, m, p, rank, u, v, w_at, apply)
        if apply:
            written += 1
    print()
    print(f"Verified: {verified}, Failed: {failed}, Written: {written}")
    if failed_fmts:
        print("FAILED:", failed_fmts)
    if not apply:
        print("Re-run with --apply to write files.")


if __name__ == "__main__":
    main()
