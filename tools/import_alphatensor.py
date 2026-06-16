#!/usr/bin/env python3
"""Import AlphaTensor matmul factorizations into the dronperminov JSON format.

Source data: github.com/google-deepmind/alphatensor/algorithms/
- factorizations_f2.npz  — 20 schemes over GF(2)
- factorizations_r.npz   — 93 schemes over standard arithmetic

Each .npz key is a string like "2,2,2" or "4,5,5" mapping to either:
- a dense 3D ndarray of shape (3, dim, rank) for cubic formats, or
- an object ndarray of length 3 holding the U, V, W matrices separately for
  non-cubic formats.

Each U/V/W is shape (slot_dim, rank); rows are flatten positions, columns are
multiplications. Convention (per the AlphaEvolve notebook from the same group,
matching dronperminov):
- U row-major A-flatten:  i·m + j
- V row-major B-flatten:  j·p + k
- W col-major C-flatten:  k·n + i

So the W from .npz can be saved as-is in the dronperminov JSON (which also
uses col-major C-flatten).

Usage:
    python3 tools/import_alphatensor.py [output_dir]
"""
from __future__ import annotations

import json
import os
import sys
import urllib.request

import numpy as np

F2_URL = "https://github.com/google-deepmind/alphatensor/raw/main/algorithms/factorizations_f2.npz"
R_URL = "https://github.com/google-deepmind/alphatensor/raw/main/algorithms/factorizations_r.npz"


def ensure_local(path: str, url: str) -> str:
    if not os.path.exists(path):
        print(f"downloading {url} → {path}", file=sys.stderr)
        urllib.request.urlretrieve(url, path)
    return path


def extract_factors(npz: np.lib.npyio.NpzFile, key: str
                    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    arr = npz[key]
    if arr.dtype == object:
        # Non-cubic: shape (3,), each element a (dim, rank) ndarray.
        return arr[0], arr[1], arr[2]
    # Cubic: shape (3, dim, rank).
    return arr[0], arr[1], arr[2]


def numpy_to_list(arr: np.ndarray) -> list[list[int]]:
    out: list[list[int]] = []
    for row in arr:
        out.append([int(v) for v in row])
    return out


def write_json(out_path: str, n: int, m: int, p: int, rank: int,
               U: np.ndarray, V: np.ndarray, W: np.ndarray, z2: bool) -> None:
    # AlphaTensor stores each factor as (dim, rank). dronperminov stores
    # them as rank × dim (per-multiplication rows). Transpose.
    u_rd = U.T
    v_rd = V.T
    w_rd = W.T  # W is already in col-major C-flatten per AlphaEvolve's convention.
    data = {
        "n": [n, m, p],
        "m": rank,
        "z2": z2,
        "field": "F_2" if z2 else "Z",
        "source": "AlphaTensor",
        "complex": False,
        "u": numpy_to_list(u_rd),
        "v": numpy_to_list(v_rd),
        "w": numpy_to_list(w_rd),
    }
    with open(out_path, "w") as f:
        json.dump(data, f, indent=2)


def import_npz(npz_path: str, z2: bool, output_dir: str) -> int:
    npz = np.load(npz_path, allow_pickle=True)
    count = 0
    for key in npz.files:
        parts = key.split(",")
        if len(parts) != 3:
            print(f"  skip unexpected key {key!r}", file=sys.stderr)
            continue
        n, m, p = int(parts[0]), int(parts[1]), int(parts[2])
        U, V, W = extract_factors(npz, key)
        rank = U.shape[1]
        # Sanity-check shapes match the format.
        expected = {(n * m, rank), (m * p, rank), (n * p, rank)}
        actual = {U.shape, V.shape, W.shape}
        if actual != expected:
            # Sometimes the convention is W as (n·p, rank) but stored as (p·n, rank)
            # — same total size, equivalent flatten. We'll accept if sizes match.
            if not (U.shape[0] == n * m and V.shape[0] == m * p and W.shape[0] == n * p):
                print(f"  WARN {key}: shapes {actual} don't match expected {expected}",
                      file=sys.stderr)
                continue

        field_slug = "F2" if z2 else "Z"
        out_name = f"alphatensor-{field_slug}_{n}x{m}x{p}_r{rank}.json"
        out_path = os.path.join(output_dir, out_name)
        write_json(out_path, n, m, p, rank, U, V, W, z2)
        count += 1
    return count


def main() -> int:
    output_dir = sys.argv[1] if len(sys.argv) > 1 else "src/main/resources/schemes"
    os.makedirs(output_dir, exist_ok=True)

    f2_path = ensure_local("/tmp/factorizations_f2.npz", F2_URL)
    r_path = ensure_local("/tmp/factorizations_r.npz", R_URL)

    nf2 = import_npz(f2_path, z2=True, output_dir=output_dir)
    print(f"  imported {nf2} F2 schemes", file=sys.stderr)
    nr = import_npz(r_path, z2=False, output_dir=output_dir)
    print(f"  imported {nr} R schemes", file=sys.stderr)
    print(f"\nDone. Total: {nf2 + nr} AlphaTensor schemes written to {output_dir}/",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
