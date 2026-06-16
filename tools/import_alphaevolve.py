#!/usr/bin/env python3
"""Import AlphaEvolve matmul decompositions into the dronperminov JSON format.

Usage:
    python3 tools/import_alphaevolve.py [notebook_path] [output_dir]

Defaults:
    notebook_path = /tmp/alphaevolve.ipynb (or downloaded into it)
    output_dir = src/main/resources/schemes

What it does:
    1. Walks the AlphaEvolve mathematical_results.ipynb notebook
       (https://github.com/google-deepmind/alphaevolve_results).
    2. For each `## Rank-R decomposition of <n,m,p> over FIELD` section:
       a. Locates the code cell that defines `decomposition_NMP = (np.array(...), ...)`.
       b. Executes that cell in a controlled namespace to obtain the U, V, W arrays.
       c. Runs the upstream `verify_tensor_decomposition(...)` from cell 4 as a
          sanity check before writing.
       d. Writes a dronperminov-format JSON file
          (alphaevolve_<n><m><p>_r<R>.json) plus a `.field` annotation
          (the AlphaEvolve fields include Z, 0.5*Z, etc., which we record but
          which the strict dronperminov format doesn't capture beyond `z2`).

Conventions:
    AlphaEvolve W uses col-major C-flatten (k*n + i, where k is C-col, i is C-row).
    That matches dronperminov's W convention — no transposition needed at this step.
    `SchemeIO.read` on the Java side then converts col-major → our internal row-major.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import urllib.request
from typing import Any

import numpy as np

NOTEBOOK_URL = "https://raw.githubusercontent.com/google-deepmind/alphaevolve_results/main/mathematical_results.ipynb"
DEFAULT_NOTEBOOK_PATH = "/tmp/alphaevolve.ipynb"

# Match: "## Rank-32 decomposition of <2,4,5> over 0.5*Z" — the notebook also has
# `##Rank-...` (no space), so allow zero-or-more whitespace after `##`.
HEADER_RE = re.compile(
    r"##\s*Rank-(?P<rank>\d+)\s+decomposition\s+of\s*<(?P<n>\d+)\s*,\s*(?P<m>\d+)\s*,\s*(?P<p>\d+)>\s+over\s+(?P<field>\S+)",
    re.IGNORECASE,
)


def ensure_notebook(path: str) -> str:
    if os.path.exists(path):
        return path
    print(f"downloading notebook to {path} ...", file=sys.stderr)
    urllib.request.urlretrieve(NOTEBOOK_URL, path)
    return path


def field_label_for_filename(field: str) -> str:
    """Turn '0.5*Z' / 'Z' into a filesystem-safe slug for the JSON filename."""
    return field.replace("*", "x").replace("/", "_").replace(" ", "_")


def is_z2(field: str) -> bool:
    return field.replace(" ", "").lower() in {"z/2", "z_2", "z2", "f2", "f_2", "gf(2)"}


def numpy_to_int_list(arr: np.ndarray) -> list[list[Any]]:
    """Encode an array of real numbers as nested int lists where possible.

    AlphaEvolve uses entries like 0.5, so half-integer coefficients survive as
    floats. We try int first; on non-integer, fall back to float.
    """
    out: list[list[Any]] = []
    for row in arr:
        encoded: list[Any] = []
        for v in row:
            f = float(v)
            if f == int(f):
                encoded.append(int(f))
            else:
                encoded.append(f)
        out.append(encoded)
    return out


def numpy_complex_to_pair_list(arr: np.ndarray) -> list[list[list[Any]]]:
    """Encode complex coefficients as nested `[re, im]` pairs.

    Real and imaginary parts each go through the int/float coercion of
    {@link numpy_to_int_list}, so e.g. 0.5+0.5j becomes [0.5, 0.5] and
    1+0j becomes [1, 0].
    """
    def coerce(v: float) -> Any:
        if v == int(v):
            return int(v)
        return v

    out: list[list[list[Any]]] = []
    for row in arr:
        encoded: list[list[Any]] = []
        for c in row:
            encoded.append([coerce(float(c.real)), coerce(float(c.imag))])
        out.append(encoded)
    return out


def is_complex(arr: np.ndarray) -> bool:
    return np.issubdtype(arr.dtype, np.complexfloating)


def write_dronperminov_json(path: str, n: int, m: int, p: int, rank: int,
                            U: np.ndarray, V: np.ndarray, W: np.ndarray,
                            field: str, source: str = "AlphaEvolve") -> None:
    """Write the scheme in the dronperminov JSON format (plus our extensions).

    U, V, W shapes from AlphaEvolve are (dim, rank). dronperminov stores them
    as rank × dim ([[m1 coefficients...], [m2 coefficients...], ...]). Transpose
    to match.

    For complex coefficients (AlphaEvolve's ⟨4,4,4⟩=48 over 0.5*C), each value
    becomes a `[re, im]` pair and the JSON sets `"complex": true` — strictly an
    extension of the dronperminov format. Java consumers that don't handle
    complex should look at the `complex` flag and skip.
    """
    # Transpose to rank × dim:
    u_rd = U.T  # shape (rank, n*m)
    v_rd = V.T  # shape (rank, m*p)
    w_rd = W.T  # shape (rank, p*n) — col-major C-flatten, same as dronperminov w
    complex_scheme = is_complex(U) or is_complex(V) or is_complex(W)
    data: dict[str, Any] = {
        "n": [n, m, p],
        "m": rank,
        "z2": is_z2(field),
        "field": field,  # extra annotation outside strict spec
        "source": source,
        "complex": complex_scheme,
    }
    if complex_scheme:
        data["u"] = numpy_complex_to_pair_list(u_rd)
        data["v"] = numpy_complex_to_pair_list(v_rd)
        data["w"] = numpy_complex_to_pair_list(w_rd)
    else:
        data["u"] = numpy_to_int_list(u_rd)
        data["v"] = numpy_to_int_list(v_rd)
        data["w"] = numpy_to_int_list(w_rd)
    with open(path, "w") as f:
        json.dump(data, f, indent=2)


def run_cell(code: str, namespace: dict) -> None:
    """Execute a notebook code cell into the given namespace.

    Filters out IPython magics (lines starting with `#@title`) and prints to
    avoid spam in the importer log.
    """
    cleaned = []
    for line in code.split("\n"):
        if line.startswith("#@title"):
            continue
        cleaned.append(line)
    exec("\n".join(cleaned), namespace)


def collect_decompositions(nb: dict) -> list[dict]:
    """Walk the notebook cells; return a list of {n, m, p, rank, field, U, V, W} dicts."""
    cells = nb["cells"]
    # Establish a baseline namespace with numpy and the verification function.
    namespace: dict[str, Any] = {}
    # Run the verification function cell so we can call it during import.
    for c in cells[:5]:
        if c["cell_type"] == "code":
            run_cell("".join(c["source"]), namespace)

    found: list[dict] = []
    for i, cell in enumerate(cells):
        if cell["cell_type"] != "markdown":
            continue
        src = "".join(cell["source"])
        m = HEADER_RE.search(src)
        if not m:
            continue
        n_ = int(m.group("n"))
        m_ = int(m.group("m"))
        p_ = int(m.group("p"))
        rank = int(m.group("rank"))
        field = m.group("field")

        # Find the next code cell (skipping any markdown notes).
        j = i + 1
        while j < len(cells) and cells[j]["cell_type"] != "code":
            j += 1
        if j == len(cells):
            print(f"  WARN: no code cell after header '{src.strip()[:80]}'", file=sys.stderr)
            continue
        code = "".join(cells[j]["source"])
        try:
            run_cell(code, namespace)
        except Exception as e:
            print(f"  ERROR executing cell {j} ({n_},{m_},{p_}/r={rank}): {e}", file=sys.stderr)
            continue

        # Look up the decomposition variable. The naming convention is
        # `decomposition_NMP` (digits concatenated).
        var_name = f"decomposition_{n_}{m_}{p_}"
        if var_name not in namespace:
            print(f"  WARN: expected variable '{var_name}' not in namespace", file=sys.stderr)
            continue
        decomp = namespace[var_name]
        if not (isinstance(decomp, tuple) and len(decomp) == 3):
            print(f"  WARN: variable '{var_name}' is not a 3-tuple", file=sys.stderr)
            continue
        U, V, W = decomp

        # Sanity-check via upstream verifier.
        try:
            namespace["verify_tensor_decomposition"](decomp, n_, m_, p_, rank)
        except AssertionError as e:
            print(f"  ERROR upstream verification failed for ⟨{n_},{m_},{p_}⟩/r={rank}: {e}",
                  file=sys.stderr)
            continue

        found.append(dict(n=n_, m=m_, p=p_, rank=rank, field=field, U=U, V=V, W=W))
    return found


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("notebook", nargs="?", default=DEFAULT_NOTEBOOK_PATH)
    ap.add_argument("output_dir", nargs="?", default="src/main/resources/schemes")
    args = ap.parse_args()

    nb_path = ensure_notebook(args.notebook)
    with open(nb_path) as f:
        nb = json.load(f)

    os.makedirs(args.output_dir, exist_ok=True)
    decompositions = collect_decompositions(nb)
    print(f"\nExtracted {len(decompositions)} decompositions.", file=sys.stderr)

    for d in decompositions:
        n, m, p, rank, field = d["n"], d["m"], d["p"], d["rank"], d["field"]
        field_slug = field_label_for_filename(field)
        fname = f"alphaevolve_{n}x{m}x{p}_r{rank}_{field_slug}.json"
        out_path = os.path.join(args.output_dir, fname)
        write_dronperminov_json(out_path, n, m, p, rank, d["U"], d["V"], d["W"], field)
        print(f"  wrote {out_path}", file=sys.stderr)

    print(f"\nDone. {len(decompositions)} schemes written to {args.output_dir}/", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
