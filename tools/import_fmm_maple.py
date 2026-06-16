#!/usr/bin/env python3
"""Import an fmm-lille .mpl.bz2 tensor file into a verified scheme JSON.

fmm-lille publishes the actual factor matrices for each catalogued
algorithm as `https://fmm.univ-lille.fr/<n>x<m>x<p>_tensor.mpl.bz2`.
The Maple file structures each rank-1 outer product as

    Tensor := TriadSet([
      Triad([Matrix(n, m, [[...]]), Matrix(m, p, [[...]]), Matrix(p, n, [[...]])]),
      ...
    ])

We parse the file (no Maple interpreter needed — the structure is
regular nested-list literals), convert U/V/W to our internal
row-major flatten, and emit a sparse dronperminov-format JSON
file under `src/main/resources/schemes/section{max_dim}/`.

Conventions per the Maple verification line:
- U_i = Matrix(n, m): coefficient of A[a, b] = U_i[a, b]
- V_i = Matrix(m, p): coefficient of B[a, b] = V_i[a, b]
- W_i = Matrix(p, n): coefficient of multiplication i in C[i_C, l_C] =
  W_i[l_C, i_C]  — note W is keyed by C-col first, matching dronperminov's
  col-major C-flatten.

Usage:
    python3 tools/import_fmm_maple.py 16x16x16            # download + import
    python3 tools/import_fmm_maple.py 6x6x6 8x8x8 9x9x9   # batch
"""
from __future__ import annotations

import argparse
import bz2
import json
import os
import re
import sys
import urllib.request

URL_PERMINOV_MIRROR = "https://raw.githubusercontent.com/dronperminov/FastMatrixMultiplication/master/schemes/known/tensor/{fmt}_tensor.mpl"
URL_FMM_UPSTREAM = "https://fmm.univ-lille.fr/{fmt}_tensor.mpl.bz2"

# Shapes for which Perminov's mirror is known to be STALE — we force-fetch
# from FMM upstream instead.
#
# Tracking issues / reasons:
#   "9x11x15": Perminov mirror has rank-981 file (and his status.json wrongly
#              reports 956 for Q, propagated from ZT). FMM upstream has the
#              fresher rank-972 decomposition. Reported as
#              https://github.com/dronperminov/FastMatrixMultiplication/issues/3
FORCE_FMM_UPSTREAM = {
    "9x11x15",
}
CACHE_DIR = "/tmp/fmm-maple-cache"
SCHEMES_DIR = "src/main/resources/schemes"


# ---------------------------------------------------------------------------
# Maple parser
# ---------------------------------------------------------------------------

import ast as _ast

# Whitelist of AST node types we accept for evaluating Maple-style
# unsimplified arithmetic tokens (e.g. `1-0`, `1*1`, `-1+(2-3)`).
_ARITH_AST_ALLOWED = (
    _ast.Expression, _ast.Constant, _ast.UnaryOp, _ast.BinOp,
    _ast.UAdd, _ast.USub, _ast.Add, _ast.Sub, _ast.Mult, _ast.Div,
)


def parse_number(tok: str) -> float:
    """Parse a Maple-style number token. Supports integers, p/q rationals,
    and unsimplified arithmetic expressions like `1-0`, `1*1`, `-1+0`
    that appear in some FMM tensor files. We compile the token into a
    restricted-AST expression containing only Constant/UnaryOp/BinOp
    nodes and evaluate it. Any other AST node (Name, Call, Attribute,
    etc.) raises — so the eval surface is bounded to pure arithmetic on
    integer/float literals.
    """
    tok = tok.strip()
    # p/q rational fast path (often appears without surrounding parens
    # in FMM Maple files; standard ast.parse handles it as Div).
    try:
        return float(tok)
    except ValueError:
        pass
    try:
        tree = _ast.parse(tok, mode="eval")
    except SyntaxError as e:
        raise ValueError(f"unparseable token: {tok!r} ({e})")
    for node in _ast.walk(tree):
        if not isinstance(node, _ARITH_AST_ALLOWED):
            raise ValueError(f"unparseable token: {tok!r} (disallowed {type(node).__name__})")
        if isinstance(node, _ast.Constant) and not isinstance(node.value, (int, float)):
            raise ValueError(f"unparseable token: {tok!r} (non-numeric constant)")
    # Safe to evaluate.
    return float(eval(compile(tree, "<tok>", "eval"), {"__builtins__": {}}, {}))


def parse_nested_list(s: str) -> list:
    """Parse a Maple-style nested list `[[a,b,c],[d,e,f],...]` into
    a list of lists of floats. Handles integer and `p/q` rational tokens.
    """
    s = s.strip()
    assert s[0] == "[" and s[-1] == "]", f"expected [...]; got {s[:30]}..."
    inner = s[1:-1].strip()
    if not inner:
        return []
    # Check whether it's a list of lists or a flat list.
    if inner[0] == "[":
        # nested
        rows = []
        depth = 0
        start = 0
        for i, c in enumerate(inner):
            if c == "[":
                if depth == 0:
                    start = i
                depth += 1
            elif c == "]":
                depth -= 1
                if depth == 0:
                    rows.append(parse_nested_list(inner[start : i + 1]))
        return rows
    # flat list of numbers
    return [parse_number(tok) for tok in split_top_level_commas(inner)]


def split_top_level_commas(s: str) -> list[str]:
    """Split on commas that aren't inside nested brackets/parens."""
    out, depth, start = [], 0, 0
    for i, c in enumerate(s):
        if c in "([":
            depth += 1
        elif c in ")]":
            depth -= 1
        elif c == "," and depth == 0:
            out.append(s[start:i])
            start = i + 1
    out.append(s[start:])
    return out


def find_matching(text: str, open_at: int, open_char: str, close_char: str) -> int:
    """Given an opening bracket at `open_at`, return the index of its matching closer."""
    assert text[open_at] == open_char
    depth = 1
    i = open_at + 1
    while i < len(text):
        c = text[i]
        if c == open_char:
            depth += 1
        elif c == close_char:
            depth -= 1
            if depth == 0:
                return i
        i += 1
    raise ValueError(f"no matching {close_char} for {open_char} at {open_at}")


def parse_matrix(text: str, start: int) -> tuple[int, int, list, int]:
    """Parse `Matrix(R, C, [[...]])` starting at `start`. Returns (rows, cols, data, end_idx)."""
    head = "Matrix("
    assert text[start:].startswith(head), f"expected Matrix( at {start}; got {text[start:start+30]!r}"
    paren_open = start + len(head) - 1  # the '('
    paren_close = find_matching(text, paren_open, "(", ")")
    content = text[paren_open + 1 : paren_close]
    # content = "R, C, [[...]]". Split on top-level commas.
    parts = split_top_level_commas(content)
    if len(parts) < 3:
        raise ValueError(f"Matrix arg count = {len(parts)}, content head = {content[:60]!r}")
    rows = int(parts[0].strip())
    cols = int(parts[1].strip())
    data_str = ",".join(parts[2:]).strip()  # rejoin in case data had commas... (it always does)
    data = parse_nested_list(data_str)
    return rows, cols, data, paren_close + 1


def parse_triads(text: str) -> list[tuple[list, list, list]]:
    """Find all `Triad([Matrix(...), Matrix(...), Matrix(...)])` blocks
    in the text and return them as a list of (U, V, W) matrix tuples.
    """
    triads = []
    i = 0
    triad_head = "Triad("
    while True:
        i = text.find(triad_head, i)
        if i == -1:
            break
        # Find the closing ) of this Triad(
        body_open = i + len(triad_head) - 1
        body_close = find_matching(text, body_open, "(", ")")
        body = text[body_open + 1 : body_close]
        # Body is `[Matrix(...), Matrix(...), Matrix(...)]`
        body_stripped = body.strip()
        assert body_stripped[0] == "[" and body_stripped[-1] == "]"
        inner = body_stripped[1:-1]
        # Find 3 Matrix( blocks
        matrices = []
        k = 0
        while len(matrices) < 3:
            k = inner.find("Matrix(", k)
            if k == -1:
                raise ValueError(f"Triad at {i}: found {len(matrices)}/3 matrices")
            rows, cols, data, k = parse_matrix(inner, k)
            matrices.append((rows, cols, data))
        triads.append(tuple(m[2] for m in matrices))
        i = body_close
    return triads


# ---------------------------------------------------------------------------
# Internal conversion + JSON writer
# ---------------------------------------------------------------------------

def convert_triads_to_sparse(triads: list[tuple[list, list, list]], n: int, m: int, p: int):
    """Build sparse u/v/w (per-multiplication position→coefficient maps) in
    dronperminov convention (W uses col-major C-flatten: position = j·n + i).
    """
    rank = len(triads)
    u_sparse = []
    v_sparse = []
    w_sparse = []
    nzU = nzV = nzW = 0
    for (umat, vmat, wmat) in triads:
        u_row, v_row, w_row = {}, {}, {}
        for i in range(n):
            for j in range(m):
                val = umat[i][j]
                if val != 0:
                    u_row[str(i * m + j)] = val
                    nzU += 1
        for j in range(m):
            for l in range(p):
                val = vmat[j][l]
                if val != 0:
                    v_row[str(j * p + l)] = val
                    nzV += 1
        # wmat is shape (p, n) per Maple; dronperminov w stores
        # position = j_C * n + i_C (col-major C-flatten), where
        # wmat[j_C][i_C] is the same entry — no extra transpose needed.
        for j_C in range(p):
            for i_C in range(n):
                val = wmat[j_C][i_C]
                if val != 0:
                    w_row[str(j_C * n + i_C)] = val
                    nzW += 1
        u_sparse.append(u_row)
        v_sparse.append(v_row)
        w_sparse.append(w_row)
    adds = (nzU - rank) + (nzV - rank) + (nzW - (n * p))
    return u_sparse, v_sparse, w_sparse, adds


def format_value(v: float):
    if v == int(v):
        return int(v)
    # Preserve simple p/q rationals exactly via numerator/denominator.
    # We don't currently have a Fraction type in our JSON; emit as float.
    return v


def write_sparse_json(out_path: str, n: int, m: int, p: int, rank: int,
                      u_sparse: list, v_sparse: list, w_sparse: list,
                      field: str, source: str) -> None:
    def coerce(d):
        return {k: format_value(v) for k, v in d.items()}

    data = {
        "n": [n, m, p],
        "m": rank,
        "z2": False,
        "field": field,
        "source": source,
        "u_sparse": [coerce(d) for d in u_sparse],
        "v_sparse": [coerce(d) for d in v_sparse],
        "w_sparse": [coerce(d) for d in w_sparse],
    }
    with open(out_path, "w") as f:
        json.dump(data, f, indent=2)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def fetch_cached(fmt: str) -> str:
    """Fetch the Maple tensor for ⟨n,m,p⟩.

    Default policy: Perminov's plain-text mirror first (he's the
    authoritative field-classifier; mirror skips bz2 decompression).
    Falls back to FMM upstream if the mirror is missing.

    Override: if {@code fmt} is in {@link FORCE_FMM_UPSTREAM}, skip
    Perminov entirely and go to FMM — used for shapes where Perminov's
    mirror is known to be stale.
    """
    os.makedirs(CACHE_DIR, exist_ok=True)
    local_mpl = os.path.join(CACHE_DIR, f"{fmt}_tensor.mpl")
    if os.path.exists(local_mpl):
        with open(local_mpl, encoding="utf-8", errors="replace") as f:
            return f.read()
    # Force-FMM exception path.
    if fmt in FORCE_FMM_UPSTREAM:
        print(f"  {fmt}: force-fetching FMM upstream (Perminov mirror known stale)",
              file=sys.stderr)
        return _fetch_fmm_bz2(fmt)
    # Try Perminov mirror (plain text, no decompression).
    perm_url = URL_PERMINOV_MIRROR.format(fmt=fmt)
    print(f"  downloading {perm_url}", file=sys.stderr)
    try:
        urllib.request.urlretrieve(perm_url, local_mpl)
        with open(local_mpl, encoding="utf-8", errors="replace") as f:
            return f.read()
    except Exception as e:
        print(f"  Perminov mirror miss ({e}); falling back to FMM bz2", file=sys.stderr)
    return _fetch_fmm_bz2(fmt)


def _fetch_fmm_bz2(fmt: str) -> str:
    local_bz2 = os.path.join(CACHE_DIR, f"{fmt}_tensor.mpl.bz2")
    if not os.path.exists(local_bz2):
        fmm_url = URL_FMM_UPSTREAM.format(fmt=fmt)
        print(f"  downloading {fmm_url}", file=sys.stderr)
        urllib.request.urlretrieve(fmm_url, local_bz2)
    with bz2.open(local_bz2, "rb") as f:
        return f.read().decode("utf-8", errors="replace")


def import_one(fmt_spec: str) -> None:
    m = re.match(r"(\d+)x(\d+)x(\d+)$", fmt_spec)
    if not m:
        raise ValueError(f"bad format spec: {fmt_spec} (expected NxMxP)")
    n, mm, p = int(m.group(1)), int(m.group(2)), int(m.group(3))
    max_dim = max(n, mm, p)

    text = fetch_cached(fmt_spec)
    triads = parse_triads(text)
    rank = len(triads)
    print(f"  parsed {rank} triads for ⟨{n},{mm},{p}⟩", file=sys.stderr)

    # Sanity check matrix shapes.
    if triads:
        ur, vr, wr = triads[0]
        assert len(ur) == n and len(ur[0]) == mm, f"U shape mismatch: {len(ur)}x{len(ur[0])} vs {n}x{mm}"
        assert len(vr) == mm and len(vr[0]) == p
        assert len(wr) == p and len(wr[0]) == n

    u_sparse, v_sparse, w_sparse, adds = convert_triads_to_sparse(triads, n, mm, p)

    # Field detection: if any coefficient is a non-integer (e.g. 1/16), call it Q.
    has_frac = any(v != int(v) for row in u_sparse + v_sparse + w_sparse for v in row.values())
    field = "Q" if has_frac else "Z"

    out_dir = os.path.join(SCHEMES_DIR, f"section{max_dim}")
    os.makedirs(out_dir, exist_ok=True)
    out_name = f"fmm-lille_{n}x{mm}x{p}_r{rank}_a{adds}.json"
    out_path = os.path.join(out_dir, out_name)
    write_sparse_json(out_path, n, mm, p, rank,
                      u_sparse, v_sparse, w_sparse,
                      field=field, source="fmm-lille")
    size_kb = os.path.getsize(out_path) // 1024
    print(f"  wrote {out_path} ({size_kb} KB, field={field})", file=sys.stderr)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("formats", nargs="+", help="format specs like '16x16x16' or '3x3x3'")
    args = ap.parse_args()
    for f in args.formats:
        try:
            import_one(f)
        except Exception as e:
            print(f"FAILED {f}: {e}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
