#!/usr/bin/env python3
"""INDEPENDENT verification of the constructed HK ⟨2,p,n⟩ schemes.

Purpose: a correctness certificate sharing NOTHING with the Java pipeline that
generated the schemes — different language, different JSON parsing, different
arithmetic (stdlib Fraction), written against the documented disk format. A
convention mistake here can only produce a false NEGATIVE (the exact bilinear
identity is convention-sensitive), so a clean pass over the full set is a
genuinely independent proof that every published scheme computes ⟨n,m,p⟩
matrix multiplication exactly, at the rank its filename and metadata claim.

Disk format (SchemeIO):
  - "n": [n, m, p]; "m": rank r.
  - dense: "u" = [r][n·m] (A row-major i·m+l), "v" = [r][m·p] (B row-major
    l·p+j), "w" = [r][n·p] with positions in dronperminov COL-MAJOR (stored
    index d = j·n + i for output cell (i,j)).
  - sparse: "u_sparse"/"v_sparse"/"w_sparse" = {"k": {"i": [pos...],
    "c": [coef...]}}, same position conventions.
  - coefficients: ints, floats, or exact "p/q" strings.

Checks per scheme:
  1. rank == r == filename token;
  2. the exact identity  Σ_k U[a][k] V[b][k] W[c][k] == T[a][b][c]  over ℚ
     for ALL (a, b, c) — via sparse accumulation;
  3. the HK formula comparison reported (informational).

Run:  python3 verify_constructed_independent.py [schemes_root]
"""
import json
import math
import re
import sys
import glob
from fractions import Fraction


def coef(x):
    if isinstance(x, str):
        num, den = x.split('/')
        return Fraction(int(num), int(den))
    if isinstance(x, float):
        f = Fraction(x).limit_denominator(10**9)
        assert abs(float(f) - x) < 1e-12, f"unrecoverable float {x}"
        return f
    return Fraction(x)


def factors(d, r, dim_u, dim_v, dim_w, n, p):
    """Return U, V, W as dicts: U[k] = {pos: Fraction}, W positions ROW-major."""
    def parse_dense(key, dim):
        rows = d[key]
        assert len(rows) == r, f"{key}: {len(rows)} != r={r}"
        out = []
        for k in range(r):
            row = rows[k]
            assert len(row) == dim
            out.append({i: coef(v) for i, v in enumerate(row) if coef(v) != 0})
        return out

    def parse_sparse(key, dim):
        sp = d[key]
        out = []
        for k in range(r):
            e = sp[str(k)]
            out.append({i: coef(c) for i, c in zip(e['i'], e['c'])})
            assert all(0 <= i < dim for i in out[-1])
        return out

    U = parse_dense('u', dim_u) if 'u' in d else parse_sparse('u_sparse', dim_u)
    V = parse_dense('v', dim_v) if 'v' in d else parse_sparse('v_sparse', dim_v)
    Wraw = parse_dense('w', dim_w) if 'w' in d else parse_sparse('w_sparse', dim_w)
    # W disk positions are col-major (d = j*n + i); convert to row-major i*p + j.
    W = []
    for k in range(r):
        W.append({(dpos % n) * p + (dpos // n): c for dpos, c in Wraw[k].items()})
    return U, V, W


def verify(path):
    d = json.load(open(path))
    n, m, p = d['n']
    r = d['m']
    fname = path.split('/')[-1]
    mm = re.match(r'(\d+)x(\d+)x(\d+)-r(\d+)-', fname)
    assert mm and [int(mm.group(i)) for i in (1, 2, 3)] == [n, m, p], f"shape/filename mismatch {fname}"
    assert int(mm.group(4)) == r, f"rank/filename mismatch {fname}"

    dim_u, dim_v, dim_w = n * m, m * p, n * p
    U, V, W = factors(d, r, dim_u, dim_v, dim_w, n, p)

    # Sparse accumulation of Σ_k U⊗V⊗W.
    acc = {}
    for k in range(r):
        for a, ua in U[k].items():
            for b, vb in V[k].items():
                uv = ua * vb
                for c, wc in W[k].items():
                    key = (a * dim_v + b) * dim_w + c
                    acc[key] = acc.get(key, Fraction(0)) + uv * wc

    # Expected tensor: T[(i*m+l), (l*p+j), (i*p+j)] = 1, else 0.
    ones = set()
    for i in range(n):
        for l in range(m):
            for j in range(p):
                ones.add(((i * m + l) * dim_v + (l * p + j)) * dim_w + (i * p + j))
    for key, v in acc.items():
        expected = 1 if key in ones else 0
        if v != expected:
            c = key % dim_w
            b = (key // dim_w) % dim_v
            a = key // dim_w // dim_v
            return False, f"cell a={a} b={b} c={c}: got {v}, expected {expected}"
    for key in ones:
        if acc.get(key, Fraction(0)) != 1:
            return False, f"missing tensor one at flat {key}"

    formula = math.ceil((3 * m * p + max(m, p)) / 2)
    return True, f"r={r} formula={formula} delta={r - formula:+d}"


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else 'src/main/resources/schemes/constructed'
    files = sorted(glob.glob(root + '/**/*.json', recursive=True))
    print(f"independent verification of {len(files)} schemes under {root}")
    bad = 0
    at_formula = 0
    for f in files:
        ok, msg = verify(f)
        if not ok:
            bad += 1
            print(f"  FAIL {f.split('/')[-1]}: {msg}")
        elif 'delta=+0' in msg:
            at_formula += 1
    print(f"RESULT: {len(files) - bad}/{len(files)} verified exactly over Q "
          f"({at_formula} at the HK formula); {bad} failures")
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
