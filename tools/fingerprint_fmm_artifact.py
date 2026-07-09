#!/usr/bin/env python3
"""Block-support fingerprint of an FMM maple tensor artifact (2026-07-09).

Usage: fingerprint_fmm_artifact.py <tensor.mpl> <aSplit> <bSplit> <cSplit>

Partitions the scheme's products by which 2x2 blocks (per operand, under the
given axis splits) their U/V/W supports touch. Used to test whether an
artifact is display-faithful (clean per-leaf clusters) and to expose the
absorbing-class device: ⟨17,17,19⟩:3266 under splits (8,8,9) shows 7 main
clusters sized BELOW their standalone leaf ranks (422 vs R(8,9,9)=430, ...)
plus 83 residue products — 84 products of leaf work replaced by 83 shared
corrections, saving exactly the 1 that separates Σleaves=3267 from 3266.
Contrast ⟨19,19,22⟩:4536, which is basis-mixed under every split (flip-graph
output; not display-faithful). See research/DISCOVERIES_PENDING_ANALYSIS.md
2026-07-09 and references/fmm-artifact-audit.md.
"""
import sys
from collections import Counter
sys.path.insert(0, __file__.rsplit('/', 1)[0])
import import_fmm_maple as m


def blocks(mat, rsplit, csplit):
    s = set()
    for i, row in enumerate(mat):
        for j, x in enumerate(row):
            if x != 0:
                s.add((0 if i < rsplit else 1, 0 if j < csplit else 1))
    return frozenset(s)


def main():
    path, ra, rb, rc = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    triads = m.parse_triads(open(path).read())
    print(len(triads), 'triads')
    sig = Counter()
    pure = 0
    for (A, B, C) in triads:
        a = blocks(A, ra, rb)
        b = blocks(B, rb, rc)
        rs, cs = (rc, ra) if len(C) != len(A) else (ra, rc)
        c = blocks(C, rs, cs)
        if len(a) == len(b) == len(c) == 1:
            pure += 1
        sig[(tuple(sorted(a)), tuple(sorted(b)), tuple(sorted(c)))] += 1
    print(f"pure single-block products: {pure}/{len(triads)}; {len(sig)} signatures")
    fmt = lambda s: '+'.join(f"{i}{j}" for i, j in s)
    for key, cnt in sorted(sig.items(), key=lambda x: -x[1]):
        a, b, c = key
        print(f"{cnt:5d}  A[{fmt(a)}] B[{fmt(b)}] C[{fmt(c)}]")


if __name__ == '__main__':
    main()
