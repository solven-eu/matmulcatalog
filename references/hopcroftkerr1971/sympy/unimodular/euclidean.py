"""Euclidean-recursive all-unimodular augmentation: B(p,m) = comb + B(m, p%m)^T tail.

Derivation: window dets of [I_p; B] reduce to (leading, trailing, sliding)
minor families of B; eliminating the comb body's permutation part shows the
tail block C (m x r, r = p % m) must satisfy exactly the SAME families
transposed, i.e. C = B(m, r)^T. Recursion terminates at r = 0 (pure comb,
proven: sliding windows are permutation matrices). Pure 0/1, no search.
"""
from comb_seam_coldfs import all_windows_unimodular

def block(p, m):
    """The m x p augmented block, recursively."""
    r = p % m
    body = p - r
    B = [[1 if (c < body and c % m == u) else 0 for c in range(p)] for u in range(m)]
    if r:
        C = block(m, r)          # r x m
        for u in range(m):
            for j in range(r):
                B[u][body + j] = C[j][u]   # transpose into the tail
    return B

def build(p, n):
    m = n - p
    rows = [[1 if c == rr else 0 for c in range(p)] for rr in range(p)]
    return rows + block(p, m) if m else rows

fails = []
for p in range(3, 33):
    line = []
    for n in range(p + 1, min(2 * p - 1, 32) + 1):
        M = build(p, n)
        ok = all_windows_unimodular(M, p, n)
        line.append('Y' if ok else '.')
        if not ok:
            fails.append((p, n))
    print(f"p={p:2d}: {''.join(line)}", flush=True)
print("FAILS:", fails if fails else "none — Euclidean recursion is all-unimodular over the entire band range", flush=True)
