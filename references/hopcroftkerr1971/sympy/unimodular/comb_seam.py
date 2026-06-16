"""Comb + seam search: all-unimodular Lemma-1 augmentations for every band (p, m).

Rows: comb of period m, phase t, on columns [0, p-r-1] (r = p mod m), plus
free 0/1 teeth on the last r columns, found by DFS with exact window checks.
"""
import sys

def det(M):
    M = [row[:] for row in M]; n = len(M); sign = 1; prev = 1
    for k in range(n-1):
        if M[k][k] == 0:
            for rr in range(k+1, n):
                if M[rr][k] != 0: M[k], M[rr] = M[rr], M[k]; sign = -sign; break
            else: return 0
        for i in range(k+1, n):
            for j in range(k+1, n):
                M[i][j] = (M[i][j]*M[k][k] - M[i][k]*M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n-1][n-1]

def all_windows_unimodular(rows, p, n):
    return all(abs(det([rows[(j+t)%n] for t in range(p)])) == 1 for j in range(n))

def build(p, m, ternary_tail=False):
    n = p + m
    r = p % m
    body = p - r
    rows0 = [[1 if c == rr else 0 for c in range(p)] for rr in range(p)]
    def comb(t):
        return [1 if (c < body and c % m == t) else 0 for c in range(p)]
    if r == 0:
        B = [comb(t) for t in range(m)]
        return B if all_windows_unimodular(rows0 + B, p, n) else None
    # DFS over tail assignments (0/1 or ternary on last r columns, per row)
    from itertools import product as iproduct
    alphabet = (0, 1, -1) if ternary_tail else (0, 1)
    sol = []
    def dfs(B):
        if sol: return
        if len(B) == m:
            if all_windows_unimodular(rows0 + B, p, n): sol.append([row[:] for row in B])
            return
        t = len(B)
        base = comb(t)
        for tail in iproduct(alphabet, repeat=r):
            row = base[:body] + list(tail)
            Bn = B + [row]
            k = len(Bn)
            # prune: leading k-minor
            if abs(det([rw[:k] for rw in Bn])) != 1: continue
            dfs(Bn)
    dfs([])
    return sol[0] if sol else None

fails = []
for p in range(3, 33):
    line = []
    for n in range(p+1, min(2*p-1, 32)+1):
        m = n - p
        B = build(p, m)
        if B is None and (p % m):
            B = build(p, m, ternary_tail=True)
        line.append('Y' if B else '.')
        if not B: fails.append((p, n))
    print(f"p={p:2d}: {''.join(line)}", flush=True)
print("FAILS:", fails if fails else "none — every band (p,n) has an all-unimodular augmentation", flush=True)
