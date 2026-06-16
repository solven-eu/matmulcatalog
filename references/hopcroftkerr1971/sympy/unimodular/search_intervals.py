import sys
from functools import lru_cache

def det(M):
    M = [row[:] for row in M]; n = len(M); sign = 1; prev = 1
    for k in range(n-1):
        if M[k][k] == 0:
            for r in range(k+1, n):
                if M[r][k] != 0: M[k], M[r] = M[r], M[k]; sign = -sign; break
            else: return 0
        for i in range(k+1, n):
            for j in range(k+1, n):
                M[i][j] = (M[i][j]*M[k][k] - M[i][k]*M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n-1][n-1]

def ok_all_windows(rows, p, n):
    return all(abs(det([rows[(j+t)%n] for t in range(p)])) == 1 for j in range(n))

def search_intervals(p, m):
    n = p + m
    cands = []
    for s in range(p):
        for l in range(2, p - s + 1):
            cands.append([1 if s <= c < s + l else 0 for c in range(p)])
    rows0 = [[1 if c == r else 0 for c in range(p)] for r in range(p)]
    sol = []
    def dfs(B):
        if sol: return
        if len(B) == m:
            if ok_all_windows(rows0 + B, p, n): sol.append([r[:] for r in B])
            return
        k = len(B) + 1
        for cand in cands:
            Bn = B + [cand]
            if abs(det([r[:k] for r in Bn])) != 1: continue
            # prune: sliding minors that are already fully determined need all rows,
            # but trailing minor of the last (m-k) cols can't be checked yet; cheap
            # extra prune: rows must be pairwise distinct
            if cand in B: continue
            dfs(Bn)
    dfs([])
    return sol

for p in range(3, 17):
    line = []
    for m in range(1, p):
        line.append('Y' if search_intervals(p, m) else '.')
    print(f"p={p:2d}: m=1..{p-1}: {''.join(line)}", flush=True)
for (p, m) in [(6,3), (9,5), (12,8)]:
    s = search_intervals(p, m)
    if s:
        print(f"example p={p},m={m}: " + " | ".join("".join(map(str,r)) for r in s[0]), flush=True)
