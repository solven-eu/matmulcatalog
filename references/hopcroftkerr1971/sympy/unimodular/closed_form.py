"""Closed-form seam test: tail C = [I_{m-r} | ones] over [I_r], across all band (p,m)."""
from comb_seam_coldfs import det, all_windows_unimodular

def build_closed(p, m):
    n = p + m
    r = p % m
    body = p - r
    B = [[1 if (c < body and c % m == u) else 0 for c in range(p)] for u in range(m)]
    if r:
        for u in range(m - r):           # first m-r rows: tooth at u, plus ones beyond col m-r-1
            for j in range(r):
                if j == u or j >= m - r:
                    B[u][body + j] = 1
        for i in range(r):               # last r rows: identity
            B[m - r + i][body + i] = 1
    rows0 = [[1 if cc == rr else 0 for cc in range(p)] for rr in range(p)]
    return B if all_windows_unimodular(rows0 + B, p, n) else None

fails = []
for p in range(3, 33):
    line = []
    for n in range(p + 1, min(2 * p - 1, 32) + 1):
        m = n - p
        ok = build_closed(p, m) is not None
        line.append('Y' if ok else '.')
        if not ok:
            fails.append((p, n))
    print(f"p={p:2d}: {''.join(line)}", flush=True)
print("FAILS:", fails if fails else "none — closed form is all-unimodular everywhere", flush=True)
