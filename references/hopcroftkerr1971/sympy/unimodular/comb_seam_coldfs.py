"""All-unimodular Lemma-1 augmentation: comb body + column-DFS seam.

Goal (task #11): an n x p matrix M = [I_p; B] (m = n-p augmented rows) with
EVERY cyclic p-row-window determinant = +-1, so the HK back-substitution
stays integer and the constructed schemes upgrade from Q to Z.

Reduction (derived 2026-06-11): window dets of [I; B] are exactly
  (a) leading k x k minors of B (cols 0..k-1),        k = 1..m
  (b) trailing w x w minors  (rows m-w.., cols p-w..), w = 1..m
  (c) sliding m x m minors over contiguous column windows.

Construction: comb body  B[u][c] = 1 iff c % m == u  on the first p-r
columns (r = p % m). PROVEN all-unimodular when r = 0 (sliding windows are
permutation matrices; leading/trailing reduce to identity blocks). For
r > 0, the last r columns ("tail") are free 0/1 bits found by DFS over
COLUMNS, pruning each completed sliding window and the trailing minors as
soon as they are determined. Exact integer arithmetic throughout.

Run: python3 comb_seam_coldfs.py [max_p]
"""
import sys


def det(M):
    M = [row[:] for row in M]
    n = len(M)
    sign = 1
    prev = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            for rr in range(k + 1, n):
                if M[rr][k] != 0:
                    M[k], M[rr] = M[rr], M[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[-1][-1]


def all_windows_unimodular(rows, p, n):
    return all(abs(det([rows[(j + t) % n] for t in range(p)])) == 1
               for j in range(n))


def build(p, m):
    """Return an all-unimodular m x p augmented block B, or None."""
    n = p + m
    r = p % m
    body = p - r

    def comb_col(c):
        return [1 if c % m == u else 0 for u in range(m)]  # column vector

    cols = [comb_col(c) for c in range(body)]  # body columns, fixed

    if r == 0:
        B = [[cols[c][u] for c in range(p)] for u in range(m)]
        assert all_windows_unimodular(
            [[1 if cc == rr else 0 for cc in range(p)] for rr in range(p)] + B,
            p, n)
        return B

    # Column-DFS over the r tail columns; candidates = all 2^m bit columns,
    # ordered comb-like-first for fast hits.
    cands = sorted(range(1 << m), key=lambda v: bin(v).count('1'))
    sol = []

    def sliding_ok(allcols, c_last):
        """Check the sliding window ENDING at column c_last (if complete)."""
        c0 = c_last - m + 1
        if c0 < 0:
            return True
        W = [[allcols[c0 + j][u] for j in range(m)] for u in range(m)]
        return abs(det(W)) == 1

    def dfs(tail):
        if sol:
            return
        allcols = cols + tail
        c_last = len(allcols) - 1
        if tail and not sliding_ok(allcols, c_last):
            return
        if len(tail) == r:
            # trailing minors (b) and leading minors (a); leading never touch
            # the tail (k <= m <= body) and are identity for the comb => +-1.
            for w in range(1, m + 1):
                W = [[allcols[p - w + j][u] for j in range(w)]
                     for u in range(m - w, m)]
                if abs(det(W)) != 1:
                    return
            B = [[allcols[c][u] for c in range(p)] for u in range(m)]
            rows0 = [[1 if cc == rr else 0 for cc in range(p)] for rr in range(p)]
            if all_windows_unimodular(rows0 + B, p, n):
                sol.append(B)
            return
        for v in cands:
            dfs(tail + [[(v >> u) & 1 for u in range(m)]])

    dfs([])
    return sol[0] if sol else None


def main():
    max_p = int(sys.argv[1]) if len(sys.argv) > 1 else 32
    fails = []
    for p in range(3, max_p + 1):
        line = []
        for nn in range(p + 1, min(2 * p - 1, 32) + 1):
            m = nn - p
            B = build(p, m)
            line.append('Y' if B else '.')
            if not B:
                fails.append((p, nn))
        print(f"p={p:2d}: {''.join(line)}", flush=True)
    print("FAILS:", fails if fails else
          "none — every band (p,n) has an all-unimodular augmentation",
          flush=True)


if __name__ == '__main__':
    main()
