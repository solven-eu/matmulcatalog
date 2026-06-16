"""(3,3,bridge-m) over ARBITRARY local rank-1 atoms — decisive block reduction.

Upgrades derive_33bridge.py from "no solution in an 80-atom candidate catalog"
(optimal-within-scope) to a structural verdict over the FULL local atom space:
the three new products may be ANY rank-1 bilinear forms u⊗v on the 6 a-vars
(rows i, b, j) x 6 x-vars (cols i, b, j) — no catalog.

Derivation (hand-checked, 2026-06-11; mirrors the working emitter pattern):
targets are the bridge sums S_ij = y_ij + y_bj, S_ji = y_ji + y_jb, and the
reusable set Sh is the 9 products of derive_33bridge.py (diagonals of i and j,
bridge diagonals, the Lemma-2 pair (i,b)'s three products). Block the 6x6
monomial space by rows RB = {a_i, a_b} vs {a_j} and cols CB = {x_i, x_b} vs
{x_j}. Sh touches only the (RB,CB) block (7 products) and the (j,j) block
(E_jj, F_jj); S_ij lives in (RB, j) with rank 2, S_ji in (j, CB) with rank 2.

Forced structure (each step is a small rank argument):
  - every atom k must have all four components u_RB, u_j, v_CB, v_j nonzero;
  - beta (the S_ji-side coefficients) has exactly two nonzero entries, say
    {1,2}, with d_1 ~ d_2 (d_k = u_k,RB (x) v_k,j); alpha likewise two, and
    NOT the same pair (else rank(S_ij block) <= 1), so wlog alpha on {1,3}
    with e_1 ~ e_3 (e_k = u_k,j (x) v_k,CB);
  - rank-2 block decompositions pin: u_1,RB, u_3,RB span U0 = the
    "virtual row" space span{a_i1+a_b1, a_i2+a_b2}; v_1,CB, v_2,CB span
    V0 = span{x_1i+x_1b, x_2i+x_2b}; v_1,j, v_3,j basis of the j-col plane;
    u_1,j, u_2,j basis of the j-row plane;
  - the (RB,CB) spill of the alpha- and beta-combos collapses to TWO dyads
        L (x) v_1,CB   with L = alpha_1 u_1,RB + alpha_3 sigma u_3,RB in U0,
        u_1,RB (x) R   with R = beta_1 v_1,CB + beta_2 lambda v_2,CB in V0,
    both provably NONZERO (L = 0 or R = 0 contradicts the rank-2 targets),
    and both must lie in span(Sh|RBCB).

Hence feasibility REQUIRES nonzero rank-1 elements in
    W = (U0 (x) V0)  intersect  span(Sh restricted to the RB x CB block).
W = 0  =>  (3,3,bridge-m) has NO 3-product completion by any local rank-1
atoms over this reusable set. W != 0 pins the dyad directions; stage 2 then
chases the remaining constraints ((j,j) diagonality + non-proportionality).

Run: python3 derive_33bridge_general.py    (exact Fraction arithmetic)
"""
from fractions import Fraction
from itertools import product as iproduct

# Coordinates: a-side e1..e4 = a_i1, a_i2, a_b1, a_b2 (the RB rows);
#              x-side f1..f4 = x_1i, x_2i, x_1b, x_2b (the CB cols).
E = [[Fraction(1 if i == k else 0) for i in range(4)] for k in range(4)]


def vec(*pairs):
    v = [Fraction(0)] * 4
    for idx, c in pairs:
        v[idx] += Fraction(c)
    return v


def dyad(a, x):
    return [a[p] * x[q] for p in range(4) for q in range(4)]


def shared_rbcb(m_bridge):
    """The 7 shared products restricted to the RB x CB block (16-dim)."""
    a_i1, a_i2, a_b1, a_b2 = 0, 1, 2, 3
    x_1i, x_2i, x_1b, x_2b = 0, 1, 2, 3
    E_ii = dyad(vec((a_i2, 1)), vec((x_2i, 1)))
    F_ii = dyad(vec((a_i1, 1)), vec((x_1i, 1)))
    if m_bridge == 1:
        diag1 = dyad(vec((a_b2, 1)), vec((x_1b, 1), (x_2b, 1)))           # A_bb
        diag2 = dyad(vec((a_b1, 1), (a_b2, -1)), vec((x_1b, 1)))          # B_bb
        # pair (i,b) methods (1,3): emits D, C on the difference + G
        d_a = vec((a_i1, 1), (a_b1, -1))
        d_a2 = vec((a_i2, 1), (a_b2, -1))
        d_x1 = vec((x_1i, 1), (x_1b, -1))
        d_x2 = vec((x_2i, 1), (x_2b, -1))
        D_ib = dyad(d_a, [p + q for p, q in zip(d_x1, d_x2)])
        C_ib = dyad([p - q for p, q in zip(d_a, d_a2)], d_x2)
        G_ib = dyad(vec((a_i1, 1), (a_b1, -1), (a_b2, 1)),
                    vec((x_1b, 1), (x_2i, -1), (x_2b, 1)))
    else:
        diag1 = dyad(vec((a_b1, 1), (a_b2, -1)), vec((x_2b, 1)))          # C_bb
        diag2 = dyad(vec((a_b1, 1)), vec((x_1b, 1), (x_2b, 1)))           # D_bb
        # pair (i,b) methods (2,3): emits A, B on the difference + G
        d_a2 = vec((a_i2, 1), (a_b2, -1))
        d_a1 = vec((a_i1, 1), (a_b1, -1))
        d_x1 = vec((x_1i, 1), (x_1b, -1))
        d_x2 = vec((x_2i, 1), (x_2b, -1))
        A_ib = dyad(d_a2, [p + q for p, q in zip(d_x1, d_x2)])
        B_ib = dyad([p - q for p, q in zip(d_a1, d_a2)], d_x1)
        G_ib = dyad(vec((a_i2, 1), (a_b2, -1), (a_b1, 1)),
                    vec((x_1i, 1), (x_1b, -1), (x_2b, -1)))
        D_ib, C_ib = A_ib, B_ib
    return [E_ii, F_ii, diag1, diag2, D_ib, C_ib, G_ib]


def rref(rows):
    rows = [r[:] for r in rows]
    mat, pivots = [], []
    for r in rows:
        for pc, pr in zip(pivots, mat):
            if r[pc] != 0:
                f = r[pc] / pr[pc]
                r = [ri - f * pi for ri, pi in zip(r, pr)]
        nz = next((i for i, v in enumerate(r) if v != 0), None)
        if nz is not None:
            pivots.append(nz)
            mat.append(r)
    return mat, pivots


def intersection(A, B):
    """Basis of span(A) intersect span(B), exact. Standard kernel trick."""
    # x in both spans <=> x = A^T s = B^T t  <=>  [A^T | -B^T] (s,t) = 0.
    na, nb = len(A), len(B)
    dim = len(A[0])
    # Build kernel of the dim x (na+nb) matrix M with columns A_i then -B_j.
    M = [[A[i][d] for i in range(na)] + [-B[j][d] for j in range(nb)]
         for d in range(dim)]
    mat, pivots = rref(M)
    free = [c for c in range(na + nb) if c not in pivots]
    basis = []
    for fc in free:
        sol = [Fraction(0)] * (na + nb)
        sol[fc] = Fraction(1)
        for pr, pc in reversed(list(zip(mat, pivots))):
            s = sum(pr[c] * sol[c] for c in range(pc + 1, na + nb))
            sol[pc] = -s / pr[pc]
        x = [sum(sol[i] * A[i][d] for i in range(na)) for d in range(dim)]
        if any(v != 0 for v in x):
            basis.append(x)
    # Independent subset.
    out, seen = [], []
    for x in basis:
        test, _ = rref(seen + [x])
        if len(test) > len(rref(seen)[0] if seen else []):
            pass
        out_mat, _ = rref(seen + [x])
        if len(out_mat) > len(seen and rref(seen)[0] or []):
            pass
        seen.append(x)
    mat2, _ = rref(basis)
    return mat2


def rank1_elements(Wbasis):
    """All rank-1 directions (as 4x4 matrices) in span(Wbasis), exact.

    dim W <= 2 in practice: solve the 2x2-minor conditions on s*W1 + t*W2.
    Returns representative (u, v) factor pairs (projective)."""
    if not Wbasis:
        return []
    out = []

    def factor(Mflat):
        M = [Mflat[4 * r:4 * r + 4] for r in range(4)]
        nzr = [r for r in range(4) if any(M[r][c] != 0 for c in range(4))]
        if not nzr:
            return None
        r0 = nzr[0]
        for r in nzr[1:]:
            # rows must be proportional
            ratio = None
            for c in range(4):
                if M[r0][c] == 0 and M[r][c] == 0:
                    continue
                if M[r0][c] == 0:
                    return None
                rr = M[r][c] / M[r0][c]
                if ratio is None:
                    ratio = rr
                elif rr != ratio:
                    return None
        u = [M[r][next(c for c in range(4) if M[r0][c] != 0)] / M[r0][
            next(c for c in range(4) if M[r0][c] != 0)] for r in range(4)]
        v = M[r0]
        return u, v

    if len(Wbasis) == 1:
        f = factor(Wbasis[0])
        return [f] if f else []
    if len(Wbasis) == 2:
        W1, W2 = Wbasis
        # rank(s W1 + t W2) <= 1: 2x2 minors are quadratics in (s, t); check
        # t = 1 with s rational roots of the gcd, plus s = 1, t = 0.
        import math

        def minors(s, t):
            M = [s * a + t * b for a, b in zip(W1, W2)]
            G = [M[4 * r:4 * r + 4] for r in range(4)]
            vals = []
            for r1 in range(4):
                for r2 in range(r1 + 1, 4):
                    for c1 in range(4):
                        for c2 in range(c1 + 1, 4):
                            vals.append(G[r1][c1] * G[r2][c2]
                                        - G[r1][c2] * G[r2][c1])
            return vals

        # Collect the minor polynomials in s (t=1) symbolically: each minor is
        # quadratic a s^2 + b s + c with a from W1-minor, c from W2-minor.
        # Brute: scan small rationals + exact root extraction per quadratic.
        from fractions import Fraction as Fr
        cands = set()
        for r1 in range(4):
            for r2 in range(r1 + 1, 4):
                for c1 in range(4):
                    for c2 in range(c1 + 1, 4):
                        def minor(M):
                            G = [M[4 * r:4 * r + 4] for r in range(4)]
                            return (G[r1][c1] * G[r2][c2]
                                    - G[r1][c2] * G[r2][c1])
                        a = minor(W1)
                        c = minor(W2)
                        Mmix = [p + q for p, q in zip(W1, W2)]
                        b = minor(Mmix) - a - c
                        if a == 0 and b == 0:
                            continue
                        if a == 0:
                            cands.add(Fr(-c, b) if b else None)
                        else:
                            disc = b * b - 4 * a * c
                            if disc >= 0:
                                root = _fr_sqrt(disc)
                                if root is not None:
                                    cands.add((-b + root) / (2 * a))
                                    cands.add((-b - root) / (2 * a))
        cands.discard(None)
        for s in cands:
            if all(v == 0 for v in minors(s, Fraction(1))):
                f = factor([s * a + b for a, b in zip(W1, W2)])
                if f:
                    out.append(f)
        if all(v == 0 for v in minors(Fraction(1), Fraction(0))):
            f = factor(W1)
            if f:
                out.append(f)
        return out
    # dim >= 3: report and bail (handle if it ever happens).
    raise NotImplementedError(f"dim W = {len(Wbasis)} — extend the analysis")


def _fr_sqrt(q):
    """Exact sqrt of a nonneg Fraction, or None if irrational."""
    import math
    n, d = q.numerator, q.denominator
    rn, rd = math.isqrt(n), math.isqrt(d)
    if rn * rn == n and rd * rd == d:
        return Fraction(rn, rd)
    return None


def main():
    U0 = [vec((0, 1), (2, 1)), vec((1, 1), (3, 1))]          # a_i + a_b rows
    V0 = [vec((0, 1), (2, 1)), vec((1, 1), (3, 1))]          # x_i + x_b cols
    P = [dyad(u, v) for u in U0 for v in V0]                  # U0 (x) V0, 4-dim
    for mb in (1, 2):
        Sh = shared_rbcb(mb)
        mat, _ = rref(Sh)
        W = intersection(Sh, P)
        print(f"(3,3,bridge-{mb}): dim Sh|RBCB = {len(mat)}, "
              f"dim W = (U0xV0 ^ Sh) = {len(W)}")
        if not W:
            print("  => W = 0: NO 3-product completion exists with ANY local "
                  "rank-1 atoms over this reusable set. IMPOSSIBLE.")
            continue
        r1 = rank1_elements(W)
        print(f"  rank-1 directions in W: {len(r1)}")
        for u, v in r1:
            print(f"    u = {u}  v = {v}")
        if not r1:
            print("  => W != 0 but contains NO rank-1 element: the two "
                  "required dyads cannot lie in it. IMPOSSIBLE.")
        else:
            print("  => stage-2 feasibility chase needed (directions pinned).")


if __name__ == '__main__':
    main()
