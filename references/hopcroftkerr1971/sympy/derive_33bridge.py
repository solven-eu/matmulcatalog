"""Derive HK same-method case (3,3,bridge-m) — needed by HK §4 Case 2 Step 3.

Setup: y_{ii} and y_{jj} both method 3 (E+F); bridge row b uses method
m_b ∈ {1, 2}. Mirrors the structure of the WORKING Java emitters
(emitSameMethodPair_*): targets are the bridge SUMS
    S_ij = y_ij + y_{b,j}        S_ji = y_ji + y_{j,b}
(the emitter then subtracts the already-computed pair (b,j) outputs), and the
reusable set is RICH: the diagonals of i, j, b plus ALL products emitted by the
Lemma-2 pair (i,b) — including the third-method virtual products on a_i − a_b.
This richer S is exactly what made the previously-"impossible" (2,2,bridge-3)
case derivable; the historical (3,3,*) enumeration that found 0 solutions used
a poorer S.

Run: python3 derive_33bridge.py     (~1 min; exact Fraction arithmetic)
"""
import itertools
from fractions import Fraction

# ---- symbols as monomial indices -------------------------------------------
# a-vars: (row, col) for row in {i, b, j} = {0, 1, 2}, col in {1, 2}
# x-vars: (row, col) for row in {1, 2}, col in {i, b, j}
AVARS = [(r, c) for r in range(3) for c in (1, 2)]      # 6
XVARS = [(r, c) for r in (1, 2) for c in range(3)]      # 6
AIDX = {v: k for k, v in enumerate(AVARS)}
XIDX = {v: k for k, v in enumerate(XVARS)}
DIM = 36  # monomial (a, x) pairs

def lin_a(coeffs):  # dict (row,col)->coef  → 6-vector
    v = [Fraction(0)] * 6
    for key, c in coeffs.items():
        v[AIDX[key]] += Fraction(c)
    return v

def lin_x(coeffs):
    v = [Fraction(0)] * 6
    for key, c in coeffs.items():
        v[XIDX[key]] += Fraction(c)
    return v

def product(la, lx):  # outer product → 36-vector
    return [la[p] * lx[q] for p in range(6) for q in range(6)]

I, B, J = 0, 1, 2

def a(row, col, s=1): return {(row, col): s}

def madd(*ds):
    out = {}
    for d in ds:
        for k, v in d.items():
            out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v != 0}

def smul(d, s): return {k: v * s for k, v in d.items()}

# x helpers: x(row∈{1,2}, col∈{I,B,J})
def x(row, col, s=1): return {(row, col): s}

# ---- HK base products (per the Java emitters) -------------------------------
def A_m(arow, xcol):   # A(a, x) = a_2 (x_1 + x_2)
    return product(lin_a(a(arow, 2)), lin_x(madd(x(1, xcol), x(2, xcol))))
def B_m(arow, xcol):   # (a_1 − a_2) x_1
    return product(lin_a(madd(a(arow, 1), a(arow, 2, -1))), lin_x(x(1, xcol)))
def C_m(arow, xcol):   # (a_1 − a_2) x_2
    return product(lin_a(madd(a(arow, 1), a(arow, 2, -1))), lin_x(x(2, xcol)))
def D_m(arow, xcol):   # a_1 (x_1 + x_2)
    return product(lin_a(a(arow, 1)), lin_x(madd(x(1, xcol), x(2, xcol))))
def E_m(arow, xcol):   # a_2 x_2
    return product(lin_a(a(arow, 2)), lin_x(x(2, xcol)))
def F_m(arow, xcol):   # a_1 x_1
    return product(lin_a(a(arow, 1)), lin_x(x(1, xcol)))

# generalized: same forms on LINEAR COMBINATIONS of rows/cols
def gen_form(kind, asig, xsig):
    # asig/xsig: dict row->sign (a-side rows; x-side cols), same structural form
    a1 = madd(*[a(r, 1, s) for r, s in asig.items()])
    a2 = madd(*[a(r, 2, s) for r, s in asig.items()])
    x1 = madd(*[x(1, c, s) for c, s in xsig.items()])
    x2 = madd(*[x(2, c, s) for c, s in xsig.items()])
    if kind == 'A': return product(lin_a(a2), lin_x(madd(x1, x2)))
    if kind == 'B': return product(lin_a(madd(a1, smul(a2, -1))), lin_x(x1))
    if kind == 'C': return product(lin_a(madd(a1, smul(a2, -1))), lin_x(x2))
    if kind == 'D': return product(lin_a(a1), lin_x(madd(x1, x2)))
    if kind == 'E': return product(lin_a(a2), lin_x(x2))
    if kind == 'F': return product(lin_a(a1), lin_x(x1))
    raise ValueError(kind)

def G_form(asig1, asig2, xsig3, xsig4):
    # G(arg1, arg2, arg3, arg4) = (arg2_1 + arg1_2) (arg3_1 − arg4_2)
    left = madd(*([a(r, 1, s) for r, s in asig2.items()] + [a(r, 2, s) for r, s in asig1.items()]))
    right_pos = [x(1, c, s) for c, s in xsig3.items()]
    right_neg = [x(2, c, -s) for c, s in xsig4.items()]
    return product(lin_a(left), lin_x(madd(*(right_pos + right_neg))))

# ---- targets ----------------------------------------------------------------
def y(r, c):  # y_{rc} = a_{r,1} x_{1,c} + a_{r,2} x_{2,c}
    return [Fraction(p) + Fraction(q) for p, q in zip(
        product(lin_a(a(r, 1)), lin_x(x(1, c))),
        product(lin_a(a(r, 2)), lin_x(x(2, c))))]

def vadd(u, v): return [p + q for p, q in zip(u, v)]

S_ij = vadd(y(I, J), y(B, J))
S_ji = vadd(y(J, I), y(J, B))

# ---- exact linear algebra ----------------------------------------------------
def in_span(basis, targets):
    """True iff every target is in the rational span of basis (exact)."""
    rows = [list(b) for b in basis]
    # Gaussian elimination to row-echelon, track pivot cols
    pivots = []
    mat = []
    for r in rows:
        r = r[:]
        for (pc, pr) in zip(pivots, mat):
            if r[pc] != 0:
                f = r[pc] / pr[pc]
                r = [ri - f * pi for ri, pi in zip(r, pr)]
        nz = next((idx for idx, v in enumerate(r) if v != 0), None)
        if nz is not None:
            pivots.append(nz)
            mat.append(r)
    for t in targets:
        r = list(t)
        for (pc, pr) in zip(pivots, mat):
            if r[pc] != 0:
                f = r[pc] / pr[pc]
                r = [ri - f * pi for ri, pi in zip(r, pr)]
        if any(v != 0 for v in r):
            return False
    return True

def run(m_bridge):
    # Shared/reusable products.
    shared = {
        'E_ii': E_m(I, I), 'F_ii': F_m(I, I),
        'E_jj': E_m(J, J), 'F_jj': F_m(J, J),
    }
    if m_bridge == 1:
        shared['A_bb'] = A_m(B, B)
        shared['B_bb'] = B_m(B, B)
        # pair (b, i) methods (1, 3) → emits D(a_i−a_b …), C(a_i−a_b …), G(a_b, a_i−a_b, x_b, x_i−x_b)
        d = {I: 1, B: -1}
        shared['D_ib'] = gen_form('D', d, d)
        shared['C_ib'] = gen_form('C', d, d)
        shared['G_ib'] = G_form({B: 1}, d, {B: 1}, d)
    else:
        shared['C_bb'] = C_m(B, B)
        shared['D_bb'] = D_m(B, B)
        # pair (b, i) methods (2, 3) → emits A(a_i−a_b …), B(a_i−a_b …), G(a_i−a_b, a_b, x_i−x_b, x_b)
        d = {I: 1, B: -1}
        shared['A_ib'] = gen_form('A', d, d)
        shared['B_ib'] = gen_form('B', d, d)
        shared['G_ib'] = G_form(d, {B: 1}, d, {B: 1})

    shared_vecs = list(shared.values())

    # Candidate new products: 6 HK forms on signed row-combinations, plus G-forms.
    sigs = [
        {I: 1, B: 1, J: -1}, {I: 1, B: -1, J: 1}, {I: -1, B: 1, J: 1}, {I: 1, B: 1, J: 1},
        {I: 1, J: -1}, {I: 1, J: 1}, {B: 1, J: -1}, {B: 1, J: 1},
    ]
    cands = {}
    for si, sig in enumerate(sigs):
        for kind in 'ABCDEF':
            cands[f'{kind}_s{si}'] = gen_form(kind, sig, sig)
        for gi, (a1, a2, x3, x4) in enumerate([
            ({J: 1}, sig, {J: 1}, sig),
            (sig, {J: 1}, sig, {J: 1}),
            ({J: 1}, sig, sig, {J: 1}),
            (sig, {J: 1}, {J: 1}, sig),
        ]):
            cands[f'G_s{si}_v{gi}'] = G_form(a1, a2, x3, x4)

    names = list(cands.keys())
    print(f"(3,3,bridge-{m_bridge}): shared={len(shared_vecs)}, candidates={len(names)}")

    # Quotient sanity: how many dims do the targets add over shared?
    base_rank_basis = shared_vecs
    need = [S_ij, S_ji]
    if in_span(base_rank_basis, need):
        print("  targets already in shared span (?!)")
        return

    hits = []
    for tri in itertools.combinations(names, 3):
        basis = shared_vecs + [cands[t] for t in tri]
        if in_span(basis, need):
            hits.append(tri)
            print(f"  SOLUTION: {tri}")
            if len(hits) >= 5:
                break
    if not hits:
        print("  no 3-product solution in this candidate catalog")

if __name__ == '__main__':
    run(1)
    run(2)
