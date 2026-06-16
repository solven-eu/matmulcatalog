"""Same-method bridge cases over the TRUE operational reusable set.

The README's (2,2,bridge-3) impossibility theorem used S = {Ci,Di,Cj,Dj,Ep,Fp}
(6 products). The EMISSION actually has 12 reusables available: the diagonals
of i, j, b (2 each) AND every product of the Lemma-2 pairs (i,b) and (b,j)
(3 each) — individually weightable in W. Targets are the raw cells y_ij, y_ji.

This searches 3-new-product completions for the four problematic cases:
    (2,2,bridge-3), (1,1,bridge-3), (3,3,bridge-1), (3,3,bridge-2)
over that true basis, with shift-form + G-form candidate atoms.

Run: python3 derive_bridge_true_reusables.py    (exact Fraction arithmetic)
"""
import itertools
from fractions import Fraction

AVARS = [(r, c) for r in range(3) for c in (1, 2)]      # rows {i,b,j} = {0,1,2}
XVARS = [(r, c) for r in (1, 2) for c in range(3)]
AIDX = {v: k for k, v in enumerate(AVARS)}
XIDX = {v: k for k, v in enumerate(XVARS)}

def lin_a(coeffs):
    v = [Fraction(0)] * 6
    for key, c in coeffs.items():
        v[AIDX[key]] += Fraction(c)
    return v

def lin_x(coeffs):
    v = [Fraction(0)] * 6
    for key, c in coeffs.items():
        v[XIDX[key]] += Fraction(c)
    return v

def product(la, lx):
    return [la[p] * lx[q] for p in range(6) for q in range(6)]

I, B, J = 0, 1, 2

def a(row, col, s=1): return {(row, col): s}
def x(row, col, s=1): return {(row, col): s}

def madd(*ds):
    out = {}
    for d in ds:
        for k, v in d.items():
            out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v != 0}

def smul(d, s): return {k: v * s for k, v in d.items()}

def gen_form(kind, asig, xsig):
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
    left = madd(*([a(r, 1, s) for r, s in asig2.items()] + [a(r, 2, s) for r, s in asig1.items()]))
    right = madd(*([x(1, c, s) for c, s in xsig3.items()] + [x(2, c, -s) for c, s in xsig4.items()]))
    return product(lin_a(left), lin_x(right))

def y(r, c):
    return [Fraction(p) + Fraction(q) for p, q in zip(
        product(lin_a(a(r, 1)), lin_x(x(1, c))),
        product(lin_a(a(r, 2)), lin_x(x(2, c))))]

def diag_products(row, m):
    s = {row: 1}
    if m == 1: return [gen_form('A', s, s), gen_form('B', s, s)]
    if m == 2: return [gen_form('C', s, s), gen_form('D', s, s)]
    return [gen_form('E', s, s), gen_form('F', s, s)]

def pair_products(r1, m1, r2, m2):
    """The 3 products the Lemma-2 (m1,m2) pair (r1,r2) emits (normalized m1<m2)."""
    if m1 > m2:
        r1, r2, m1, m2 = r2, r1, m2, m1
    if (m1, m2) == (1, 2):
        ssum = {r1: 1, r2: 1}
        return [gen_form('E', ssum, ssum), gen_form('F', ssum, ssum),
                G_form({r1: 1}, {r2: 1}, {r1: 1}, {r2: 1})]
    if (m1, m2) == (1, 3):
        d = {r2: 1, r1: -1}
        return [gen_form('D', d, d), gen_form('C', d, d),
                G_form({r1: 1}, d, {r1: 1}, d)]
    if (m1, m2) == (2, 3):
        d = {r2: 1, r1: -1}
        return [gen_form('A', d, d), gen_form('B', d, d),
                G_form(d, {r1: 1}, d, {r1: 1})]
    raise ValueError((m1, m2))

def in_span(basis, targets):
    pivots, mat = [], []
    for row in basis:
        r = list(row)
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

def run(m_pair, m_bridge):
    shared = []
    shared += diag_products(I, m_pair)
    shared += diag_products(J, m_pair)
    shared += diag_products(B, m_bridge)
    shared += pair_products(I, m_pair, B, m_bridge)
    shared += pair_products(B, m_bridge, J, m_pair)
    targets = [y(I, J), y(J, I)]

    sigs = [
        {I: 1, B: 1, J: -1}, {I: 1, B: -1, J: 1}, {I: -1, B: 1, J: 1}, {I: 1, B: 1, J: 1},
        {I: 1, J: -1}, {I: 1, J: 1}, {B: 1, J: -1}, {B: 1, J: 1}, {I: 1, B: -1}, {I: 1, B: 1},
    ]
    cands = {}
    for si, sig in enumerate(sigs):
        for kind in 'ABCDEF':
            cands[f'{kind}_s{si}'] = gen_form(kind, sig, sig)
        for gi, (a1, a2, x3, x4) in enumerate([
            ({J: 1}, sig, {J: 1}, sig), (sig, {J: 1}, sig, {J: 1}),
            ({J: 1}, sig, sig, {J: 1}), (sig, {J: 1}, {J: 1}, sig),
            ({I: 1}, sig, {I: 1}, sig), (sig, {I: 1}, sig, {I: 1}),
        ]):
            cands[f'G_s{si}_v{gi}'] = G_form(a1, a2, x3, x4)

    names = list(cands.keys())
    print(f"({m_pair},{m_pair},bridge-{m_bridge}): shared={len(shared)}, candidates={len(names)}", flush=True)

    # Pre-reduce: eliminate the shared basis ONCE; reduce candidates + targets
    # modulo it. Per-triple work is then a tiny rank check in the quotient.
    pivots, mat = [], []
    def reduce_mod(vec):
        r = list(vec)
        for (pc, pr) in zip(pivots, mat):
            if r[pc] != 0:
                f = r[pc] / pr[pc]
                r = [ri - f * pi for ri, pi in zip(r, pr)]
        return r
    for row in shared:
        r = reduce_mod(row)
        nz = next((idx for idx, v in enumerate(r) if v != 0), None)
        if nz is not None:
            pivots.append(nz)
            mat.append(r)
    tred = [reduce_mod(t) for t in targets]
    if all(all(v == 0 for v in t) for t in tred):
        print("  targets ALREADY in shared span — 0 new products needed (!)", flush=True)
        return
    cred = {}
    for nm in names:
        r = reduce_mod(cands[nm])
        if any(v != 0 for v in r):
            cred[nm] = r
    rnames = list(cred.keys())
    print(f"  quotient: {len(rnames)} candidates survive reduction", flush=True)

    def covers(tri):
        piv2, mat2 = [], []
        for nm in tri:
            r = list(cred[nm])
            for (pc, pr) in zip(piv2, mat2):
                if r[pc] != 0:
                    f = r[pc] / pr[pc]
                    r = [ri - f * pi for ri, pi in zip(r, pr)]
            nz = next((idx for idx, v in enumerate(r) if v != 0), None)
            if nz is not None:
                piv2.append(nz)
                mat2.append(r)
        for t in tred:
            r = list(t)
            for (pc, pr) in zip(piv2, mat2):
                if r[pc] != 0:
                    f = r[pc] / pr[pc]
                    r = [ri - f * pi for ri, pi in zip(r, pr)]
            if any(v != 0 for v in r):
                return False
        return True

    hits = 0
    for tri in itertools.combinations(rnames, 3):
        if covers(tri):
            print(f"  SOLUTION: {tri}", flush=True)
            hits += 1
            if hits >= 4:
                break
    if hits == 0:
        print("  no 3-product solution in this candidate catalog", flush=True)

if __name__ == '__main__':
    run(2, 3)
    run(1, 3)
    run(3, 1)
    run(3, 2)
