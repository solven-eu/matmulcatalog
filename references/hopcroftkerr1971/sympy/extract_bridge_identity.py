"""Extract the explicit W-weights for the (m,m,bridge-b) solutions found by
derive_bridge_true_reusables.py — produces the exact identity to implement in
the Java emitter."""
from fractions import Fraction
import importlib
m = importlib.import_module('derive_bridge_true_reusables')

def solve_weights(m_pair, m_bridge, tri_sigs):
    shared = []
    labels = []
    for row, lbl in [(m.I, 'i'), (m.J, 'j'), (m.B, 'b')]:
        mm = m_pair if lbl in ('i', 'j') else m_bridge
        prods = m.diag_products(row, mm)
        kinds = {1: 'AB', 2: 'CD', 3: 'EF'}[mm]
        for q, P in enumerate(prods):
            shared.append(P); labels.append(f"{kinds[q]}({lbl})")
    for (r1, m1, r2, m2, tag) in [(m.I, m_pair, m.B, m_bridge, 'ib'), (m.B, m_bridge, m.J, m_pair, 'bj')]:
        prods = m.pair_products(r1, m1, r2, m2)
        for q, P in enumerate(prods):
            shared.append(P); labels.append(f"pair{tag}#{q}")
    cands = {}
    sigs = [
        {m.I: 1, m.B: 1, m.J: -1}, {m.I: 1, m.B: -1, m.J: 1}, {m.I: -1, m.B: 1, m.J: 1},
        {m.I: 1, m.B: 1, m.J: 1},
        {m.I: 1, m.J: -1}, {m.I: 1, m.J: 1}, {m.B: 1, m.J: -1}, {m.B: 1, m.J: 1},
        {m.I: 1, m.B: -1}, {m.I: 1, m.B: 1},
    ]
    for si, sig in enumerate(sigs):
        for kind in 'ABCDEF':
            cands[f'{kind}_s{si}'] = m.gen_form(kind, sig, sig)
        for gi, (a1, a2, x3, x4) in enumerate([
            ({m.J: 1}, sig, {m.J: 1}, sig), (sig, {m.J: 1}, sig, {m.J: 1}),
            ({m.J: 1}, sig, sig, {m.J: 1}), (sig, {m.J: 1}, {m.J: 1}, sig),
            ({m.I: 1}, sig, {m.I: 1}, sig), (sig, {m.I: 1}, sig, {m.I: 1}),
        ]):
            cands[f'G_s{si}_v{gi}'] = m.G_form(a1, a2, x3, x4)
    basis = shared + [cands[t] for t in tri_sigs]
    names = labels + list(tri_sigs)
    # Solve basis^T w = target for each target via Gaussian elimination with
    # augmented tracking (least-structure: use sympy-free exact solve).
    for tname, target in [('y_ij', m.y(m.I, m.J)), ('y_ji', m.y(m.J, m.I))]:
        # Build matrix: rows = 36 monomials, cols = basis vectors; solve M w = t.
        rows = 36; cols = len(basis)
        M_ = [[basis[c][r] for c in range(cols)] for r in range(rows)]
        t = [target[r] for r in range(rows)]
        # Gaussian elimination on augmented [M | t]
        aug = [M_[r] + [t[r]] for r in range(rows)]
        piv_cols = []
        rr = 0
        for c in range(cols):
            pr = None
            for r2 in range(rr, rows):
                if aug[r2][c] != 0: pr = r2; break
            if pr is None: continue
            aug[rr], aug[pr] = aug[pr], aug[rr]
            pv = aug[rr][c]
            aug[rr] = [v / pv for v in aug[rr]]
            for r2 in range(rows):
                if r2 != rr and aug[r2][c] != 0:
                    f = aug[r2][c]
                    aug[r2] = [v2 - f * v1 for v1, v2 in zip(aug[rr], aug[r2])]
            piv_cols.append(c)
            rr += 1
        # Check consistency
        for r2 in range(rr, rows):
            if aug[r2][cols] != 0:
                print(f"  {tname}: INCONSISTENT?!"); return
        w = [Fraction(0)] * cols
        for irow, c in enumerate(piv_cols):
            w[c] = aug[irow][cols]
        terms = [f"{'+' if w[c] > 0 else ''}{w[c]}·{names[c]}" for c in range(cols) if w[c] != 0]
        print(f"  {tname} = " + " ".join(terms))

if __name__ == '__main__':
    print("(2,2,bridge-3) solution 1: E_s0, F_s0, G_s6_v5  [sub = a_i+a_b-a_j; G sig (B,-J) v5=(sig,I,sig,I)]")
    solve_weights(2, 3, ('E_s0', 'F_s0', 'G_s6_v5'))
    print("(1,1,bridge-3) solution 1: E_s0, F_s0, G_s6_v4  [v4=(I,sig,I,sig)]")
    solve_weights(1, 3, ('E_s0', 'F_s0', 'G_s6_v4'))
