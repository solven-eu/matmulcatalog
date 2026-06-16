"""Generic derivation for all 6 same-method HK page-10 cases.
For each (mij, mBridge), produce formulas for y_{ij}+y_{i+1,j} and
y_{ji}+y_{j,i+1}, identifying 3 NEW products + 4 SHARED.
"""
import sympy as sp
from itertools import product as iprod

i, ipp, j = 1, 2, 3
A = {(r, c): sp.symbols(f'a{r}{c}', commutative=True) for r in [i, ipp, j] for c in [1, 2]}
X = {(r, c): sp.symbols(f'x{r}{c}', commutative=True) for r in [1, 2] for c in [i, ipp, j]}
def target(r, c): return A[(r,1)]*X[(1,c)] + A[(r,2)]*X[(2,c)]
def Ap(ai2, x1l, x2l): return ai2 * (x1l + x2l)
def Bp(ai1, ai2, x1i): return (ai1 - ai2) * x1i
def Cp(aj1, aj2, x2j): return (aj1 - aj2) * x2j
def Dp(aj1, x1l, x2l): return aj1 * (x1l + x2l)
def Ep(ai2, aj2, x2i, x2j): return (ai2 + aj2) * (x2i + x2j)
def Fp(ai1, aj1, x1i, x1j): return (ai1 + aj1) * (x1i + x1j)
def Gp(s2, s1, s1k, s2l): return (s1 + s2) * (s1k - s2l)

sub_a_1 = A[(i,1)] + A[(ipp,1)] - A[(j,1)]
sub_a_2 = A[(i,2)] + A[(ipp,2)] - A[(j,2)]
sub_x_1 = X[(1,i)] + X[(1,ipp)] - X[(1,j)]
sub_x_2 = X[(2,i)] + X[(2,ipp)] - X[(2,j)]

# Single-arg E and F (used when method 3 diagonal: E(a,x)=a_2*x_2, F(a,x)=a_1*x_1)
def E1(ai2, x2l): return ai2 * x2l
def F1(ai1, x1l): return ai1 * x1l

def shared_atoms(mij, mBridge):
    """Products already emitted (as shared) given the methods."""
    a = {}
    # y_ii method mij at (a_i, x_i)
    if mij == 1:
        a['A_ii'] = Ap(A[(i,2)], X[(1,i)], X[(2,i)])
        a['B_ii'] = Bp(A[(i,1)], A[(i,2)], X[(1,i)])
    elif mij == 2:
        a['C_ii'] = Cp(A[(i,1)], A[(i,2)], X[(2,i)])
        a['D_ii'] = Dp(A[(i,1)], X[(1,i)], X[(2,i)])
    elif mij == 3:
        a['E_ii'] = E1(A[(i,2)], X[(2,i)])
        a['F_ii'] = F1(A[(i,1)], X[(1,i)])
    # y_jj method mij at (a_j, x_j)
    if mij == 1:
        a['A_jj'] = Ap(A[(j,2)], X[(1,j)], X[(2,j)])
        a['B_jj'] = Bp(A[(j,1)], A[(j,2)], X[(1,j)])
    elif mij == 2:
        a['C_jj'] = Cp(A[(j,1)], A[(j,2)], X[(2,j)])
        a['D_jj'] = Dp(A[(j,1)], X[(1,j)], X[(2,j)])
    elif mij == 3:
        a['E_jj'] = E1(A[(j,2)], X[(2,j)])
        a['F_jj'] = F1(A[(j,1)], X[(1,j)])
    # y_{i+1,i+1} method mBridge at (a_{i+1}, x_{i+1})
    if mBridge == 1:
        a['A_ipp'] = Ap(A[(ipp,2)], X[(1,ipp)], X[(2,ipp)])
        a['B_ipp'] = Bp(A[(ipp,1)], A[(ipp,2)], X[(1,ipp)])
    elif mBridge == 2:
        a['C_ipp'] = Cp(A[(ipp,1)], A[(ipp,2)], X[(2,ipp)])
        a['D_ipp'] = Dp(A[(ipp,1)], X[(1,ipp)], X[(2,ipp)])
    elif mBridge == 3:
        a['E_ipp'] = E1(A[(ipp,2)], X[(2,ipp)])
        a['F_ipp'] = F1(A[(ipp,1)], X[(1,ipp)])
    # Pair (i, i+1) Lemma 2 — emits E and F applied to (a_i+a_{i+1}, x_i+x_{i+1}) when methods (mij, mBridge) involve {1,2}.
    # For all mij≠mBridge, Lemma 2 emits the "sum" products E_adj and F_adj.
    a['E_adj'] = Ep(A[(i,2)], A[(ipp,2)], X[(2,i)], X[(2,ipp)])
    a['F_adj'] = Fp(A[(i,1)], A[(ipp,1)], X[(1,i)], X[(1,ipp)])
    return a

def candidate_new():
    """Generic catalog of plausible NEW products (subbed-args variants)."""
    return {
        'A_sub': Ap(sub_a_2, sub_x_1, sub_x_2),
        'B_sub': Bp(sub_a_1, sub_a_2, sub_x_1),
        'C_sub': Cp(sub_a_1, sub_a_2, sub_x_2),
        'D_sub': Dp(sub_a_1, sub_x_1, sub_x_2),
        # G's 4-arg variants (only the substituted-args versions are relevant for "new"):
        'G_v1': Gp(A[(j,2)], sub_a_1, X[(1,j)], sub_x_2),     # (1,1,b2)-style
        'G_v2': Gp(sub_a_2, A[(j,1)], sub_x_1, X[(2,j)]),      # (2,2,b1)-style
        'G_v3': Gp(A[(j,2)], sub_a_1, sub_x_1, X[(2,j)]),
        'G_v4': Gp(sub_a_2, A[(j,1)], X[(1,j)], sub_x_2),
    }

def search(mij, mBridge):
    shared = shared_atoms(mij, mBridge)
    new = candidate_new()
    all_atoms = {**shared, **new}
    names = list(all_atoms.keys())
    expanded = {n: sp.expand(v) for n, v in all_atoms.items()}

    targets = {
        'F1': sp.expand(target(i,j) + target(ipp,j)),         # y_ij + y_(i+1)j
        'F2': sp.expand(target(j,i) + target(j,ipp)),         # y_ji + y_j(i+1)
    }

    def find_formula(tgt):
        results = []
        for coefs in iprod([-1, 0, 1], repeat=len(names)):
            nz = sum(1 for c in coefs if c != 0)
            if nz == 0 or nz > 4: continue
            s = sum(c * expanded[n] for c, n in zip(coefs, names) if c != 0)
            if sp.expand(s - tgt) == 0:
                terms = [(n, c) for c, n in zip(coefs, names) if c != 0]
                new_used = sum(1 for n, _ in terms if n in new)
                results.append((new_used, terms))
        results.sort()
        return results

    f1 = find_formula(targets['F1'])
    f2 = find_formula(targets['F2'])
    return f1, f2

for mij, mBridge in [(1,3), (2,3), (3,1), (3,2)]:
    print(f'\n=== ({mij},{mij},bridge{mBridge}) ===')
    f1, f2 = search(mij, mBridge)
    print(f'  F1 (y_ij+y_(i+1)j) solutions: {len(f1)}')
    if f1:
        nu, terms = f1[0]
        print(f'    BEST ({nu} new): ', '  '.join(f'{c:+d}·{n}' for n, c in terms))
    print(f'  F2 (y_ji+y_(j,i+1)) solutions: {len(f2)}')
    if f2:
        nu, terms = f2[0]
        print(f'    BEST ({nu} new): ', '  '.join(f'{c:+d}·{n}' for n, c in terms))
