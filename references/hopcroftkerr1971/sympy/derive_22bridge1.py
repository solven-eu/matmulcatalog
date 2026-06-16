"""Derive HK same-method case (2,2,bridge1) by symbolic enumeration.

Setup: y_{ii} method 2 (-C+D), y_{jj} method 2, y_{i+1,i+1} method 1 (A+B).
By analogy with (1,1,bridge2), the 3 new products should involve substituted
arguments. Enumerate candidate combinations and verify symbolically.
"""
import sympy as sp
import itertools

i, ipp, j = 1, 2, 3
A = {(r, c): sp.symbols(f'a{r}{c}', commutative=True) for r in [i, ipp, j] for c in [1, 2]}
X = {(r, c): sp.symbols(f'x{r}{c}', commutative=True) for r in [1, 2] for c in [i, ipp, j]}

def target(r, c): return A[(r,1)]*X[(1,c)] + A[(r,2)]*X[(2,c)]
def A_p(ai2, x1l, x2l): return ai2 * (x1l + x2l)
def B_p(ai1, ai2, x1i): return (ai1 - ai2) * x1i
def C_p(aj1, aj2, x2j): return (aj1 - aj2) * x2j
def D_p(aj1, x1l, x2l): return aj1 * (x1l + x2l)
def E_p(ai2, aj2, x2i, x2j): return (ai2 + aj2) * (x2i + x2j)
def F_p(ai1, aj1, x1i, x1j): return (ai1 + aj1) * (x1i + x1j)
def G_p(s2, s1, s1k, s2l): return (s1 + s2) * (s1k - s2l)

# Substituted sums for (2,2,bridge1) — analogous to (1,1,bridge2) but
# possibly different sub-expression structure:
sub_a_1 = A[(i,1)] + A[(ipp,1)] - A[(j,1)]   # -a_j + a_i + a_{i+1}
sub_a_2 = A[(i,2)] + A[(ipp,2)] - A[(j,2)]
sub_x_1 = X[(1,i)] + X[(1,ipp)] - X[(1,j)]
sub_x_2 = X[(2,i)] + X[(2,ipp)] - X[(2,j)]

# Build catalog of available product instances.
# For the (2,2,bridge1) case, available SHARED products:
# - From y_{ii} method 2:  C(a_i,x_i), D(a_i,x_i)
# - From y_{jj} method 2:  C(a_j,x_j), D(a_j,x_j)
# - From y_{i+1,i+1} m1:   A(a_{i+1},x_{i+1}), B(a_{i+1},x_{i+1})
# - From pair (i, i+1) Lemma 2 case (2,1): E(a_i+a_{i+1}, x_i+x_{i+1}),
#   F(a_i+a_{i+1}, x_i+x_{i+1}), G(a_{i+1}, a_i, x_{i+1}, x_i)

shared = {
    'C_ii': C_p(A[(i,1)], A[(i,2)], X[(2,i)]),
    'D_ii': D_p(A[(i,1)], X[(1,i)], X[(2,i)]),
    'C_jj': C_p(A[(j,1)], A[(j,2)], X[(2,j)]),
    'D_jj': D_p(A[(j,1)], X[(1,j)], X[(2,j)]),
    'A_ipp': A_p(A[(ipp,2)], X[(1,ipp)], X[(2,ipp)]),
    'B_ipp': B_p(A[(ipp,1)], A[(ipp,2)], X[(1,ipp)]),
    'E_iIpp': E_p(A[(i,2)], A[(ipp,2)], X[(2,i)], X[(2,ipp)]),
    'F_iIpp': F_p(A[(i,1)], A[(ipp,1)], X[(1,i)], X[(1,ipp)]),
}

# Candidate NEW products (subbed-args versions of A/B/C/D/E/F/G):
new_prod_candidates = {
    'A_sub': A_p(sub_a_2, sub_x_1, sub_x_2),
    'B_sub': B_p(sub_a_1, sub_a_2, sub_x_1),
    'C_sub': C_p(sub_a_1, sub_a_2, sub_x_2),
    'D_sub': D_p(sub_a_1, sub_x_1, sub_x_2),
    'G_v1':  G_p(A[(j,2)], sub_a_1, X[(1,j)], sub_x_2),
    'G_v2':  G_p(sub_a_2, A[(j,1)], sub_x_1, X[(2,j)]),
    'G_v3':  G_p(A[(j,2)], sub_a_1, sub_x_1, X[(2,j)]),
    'G_v4':  G_p(sub_a_2, A[(j,1)], X[(1,j)], sub_x_2),
}

tgt = target(i, j) + target(ipp, j)
all_atoms = {**shared, **new_prod_candidates}
names = list(all_atoms.keys())

# Search small linear combinations: coefficient ∈ {-1, 0, 1} for each atom.
# Limit to ≤ 5 nonzero atoms (since 4 shared + 3 new = 7 atoms in (1,1,bridge2)).
from itertools import product as iprod
import random
random.seed(0xC0FFEE)

print(f'Total atoms: {len(names)}')
found = []
expanded_atoms = {n: sp.expand(v) for n, v in all_atoms.items()}
tgt_expand = sp.expand(tgt)

# Smart search: subtract target as a fixed polynomial; find subset of atoms
# whose signed sum equals tgt.
# Try {-1, 0, +1} coefficients exhaustively over 8 atoms = 3^8 = 6561 combos.
combos_checked = 0
for coefs in iprod([-1, 0, 1], repeat=len(names)):
    combos_checked += 1
    nonzero = sum(1 for c in coefs if c != 0)
    if nonzero == 0 or nonzero > 6: continue
    s = sum(c * expanded_atoms[n] for c, n in zip(coefs, names) if c != 0)
    if sp.expand(s - tgt_expand) == 0:
        terms = [(n, c) for c, n in zip(coefs, names) if c != 0]
        found.append(terms)
        if len(found) >= 10: break
print(f'Combos checked: {combos_checked}')
print(f'Candidate formulas found: {len(found)}')
for terms in found[:5]:
    print(' ', '  '.join(f'{c:+d}·{n}' for n, c in terms))
