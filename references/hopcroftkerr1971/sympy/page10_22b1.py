import sympy as sp
i, ipp, j = 1, 2, 3
A = {(r, c): sp.symbols(f'a{r}{c}', commutative=True) for r in [i, ipp, j] for c in [1, 2]}
X = {(r, c): sp.symbols(f'x{r}{c}', commutative=True) for r in [1, 2] for c in [i, ipp, j]}
def target(r, c): return A[(r,1)]*X[(1,c)] + A[(r,2)]*X[(2,c)]
def A_(ai2, x1l, x2l): return ai2 * (x1l + x2l)
def B_(ai1, ai2, x1i): return (ai1 - ai2) * x1i
def C_(aj1, aj2, x2j): return (aj1 - aj2) * x2j
def D_(aj1, x1l, x2l): return aj1 * (x1l + x2l)
def E_(ai2, aj2, x2i, x2j): return (ai2 + aj2) * (x2i + x2j)
def F_(ai1, aj1, x1i, x1j): return (ai1 + aj1) * (x1i + x1j)
def G_(s2, s1, s1k, s2l): return (s1 + s2) * (s1k - s2l)

# (2,2,bridge1): y_ii method 2 (-C+D), y_jj method 2, y_{i+1,i+1} method 1 (A+B).
# By symmetry with (1,1,bridge2), guess: swap A↔C, B↔D throughout.
# Original formula 1: B(a_j,x_j) + C(sub) + E(a_i+a_{i+1}, x_i+x_{i+1}) + G(...)
# Guess: D(a_j,x_j) + A(sub) + ??? + G(?)
# Hmm not clean - A and D have different structures. Let's enumerate.

# Available products in (2,2,bridge1):
# - C(a_i,x_i), D(a_i,x_i), C(a_j,x_j), D(a_j,x_j), A(a_{i+1},x_{i+1}), B(a_{i+1},x_{i+1})
# - Adjacent pair (i, i+1) methods (2,1) Lemma 2: E(a_i+a_{i+1}, x_i+x_{i+1}),
#   F(a_i+a_{i+1}, x_i+x_{i+1}), G(...)
# - 3 NEW: ?, ?, ?

# Try by analogy: swap (a, x) interpretation in C, D vs A, B definitions:
#   A(a, x) = a_2 (x_1 + x_2)              vs    D(a, x) = a_1 (x_1 + x_2)
#   B(a, x) = (a_1 - a_2) x_1              vs    C(a, x) = (a_1 - a_2) x_2
# So {A, C} differ in (a_2/a_1) and (x_sum/x_2). {B, D} similar.
# The (1,1,bridge2) NEW products were C(sub), D(sub), G(sub).
# By the symmetry: (2,2,bridge1) NEW should be A(sub), B(sub), G(sub_diff_args)?

# Try formula candidates and let sympy verify.
sub_a_1 = A[(i,1)] + A[(ipp,1)] - A[(j,1)]
sub_a_2 = A[(i,2)] + A[(ipp,2)] - A[(j,2)]
sub_x_1 = X[(1,i)] + X[(1,ipp)] - X[(1,j)]
sub_x_2 = X[(2,i)] + X[(2,ipp)] - X[(2,j)]

# Guess by swap (B,C)↔(D,A), keep G, keep E↔F:
candidates = [
    ('+D(a_j) + A(sub) + F + G', D_(A[(j,1)], X[(1,j)], X[(2,j)]) + A_(sub_a_2, sub_x_1, sub_x_2)
      + F_(A[(i,1)], A[(ipp,1)], X[(1,i)], X[(1,ipp)]) + G_(A[(j,2)], sub_a_1, X[(1,j)], sub_x_2)),
    ('-C(a_j) - B(sub) + E + G', -C_(A[(j,1)], A[(j,2)], X[(2,j)]) - B_(sub_a_1, sub_a_2, sub_x_1)
      + E_(A[(i,2)], A[(ipp,2)], X[(2,i)], X[(2,ipp)]) + G_(A[(j,2)], sub_a_1, X[(1,j)], sub_x_2)),
    ('+D(a_j) - B(sub) + F - G', D_(A[(j,1)], X[(1,j)], X[(2,j)]) - B_(sub_a_1, sub_a_2, sub_x_1)
      + F_(A[(i,1)], A[(ipp,1)], X[(1,i)], X[(1,ipp)]) - G_(A[(j,2)], sub_a_1, X[(1,j)], sub_x_2)),
    ('+D(a_j) + A(sub) - E + G', D_(A[(j,1)], X[(1,j)], X[(2,j)]) + A_(sub_a_2, sub_x_1, sub_x_2)
      - E_(A[(i,2)], A[(ipp,2)], X[(2,i)], X[(2,ipp)]) + G_(A[(j,2)], sub_a_1, X[(1,j)], sub_x_2)),
]
tgt = target(i, j) + target(ipp, j)
for label, formula in candidates:
    diff = sp.expand(formula - tgt)
    print(f'  {label}: {"PASS" if diff == 0 else "FAIL"}')
