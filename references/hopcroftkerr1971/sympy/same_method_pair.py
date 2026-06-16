"""Symbolic derivation of page-10 same-method formulas.

Paper gives the (m1, m1, bridge m2) case explicitly. By symmetry
we want (m2, m2, bridge m1) and (m3, m3, bridge m1 or m2) too.
Verify each derivation symbolically against the matmul tensor.
"""
import sympy as sp

# Three rows i, i+1, j (where j = i+2 or some other "same-method" pair).
# For the (m1, m2, m1) bridge case:
#   y_{ii}     uses method 1 (A+B)
#   y_{i+1,i+1} uses method 2 (-C+D)
#   y_{jj}     uses method 1 (A+B)
# Pair (i+1, j) uses Lemma 2 case (2, 1) (or equivalently (1, 2) swapped).
# Pair (i, j) needs the same-method fallback.

i, ipp, j = 1, 2, 3
A = {(r, c): sp.symbols(f'a{r}{c}', commutative=True) for r in [i, ipp, j] for c in [1, 2]}
X = {(r, c): sp.symbols(f'x{r}{c}', commutative=True) for r in [1, 2] for c in [i, ipp, j]}

def matmul_target(r, c):
    return A[(r, 1)] * X[(1, c)] + A[(r, 2)] * X[(2, c)]

# Page-10 formula: y_{ij} + y_{i+1, j} for the (m1, m1, bridge m2) case.
# y_{ij}+y_{i+1,j} = B(a_j, x_j) + C(-a_j+a_i+a_{i+1}, -x_j+x_i+x_{i+1})
#                  + E(a_j+a_{i+1}, x_j+x_{i+1}) + G(a_j, -a_j+a_{i+1}+a_i, x_j, -x_j+x_i+x_{i+1})
#
# Wait — paper formula reads: G(a_j, -a_j+a_i+a_{i+1}, x_j, -x_j+x_i+x_{i+1})
# So G args: arg1=a_j, arg2=-a_j+a_i+a_{i+1}, arg3=x_j, arg4=-x_j+x_i+x_{i+1}
#   arg1.sub2 = a_{j,2}, arg2.sub1 = -a_{j,1}+a_{i,1}+a_{i+1,1}
#   arg3.sub1 = x_{1,j}, arg4.sub2 = -x_{2,j}+x_{2,i}+x_{2,i+1}
#   G = (arg2.sub1 + arg1.sub2) * (arg3.sub1 - arg4.sub2)

def B_func(ai1, ai2, x1i): return (ai1 - ai2) * x1i
def C_func(aj1, aj2, x2j): return (aj1 - aj2) * x2j
def E_func(ai2, aj2, x2i, x2j): return (ai2 + aj2) * (x2i + x2j)
def G_func(arg1_sub2, arg2_sub1, arg3_sub1, arg4_sub2):
    return (arg2_sub1 + arg1_sub2) * (arg3_sub1 - arg4_sub2)

# Compute claim_sum_ij = B(a_j, x_j) + C(...) + E(a_j+a_{i+1}, x_j+x_{i+1}) + G(...)
B_term = B_func(A[(j, 1)], A[(j, 2)], X[(1, j)])
C_term = C_func(-A[(j, 1)] + A[(i, 1)] + A[(ipp, 1)],
                -A[(j, 2)] + A[(i, 2)] + A[(ipp, 2)],
                -X[(2, j)] + X[(2, i)] + X[(2, ipp)])
E_term = E_func(A[(j, 2)], A[(ipp, 2)], X[(2, j)], X[(2, ipp)])
G_term = G_func(A[(j, 2)],
                -A[(j, 1)] + A[(i, 1)] + A[(ipp, 1)],
                X[(1, j)],
                -X[(2, j)] + X[(2, i)] + X[(2, ipp)])

claim_sum_ij = B_term + C_term + E_term + G_term
target_sum_ij = matmul_target(i, j) + matmul_target(ipp, j)
diff = sp.expand(claim_sum_ij - target_sum_ij)
print(f'(y_ij + y_(i+1)j) = sum_formula?  {"PASS" if diff == 0 else "FAIL " + str(diff)}')

# Page-10 second formula: y_{ji} + y_{j, i+1} = A(a_j, x_j) - D(...) + F(a_j+a_{i+1}, x_j+x_{i+1}) - G(...)
def A_func(ai2, x1l, x2l): return ai2 * (x1l + x2l)
def D_func(aj1, x1l, x2l): return aj1 * (x1l + x2l)
def F_func(ai1, aj1, x1i, x1j): return (ai1 + aj1) * (x1i + x1j)

A_term = A_func(A[(j, 2)], X[(1, j)], X[(2, j)])
D_term = D_func(-A[(j, 1)] + A[(i, 1)] + A[(ipp, 1)],
                -X[(1, j)] + X[(1, i)] + X[(1, ipp)],
                -X[(2, j)] + X[(2, i)] + X[(2, ipp)])
F_term = F_func(A[(j, 1)], A[(ipp, 1)], X[(1, j)], X[(1, ipp)])
# Same G as before
claim_sum_ji = A_term - D_term + F_term - G_term
target_sum_ji = matmul_target(j, i) + matmul_target(j, ipp)
diff = sp.expand(claim_sum_ji - target_sum_ji)
print(f'(y_ji + y_j(i+1)) = sum_formula?  {"PASS" if diff == 0 else "FAIL " + str(diff)}')
