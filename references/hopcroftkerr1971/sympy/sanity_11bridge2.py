"""Sanity: bilinear rank of (1,1,bridge-2) residuals — should be ≤ 3 (page-10 uses 3 new products)."""
from sympy import symbols, expand, Matrix, solve

a_i_1, a_i_2 = symbols('a_i_1 a_i_2')
a_p_1, a_p_2 = symbols('a_p_1 a_p_2')
a_j_1, a_j_2 = symbols('a_j_1 a_j_2')
x_i_1, x_i_2 = symbols('x_i_1 x_i_2')
x_p_1, x_p_2 = symbols('x_p_1 x_p_2')
x_j_1, x_j_2 = symbols('x_j_1 x_j_2')

a_vars = [a_i_1, a_i_2, a_p_1, a_p_2, a_j_1, a_j_2]
x_vars = [x_i_1, x_i_2, x_p_1, x_p_2, x_j_1, x_j_2]
mono = [av * xv for av in a_vars for xv in x_vars]
N = 36

def to_vec(e):
    e = expand(e)
    return [int(e.coeff(m)) for m in mono]

def A(a1,a2,x1,x2): return a2*(x1+x2)
def B(a1,a2,x1,x2): return (a1-a2)*x1
def C(a1,a2,x1,x2): return (a1-a2)*x2
def D(a1,a2,x1,x2): return a1*(x1+x2)
def E(a1,a2,x1,x2): return a2*x2
def F(a1,a2,x1,x2): return a1*x1

# (1,1,bridge-2): y_ii method 1, y_jj method 1, bridge y_{p,p} method 2.
shared = [
    A(a_i_1, a_i_2, x_i_1, x_i_2),  # A_ii
    B(a_i_1, a_i_2, x_i_1, x_i_2),  # B_ii
    A(a_j_1, a_j_2, x_j_1, x_j_2),  # A_jj
    B(a_j_1, a_j_2, x_j_1, x_j_2),  # B_jj
    C(a_p_1, a_p_2, x_p_1, x_p_2),  # C_pp (bridge method 2)
    D(a_p_1, a_p_2, x_p_1, x_p_2),  # D_pp
    # Bridge pair (i, p) Case (1, 2) emits E_combined, F_combined, G
    E(a_i_1 + a_p_1, a_i_2 + a_p_2, x_i_1 + x_p_1, x_i_2 + x_p_2),
    F(a_i_1 + a_p_1, a_i_2 + a_p_2, x_i_1 + x_p_1, x_i_2 + x_p_2),
    (a_p_1 + a_i_2) * (x_i_1 - x_p_2),  # G_pair
]
S_mat = Matrix([to_vec(s) for s in shared]).T

T1 = (a_i_1 + a_p_1)*x_j_1 + (a_i_2 + a_p_2)*x_j_2
T2 = a_j_1*(x_i_1 + x_p_1) + a_j_2*(x_i_2 + x_p_2)
T1_col = Matrix(to_vec(T1))
T2_col = Matrix(to_vec(T2))

def proj_resid(T, S):
    StS = S.T * S
    StT = S.T * T
    c = symbols(f'c_:{S.cols}')
    eqs = [sum(StS[i,k]*c[k] for k in range(S.cols)) - StT[i,0] for i in range(S.cols)]
    sol = solve(eqs, c, dict=True)
    if not sol: return T
    p = Matrix.zeros(36, 1)
    for k in range(S.cols):
        v = sol[0].get(c[k], 0)
        for r in range(36): p[r,0] += S[r,k] * v
    return T - p

R1 = proj_resid(T1_col, S_mat)
R2 = proj_resid(T2_col, S_mat)

def to_mat(R):
    M = Matrix.zeros(6, 6)
    R_expr = expand(sum(R[i, 0] * mono[i] for i in range(36)))
    for ai, av in enumerate(a_vars):
        for xj, xv in enumerate(x_vars):
            M[ai, xj] = R_expr.coeff(av * xv)
    return M

M1 = to_mat(R1)
M2 = to_mat(R2)
joint = Matrix.hstack(M1, M2)

print(f"(1,1,bridge-2) sanity:")
print(f"  rank(M1) = {M1.rank()}   (R1 alone)")
print(f"  rank(M2) = {M2.rank()}   (R2 alone)")
print(f"  rank([M1 | M2]) = {joint.rank()}   ← joint bilinear rank (need this many new products)")
print()
print(f"  page-10 uses 3 new products. So joint rank ≤ 3 must hold.")
