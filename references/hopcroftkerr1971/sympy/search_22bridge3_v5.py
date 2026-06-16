"""
Math-driven (2,2,bridge-3) analysis.

Instead of enumerating atoms and checking, compute the residual vectors
R1, R2 explicitly and examine their monomial structure. Then try to
factorize them as bilinear rank-1 atoms directly.
"""
from sympy import symbols, expand, Matrix, solve, Symbol, sqrt, Rational, Poly

a_i_1, a_i_2 = symbols('a_i_1 a_i_2')
a_p_1, a_p_2 = symbols('a_p_1 a_p_2')
a_j_1, a_j_2 = symbols('a_j_1 a_j_2')
x_i_1, x_i_2 = symbols('x_i_1 x_i_2')
x_p_1, x_p_2 = symbols('x_p_1 x_p_2')
x_j_1, x_j_2 = symbols('x_j_1 x_j_2')

a_vars = [a_i_1, a_i_2, a_p_1, a_p_2, a_j_1, a_j_2]
x_vars = [x_i_1, x_i_2, x_p_1, x_p_2, x_j_1, x_j_2]
mono = [av * xv for av in a_vars for xv in x_vars]
N = len(mono)

def to_vec(expr):
    expr = expand(expr)
    return [int(expr.coeff(m)) for m in mono]

def vec_to_expr(vec):
    return sum(c * m for c, m in zip(vec, mono) if c != 0)

def A(a1, a2, x1, x2): return a2 * (x1 + x2)
def B(a1, a2, x1, x2): return (a1 - a2) * x1
def C(a1, a2, x1, x2): return (a1 - a2) * x2
def D(a1, a2, x1, x2): return a1 * (x1 + x2)
def E(a1, a2, x1, x2): return a2 * x2
def F(a1, a2, x1, x2): return a1 * x1

d_a1, d_a2 = a_p_1 - a_i_1, a_p_2 - a_i_2
d_x1, d_x2 = x_p_1 - x_i_1, x_p_2 - x_i_2
shared = [
    C(a_i_1, a_i_2, x_i_1, x_i_2),
    D(a_i_1, a_i_2, x_i_1, x_i_2),
    C(a_j_1, a_j_2, x_j_1, x_j_2),
    D(a_j_1, a_j_2, x_j_1, x_j_2),
    E(a_p_1, a_p_2, x_p_1, x_p_2),
    F(a_p_1, a_p_2, x_p_1, x_p_2),
    A(d_a1, d_a2, d_x1, d_x2),
    B(d_a1, d_a2, d_x1, d_x2),
    (d_a1 + a_i_2) * (d_x1 - x_i_2),
]
shared_names = ["C_ii", "D_ii", "C_jj", "D_jj", "E_pp", "F_pp", "A_diff", "B_diff", "G_pair"]
S_mat = Matrix([to_vec(s) for s in shared]).T

T1 = (a_i_1 + a_p_1) * x_j_1 + (a_i_2 + a_p_2) * x_j_2
T2 = a_j_1 * (x_i_1 + x_p_1) + a_j_2 * (x_i_2 + x_p_2)
T1_col = Matrix(to_vec(T1))
T2_col = Matrix(to_vec(T2))

def proj_resid(T_col, S):
    """T - proj(T onto span(S)) via Gram-system solve over Q."""
    StS = S.T * S
    StT = S.T * T_col
    csym = symbols(f'c_:{S.cols}')
    eqs = [sum(StS[i, k] * csym[k] for k in range(S.cols)) - StT[i, 0]
           for i in range(S.cols)]
    sol = solve(eqs, csym, dict=True)
    if not sol:
        return T_col
    proj = Matrix.zeros(N, 1)
    for c in range(S.cols):
        v = sol[0].get(csym[c], 0)
        for r in range(N):
            proj[r, 0] += S[r, c] * v
    return T_col - proj

print("=== Computing residual vectors ===")
R1 = proj_resid(T1_col, S_mat)
R2 = proj_resid(T2_col, S_mat)

R1_expr = expand(vec_to_expr([R1[i, 0] for i in range(N)]))
R2_expr = expand(vec_to_expr([R2[i, 0] for i in range(N)]))

print("\nR1 (residual of T1 from shared span):")
print(f"  = {R1_expr}")
print("\nR2 (residual of T2 from shared span):")
print(f"  = {R2_expr}")

# ── Examine R1 and R2 structure: which variables appear? ──
print("\n=== Variable involvement ===")
print(f"  R1 involves a_vars: {[v for v in a_vars if R1_expr.has(v)]}")
print(f"  R1 involves x_vars: {[v for v in x_vars if R1_expr.has(v)]}")
print(f"  R2 involves a_vars: {[v for v in a_vars if R2_expr.has(v)]}")
print(f"  R2 involves x_vars: {[v for v in x_vars if R2_expr.has(v)]}")

# ── Try to factorize R1 as a sum of bilinear rank-1 atoms ──
print("\n=== Bilinear-rank analysis of R1 ===")
# Represent R1 as a 6×6 matrix M[a-row][x-col] of coefficients.
# Rank of this matrix = bilinear rank of R1.
M1 = Matrix.zeros(6, 6)
for ai, av in enumerate(a_vars):
    for xj, xv in enumerate(x_vars):
        c = R1_expr.coeff(av * xv)
        M1[ai, xj] = c
print(f"  R1 coefficient matrix M1:")
for i in range(6):
    print(f"    {a_vars[i].name}: {[M1[i, j] for j in range(6)]}")
print(f"  rank(M1) = {M1.rank()}")

M2 = Matrix.zeros(6, 6)
for ai, av in enumerate(a_vars):
    for xj, xv in enumerate(x_vars):
        c = R2_expr.coeff(av * xv)
        M2[ai, xj] = c
print(f"\n  R2 coefficient matrix M2:")
for i in range(6):
    print(f"    {a_vars[i].name}: {[M2[i, j] for j in range(6)]}")
print(f"  rank(M2) = {M2.rank()}")

# ── A residual of bilinear rank k can be written as the sum of k
# rank-1 atoms via singular value decomposition or, exactly over Q,
# via LU decomposition.  Each rank-1 atom is α(a)·β(x) where α, β are
# linear forms.
# ── If rank(M1) = 1, then R1 IS a single bilinear product.
# ── If rank(M1) = 2 and rank(M2) = 2 and rank([M1 | M2]) = 2, then
# the SAME pair of (α, β) can express both as linear combinations.
print("\n=== Joint rank ===")
M_joint = Matrix.hstack(M1, M2)
print(f"  rank([M1 | M2]) = {M_joint.rank()}")

# Joint rank-2: 2 row spaces and 2 column spaces span both residuals.
# We can pick the basis directly.

# ── Find a basis of the joint row space ──
print("\n=== Row-space basis for both residuals ===")
# rref the combined matrix; pivot rows give the row basis
rref_joint, pivots = M_joint.rref()
print(f"  rref pivots: {pivots}")

# ── Display M1, M2 in factored form: M_i = sum_k α_k ⊗ β_k where α_k spans row space ──
# Use M1's rank-r decomposition: extract pivots.
print("\n=== Rank-1 decomposition of R1 ===")
U1, P1 = M1.T.rref()
print(f"  M1 column-space pivots: {P1}")

print("\n=== Decomposition strategy ===")
print(f"  If joint row-space (in a-vars) has basis (α_1, α_2)")
print(f"  and joint column-space (in x-vars) has basis (β_1, β_2),")
print(f"  then both R1 and R2 are 2x2 combinations: R_t = sum c^t_{{ij}} α_i β_j.")
print(f"  Need to express R1, R2 as a sum of ≤ 3 bilinear products α_k(a)·β_k(x).")
print(f"  Schmidt / outer-product decomposition over Q gives the minimal k.")

# Try directly: factor M1 as α₁ β₁^T + α₂ β₂^T using rank decomposition.
# sympy's Matrix.LDLdecomposition gives one such basis.
print("\n=== Direct rank-r decomposition of M1 ===")
r1 = M1.rank()
print(f"  rank(M1) = {r1}")
# Get any basis: use rref to find α-basis (rows of pivot columns of M1^T rref)
# and β-basis (columns derived from the reduced form).
M1_rref, M1_pivots = M1.rref()
print(f"  rref(M1) pivot columns (β-basis indices): {M1_pivots}")

# Take β as the rows of M1[pivot_rows, :], then α from the corresponding decomposition.
M1T_rref, M1T_pivots = M1.T.rref()
print(f"  rref(M1^T) pivot columns (α-basis indices): {M1T_pivots}")

if r1 > 0:
    # Each rank-1 atom: α_k (linear combo in a_vars), β_k (linear combo in x_vars)
    # Form α as combinations of basis rows. Use:
    # alpha_k = some linear combo of a_vars that picks out the k-th pivot row of M1
    # beta_k = the corresponding row in M1 normalized for that pivot
    alphas = []
    betas = []
    Mr = M1.copy()  # we'll reduce in place
    for k in range(r1):
        # Find a non-zero row in Mr
        row = None
        for i in range(6):
            if any(Mr[i, j] != 0 for j in range(6)):
                row = i; break
        if row is None:
            break
        # Find a non-zero column in that row
        col = None
        for j in range(6):
            if Mr[row, j] != 0:
                col = j; break
        # Pivot value
        piv = Mr[row, col]
        # Build α_k from the column-vector: α_k = sum_i (Mr[i, col]/piv) * a_vars[i]
        alpha_coefs = [Mr[i, col] / piv for i in range(6)]
        # Build β_k from the row-vector
        beta_coefs = [Mr[row, j] for j in range(6)]
        # Append
        alpha_expr = sum(c * av for c, av in zip(alpha_coefs, a_vars))
        beta_expr = sum(c * xv for c, xv in zip(beta_coefs, x_vars))
        alphas.append(alpha_expr)
        betas.append(beta_expr)
        # Subtract α_k β_k^T from Mr
        for i in range(6):
            for j in range(6):
                Mr[i, j] -= alpha_coefs[i] * beta_coefs[j]
    print(f"\n  R1 = sum of {len(alphas)} rank-1 atoms:")
    for k in range(len(alphas)):
        print(f"    atom {k+1}: ({alphas[k]}) · ({betas[k]})")
    # Verify
    R1_check = sum(alphas[k] * betas[k] for k in range(len(alphas)))
    print(f"\n  Verify R1 = {expand(R1_expr - R1_check) == 0}")

print("\n=== Done ===")
