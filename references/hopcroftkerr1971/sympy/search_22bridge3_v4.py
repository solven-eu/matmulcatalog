"""
Residual-aware search for (2,2,bridge-3).

Key insight from v3: the residual subspace (what new products MUST cover)
has rank only 2. So an atom is useful only if its projection onto the
residual subspace is non-zero. We can prune the candidate catalog
massively by computing each atom's residual-projection up front and
keeping only the contributing ones.

Then for k=2: find pairs whose residual-projections span the 2-dim
residual subspace.

Search expansion vs v3:
 - asymmetric shifts: shift_a-axis-1 ≠ shift_a-axis-2 (and same for x)
 - coefficient range {-2, -1, 0, 1, 2}
 - G-style atoms with mixed row/col forms
"""
from sympy import symbols, expand, Matrix, solve
import itertools, time

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

def A(a1, a2, x1, x2): return a2 * (x1 + x2)
def B(a1, a2, x1, x2): return (a1 - a2) * x1
def C(a1, a2, x1, x2): return (a1 - a2) * x2
def D(a1, a2, x1, x2): return a1 * (x1 + x2)
def E(a1, a2, x1, x2): return a2 * x2
def F(a1, a2, x1, x2): return a1 * x1

# Shared & targets
d_a1, d_a2 = a_p_1 - a_i_1, a_p_2 - a_i_2
d_x1, d_x2 = x_p_1 - x_i_1, x_p_2 - x_i_2
shared = [
    ("C_ii", C(a_i_1, a_i_2, x_i_1, x_i_2)),
    ("D_ii", D(a_i_1, a_i_2, x_i_1, x_i_2)),
    ("C_jj", C(a_j_1, a_j_2, x_j_1, x_j_2)),
    ("D_jj", D(a_j_1, a_j_2, x_j_1, x_j_2)),
    ("E_pp", E(a_p_1, a_p_2, x_p_1, x_p_2)),
    ("F_pp", F(a_p_1, a_p_2, x_p_1, x_p_2)),
    ("A_diff", A(d_a1, d_a2, d_x1, d_x2)),
    ("B_diff", B(d_a1, d_a2, d_x1, d_x2)),
    ("G_pair", (d_a1 + a_i_2) * (d_x1 - x_i_2)),
]
S_mat = Matrix([to_vec(s[1]) for s in shared]).T  # 36 x 9
T1 = (a_i_1 + a_p_1) * x_j_1 + (a_i_2 + a_p_2) * x_j_2
T2 = a_j_1 * (x_i_1 + x_p_1) + a_j_2 * (x_i_2 + x_p_2)
T1_col = Matrix(to_vec(T1))
T2_col = Matrix(to_vec(T2))

# Compute residual subspace basis
def proj_resid(T_col, S):
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

print("Computing residual subspace ...")
t0 = time.time()
R1 = proj_resid(T1_col, S_mat)
R2 = proj_resid(T2_col, S_mat)
R_basis = R1.row_join(R2)
print(f"  residual rank = {R_basis.rank()}  ({time.time()-t0:.1f}s)")
print(f"  R1 norm-1 = {sum(abs(R1[i,0]) for i in range(N))}")
print(f"  R2 norm-1 = {sum(abs(R2[i,0]) for i in range(N))}")

# Build a wider atom catalog WITH asymmetric shifts.
print("\nBuilding candidate atoms (asymmetric shifts, coefs ∈ {-1,0,1}) ...")
t0 = time.time()
shifts = [c for c in itertools.product([-1, 0, 1], repeat=3) if any(c)]
print(f"  {len(shifts)} non-zero shift vectors per axis")

candidates = []
# A..F atoms with INDEPENDENT shifts on axis 1 vs axis 2
# Build sa1 = shift over a_*_1, sa2 = shift over a_*_2 (independent), sx1, sx2 similar
for ca1 in shifts:
    sa1 = ca1[0]*a_i_1 + ca1[1]*a_p_1 + ca1[2]*a_j_1
    for ca2 in shifts:
        sa2 = ca2[0]*a_i_2 + ca2[1]*a_p_2 + ca2[2]*a_j_2
        for cx1 in shifts:
            sx1 = cx1[0]*x_i_1 + cx1[1]*x_p_1 + cx1[2]*x_j_1
            for cx2 in shifts:
                sx2 = cx2[0]*x_i_2 + cx2[1]*x_p_2 + cx2[2]*x_j_2
                # 6 methods
                for mname, fn in [("A", A), ("B", B), ("C", C), ("D", D), ("E", E), ("F", F)]:
                    atom = fn(sa1, sa2, sx1, sx2)
                    v = to_vec(atom)
                    if any(v):
                        candidates.append((f"{mname}(a:{ca1}/{ca2},x:{cx1}/{cx2})", v))
print(f"  {len(candidates)} method atoms  ({time.time()-t0:.1f}s)")

# Project each onto residual; keep only those with non-trivial projection
print("\nResidual-projection pruning ...")
t0 = time.time()
# Residual subspace basis: column space of R_basis (rank 2)
# Project atom v onto R_basis^⊥: subtract projection onto span(R_basis^⊥)
# Simpler: check if v has a non-zero component in span(R_basis).
# Atom v contributes iff <S; v> rank > rank(S), i.e., v not in span(S).
# But we want v to also have a component in span(R_basis).
# Cleanest: rank of [S | R_basis] vs rank of [S | R_basis | v].
SR = S_mat.row_join(R_basis)
SR_rank = SR.rank()
contributing = []
for lbl, v in candidates:
    vc = Matrix(v)
    full_rank = SR.row_join(vc).rank()
    # An atom contributes iff it's outside span(S) AND brings rank progress.
    # We want atoms in span(SR) ∖ span(S) ideally — those project to the residual.
    if S_mat.row_join(vc).rank() > S_mat.rank():
        # Check projection onto residual subspace
        # Atom v ∈ span(SR) iff rank([SR|v]) == SR_rank
        if SR.row_join(vc).rank() == SR_rank:
            contributing.append((lbl, v))
print(f"  {len(contributing)} atoms have non-trivial residual-projection  ({time.time()-t0:.1f}s)")

# Dedup: keep one representative per residual-projection equivalence class
print("\nDeduping residual-projection equivalence classes ...")
t0 = time.time()
# Compute each atom's projection onto residual subspace.
# Project v onto span(R_basis) by solving (R_basis^T R_basis) c = R_basis^T v
def project_onto_R(v):
    R = R_basis
    RtR = R.T * R
    Rtv = R.T * v
    csym = symbols(f'p_:{R.cols}')
    eqs = [sum(RtR[i, k] * csym[k] for k in range(R.cols)) - Rtv[i, 0]
           for i in range(R.cols)]
    sol = solve(eqs, csym, dict=True)
    if not sol:
        return Matrix.zeros(N, 1)
    p = Matrix.zeros(N, 1)
    for c in range(R.cols):
        val = sol[0].get(csym[c], 0)
        for r in range(N):
            p[r, 0] += R[r, c] * val
    return p

# This is O(N * eqs) per atom. Limit to fewer atoms.
sampled = contributing[:200]
projs = [(lbl, v, project_onto_R(Matrix(v))) for lbl, v in sampled]
print(f"  computed {len(projs)} residual-projections  ({time.time()-t0:.1f}s)")

# Now search pairs (k=2)
print("\nSearching for k=2 solutions (pair of atoms spanning residual)...")
t0 = time.time()
hits = 0
for i in range(len(projs)):
    if time.time() - t0 > 120:
        print(f"  timeout after 2 min, i={i}")
        break
    lbl1, v1, p1 = projs[i]
    if all(p1[r,0] == 0 for r in range(N)):
        continue
    for j in range(i+1, len(projs)):
        lbl2, v2, p2 = projs[j]
        # Pair (p1, p2) must span the residual subspace (= R_basis)
        pair_mat = Matrix([list(p1.T)[0], list(p2.T)[0]]).T
        # rank(pair_mat) must equal 2 (= residual rank), AND R_basis ⊂ span(pair_mat)
        if pair_mat.rank() != 2:
            continue
        # Check R_basis cols are in span(pair_mat)
        if (pair_mat.row_join(R_basis)).rank() == 2:
            hits += 1
            print(f"  k=2 HIT #{hits}: {lbl1} ∧ {lbl2}")
            # Compute T1, T2 coefficients
            M = S_mat.row_join(Matrix(v1)).row_join(Matrix(v2))
            all_lbls = [s[0] for s in shared] + [lbl1, lbl2]
            for tname, tcol in [("T1", T1_col), ("T2", T2_col)]:
                qs = symbols(f'q_:{M.cols}')
                eqs = [sum(M[r, c] * qs[c] for c in range(M.cols)) - tcol[r, 0]
                       for r in range(N)]
                sol = solve(eqs, qs, dict=True)
                if sol:
                    print(f"    {tname} =")
                    for ci, lbl in enumerate(all_lbls):
                        v = sol[0].get(qs[ci], 0)
                        if v != 0:
                            print(f"      {'+' if v > 0 else '-'} {abs(v) if abs(v) != 1 else ''}{lbl}")
            if hits >= 3:
                break
    if hits >= 3:
        break

if hits == 0:
    print("\n  no k=2 solutions found in this 200-atom sample.")
    print("  → either need wider atom catalog or k=3+ products.")
print(f"\nDone. {hits} hits, {time.time()-t0:.1f}s total search.")
