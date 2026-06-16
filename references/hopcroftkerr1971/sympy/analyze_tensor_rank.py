"""
Tensor-rank analysis of the (2,2,bridge-3) residual.

The right framing: find the minimum k such that there exist
 - linear combos s1, s2 ∈ span(shared) [free; don't count]
 - k bilinear rank-1 atoms (α_i ⊗ β_i)
 - scalar pairs (c1_i, c2_i)
satisfying
   T1 + s1 = Σ_i c1_i · α_i ⊗ β_i
   T2 + s2 = Σ_i c2_i · α_i ⊗ β_i

This is equivalent to asking the tensor rank of (T1, T2) viewed as
a 6×6×2 tensor MODULO span(shared) ⊗ {e_1, e_2}.

Strategy: for each k = 3, 4, …, attempt to find a solution using a
combination of:
 (a) symbolic constraint solving (sympy.solve / nonlinsolve), and
 (b) heuristic "subtract a shared atom and recompute rank" search.

This script focuses on (b) — sympy nonlinsolve over 40+ variables
gets exponentially slow. The heuristic still provides clean results.
"""
from sympy import symbols, expand, Matrix, solve, Rational
import itertools, time

a_i_1, a_i_2 = symbols('a_i_1 a_i_2')
a_p_1, a_p_2 = symbols('a_p_1 a_p_2')
a_j_1, a_j_2 = symbols('a_j_1 a_j_2')
x_i_1, x_i_2 = symbols('x_i_1 x_i_2')
x_p_1, x_p_2 = symbols('x_p_1 x_p_2')
x_j_1, x_j_2 = symbols('x_j_1 x_j_2')

a_vars = [a_i_1, a_i_2, a_p_1, a_p_2, a_j_1, a_j_2]
x_vars = [x_i_1, x_i_2, x_p_1, x_p_2, x_j_1, x_j_2]

def bilinear_matrix(expr):
    """Coefficient matrix M ∈ Q^{6×6} where M[i,j] = coeff of a_vars[i] · x_vars[j]."""
    expr = expand(expr)
    M = Matrix.zeros(6, 6)
    for i, av in enumerate(a_vars):
        for j, xv in enumerate(x_vars):
            M[i, j] = Rational(expr.coeff(av * xv))
    return M

def A(a1,a2,x1,x2): return a2*(x1+x2)
def B(a1,a2,x1,x2): return (a1-a2)*x1
def C(a1,a2,x1,x2): return (a1-a2)*x2
def D(a1,a2,x1,x2): return a1*(x1+x2)
def E(a1,a2,x1,x2): return a2*x2
def F(a1,a2,x1,x2): return a1*x1

d_a1, d_a2 = a_p_1 - a_i_1, a_p_2 - a_i_2
d_x1, d_x2 = x_p_1 - x_i_1, x_p_2 - x_i_2
shared_exprs = [
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
shared_mats = [(n, bilinear_matrix(e)) for n, e in shared_exprs]

T1_mat = bilinear_matrix((a_i_1 + a_p_1)*x_j_1 + (a_i_2 + a_p_2)*x_j_2)
T2_mat = bilinear_matrix(a_j_1*(x_i_1 + x_p_1) + a_j_2*(x_i_2 + x_p_2))

print("=== Baseline: rank of (T1, T2) before modification ===")
joint_baseline = Matrix.hstack(T1_mat, T2_mat)
print(f"  rank(T1) = {T1_mat.rank()}, rank(T2) = {T2_mat.rank()}, joint = {joint_baseline.rank()}")

# Strategy: try all integer combinations of shared atoms (limited range) to
# add to T1 and T2, and compute the joint matrix rank of the modified pair.
# Minimum joint rank over modifications ≈ tensor rank of T modulo shared.
# (Strictly: joint matrix rank is a LOWER bound on tensor rank, but tight
# in many cases.)

print("\n=== Searching modifications to reduce joint matrix rank ===")
print("    (lower joint rank → fewer new products needed)")
t0 = time.time()
best_joint_rank = joint_baseline.rank()
best_mod = None

# Try modifications: T1 += sum(c1_i * S_i), T2 += sum(c2_i * S_i) with c ∈ {-1, 0, 1}.
# 9 shared × 2 targets = 18 coefficients, 3^18 ≈ 387M — too many.
# Restrict to varying 4 most "central" shared atoms (the ones with mixed-var support).
# G_pair and A_diff, B_diff have cross-monomials; the diagonals are local to single rows.
relevant = [6, 7, 8]  # indices in shared_mats: A_diff, B_diff, G_pair

trials = 0
for c_combo in itertools.product([-1, 0, 1], repeat=len(relevant)*2):
    c1_vec = c_combo[:len(relevant)]
    c2_vec = c_combo[len(relevant):]
    if all(c == 0 for c in c_combo):
        continue
    M1 = T1_mat.copy()
    M2 = T2_mat.copy()
    for k, idx in enumerate(relevant):
        if c1_vec[k] != 0: M1 = M1 + c1_vec[k] * shared_mats[idx][1]
        if c2_vec[k] != 0: M2 = M2 + c2_vec[k] * shared_mats[idx][1]
    joint = Matrix.hstack(M1, M2)
    jr = joint.rank()
    trials += 1
    if jr < best_joint_rank:
        best_joint_rank = jr
        best_mod = (c1_vec, c2_vec)
        print(f"  trial {trials}: joint rank dropped to {jr} via mod {c1_vec}/{c2_vec}")
        if jr <= 3:
            break

print(f"\n  Searched {trials} modifications in {time.time()-t0:.1f}s")
print(f"  Best modified joint rank: {best_joint_rank}")
if best_mod:
    c1_vec, c2_vec = best_mod
    print(f"  Modifications: T1 += {dict((shared_mats[relevant[k]][0], c1_vec[k]) for k in range(len(relevant)) if c1_vec[k])}")
    print(f"                 T2 += {dict((shared_mats[relevant[k]][0], c2_vec[k]) for k in range(len(relevant)) if c2_vec[k])}")

# If joint rank is 4 even after modification, that's a strong indication
# that 4 new products are needed (modulo shared) and impossibility for k=3.
# If joint rank drops to 3, a 3-product solution exists.

if best_joint_rank <= 3:
    print("\n*** A 3-or-fewer-product solution exists modulo shared atoms. ***")
    print("    Tensor rank ≤ joint matrix rank = 3.")
elif best_joint_rank == 4:
    print("\n*** Joint matrix rank cannot drop below 4 across simple modifications. ***")
    print("    Strong evidence for 4-product requirement (HK formula +1 per pair).")
    # But search was limited — extend below if needed.

# Extended search: 5-shared-atom modification (limit cube to ±1 still).
if best_joint_rank > 3:
    print("\n=== Extended search: vary all 9 shared atoms ===")
    print("    (3^18 = 387M too large for full enumeration; sample randomly)")
    import random
    random.seed(42)
    t0 = time.time()
    n_trials = 200000
    for trial in range(n_trials):
        c1_vec = [random.choice([-1, 0, 0, 1]) for _ in range(9)]
        c2_vec = [random.choice([-1, 0, 0, 1]) for _ in range(9)]
        if all(c == 0 for c in c1_vec + c2_vec):
            continue
        M1 = T1_mat + sum(c1_vec[k] * shared_mats[k][1] for k in range(9) if c1_vec[k])
        M2 = T2_mat + sum(c2_vec[k] * shared_mats[k][1] for k in range(9) if c2_vec[k])
        joint = Matrix.hstack(M1, M2)
        jr = joint.rank()
        if jr < best_joint_rank:
            best_joint_rank = jr
            print(f"  trial {trial}: joint rank → {jr}")
            print(f"    T1 mod: {dict((shared_mats[k][0], c1_vec[k]) for k in range(9) if c1_vec[k])}")
            print(f"    T2 mod: {dict((shared_mats[k][0], c2_vec[k]) for k in range(9) if c2_vec[k])}")
            if jr <= 3:
                break
        if time.time() - t0 > 60:
            print(f"  (1-min timeout after {trial} trials)")
            break

print(f"\nFinal best modified joint matrix rank: {best_joint_rank}")
print(f"  → conservative tensor-rank estimate: ≥ {best_joint_rank} new products")
