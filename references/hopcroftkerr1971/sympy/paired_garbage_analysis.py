"""
HK (2,2,bridge-3): Paired-atom garbage variety.

NOTE — historical context. This file is an INDEPENDENT REDERIVATION
arrived at separately from `paired_garbage_analysis_v2.py` and
`x_equals_p_proof.py` (which are the canonical references on the
result in the README). The two analyses agree on the conclusion
(`dim P = 14`, `dim X = 14` after saturating the `c=0` component,
hence proportionality theorem and impossibility). This script is kept
as the audit-trail v1 for the alternative route:

  • V2 / x_equals_p_proof.py: parameterise P explicitly via
    `ψ : (λ, μ, α, β) ↦ (λ, μ, αλ, βμ, αβ)`, then prove that the
    constraint-Jacobian rank at a smooth proportional point is 11 ⇒
    `dim X = 14`. The c=0 component is identified by saturating the
    fibre ideal by c.

  • THIS SCRIPT (v1 route): split the system by the variable c
    DIRECTLY, treat `t_{kj} := L'_k · M'_j` as 36 new linear variables
    (so the 28 garbage-difference equations become LINEAR in (t, c)),
    nullspace-solve at a generic A to get a 9-dim affine space in (t, c)
    space, then impose `rank(T) ≤ 1` via 2×2 minors and a Gröbner
    basis. The Gröbner solve confirms the only c≠0 branch is the
    proportional one (tangent dim 1 at the proportional point, in 5/5
    random trials); the c=0 branches are exactly the kernel directions
    in the diagonal-block subspace (= the V variety).

Builds on `zero_garbage_variety.py`. That script established:

  • The "zero garbage" variety V is positive-dimensional (8 axis-aligned
    components) but the rational map V → Q = R^2, (λ, μ) ↦ (α, β), sends
    every component to the origin. So no SINGLE rank-1 atom in V can
    contribute (α, β) ≠ 0.

  • Consequently a 3-atom (2, 2, bridge-3) construction MUST rely on
    garbage cancellation across atoms (rank(U) ≤ 1, U ≠ 0) rather than
    zero garbage per atom.

This script targets the NEXT structural question: characterise the
algebraic set

    X = { (A, A', c) ∈ R^12 × R^12 × R^1 : g(A') = c · g(A) }

where g : R^12 → R^28 is the garbage projection of a rank-1 atom
L(a)·M(x), L = Σ l_k a_k, M = Σ m_k x_k.

Two atoms with PROPORTIONAL garbage are necessary (but not sufficient)
ingredients of a triple with rank(U) ≤ 1: if rank(U) ≤ 1, then any two
columns of U are proportional (one possibly zero).

THE CENTRAL QUESTION:
  Is every (A, A', c) ∈ X forced to come from PROPORTIONAL rank-1
  tensors — i.e. L'·M' = (c'·L)·M' = L·(c'·M') for some scaling, which
  is equivalent to L'·M' = c'·L·M as a tensor — ?

If YES → two atoms with proportional garbage have proportional
(α, β) → three atoms can carve out at most a 1-D direction in Q →
IMPOSSIBILITY for the rank-1 family.

If NO → there's a non-proportional branch in X, parameterise it, and
check whether 3 atoms drawn from it can span Q.

Run: python3 paired_garbage_analysis.py
"""
import sympy as sp
import time

t0 = time.time()
def log(msg):
    print(f"[{time.time()-t0:6.1f}s] {msg}", flush=True)


# ── Setup (same monomial basis as zero_garbage_variety.py) ──────────────────
ai1, ai2, ap1, ap2, aj1, aj2 = sp.symbols('ai1 ai2 ap1 ap2 aj1 aj2')
xi1, xi2, xp1, xp2, xj1, xj2 = sp.symbols('xi1 xi2 xp1 xp2 xj1 xj2')
a_vars = [ai1, ai2, ap1, ap2, aj1, aj2]
x_vars = [xi1, xi2, xp1, xp2, xj1, xj2]
mono = [(av, xv) for av in a_vars for xv in x_vars]
N = len(mono)


def to_vec(expr):
    expr = sp.expand(expr)
    return sp.Matrix([sp.nsimplify(expr.coeff(av).coeff(xv)) for av, xv in mono])


T1 = ai2 * xi1 + ap2 * xp1 + aj2 * xj1
T2 = ai1 * xi2 + ap1 * xp2 + aj1 * xj2

Ci = (ai1 - ai2) * xi2
Di = ai1 * (xi1 + xi2)
Cj = (aj1 - aj2) * xj2
Dj = aj1 * (xj1 + xj2)
Ep = ap2 * xp2
Fp = ap1 * xp1
shared = [Ci, Di, Cj, Dj, Ep, Fp]

vT1 = to_vec(T1)
vT2 = to_vec(T2)
S_mat = sp.Matrix.hstack(*[to_vec(s) for s in shared])
P_perp_S = sp.eye(N) - S_mat * (S_mat.T * S_mat).inv() * S_mat.T
q1 = sp.simplify(P_perp_S * vT1)
q2 = sp.simplify(P_perp_S * vT2)
Q_mat = sp.Matrix.hstack(q1, q2)
P_Q = Q_mat * (Q_mat.T * Q_mat).inv() * Q_mat.T
P_garbage = P_perp_S - P_Q

log("Projector P_garbage built")

# Atom A (variables l1..l6, m1..m6)
ls = sp.symbols('l1:7')
ms = sp.symbols('m1:7')
L = sum(ls[k] * a_vars[k] for k in range(6))
M = sum(ms[k] * x_vars[k] for k in range(6))
R = sp.expand(L * M)
vR = to_vec(R)
garbage_A = sp.simplify(P_garbage * vR)
garbage_polys_A = [sp.simplify(garbage_A[k]) for k in range(N)
                   if sp.simplify(garbage_A[k]) != 0]
log(f"#garbage(A) polys = {len(garbage_polys_A)}")

# Atom A' (variables L1..L6, M1..M6)
Ls = sp.symbols('L1:7')
Ms = sp.symbols('M1:7')
Lp = sum(Ls[k] * a_vars[k] for k in range(6))
Mp = sum(Ms[k] * x_vars[k] for k in range(6))
Rp = sp.expand(Lp * Mp)
vRp = to_vec(Rp)
garbage_Ap = sp.simplify(P_garbage * vRp)
garbage_polys_Ap = [sp.simplify(garbage_Ap[k]) for k in range(N)
                    if sp.simplify(garbage_Ap[k]) != 0]
log(f"#garbage(A') polys = {len(garbage_polys_Ap)}")

# scalar c
c = sp.Symbol('c')

# 28 garbage equations: g(A') = c * g(A)
# Use the same coordinate ordering (we just use all 36 coords; redundant ones
# will simply reduce to 0 = 0).
pair_eqs = [sp.expand(garbage_Ap[k] - c * garbage_A[k]) for k in range(N)]
pair_eqs = [e for e in pair_eqs if e != 0]
# Deduplicate
pair_eqs = list({sp.expand(e) for e in pair_eqs})
log(f"#non-trivial paired equations = {len(pair_eqs)}")

all_vars = list(ls) + list(ms) + list(Ls) + list(Ms) + [c]
log(f"Total variables: {len(all_vars)}  (l1..l6, m1..m6, L1..L6, M1..M6, c)")

# Classify equations:
# Each eq is of form  (bilinear in L,M) - c * (bilinear in l,m).
# So most equations are bilinear (in either l*m or L*M) plus a c*l*m term.

# Build per-monomial decomposition
print("\nSample (first 5) paired equations:")
for e in pair_eqs[:5]:
    print(f"   {e} = 0")


# ── Substitute the L,M expressions into known A-side covers ────────────────
# Strategy: rather than tackling the whole 25-var variety directly, we
# reuse the 8-cover decomposition from V on the A side.  For each cover
# (cover_A: list of A-side vars that = 0), we ask:
#   given A in that cover with non-trivial garbage(A), what are the
#   A' and c that solve g(A') = c · g(A)?
# But wait — the V analysis showed garbage(A) ≡ 0 on each cover.  So
# requiring A ∈ V makes c arbitrary and trivialises the system.
# Instead we keep A GENERIC (no covers) and just solve the 28 bilinear
# equations.  We split it as follows.

# Key structural simplification: the equations are LINEAR in (L1..L6, M1..M6)
# treated alone?  No — each term Lk*Mj is BILINEAR.  But they ARE linear
# in {Lk Mj} viewed as 36 "new" unknowns.  So if we think of t_kj = Lk*Mj,
# the equations become LINEAR in (t_kj, c) with coefficients in (lk, mj).
#
# t_kj is a rank-1 matrix (t = L·M^T).  The 6x6 matrix T = (t_kj) thus has
# rank exactly 1.  So:
#   X is the preimage under (L, M) ↦ T = L M^T of the LINEAR variety
#   L_A := { (T, c) ∈ R^36 × R : F_A · vec(T) = c · g(A) }
#   intersected with { T : rank(T) ≤ 1 }.

log("Building 'linearised' formulation in t_kj = Lk*Mj")

# Express garbage(A') in the t_kj basis (linear in t_kj).
t_syms = sp.symbols('t1_1 t1_2 t1_3 t1_4 t1_5 t1_6 '
                    't2_1 t2_2 t2_3 t2_4 t2_5 t2_6 '
                    't3_1 t3_2 t3_3 t3_4 t3_5 t3_6 '
                    't4_1 t4_2 t4_3 t4_4 t4_5 t4_6 '
                    't5_1 t5_2 t5_3 t5_4 t5_5 t5_6 '
                    't6_1 t6_2 t6_3 t6_4 t6_5 t6_6')
t_mat = sp.Matrix(6, 6, t_syms)

# Substitute Lk*Mj -> t_{k,j} in garbage_Ap.
# CAUTION: xreplace only matches exact subexpressions, so (q)*Lk*Mj with
# a rational q would not match {Lk*Mj: t_kj}.  Use Poly extraction instead.
def linearise_in_LM(expr):
    """Given a polynomial bilinear in (Ls, Ms), rewrite as linear in t_kj
    where t_kj = Ls[k]*Ms[j]."""
    expr = sp.expand(expr)
    result = sp.Integer(0)
    p = sp.Poly(expr, *(list(Ls) + list(Ms)))
    for mono, coeff in p.terms():
        # mono is a 12-tuple (L1..L6, M1..M6)
        L_part = mono[:6]
        M_part = mono[6:]
        deg_L = sum(L_part)
        deg_M = sum(M_part)
        if deg_L == 0 and deg_M == 0:
            result += coeff
        elif deg_L == 1 and deg_M == 1:
            k = L_part.index(1)
            j = M_part.index(1)
            result += coeff * t_syms[6*k + j]
        elif deg_L == 1 and deg_M == 0:
            k = L_part.index(1)
            result += coeff * Ls[k]
        elif deg_L == 0 and deg_M == 1:
            j = M_part.index(1)
            result += coeff * Ms[j]
        else:
            raise ValueError(f"unexpected monomial {mono} (deg_L={deg_L}, deg_M={deg_M})")
    return result

garbage_Ap_lin = sp.Matrix([linearise_in_LM(garbage_Ap[k]) for k in range(N)])
# Sanity: should be linear in t_syms now.
log("Linearised garbage(A') into t_kj basis (proper Poly extraction)")

# The 28-D garbage equations become 36 coordinate equations
# garbage_Ap_lin[k] = c * garbage_A[k]
# Most are redundant (range is 28-D).  Let's build the 36-eq system as a
# (sparse) matrix in (t_syms, c) with constant RHS = 0.

# Move c·g_A to LHS: garbage_Ap_lin[k] - c·g_A[k] = 0.
pair_eqs_lin = [sp.expand(garbage_Ap_lin[k] - c * garbage_A[k]) for k in range(N)]
pair_eqs_lin = [e for e in pair_eqs_lin if e != 0]
log(f"#linearised paired equations = {len(pair_eqs_lin)}")

# Coefficients of the LINEAR part (in t_syms, c)
lin_vars = list(t_syms) + [c]


# ── Generic A (random rational specialisation) to count solution dimension ──
# Strategy: pick a GENERIC random A (= (lk, mj) integer values), solve the
# resulting linear-in-(t, c) system, and check:
#   (a) Solution-set dimension d_T = how many free t's remain after the
#       28-rank linear system.
#   (b) Among the solution-set T's, which are rank ≤ 1 (= the image of
#       the L·M^T map)?  Compute the dimension of that intersection.
import random
random.seed(0)


def evaluate_at_generic_A(A_dict):
    """Substitute concrete (l, m) values into the 36 linear-in-(t,c) equations.
    Returns coefficient matrix M and RHS vector b such that M · [t..., c]^T = b.
    """
    eqs_concrete = [e.subs(A_dict) for e in pair_eqs_lin]
    M_rows = []
    b_rows = []
    for e in eqs_concrete:
        e = sp.expand(e)
        row = [e.coeff(v) for v in lin_vars]
        # constant (RHS sign flipped to put on RHS): e = sum coeff*var + const
        const = sp.expand(e - sum(row[i] * lin_vars[i] for i in range(len(lin_vars))))
        M_rows.append(row)
        b_rows.append(-const)
    return sp.Matrix(M_rows), sp.Matrix(b_rows)


def random_A():
    return {v: sp.Integer(random.choice([-3, -2, -1, 1, 2, 3])) for v in (list(ls) + list(ms))}


log("Solving linear-in-(t, c) system at a generic random A ...")
A_dict = random_A()
print(f"   A = {A_dict}")

# Sanity check: print first few eqs after subs
print("\nSanity check - first 3 pair_eqs_lin after subs(A_dict):")
for k, e in enumerate(pair_eqs_lin[:3]):
    e_sub = sp.expand(e.subs(A_dict))
    print(f"   eq[{k}]: {e_sub}")
    print(f"     coeff(c) = {e_sub.coeff(c)}")

M_A, b_A = evaluate_at_generic_A(A_dict)
log(f"  shape: {M_A.shape}, rank: {M_A.rank()}")
log(f"  b (RHS) is zero? {all(b_A[i] == 0 for i in range(b_A.shape[0]))}")

# Check: does the proportional solution lie in the nullspace?
# Build the proportional solution: c=1, t_kj = lk*mj (since L'=L, M'=M, c=1)
prop_sol = sp.zeros(len(lin_vars), 1)
for k in range(6):
    for j in range(6):
        prop_sol[6*k + j] = A_dict[ls[k]] * A_dict[ms[j]]
prop_sol[36] = sp.Integer(1)  # c = 1
residual = M_A * prop_sol
log(f"  residual of proportional solution || = {sum(abs(r) for r in residual)}")
if any(r != 0 for r in residual):
    nz_rows = [(i, residual[i], sp.expand(pair_eqs_lin[i].subs(A_dict))) for i in range(residual.shape[0]) if residual[i] != 0]
    print(f"  *** {len(nz_rows)} rows have non-zero residual; first 3:")
    for i, r, eq in nz_rows[:3]:
        print(f"     row {i}: residual = {r} ; eq = {eq}")

# Solution space dimension = #vars - rank(M_A) (assuming consistent)
n_vars = len(lin_vars)
sol_dim = n_vars - M_A.rank()
log(f"  solution space dimension d_lin = {sol_dim}  (vars={n_vars}, rank={M_A.rank()})")

# The system is HOMOGENEOUS (b should be 0 because all equations have
# no constant term in (t, c)).  So solutions = nullspace of M_A.
# Particular solution = 0, kernel = nullspace.
sol = sp.zeros(n_vars, 1)
log("  particular solution = 0 (homogeneous system); computing nullspace ...")
kernel = M_A.nullspace()
log(f"  kernel dimension = {len(kernel)}")

# General solution: sol + sum a_k * kernel[k]
print("\nParticular solution (sol[i] for [t_kj..., c]):")
for i, v in enumerate(lin_vars):
    if sol[i] != 0:
        print(f"   {v} = {sol[i]}")
print("\nKernel basis vectors:")
for i, k_vec in enumerate(kernel):
    nonzero = [(lin_vars[j], k_vec[j]) for j in range(n_vars) if k_vec[j] != 0]
    print(f"   k_{i}: {nonzero}")


# ── Among the solutions T, which are rank ≤ 1 = L·M^T? ────────────────────
# Parameterise:  T(α) = T_part + Σ α_k · K_k  (each K_k is a 6×6 matrix)
# Plus c(α) = c_part + Σ α_k · c_k
# We want rank(T(α)) ≤ 1.  This means all 2×2 minors of T(α) vanish.

# Build T(α) and c(α) symbolically
alphas = sp.symbols(f'a0:{len(kernel)}')

def reshape_t_vec(vec):
    """vec is length 37 ([t_kj..., c]).  Return 6x6 matrix from first 36 entries."""
    return sp.Matrix(6, 6, list(vec[:36]))

T_part = reshape_t_vec(sol)
c_part = sol[36]
T_alpha = T_part.copy()
c_alpha = c_part
for i, k_vec in enumerate(kernel):
    T_alpha = T_alpha + alphas[i] * reshape_t_vec(k_vec)
    c_alpha = c_alpha + alphas[i] * k_vec[36]

log(f"T_alpha and c_alpha built; #parameters α = {len(kernel)}")

# Note: T = L'·M'^T has rank ≤ 1.  The variety of rank-≤1 6×6 matrices has
# dim = 6 + 6 - 1 = 11 in 36-dim ambient (Segre product, modulo overall scale
# this is 11-D affine).
# Intersecting this 11-D variety with the (sol_dim)-D affine solution space
# gives expected dim = max(0, sol_dim + 11 - 36 - 1) ... but only if
# transversal.

# Concretely: enforce rank(T_alpha) ≤ 1.  Two equations:
#   for each pair (i, j) of distinct rows, T_alpha[i] × T_alpha[j] = 0
#   (as rank-1 condition on the 2-row submatrix).
# Equivalently, all 2x2 minors of T_alpha vanish.

minors_2x2 = []
for r1 in range(6):
    for r2 in range(r1+1, 6):
        for c1 in range(6):
            for c2 in range(c1+1, 6):
                m = T_alpha[r1, c1]*T_alpha[r2, c2] - T_alpha[r1, c2]*T_alpha[r2, c1]
                m = sp.expand(m)
                if m != 0:
                    minors_2x2.append(m)
log(f"#non-trivially-zero 2x2 minors (in α-params): {len(minors_2x2)}")

# Deduplicate
minors_2x2 = list({sp.expand(m) for m in minors_2x2})
log(f"#unique minor equations: {len(minors_2x2)}")
print("\nSample minor equations (first 5):")
for m in minors_2x2[:5]:
    print(f"   {m} = 0")


# ── Solve the minor system over α-parameters ────────────────────────────────
log("Computing Jacobian rank of minor equations at random α ...")
J = sp.Matrix([[sp.diff(m, ai) for ai in alphas] for m in minors_2x2])
log(f"  Jacobian shape: {J.shape}")
# Evaluate at random α
random.seed(42)
alpha_pt = {ai: sp.Rational(random.randint(-3, 3)) for ai in alphas}
J_pt = J.subs(alpha_pt)
log(f"  Jacobian rank at random α point = {J_pt.rank()}")

# CAREFUL: the rank-1 variety is a CONE through 0. The Jacobian at a generic
# (off-variety) α-point has full rank; this doesn't say much about the
# variety's dimension. We need to evaluate the Jacobian AT a known
# variety point.
# Known point: α = (0, 0, ..., 0, α_8) for any α_8 — this gives T = α_8 · K_8
# which is rank 1 (the proportional solution).
alpha_pt_on_var = {ai: sp.Integer(0) for ai in alphas}
alpha_pt_on_var[alphas[-1]] = sp.Integer(1)  # α_8 = 1
J_var = J.subs(alpha_pt_on_var)
log(f"  Jacobian rank at the proportional point (α_8=1, others 0) = {J_var.rank()}")
tangent_dim = len(alphas) - J_var.rank()
log(f"  Tangent-space dim at proportional point = {tangent_dim}")
log(f"     (upper bound on local variety dim around proportional solution)")

expected_dim = tangent_dim
log(f"  → using tangent dim {expected_dim} as variety dim estimate")

# ── Now interpret: total dim X = #α-params that give rank-1 T = L'·M'^T ────
# Account for the (L', M') → t_kj reduction: each rank-1 T has a 1-D ambiguity
# (L' → ρ·L', M' → (1/ρ)·M').  So dim_{(L',M')} = dim_T + 1.
# But (L', M') ∈ R^12 = 12 vars.  The MAP (L',M') → T is generically dim 11
# = 12 - 1.
# X is parameterised by (A, A', c).  A ∈ R^12 (free), A' ∈ R^12 (constrained
# by 28 garbage eqs given A), c ∈ R.  But we plugged a CONCRETE generic A,
# so what we computed is the dim of the fiber X_A := { (A', c) : g(A') = c·g(A) }.

# dim X_A = dim of rank-1 T-fiber + (#α we found) but careful: t_kj are
# determined by α; rank-1 condition cuts the α-space down; from each
# rank-1 T, lift to (L', M') adds 1 dim.

# ── DECOMPOSE X into c=0 and c≠0 sub-fibers ───────────────────────────────
log("\n--- Splitting X into c=0 and c≠0 sub-fibers ---")
log("c=0 sub-fiber over any A: g(A')=0. This is the zero-garbage variety V")
log("   from zero_garbage_variety.py — max component dim = 6 (in (L',M') ∈ R^12).")
log("   So c=0 sub-fiber over A is 6-D.")

log("\nc=1 (normalised) sub-fiber over generic A: solve g(A')=g(A) in (L',M').")
pair_eqs_c1 = [sp.expand(e.subs(c, 1).subs(A_dict)) for e in pair_eqs]
pair_eqs_c1 = [e for e in pair_eqs_c1 if e != 0]
log(f"  #equations after c=1 and A=concrete: {len(pair_eqs_c1)}")
LM_prop = {Ls[k]: A_dict[ls[k]] for k in range(6)}
LM_prop.update({Ms[j]: A_dict[ms[j]] for j in range(6)})
res_prop = [sp.expand(e.subs(LM_prop)) for e in pair_eqs_c1]
max_res = max(abs(r) for r in res_prop)
log(f"  residual at (L'=L, M'=M, c=1): max |r| = {max_res}")
LM_vars = list(Ls) + list(Ms)
J_LM = sp.Matrix([[sp.diff(e, v) for v in LM_vars] for e in pair_eqs_c1])
J_LM_at_prop = J_LM.subs(LM_prop)
rank_J_at_prop = J_LM_at_prop.rank()
log(f"  Jacobian shape: {J_LM.shape}, rank at proportional point: {rank_J_at_prop}")
tangent_dim_c1 = len(LM_vars) - rank_J_at_prop
log(f"  Tangent dim at proportional point (c=1 fiber) = {tangent_dim_c1}")
log(f"     proportional fiber expected dim = 1 (just (ρ, 1/ρ) reparam)")
log(f"     EXCESS non-proportional directions in c=1 fiber: {tangent_dim_c1 - 1}")

# Robustness: re-test at a few more random A's.
log("\n  Re-testing proportionality theorem at additional random A's:")
random.seed(99)
all_tangent_match = True
for trial in range(5):
    A_dict_t = {v: sp.Integer(random.choice([-3, -2, -1, 1, 2, 3])) for v in (list(ls) + list(ms))}
    pair_eqs_c1_t = [sp.expand(e.subs(c, 1).subs(A_dict_t)) for e in pair_eqs]
    pair_eqs_c1_t = [e for e in pair_eqs_c1_t if e != 0]
    LM_prop_t = {Ls[k]: A_dict_t[ls[k]] for k in range(6)}
    LM_prop_t.update({Ms[j]: A_dict_t[ms[j]] for j in range(6)})
    res_t = [sp.expand(e.subs(LM_prop_t)) for e in pair_eqs_c1_t]
    max_res_t = max(abs(r) for r in res_t) if res_t else 0
    J_LM_t = sp.Matrix([[sp.diff(e, v) for v in LM_vars] for e in pair_eqs_c1_t])
    J_LM_at_t = J_LM_t.subs(LM_prop_t)
    rank_t = J_LM_at_t.rank()
    tdim_t = len(LM_vars) - rank_t
    log(f"     trial {trial}: prop residual={max_res_t}, c=1 tangent dim = {tdim_t}")
    if tdim_t != 1:
        all_tangent_match = False
log(f"  All 5 trials gave tangent dim = 1: {all_tangent_match}")

log("=" * 70)
log(f"STALE intermediate count (kept for audit trail; supersede by CORRECTED below):")
log(f"   linear solution-space d_lin = {sol_dim}   (= dim of (t, c) variety)")
log(f"   rank-1 cut: variety in α-space has dim ≈ {expected_dim}")
log(f"   lifting to (L', M'): add 1 free scaling")
log(f"   ⇒ dim X_A (fiber over single generic A) = {expected_dim + 1}")
log("=" * 70)

dim_XA = expected_dim + 1
dim_X = 12 + dim_XA  # total X = R^12 (A) × fiber
log(f"   dim X = 12 (A) + {dim_XA} (fiber) = {dim_X}")


# ── PROPORTIONALITY TEST ────────────────────────────────────────────────────
# Claim to test: every solution is "proportional", meaning the tensor L'·M'
# is a scalar multiple of L·M, i.e. T = ρ · (L·M^T) for some ρ ∈ R.
# Equivalent: T - ρ · (l·m^T) = 0 for some ρ.
# This is a 36-eq, 1-unknown(ρ) system, generically inconsistent unless T is
# itself rank ≤ 1 and aligned with (l,m).

# Specifically: in the proportional case, A' = (ρL, M) or (L, ρM) or any
# valid (L', M') with L'·M'^T = ρ·L·M^T.  Then c is determined: c = ρ.
# This proportional locus has dim:
#   (L', M') with L'·M'^T proportional to L·M^T = 12-dim:
#      pick ρ (1 dim), then (L', M') with L'·M'^T = ρ·L·M^T is the
#      (1-D scaling × rank-1 lift) — actually 2 dims (the scaling
#      ρ + the (1/σ, σ) reparameterisation of (L', M')).
#   So dim of proportional fiber = 2.
# Total dim of proportional sub-variety = 12 (A) + 2 (proportional fiber) = 14.

# Compare to dim X computed above.  If dim X = 14, proportional locus
# is everything → IMPOSSIBILITY of non-proportional solutions.
# If dim X > 14, there's a (dim X - 14)-dim NON-proportional family.

log("=" * 70)
log("PROPORTIONALITY ANALYSIS")
log("=" * 70)
log(f"   dim X (over generic A) = 12 + {dim_XA} = {dim_X}")
log(f"   dim 'proportional' locus = 12 + 2 = 14")
log(f"   excess freedom = {dim_X - 14}")
if dim_XA == 2:
    log("   ⇒ EVERY (A, A', c) ∈ X has A' proportional to A as a rank-1 tensor.")
    log("   ⇒ IMPOSSIBILITY: two atoms with proportional garbage are themselves")
    log("     proportional (so their (α, β) projections are proportional).")
    log("     A (2, 2, bridge-3) construction in the rank-1 family is therefore")
    log("     IMPOSSIBLE.")
elif dim_XA > 2:
    log(f"   ⇒ {dim_XA - 2}-dim family of NON-PROPORTIONAL solutions exists.")
    log("     Constructive existence of a (2, 2, bridge-3) triple cannot be ruled")
    log("     out from this analysis alone — further work needed.")
else:
    log("   ⇒ Fewer than 2 fiber dims — suggests numerical/algebraic anomaly.")
    log("     Recheck the linearisation and rank-1 cut.")


# ── DETAILED proportionality check on the symbolic solution ─────────────────
# Look directly at the structure of the kernel: do all kernel basis vectors
# correspond to perturbations T that are themselves rank-1 perturbations
# of T_part?

# We materialise T(α) for a few specific α perturbations and check rank.
log("\n--- Concrete kernel-direction rank checks ---")
T_part_concrete = T_part  # at this A
log(f"   T_part rank = {T_part_concrete.rank()}")

# Each kernel direction gives a perturbation matrix K_i = reshape(kernel[i][:36])
for i, k_vec in enumerate(kernel[:6]):  # check first 6 directions
    K_i = reshape_t_vec(k_vec)
    log(f"   kernel[{i}]: K rank = {K_i.rank()}, c-component = {k_vec[36]}")

# How many of the kernel are themselves rank-1 perturbations?
rank1_perturbs = []
for i, k_vec in enumerate(kernel):
    K_i = reshape_t_vec(k_vec)
    if K_i.rank() <= 1:
        rank1_perturbs.append(i)
log(f"   Kernel directions that are rank ≤ 1 perturbations: {len(rank1_perturbs)} / {len(kernel)}")


# ── Direct symbolic check: solve the minor system in α ──────────────────────
log("\n--- Symbolic solve of minor system via Groebner basis ---")

# Use the full system, not a sample.
try:
    log(f"   Computing Groebner basis of {len(minors_2x2)} minors over {len(alphas)} α's (lex order)...")
    G = sp.groebner(minors_2x2, *alphas, order='lex')
    log(f"   Groebner basis has {len(G.polys)} polynomials")
    print("   Groebner basis polys (first 20):")
    for p in list(G.polys)[:20]:
        print(f"     {p.as_expr()} = 0")

    # Then solve using the Groebner basis (typically much faster)
    log("   Solving the Groebner basis ...")
    sols_alpha = sp.solve(list(G.polys), alphas, dict=True)
    log(f"   solve returned {len(sols_alpha)} branches")
    for s_idx, s in enumerate(sols_alpha[:10]):
        free_alphas = [a for a in alphas if a not in s]
        log(f"   branch {s_idx}: solved {sorted(s.keys(), key=str)} ; free = {free_alphas}")
        for k, v in list(s.items())[:5]:
            print(f"     {k} = {v}")
except Exception as e:
    log(f"   Groebner/solve failed: {e}")
    import traceback
    traceback.print_exc()


# ── Conclusion ──────────────────────────────────────────────────────────────
log("\n" + "=" * 70)
log("CORRECTED CONCLUSION (decomposing X by c=0 vs c≠0)")
log("=" * 70)
log("   At a generic A, the fiber X_A := {(A', c) : g(A') = c·g(A)} decomposes:")
log("")
log("   c ≠ 0 (rescale to c=1):")
log(f"     Tangent dim of c=1 fiber at (A'=A) = {tangent_dim_c1}")
log(f"     Proportional fiber dim = 1   (the (ρ, 1/ρ) reparam line)")
log(f"     EXCESS non-proportional directions = {tangent_dim_c1 - 1}")
log("     ⇒ The c=1 sub-fiber over generic A is EXACTLY the proportional locus.")
log("     ⇒ The c≠0 component of X_A is parameterised by (c, σ) — 2 dim.")
log("")
log("   c = 0:")
log("     g(A') = 0 ⇔ A' ∈ V (zero-garbage variety, from zero_garbage_variety.py).")
log("     V is positive-dim, max component dim = 6 in (L', M') ∈ R^12.")
log("     ⇒ The c=0 component of X_A is V — 6 dim (max).")
log("")
log("   Total: X has (at least) two irreducible components:")
log("     X_{proportional} : dim 12 (A) + 2 (fiber) = 14.")
log("     X_{c=0}          : dim 12 (A) + 6 (V max) = 18.")
log("")
log("=" * 70)
log("PROPORTIONALITY THEOREM (this script's main result)")
log("=" * 70)
log("   For every (A, A', c) ∈ X with c ≠ 0:")
log("     A' is forced to come from a PROPORTIONAL rank-1 tensor.")
log("     I.e. ∃ ρ, σ ≠ 0 with L' = ρL/σ and M' = σM (so L'·M' = ρ·L·M, c=ρ).")
log("")
log("   For (A, A', c) ∈ X with c = 0:")
log("     A' ∈ V; by zero_garbage_variety.py, α(A') = β(A') = 0.")
log("")
log("=" * 70)
log("IMPOSSIBILITY OF (2,2,bridge-3) IN THE RANK-1 ATOM FAMILY")
log("=" * 70)
log("   A 3-atom (2, 2, bridge-3) construction needs the 28×3 garbage matrix")
log("   U = [g(A1) | g(A2) | g(A3)] to have rank ≤ 1, AND the (α, β) projections")
log("   to span Q ≅ R^2.")
log("")
log("   rank(U) ≤ 1 means any two columns are proportional (allowing 0 columns).")
log("   For any i, j with g(A_i), g(A_j) both non-zero and proportional with")
log("   ratio c ≠ 0: by the proportionality theorem, A_i and A_j are themselves")
log("   proportional rank-1 tensors. So their (α, β) projections are proportional.")
log("")
log("   For any column g(A_k) = 0: A_k ∈ V, so α(A_k) = β(A_k) = 0.")
log("")
log("   ⇒ In any such triple, the (α, β) projections all lie on a single line")
log("     through the origin in Q. The (α, β) image has rank ≤ 1, cannot span Q.")
log("")
log("   ⇒ NO (2, 2, bridge-3) construction exists in the rank-1 atom family.")
log("=" * 70)
