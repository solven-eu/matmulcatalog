"""
HK (2,2,bridge-3) quotient analysis — CORRECTED.

The first version (quotient_analysis.py) computed only the (α, β) component
of the projection of a rank-1 tensor onto Q. That answer is necessary but
NOT sufficient: each rank-1 atom L·M, when projected to S^⊥, generically
has 30 non-zero coordinates — only 2 of them lie in Q = span(q1, q2);
the other 28 are "garbage" in (S^⊥ \\ Q).

For three atoms P_1, P_2, P_3 to combine via output weights W into
exactly (T1, T2) mod S, we need BOTH:
  (1) Their Q-projections to span Q (the (α, β) condition).
  (2) Their 28-D garbage components to lie in a subspace small enough
      that W has 2 linearly independent left-null vectors that ALSO
      satisfy the Q condition.

Specifically: if u_k ∈ R^28 is atom k's garbage, we need W_1, W_2 ∈ R^3 to
satisfy  W_1·U = 0, W_2·U = 0  (where U = [u_1 u_2 u_3]) AND
W_1·(α_k) = 1, W_1·(β_k) = 0, W_2·(α_k) = 0, W_2·(β_k) = 1.

The W left-null-space of U has dimension 3 - rank(U). We need it ≥ 2,
i.e. rank(U) ≤ 1.

This script checks:
  • Geometric condition: does a single rank-1 atom L·M have ZERO garbage?
    (i.e. is proj_{S^⊥ \\ Q}(L·M) = 0 ever achievable for L, M ≠ 0?)
  • If atoms with zero garbage exist as a positive-dimensional variety,
    then 3 atoms can be picked from it whose (α, β) span Q → solution exists.

Run: python3 quotient_analysis_v2.py
"""
import sympy as sp


# ── Setup (same as v1) ──────────────────────────────────────────────────────
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
P_perp_S = sp.eye(N) - S_mat * (S_mat.T * S_mat).inv() * S_mat.T  # 36→36, proj onto S^⊥
q1 = sp.simplify(P_perp_S * vT1)
q2 = sp.simplify(P_perp_S * vT2)
Q_mat = sp.Matrix.hstack(q1, q2)

# ── Build the projector onto Q (a 2-D subspace of S^⊥) ─────────────────────
# Combined orthonormal-ish: in the inherited inner product on R^36.
P_Q = Q_mat * (Q_mat.T * Q_mat).inv() * Q_mat.T  # 36→36, proj onto span(q1,q2)

# Projector onto (S^⊥ \ Q) = S^⊥ ∩ Q^⊥. In the ambient: P_perp_S - P_Q.
P_garbage = P_perp_S - P_Q  # 36→36, proj onto the 28-dim "garbage" subspace
print(f"rank(P_garbage) = {P_garbage.rank()}  (expect 36 - 6 - 2 = 28)")


# ── Generic rank-1 atom L · M ──────────────────────────────────────────────
ls = sp.symbols('l1:7')
ms = sp.symbols('m1:7')
L = sum(ls[k] * a_vars[k] for k in range(6))
M = sum(ms[k] * x_vars[k] for k in range(6))
R = sp.expand(L * M)
vR = to_vec(R)


# ── The garbage component of vR ────────────────────────────────────────────
garbage = sp.simplify(P_garbage * vR)
# Pick a *basis* of the 28-D garbage subspace — any 28 linearly independent
# columns of P_garbage will do. We just need to know the image dimension as
# a function of (λ, μ).
# Easier: count the number of independent non-zero coordinates of garbage
# as polynomials in (λ, μ).
print(f"\ngarbage vector has {sum(1 for k in range(N) if sp.simplify(garbage[k]) != 0)}"
      f" non-zero coordinates (as polynomials in λ, μ)")


# ── Variety of "zero garbage" rank-1 atoms ─────────────────────────────────
# We want to know: does the system { garbage[k] = 0 for k in 0..N-1 } have
# non-trivial solutions in (λ, μ)? Each garbage[k] is bilinear in (λ, μ).
# The variety is cut out by these bilinear equations.
all_vars = list(ls) + list(ms)
garbage_polys = [sp.simplify(garbage[k]) for k in range(N) if sp.simplify(garbage[k]) != 0]
print(f"\n{len(garbage_polys)} independent garbage-vanishing equations to solve.")
print("First 5 equations:")
for p in garbage_polys[:5]:
    print(f"   {p} = 0")


# Try to solve symbolically — bilinear in (λ_1..λ_6, μ_1..μ_6).
# Sympy's solve may be slow on 28+ bilinear equations. Strategy: try a small
# number, see if the solution set is non-trivial.
print("\nAttempting partial symbolic solve …")
try:
    # Solve first ~10 equations to see structure
    sols = sp.solve(garbage_polys[:10], all_vars, dict=True)
    print(f"  partial solve over first 10 eqs returned {len(sols)} solution branches")
    for i, s in enumerate(sols[:3]):
        print(f"    branch {i}: {s}")
except Exception as e:
    print(f"  partial solve failed: {e}")


# ── Numerical check: evaluate garbage at many random points ────────────────
# If the garbage is generically non-zero (it should be), then no single
# rank-1 atom has zero garbage. The question becomes: can THREE atoms be
# chosen such that their garbages span a 1-D subspace and (α, β) span Q?
import random
random.seed(0)

garbage_norm_samples = []
for _ in range(50):
    subs = {v: sp.Rational(random.randint(-3, 3)) for v in all_vars}
    g = garbage.subs(subs)
    nrm = sum(abs(g[k]) for k in range(N))
    garbage_norm_samples.append(float(nrm))

zero_count = sum(1 for n in garbage_norm_samples if n == 0)
nonzero_count = len(garbage_norm_samples) - zero_count
print(f"\nNumerical: out of 50 random integer-{{-3..3}}-valued rank-1 atoms,")
print(f"  {zero_count} have ZERO garbage (= live entirely in span(S, q1, q2))")
print(f"  {nonzero_count} have non-zero garbage")


# ── Rank of garbage-map as a function of (λ, μ) ────────────────────────────
# As (λ, μ) varies, the garbage vector spans a subset of R^28. Its rank
# (as a polynomial map) tells us:
#   • rank = 0: garbage is identically zero ⇒ ALL atoms are valid ⇒ trivial existence.
#   • rank = k < 28: garbage hits a k-dim subspace; need k+1 atoms to cancel.
#   • rank = 28: garbage spans all of R^28; cancellation requires very specific tuning.
print("\nComputing Jacobian rank of (λ,μ) → garbage(λ,μ) …")
J_garbage = sp.Matrix([[sp.diff(g, v) for v in all_vars] for g in garbage_polys])
print(f"  Jacobian is {J_garbage.shape[0]} × {J_garbage.shape[1]} (= garbage eqs × λ,μ)")
J_at_pt = J_garbage.subs({v: sp.Rational(random.randint(1, 5)) for v in all_vars})
print(f"  rank at random point = {J_at_pt.rank()}")


# ── Final analysis ─────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("INTERPRETATION:")
print("=" * 70)
print("""
The (α, β) Q-projection alone is necessary but NOT sufficient for a
3-atom (2,2,bridge-3) construction. Each atom L·M, projected to S^⊥,
carries a 28-D 'garbage' component in (S^⊥ \\ Q) that must cancel in
the W-weighted sum of the three atoms.

Specifically: let U = [u_1 u_2 u_3] be the 28×3 matrix of the three
atoms' garbage components. The output weights W_1, W_2 ∈ R^3 must lie
in the left null space of U (a 3 - rank(U)-dimensional subspace), AND
must produce (1, 0) and (0, 1) on the (α, β) projections.

For a solution to exist:
  rank(U) ≤ 1   AND   the (α, β) 2×3 matrix has rank 2 with W ∈ ker(U^T).

Conclusions hinge on:
  • Whether atoms with zero garbage form a positive-dimensional variety.
  • If not, whether the 'garbage' image is low-dimensional enough that
    three atoms can be picked with rank-1 garbage AND linearly
    independent (α, β).

A definitive impossibility result requires showing rank(U) ≥ 2 for ALL
choices of three rank-1 atoms with (α, β)-projections spanning Q — a
much stronger statement than "the (α, β) Jacobian has rank 2".
""")
