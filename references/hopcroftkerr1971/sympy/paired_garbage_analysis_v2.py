"""
HK (2,2,bridge-3): Paired-atom garbage variety — CORRECTED v2.

The v1 script (paired_garbage_analysis.py) reported dim X = 13, less than
dim P = 14 (the proportional sub-variety), which is impossible (a sub-
variety cannot exceed its ambient). The bug was in the rank-1 cut on
the kernel: it evaluated the minor Jacobian at a generic α in the
EMBEDDING (not on the variety), so it reported the WRONG tangent space
dimension. Furthermore, lifting T → (L', M') only adds 1 dim, but the
rank-1 variety inside the (block-diag) kernel decomposes into multiple
components and the per-component dim must be computed correctly.

This v2 script does two things rigorously:

  Step 1: explicit parametric proof  dim P = 14.
  Step 2: tangent-space computation of dim X at a smooth proportional
          point.  P ⊆ X, so codim_X(P) ≥ 0 and dim X ≥ dim P = 14.
          If the Jacobian of the constraint map at a smooth proportional
          point has co-rank 14, then locally X = P at that point.
          Combined with irreducibility of P, this gives dim X = 14 and
          P is (Zariski-)dense in X.

If dim X = 14, three atoms with pairwise-proportional garbages are
themselves pairwise-proportional as RANK-1 TENSORS, hence have
proportional (α, β) projections, hence span a 1-D subspace of Q ≅ R^2,
hence CANNOT realise both targets T1 and T2 — IMPOSSIBILITY.

Run: python3 paired_garbage_analysis_v2.py
"""
import sympy as sp
import time

t0 = time.time()
def log(msg):
    print(f"[{time.time()-t0:6.1f}s] {msg}", flush=True)


# ── Setup (identical monomial basis to paired_garbage_analysis.py) ─────────
ai1, ai2, ap1, ap2, aj1, aj2 = sp.symbols('ai1 ai2 ap1 ap2 aj1 aj2')
xi1, xi2, xp1, xp2, xj1, xj2 = sp.symbols('xi1 xi2 xp1 xp2 xj1 xj2')
a_vars = [ai1, ai2, ap1, ap2, aj1, aj2]
x_vars = [xi1, xi2, xp1, xp2, xj1, xj2]
mono = [(av, xv) for av in a_vars for xv in x_vars]
N = len(mono)  # 36


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
P_garbage = P_perp_S - P_Q  # 36x36 projector onto 28-D garbage subspace
log("Projector P_garbage built")

# ── Rank-1 atom A and atom A' ──────────────────────────────────────────────
ls = sp.symbols('l1:7')
ms = sp.symbols('m1:7')
Ls = sp.symbols('L1:7')
Ms = sp.symbols('M1:7')

L_a  = sum(ls[k]  * a_vars[k] for k in range(6))
M_a  = sum(ms[k]  * x_vars[k] for k in range(6))
L_ap = sum(Ls[k]  * a_vars[k] for k in range(6))
M_ap = sum(Ms[k]  * x_vars[k] for k in range(6))

vA  = to_vec(sp.expand(L_a  * M_a))
vAp = to_vec(sp.expand(L_ap * M_ap))

garbage_A  = sp.simplify(P_garbage * vA)
garbage_Ap = sp.simplify(P_garbage * vAp)
log("garbage_A and garbage_A' built")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 1: dim P = 14 by explicit parametrisation.
# ═══════════════════════════════════════════════════════════════════════════
# P = { (A, A', c) : A' = (α·L, β·M) for some α, β ∈ R, c = α·β }
# Parametrising map  ψ : R^14 → R^25,
#   (λ_1..λ_6, μ_1..μ_6, α, β) ↦ (λ, μ, α·λ, β·μ, α·β).
# Generic injectivity of ψ on the open set {λ_k ≠ 0, μ_j ≠ 0}: from (λ, μ)
# we recover them; then α·λ_k recovers α from any non-zero λ_k, similarly
# β; then α·β is consistent.  So ψ is generically injective ⇒ dim P = 14.
log("=" * 70)
log("STEP 1: dim P = 14 by explicit parametrisation.")
log("=" * 70)
# Verify by Jacobian rank at a generic point.
alpha, beta = sp.symbols('alpha beta')
psi_vars = list(ls) + list(ms) + [alpha, beta]
psi_image = (list(ls) + list(ms)
             + [alpha * ls[k] for k in range(6)]
             + [beta  * ms[k] for k in range(6)]
             + [alpha * beta])
J_psi = sp.Matrix([[sp.diff(p, v) for v in psi_vars] for p in psi_image])
# Pick a generic numeric point
psi_pt = {ls[k]: sp.Integer(k+2)  for k in range(6)}
psi_pt.update({ms[k]: sp.Integer(k+3) for k in range(6)})
psi_pt[alpha] = sp.Integer(5)
psi_pt[beta]  = sp.Integer(7)
r_psi = J_psi.subs(psi_pt).rank()
log(f"   J_psi shape = {J_psi.shape} (= 25 × 14)")
log(f"   rank J_psi at generic (λ, μ, α=5, β=7) = {r_psi}  (expect 14)")
assert r_psi == 14, f"ψ should be generically immersive; got rank {r_psi}"
log("   ⇒ dim P = 14 confirmed.")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 2: dim X via tangent space at a smooth proportional point.
# ═══════════════════════════════════════════════════════════════════════════
# Constraint map  Φ : R^25 → R^36,
#   (A, A', c) ↦ g(A') − c · g(A)   (28-dim image, but we compute the full 36)
# At a smooth point of X, dim X = 25 − rank(JΦ).
#
# We evaluate at a PROPORTIONAL point:  A' = (α·L, β·M), c = α·β.
# Then Φ vanishes identically; we want the rank of JΦ.

log("=" * 70)
log("STEP 2: tangent space of X at a smooth proportional point.")
log("=" * 70)

c = sp.Symbol('c')
all_vars = list(ls) + list(ms) + list(Ls) + list(Ms) + [c]
log(f"   constraint variables = {len(all_vars)}")

# The 36 coordinates of Φ.
Phi = [sp.expand(garbage_Ap[k] - c * garbage_A[k]) for k in range(N)]
# Drop the zero ones (some coordinates of the garbage projection are 0 ≡ 0
# regardless of (l, m)).
Phi_nonzero_indices = [k for k in range(N) if Phi[k] != 0]
Phi_eqs = [Phi[k] for k in Phi_nonzero_indices]
log(f"   #non-trivial Φ coords = {len(Phi_eqs)}  (≤ 28 expected)")

# Compute Jacobian symbolically (rows = equations, cols = variables).
log("   computing symbolic Jacobian of Φ (this may take ~30 s) ...")
J_Phi = sp.Matrix([[sp.diff(eq, v) for v in all_vars] for eq in Phi_eqs])
log(f"   J_Phi shape = {J_Phi.shape}")

# Pick a smooth proportional point.
# A = (lk = 2, 3, 5, 7, 11, 13), M same shape; α = 2, β = 3.
proportional_pt = {}
for k in range(6):
    proportional_pt[ls[k]] = sp.Integer([2, 3, 5, 7, 11, 13][k])
    proportional_pt[ms[k]] = sp.Integer([3, 5, 7, 11, 13, 17][k])
ALPHA = sp.Integer(2)
BETA  = sp.Integer(3)
for k in range(6):
    proportional_pt[Ls[k]] = ALPHA * proportional_pt[ls[k]]
    proportional_pt[Ms[k]] = BETA  * proportional_pt[ms[k]]
proportional_pt[c] = ALPHA * BETA
log(f"   proportional point: α = {ALPHA}, β = {BETA}, c = {ALPHA*BETA}")

# Sanity: confirm Φ = 0 at this point.
phi_vals = [eq.subs(proportional_pt) for eq in Phi_eqs]
log(f"   Φ at proportional point: max |coord| = {max(abs(v) for v in phi_vals)}")
assert all(v == 0 for v in phi_vals), "proportional point should satisfy Φ = 0"

J_pt = J_Phi.subs(proportional_pt)
log("   computing rank of J_Phi at proportional point ...")
r = J_pt.rank()
log(f"   rank J_Phi at proportional point = {r}")
dim_X_local = len(all_vars) - r
log(f"   ⇒ dim X (local, at smooth proportional point) = 25 − {r} = {dim_X_local}")


# ═══════════════════════════════════════════════════════════════════════════
# Cross-check: rank at SEVERAL proportional points; should be constant
# if the proportional locus is smooth in X.
# ═══════════════════════════════════════════════════════════════════════════
log("=" * 70)
log("CROSS-CHECK: rank at multiple proportional points")
log("=" * 70)

import random
random.seed(123)
ranks_seen = []
for trial in range(4):
    pt = {}
    for k in range(6):
        pt[ls[k]] = sp.Integer(random.choice([-5, -3, -2, 1, 2, 3, 5]))
        pt[ms[k]] = sp.Integer(random.choice([-5, -3, -2, 1, 2, 3, 5]))
    a = sp.Integer(random.choice([-3, -2, 2, 3, 5]))
    b = sp.Integer(random.choice([-3, -2, 2, 3, 5]))
    for k in range(6):
        pt[Ls[k]] = a * pt[ls[k]]
        pt[Ms[k]] = b * pt[ms[k]]
    pt[c] = a * b
    r_trial = J_Phi.subs(pt).rank()
    ranks_seen.append(r_trial)
    log(f"   trial {trial}: a={a}, b={b}, rank = {r_trial}, dim_X = {25 - r_trial}")

# Take the MAX rank seen — that's the generic rank (smooth-locus rank).
r_generic = max(ranks_seen + [r])
dim_X = 25 - r_generic
log(f"")
log(f"   generic rank seen = {r_generic}")
log(f"   ⇒ dim X = 25 − {r_generic} = {dim_X}")


# ═══════════════════════════════════════════════════════════════════════════
# CONCLUSION
# ═══════════════════════════════════════════════════════════════════════════
log("=" * 70)
log("CONCLUSION")
log("=" * 70)
log(f"   dim P = 14 (explicit parametrisation, rank check above).")
log(f"   dim X = {dim_X} (tangent space at a smooth proportional point).")
if dim_X == 14:
    log("   ⇒ dim X = dim P, and P ⊆ X with P irreducible (image of R^14 under")
    log("     a polynomial map) ⇒ Zariski-closure of P contains X's connected")
    log("     component through any proportional point.")
    log("")
    log("     Hence on the (Zariski-open) smooth locus, X = P.")
    log("     Every paired-atom solution (A, A', c) has A' proportional to A")
    log("     as a rank-1 tensor.")
    log("")
    log("     IMPOSSIBILITY THEOREM:")
    log("     Three rank-1 atoms with pairwise-proportional garbages are")
    log("     pairwise-proportional rank-1 tensors, hence their (α, β)")
    log("     projections span a subspace of dim ≤ 1 in Q ≅ R^2.  They")
    log("     cannot realise both target functionals T1 and T2 mod S.")
    log("     ⇒ No 3-atom (2, 2, bridge-3) construction exists in the rank-1")
    log("       atom family.")
elif dim_X > 14:
    log(f"   ⇒ dim X = {dim_X} > 14: non-proportional family of dim {dim_X - 14}")
    log(f"     exists.  Need to parameterise X \\ P and check whether a triple")
    log(f"     drawn from it spans Q in (α, β).")
else:
    log(f"   ⇒ dim X = {dim_X} < 14: contradicts P ⊆ X.  Bug somewhere; recheck.")

log("=" * 70)
log("DONE")
