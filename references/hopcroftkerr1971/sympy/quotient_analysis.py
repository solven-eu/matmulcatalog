"""
HK (2,2,bridge-3) quotient analysis.

Tasks (from user spec):
  1. Explicit basis (q1, q2) for Q = span(T1, T2, S) / S; express T1, T2 in it.
  2. Projection of a general rank-1 L·M into Q as α(λ,μ) q1 + β(λ,μ) q2.
  3. Image geometry of (λ,μ) → (α, β): all of Q? Zariski-dense? on a curve?
     Generic Jacobian rank?
  4. Consequence: impossibility result for a 3-product (2,2,bridge-3)
     derivation in any rank-1 atom family, or constructive recipe.

Run: python3 quotient_analysis.py
"""
import sympy as sp


# ── Setup ───────────────────────────────────────────────────────────────────
ai1, ai2, ap1, ap2, aj1, aj2 = sp.symbols('ai1 ai2 ap1 ap2 aj1 aj2')
xi1, xi2, xp1, xp2, xj1, xj2 = sp.symbols('xi1 xi2 xp1 xp2 xj1 xj2')

a_vars = [ai1, ai2, ap1, ap2, aj1, aj2]
x_vars = [xi1, xi2, xp1, xp2, xj1, xj2]

# 36 monomials a_α_β · x_α'_β' as the basis of our ambient tensor space
mono = [(av, xv) for av in a_vars for xv in x_vars]
N = len(mono)


def to_vec(expr):
    """Map a polynomial in (a_*, x_*) to its coordinate vector in R^36."""
    expr = sp.expand(expr)
    out = []
    for av, xv in mono:
        # coefficient of av*xv in `expr`
        c = expr.coeff(av).coeff(xv)
        out.append(sp.nsimplify(c))
    return sp.Matrix(out)


# Targets
T1 = ai2 * xi1 + ap2 * xp1 + aj2 * xj1
T2 = ai1 * xi2 + ap1 * xp2 + aj1 * xj2

# Reusable products (the "shared" set S)
Ci = (ai1 - ai2) * xi2
Di = ai1 * (xi1 + xi2)
Cj = (aj1 - aj2) * xj2
Dj = aj1 * (xj1 + xj2)
Ep = ap2 * xp2
Fp = ap1 * xp1
shared = [Ci, Di, Cj, Dj, Ep, Fp]
shared_names = ['Ci', 'Di', 'Cj', 'Dj', 'Ep', 'Fp']

vT1 = to_vec(T1)
vT2 = to_vec(T2)
S_mat = sp.Matrix.hstack(*[to_vec(s) for s in shared])

print(f"Ambient dim = {N}, |S| = {S_mat.shape[1]}, rank(S) = {S_mat.rank()}")

# Confirm T1, T2 contribute 2 new dimensions modulo S
full = sp.Matrix.hstack(S_mat, vT1, vT2)
print(f"rank(S, T1, T2) = {full.rank()}  (expect rank(S) + 2 = {S_mat.rank() + 2})")


# ── Task 1: explicit basis (q1, q2) for Q = span(T1, T2, S) / S ─────────────
# Strategy: orthogonal-projection onto S^⊥ (Euclidean inner product on R^36).
# Then q1 := P_⊥ vT1, q2 := P_⊥ vT2 give explicit representatives in S^⊥
# of the cosets T1+S, T2+S — these are basis vectors of Q ≅ S^⊥ ∩ span(T1,T2,S).

SS_inv = (S_mat.T * S_mat).inv()
P_S = S_mat * SS_inv * S_mat.T            # projector onto span(S)
P_perp = sp.eye(N) - P_S                  # projector onto S^⊥
q1 = sp.simplify(P_perp * vT1)
q2 = sp.simplify(P_perp * vT2)
Q_mat = sp.Matrix.hstack(q1, q2)
print(f"\n[Task 1] dim Q = rank([q1,q2]) = {Q_mat.rank()}  (must be 2)")


def render_in_mono(vec):
    """Pretty-print a 36-vector as a sum of a*x monomials with non-zero coefs."""
    terms = []
    for k, (av, xv) in enumerate(mono):
        c = sp.simplify(vec[k])
        if c == 0:
            continue
        if c == 1:
            terms.append(f"{av}*{xv}")
        elif c == -1:
            terms.append(f"-{av}*{xv}")
        else:
            terms.append(f"({c})*{av}*{xv}")
    return " + ".join(terms) if terms else "0"


print(f"\nq1 (= P_⊥ T1):  {render_in_mono(q1)}")
print(f"q2 (= P_⊥ T2):  {render_in_mono(q2)}")
print("\nIn this basis, by construction:")
print("  T1 ≡ 1·q1 + 0·q2  (mod S)")
print("  T2 ≡ 0·q1 + 1·q2  (mod S)")


# ── Task 2: project a general rank-1 L·M onto Q ─────────────────────────────
ls = sp.symbols('l1:7')
ms = sp.symbols('m1:7')

L = sum(ls[k] * a_vars[k] for k in range(6))
M = sum(ms[k] * x_vars[k] for k in range(6))
R = sp.expand(L * M)
vR = to_vec(R)

# Project to Q (= S^⊥ ∩ span(T1,T2,S)). After P_⊥, the part outside
# span(T1,T2,S) is irrelevant: we want the [q1, q2]-coordinates only.
vR_perp = P_perp * vR

# Solve [q1 q2] [α; β] = vR_perp in least-squares sense
# (vR_perp may have a component outside span(q1,q2); we project to span(q1,q2))
QQ = Q_mat.T * Q_mat
alpha_beta = sp.simplify(QQ.inv() * Q_mat.T * vR_perp)
alpha = sp.simplify(sp.expand(alpha_beta[0]))
beta = sp.simplify(sp.expand(alpha_beta[1]))

print(f"\n[Task 2] α(λ,μ) =\n  {alpha}")
print(f"\n[Task 2] β(λ,μ) =\n  {beta}")


# ── Task 3: image geometry — Jacobian rank ──────────────────────────────────
all_vars = list(ls) + list(ms)
J = sp.Matrix([[sp.diff(alpha, v) for v in all_vars],
               [sp.diff(beta,  v) for v in all_vars]])

# Symbolic rank is expensive; evaluate at a random rational point first
import random
random.seed(0)
subs_pt = {v: sp.Rational(random.randint(-5, 5)) for v in all_vars}
# Avoid the all-zero / degenerate point
for v in all_vars:
    if subs_pt[v] == 0:
        subs_pt[v] = sp.Rational(1)

J_at_pt = J.subs(subs_pt)
print(f"\n[Task 3] Jacobian J = ∂(α,β)/∂(λ,μ) is a 2×{len(all_vars)} matrix.")
print(f"         rank at random integer point {dict((str(v), subs_pt[v]) for v in all_vars)}:")
print(f"         rank = {J_at_pt.rank()}")

# Try a SECOND random point to be sure
random.seed(42)
subs_pt2 = {v: sp.Rational(random.randint(-5, 5)) for v in all_vars}
for v in all_vars:
    if subs_pt2[v] == 0:
        subs_pt2[v] = sp.Rational(1)
J_at_pt2 = J.subs(subs_pt2)
print(f"         rank at a second random point: {J_at_pt2.rank()}")

# Generic rank: try symbolic rank of J (may be slow)
print("\n         Attempting symbolic rank of J (may take a moment) …")
try:
    sym_rank = J.rank()
    print(f"         symbolic rank(J) = {sym_rank}")
except Exception as e:
    print(f"         symbolic rank failed: {e}")


# ── Task 4: consequence ─────────────────────────────────────────────────────
print("\n[Task 4] Consequences:")
print("  • If generic rank(J) = 2, the map (λ,μ) → (α,β) is locally surjective")
print("    on Q ≅ R², so its image is Zariski-dense in Q. There exist rank-1")
print("    atoms whose projections span Q ⇒ a 3-atom solution exists (in fact 2 suffice).")
print("  • If generic rank(J) = 1, all rank-1 atom projections lie on a 1-D curve")
print("    inside Q ⇒ NO finite number of rank-1 atoms can span Q ⇒")
print("    (2,2,bridge-3) admits NO 3-product solution in any rank-1 family.")
print("  • If generic rank(J) = 0, α and β are constant on the rank-1 variety —")
print("    extremely unlikely given the explicit formulas above.")
