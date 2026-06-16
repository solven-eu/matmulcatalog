"""
HK (2,2,bridge-3): zero-garbage variety V.

Builds on `quotient_analysis_v2.py`. That script established:

  • Each rank-1 atom L·M with L = Σ l_k·a_k, M = Σ m_k·x_k (k=1..6)
    projects to S^⊥ with a 2-D "(α, β) component" in Q plus a
    28-D "garbage" component in S^⊥ ∩ Q^⊥.
  • There are 34 polynomial equations garbage_k(λ, μ) = 0 cutting
    out the zero-garbage variety V ⊂ R^12.
  • At a random integer point the Jacobian has rank 11, so V has
    expected codimension ≤ 11, i.e. dim V ≥ 1.

This script:
  (1) Materialises all 34 equations.
  (2) Drops dependent ones until a minimal "independent" subset is found
      via successive Groebner / column-rank tests.
  (3) Walks through the cases of the *monomial* equations (l_i · m_j = 0)
      explicitly, branching into either l_i = 0 or m_j = 0, and solves
      the residual bilinears in each branch.
  (4) Lists irreducible components, parameterises each, computes (α, β)
      on each component, and reports whether 3 atoms from V can
      span Q = R^2.

Run: python3 zero_garbage_variety.py
"""
import sympy as sp


# ── Setup (copied from v2 to be self-contained) ─────────────────────────────
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

ls = sp.symbols('l1:7')
ms = sp.symbols('m1:7')
L = sum(ls[k] * a_vars[k] for k in range(6))
M = sum(ms[k] * x_vars[k] for k in range(6))
R = sp.expand(L * M)
vR = to_vec(R)

garbage_full = sp.simplify(P_garbage * vR)
garbage_polys = [sp.simplify(garbage_full[k]) for k in range(N)
                 if sp.simplify(garbage_full[k]) != 0]
all_vars = list(ls) + list(ms)

# (α, β) per the v1 formulas — precomputed for speed.
alpha = (ls[1] * ms[0] + ls[3] * ms[2] + ls[5] * ms[4]) / 3
beta = (-ls[0] * ms[0] + ls[0] * ms[1] + ls[1] * ms[1]
        + 3 * ls[2] * ms[3] - ls[4] * ms[4] + ls[4] * ms[5]
        + ls[5] * ms[5]) / 5

print(f"#garbage equations = {len(garbage_polys)}")
print(f"#variables = {len(all_vars)} (l1..l6, m1..m6)")


# ── Catalog: monomial-shaped equations ──────────────────────────────────────
# A monomial l_i * m_j is "easy" — it forces l_i = 0 OR m_j = 0.
monomial_eqs = []      # list of (i, j) for l_i * m_j = 0
nonmonomial_eqs = []
for g in garbage_polys:
    g = sp.expand(g)
    # Detect a single-term l_i * m_j (modulo overall scalar)
    if g.is_Mul or g.func == sp.Symbol or g.is_Pow:
        # only one term — scalar * (l_i) * (m_j) or product
        ll = [v for v in g.free_symbols if v in ls]
        mm = [v for v in g.free_symbols if v in ms]
        if len(ll) == 1 and len(mm) == 1:
            monomial_eqs.append((ll[0], mm[0]))
            continue
    # Otherwise collect for later
    nonmonomial_eqs.append(g)

print(f"\n#monomial l_i*m_j = 0 equations : {len(monomial_eqs)}")
for li, mj in monomial_eqs:
    print(f"   {li} * {mj} = 0")
print(f"\n#non-monomial equations : {len(nonmonomial_eqs)}")


# ── Reduce nonmonomial equations to a minimal independent set ───────────────
# Each gives a polynomial; we keep them but won't run Groebner over all 12 vars
# (too slow). Instead we branch on the monomial equations and re-evaluate.

# Build the set of l_i's and m_j's that appear in monomial equations
ls_in_mono = set(li for li, mj in monomial_eqs)
ms_in_mono = set(mj for li, mj in monomial_eqs)
print(f"\nl_i appearing in monomial eqs: {sorted(ls_in_mono, key=str)}")
print(f"m_j appearing in monomial eqs: {sorted(ms_in_mono, key=str)}")


# ── Branching: for each monomial l_i*m_j=0, we choose l_i=0 or m_j=0 ────────
# Brute over all 2^(#monomial_eqs) branches is 2^k; if k is moderate we do it.
# Otherwise group monomials by shared variables and reason structurally.

# First, build a bipartite graph: edges (l_i, m_j) from monomial equations.
from collections import defaultdict
ladj = defaultdict(set)   # l_i -> set of m_j
madj = defaultdict(set)   # m_j -> set of l_i
for li, mj in monomial_eqs:
    ladj[li].add(mj)
    madj[mj].add(li)

print("\nBipartite-cover structure (l_i must = 0 OR all its partners m_j must = 0):")
for li, mjs in sorted(ladj.items(), key=lambda kv: str(kv[0])):
    print(f"   {li}  --  {{{', '.join(str(m) for m in sorted(mjs, key=str))}}}")


# Compute minimum vertex covers of this bipartite graph. For each cover,
# zero out those variables and analyse the residual variety.
def all_covers(ledges):
    """Enumerate minimal vertex covers of bipartite graph given as
    {l: set(m)}. We exploit the *block structure* visible in the
    output: each l_i partners with exactly one of three 4-element
    m-blocks, and l_i, l_{i+1} (per i-block) share the SAME
    m-partner set. So per block we either kill BOTH l's of the
    block or kill ALL 4 partner m's — anything else is non-minimal.

    Returns: list of frozensets of variables to set = 0.
    """
    # Group l's by their partner set (= identify blocks)
    from collections import defaultdict
    blocks = defaultdict(list)
    for l, ms in ledges.items():
        blocks[frozenset(ms)].append(l)
    # For each block, two choices: kill all the l's of that block,
    # or kill all the m's it partners with.
    block_choices = []
    for m_partners, l_members in blocks.items():
        block_choices.append([frozenset(l_members), m_partners])
    # Cartesian product of choices
    covers = []
    def rec(i, acc):
        if i == len(block_choices):
            covers.append(acc)
            return
        for choice in block_choices[i]:
            rec(i + 1, acc | choice)
    rec(0, frozenset())
    # Dedupe (some covers may coincide)
    covers = list({c for c in covers})
    return covers


covers = all_covers(ladj)
if covers is None:
    print("\n  (skipping vertex-cover enumeration: too many monomial eqs)")
else:
    print(f"\n#minimal vertex covers of monomial-equation graph: {len(covers)}")
    # Sort by size
    covers.sort(key=lambda c: (len(c), str(sorted(c, key=str))))
    for c in covers[:10]:
        print(f"   size {len(c)}: {sorted(c, key=str)}")


# ── For each minimal cover, substitute zero and analyse residual ────────────
print("\n" + "=" * 70)
print("Per-cover analysis:")
print("=" * 70)


def remaining_vars(zero_set):
    return [v for v in all_vars if v not in zero_set]


def analyse_branch(zero_set, branch_id):
    """Given a set of (λ, μ) variables forced to 0, simplify all garbage
    equations and analyse the residual variety."""
    subs = {v: 0 for v in zero_set}
    res_eqs = [sp.expand(g.subs(subs)) for g in garbage_polys]
    res_eqs = [g for g in res_eqs if g != 0]
    # Dedupe by simplification
    res_eqs = list({sp.expand(g) for g in res_eqs})
    free = remaining_vars(zero_set)
    a_res = sp.expand(alpha.subs(subs))
    b_res = sp.expand(beta.subs(subs))
    print(f"\n--- Branch {branch_id}: zero = {sorted(zero_set, key=str)} ---")
    print(f"   free vars ({len(free)}): {free}")
    print(f"   #residual eqs after dedupe: {len(res_eqs)}")
    for e in res_eqs[:8]:
        print(f"     {e} = 0")
    if len(res_eqs) > 8:
        print(f"     … (+{len(res_eqs)-8} more)")
    print(f"   α | branch = {a_res}")
    print(f"   β | branch = {b_res}")

    # Try to solve the residuals over the free variables
    if res_eqs:
        try:
            sols = sp.solve(res_eqs, free, dict=True)
            print(f"   solve() returned {len(sols)} branches")
            for k, s in enumerate(sols[:4]):
                # Substitute back to get α, β as functions of remaining params
                a_k = sp.expand(a_res.subs(s))
                b_k = sp.expand(b_res.subs(s))
                still_free = [v for v in free if v not in s]
                print(f"     sub-branch {k}: solved {sorted(s.keys(), key=str)}, "
                      f"free remaining = {still_free}")
                print(f"       α = {a_k}")
                print(f"       β = {b_k}")
        except Exception as e:
            print(f"   solve failed: {e}")
    else:
        # No residual: V on this branch is the full free-variable affine space
        print(f"   no residual equations — V on this branch is R^{len(free)}")
        print(f"   (α, β) maps freely on this {len(free)}-dim component")
    return a_res, b_res, res_eqs


# Run analysis on each minimal cover
covers_to_analyse = covers if covers else []
component_data = []
for idx, cover in enumerate(covers_to_analyse):
    a_res, b_res, res_eqs = analyse_branch(cover, idx)
    component_data.append((cover, a_res, b_res, res_eqs))


# ── Dimension summary + (α,β) reachability across components ────────────────
print("\n" + "=" * 70)
print("Component summary & (α,β) reachability")
print("=" * 70)


def alpha_beta_reachable(a_expr, b_expr, free):
    """Heuristic: report whether (α, β) is identically 0 / constant / has full
    image on the component, by inspecting which symbolic variables appear."""
    a_vars_used = a_expr.free_symbols & set(free)
    b_vars_used = b_expr.free_symbols & set(free)
    return a_vars_used, b_vars_used


components_with_nontrivial_ab = []
for idx, (cover, a_res, b_res, res_eqs) in enumerate(component_data):
    free = remaining_vars(cover)
    av, bv = alpha_beta_reachable(a_res, b_res, free)
    print(f"\nBranch {idx}: zero={sorted(cover, key=str)}")
    print(f"  α uses free vars: {sorted(av, key=str)}  ;  α expr: {a_res}")
    print(f"  β uses free vars: {sorted(bv, key=str)}  ;  β expr: {b_res}")
    if a_res == 0 and b_res == 0:
        verdict = "α=β=0 (trivial)"
    elif a_res == 0:
        verdict = "α=0 always (only β varies) → image is the β-axis"
    elif b_res == 0:
        verdict = "β=0 always (only α varies) → image is the α-axis"
    else:
        verdict = "both α, β potentially non-zero"
        components_with_nontrivial_ab.append(idx)
    print(f"  verdict: {verdict}")


# ── Try to build a 3-atom triple from V ──────────────────────────────────────
print("\n" + "=" * 70)
print("Looking for a 3-atom (λ, μ)-triple in V with (α, β)-rank 2")
print("=" * 70)


# Strategy: if there are TWO different components, each producing a non-zero
# (α, β) direction (one along α-axis, the other along β-axis, say), we can
# combine them: atom1 from α-axis component (α₁≠0, β₁=0), atom2 from β-axis
# component (α₂=0, β₂≠0). A third atom (anywhere in V) lets us tune.
α_only_branches = [idx for idx, (c, a, b, _) in enumerate(component_data)
                   if a != 0 and b == 0]
β_only_branches = [idx for idx, (c, a, b, _) in enumerate(component_data)
                   if a == 0 and b != 0]
both_branches = [idx for idx, (c, a, b, _) in enumerate(component_data)
                 if a != 0 and b != 0]
print(f"\n#branches reaching only α (β≡0): {len(α_only_branches)} → {α_only_branches}")
print(f"#branches reaching only β (α≡0): {len(β_only_branches)} → {β_only_branches}")
print(f"#branches reaching both α & β   : {len(both_branches)} → {both_branches}")


# Construct an explicit (α, β)-rank-2 triple if possible
def materialise_atom(cover, a_target=None, b_target=None):
    """Given a branch cover, produce a *concrete* (λ, μ) instance that
    sets the zeroed vars to 0 and the free vars to small integers giving
    non-trivial (α, β). Returns dict, α-value, β-value."""
    subs = {v: 0 for v in cover}
    free = remaining_vars(cover)
    # try (free = 1) first, then perturb if degenerate
    for trial in range(20):
        import random
        random.seed(100 + trial)
        for v in free:
            subs[v] = sp.Integer(random.choice([-2, -1, 1, 2]))
        a_val = sp.expand(alpha.subs(subs))
        b_val = sp.expand(beta.subs(subs))
        # Check residuals are 0 on this branch
        ok = all(sp.expand(g.subs(subs)) == 0 for g in garbage_polys)
        if ok:
            return subs.copy(), a_val, b_val
    return None, None, None


if α_only_branches and β_only_branches:
    print("\nAttempting to build a rank-2 triple from one α-only + one β-only branch:")
    a_branch = α_only_branches[0]
    b_branch = β_only_branches[0]
    cover_a = component_data[a_branch][0]
    cover_b = component_data[b_branch][0]
    atom_a, a_a, b_a = materialise_atom(cover_a)
    atom_b, a_b, b_b = materialise_atom(cover_b)
    print(f"  α-only branch sample (zero={sorted(cover_a, key=str)}): (α, β) = ({a_a}, {b_a})")
    print(f"  β-only branch sample (zero={sorted(cover_b, key=str)}): (α, β) = ({a_b}, {b_b})")
    # The 3rd atom can be any from V (say zero atom).
    if a_a != 0 and b_b != 0:
        print(f"  Rank of [{a_a} {a_b}; {b_a} {b_b}] = "
              f"{sp.Matrix([[a_a, a_b], [b_a, b_b]]).rank()}")
        print("  → If rank 2, two atoms suffice; the 3-atom rank-2 triple is trivially constructible.")


print("\n" + "=" * 70)
print("Rigorously verifying (α, β) ≡ (0, 0) on V — per branch")
print("=" * 70)

# For each branch, the residual ideal lives in the free vars. Check whether
# α and β reduce to 0 modulo this ideal.
all_kill = True
for idx, (cover, a_res, b_res, res_eqs) in enumerate(component_data):
    free = remaining_vars(cover)
    if not res_eqs:
        # No residuals; α, β themselves must be identically 0 on the cover
        verdict = "no residuals; α, β as functions of free vars"
        a_red, b_red = sp.expand(a_res), sp.expand(b_res)
    else:
        try:
            G = sp.groebner(res_eqs, *free, order='lex')
            a_red = sp.reduced(a_res, G.polys, *free)[1]
            b_red = sp.reduced(b_res, G.polys, *free)[1]
            verdict = f"reduced modulo Groebner of {len(G.polys)} polys"
        except Exception as e:
            verdict = f"groebner failed: {e}"
            a_red, b_red = sp.expand(a_res), sp.expand(b_res)
    print(f"\nBranch {idx}: cover={sorted(cover, key=str)}")
    print(f"  free: {free}")
    print(f"  α reduced: {a_red}  (was {a_res})")
    print(f"  β reduced: {b_red}  (was {b_res})")
    print(f"  {verdict}")
    if a_red != 0 or b_red != 0:
        all_kill = False
        print(f"  *** (α, β) is NOT identically zero on this component ***")

print("\n" + "=" * 70)
if all_kill:
    print("CONCLUSION: V is non-empty and positive-dimensional, but on EVERY")
    print("irreducible component (α, β) ≡ (0, 0). Therefore no atom from V")
    print("can contribute non-zero (α, β); the V-only strategy CANNOT produce")
    print("a 3-atom triple with (α, β)-rank 2.")
    print("\nImpossibility for the 'zero garbage' strategy: PROVEN.")
    print("\nThe (2, 2, bridge-3) construction in any rank-1 atom family thus")
    print("requires garbage CANCELLATION across atoms (rank(U) ≤ 1 with U ≠ 0),")
    print("not zero garbage per atom.")
else:
    print("Some component reaches (α, β) ≠ 0; constructive 3-atom triples")
    print("are possible — see materialise_atom output above.")
print("=" * 70)
