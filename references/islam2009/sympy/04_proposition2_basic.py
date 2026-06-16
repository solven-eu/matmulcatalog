"""
Islam 2009 §5.4 — Phase A, step 4: verify Eq 5.3 (page 46).

After applying Lemma 1 to each of the 4 "T-products" of the 2×2 block
schedule (and Lemma 2 to eliminate the R-correction terms because they
cancel pairwise), the joint product T(Ã, B̃, C̃) — equivalently the
output C̃ = Ã·B̃ in polynomial form — can be expressed as

    T(Ã, B̃, C̃) = (1/3) S(Ã^{1,1}, B̃^{1,1}, C̃^{1,1}, Ã^{1,1}, B̃^{1,1}, C̃^{1,1}, Ã^{1,1}, B̃^{1,1}, C̃^{1,1})
                  + S(Ã^{1,2}, B̃^{2,1}, C̃^{1,1}, -Ã^{1,1}, B̃^{1,2}, -C̃^{1,2}, Ã^{2,1}, B̃^{1,1}, C̃^{2,1})
                  + S(Ã^{1,2}, B̃^{2,2}, C̃^{1,2}, Ã^{2,1}, B̃^{1,2}, C̃^{2,2}, -Ã^{2,2}, B̃^{2,1}, -C̃^{2,1})
                  + (1/3) S(Ã^{2,2}, B̃^{2,2}, C̃^{2,2}, Ã^{2,2}, B̃^{2,2}, C̃^{2,2}, Ã^{2,2}, B̃^{2,2}, C̃^{2,2})
                  − Σ_{8 triples} U(·, ·, ·, ·, ·)

where the 8 U-correction terms are listed in Eq 5.3. This formula is the
heart of the TA construction; once it holds, the improved version (§5.5,
§5.6) just refines the count via Lemmas 3, 4 + index-m savings.

This step verifies Eq 5.3 SYMBOLICALLY for the smallest non-trivial case
n_b = 2 (block size 2, full input matrices 4×4, m = n_b + 1 = 3).

Run: python3 04_proposition2_basic.py
"""
import sympy as sp


def build_padding(n):
    I_n = sp.eye(n)
    u = sp.ones(n, 1)
    L = sp.Matrix.vstack(I_n, -u.T)
    Rblock = I_n - sp.Rational(1, n + 1) * u * u.T
    Rcol = -sp.Rational(1, n + 1) * u
    R = sp.Matrix.hstack(Rblock, Rcol)
    return L, R


def symbolic_block(prefix, n, i, j):
    return sp.Matrix(n, n, lambda r, c:
                     sp.Symbol(f"{prefix}_{i}{j}_{r+1}{c+1}"))


def symbolic_matrix(prefix, m):
    return sp.Matrix(m, m, lambda r, c: sp.Symbol(f"{prefix}_{r+1}{c+1}"))


def T_polynomial(A, B, C):
    """Σ_{i,j,k} a_{i,j} b_{j,k} c_{k,i}."""
    m = A.rows
    out = sp.Integer(0)
    for i in range(1, m + 1):
        for j in range(1, m + 1):
            for k in range(1, m + 1):
                out += A[i-1, j-1] * B[j-1, k-1] * C[k-1, i-1]
    return out


def S_full(A, B, C, U, V, W, X, Y, Z, m):
    out = sp.Integer(0)
    for i in range(1, m + 1):
        for j in range(1, m + 1):
            for k in range(1, m + 1):
                f1 = A[i-1, j-1] + U[j-1, k-1] + X[k-1, i-1]
                f2 = B[j-1, k-1] + V[k-1, i-1] + Y[i-1, j-1]
                f3 = C[k-1, i-1] + W[i-1, j-1] + Z[j-1, k-1]
                out += f1 * f2 * f3
    return out


def U_term(A, Y, C, W, Z):
    """U(A, Y, C, W, Z) = Σ_{i,j,k} a_{i,j} y_{i,j} (c_{k,i} + w_{i,j} + z_{j,k})."""
    m = A.rows
    out = sp.Integer(0)
    for i in range(1, m + 1):
        for j in range(1, m + 1):
            for k in range(1, m + 1):
                out += A[i-1, j-1] * Y[i-1, j-1] * (
                    C[k-1, i-1] + W[i-1, j-1] + Z[j-1, k-1])
    return out


def test(n_b):
    m = n_b + 1
    print(f"\n── n_b = {n_b}  (full size {2*n_b}×{2*n_b},  m = n_b + 1 = {m}) ──")

    L, R = build_padding(n_b)
    A_blocks = {(i, j): symbolic_block("A", n_b, i, j) for i in (1, 2) for j in (1, 2)}
    B_blocks = {(i, j): symbolic_block("B", n_b, i, j) for i in (1, 2) for j in (1, 2)}

    At = {(i, j): L * A_blocks[(i, j)] * R       for (i, j) in A_blocks}
    Bt = {(i, j): L * B_blocks[(i, j)] * L.T     for (i, j) in B_blocks}

    # Output blocks C̃^{i,j}: symbolic m × m matrices (placeholder for output
    # variables in the polynomial T(Ã, B̃, C̃)).
    Ct = {(i, j): sp.Matrix(m, m, lambda r, c, _i=i, _j=j:
                            sp.Symbol(f"ct_{_i}{_j}_{r+1}{c+1}"))
          for i in (1, 2) for j in (1, 2)}

    # Build Ã, B̃, C̃ as 2m × 2m matrices.
    At_full = sp.Matrix.vstack(
        sp.Matrix.hstack(At[(1, 1)], At[(1, 2)]),
        sp.Matrix.hstack(At[(2, 1)], At[(2, 2)]),
    )
    Bt_full = sp.Matrix.vstack(
        sp.Matrix.hstack(Bt[(1, 1)], Bt[(1, 2)]),
        sp.Matrix.hstack(Bt[(2, 1)], Bt[(2, 2)]),
    )
    Ct_full = sp.Matrix.vstack(
        sp.Matrix.hstack(Ct[(1, 1)], Ct[(1, 2)]),
        sp.Matrix.hstack(Ct[(2, 1)], Ct[(2, 2)]),
    )

    print("   computing LHS = T(Ã, B̃, C̃) directly ...")
    lhs = sp.expand(T_polynomial(At_full, Bt_full, Ct_full))

    print("   computing RHS via Eq 5.3 (4 S-terms + 8 U-terms) ...")
    # 4 S-terms (page 46 Eq 5.3 — note the 1/3 factor on diagonal terms,
    # and the sign changes in the disjoint triples).
    A11, A12, A21, A22 = At[(1, 1)], At[(1, 2)], At[(2, 1)], At[(2, 2)]
    B11, B12, B21, B22 = Bt[(1, 1)], Bt[(1, 2)], Bt[(2, 1)], Bt[(2, 2)]
    C11, C12, C21, C22 = Ct[(1, 1)], Ct[(1, 2)], Ct[(2, 1)], Ct[(2, 2)]

    s_diag1 = sp.Rational(1, 3) * S_full(A11, B11, C11, A11, B11, C11, A11, B11, C11, m)
    s_off1  = S_full(A12, B21, C11, -A11, B12, -C12, A21, B11, C21, m)
    s_off2  = S_full(A12, B22, C12, A21, B12, C22, -A22, B21, -C21, m)
    s_diag2 = sp.Rational(1, 3) * S_full(A22, B22, C22, A22, B22, C22, A22, B22, C22, m)

    # 8 U-terms — page 46 / page 50 Eq 5.3.
    u_terms = (
        + U_term(A11, B11, C11, C11, C11)
        + U_term(A12, B11, C11, -C12, C21)
        + U_term(-A11, B21, -C12, C21, C11)
        + U_term(A21, B12, C21, C11, -C12)
        + U_term(A12, B21, C12, C22, -C21)
        + U_term(A21, B22, C22, -C21, C12)
        + U_term(-A22, B12, -C21, C12, C22)
        + U_term(A22, B22, C22, C22, C22)
    )

    rhs = sp.expand(s_diag1 + s_off1 + s_off2 + s_diag2 - u_terms)

    diff = sp.expand(lhs - rhs)
    if diff == 0:
        print(f"   OK — Eq 5.3 holds for n_b = {n_b}, m = {m}")
        return True
    print(f"   FAIL — residual has {sum(1 for _ in diff.as_ordered_terms())} nonzero terms")
    for t in diff.as_ordered_terms()[:5]:
        print(f"     {t}")
    return False


if __name__ == "__main__":
    print("Islam 2009 §5.4: Eq 5.3 (basic Proposition 2)")
    print("=" * 70)
    ok = test(n_b=2)
    print("\n" + "=" * 70)
    if ok:
        print("[Phase A step 4] PASS — basic block-product formula (Eq 5.3) holds.")
        print()
        print("Note: this is the BASIC version. The improvements in §5.5 (drop")
        print("index-m terms in s_0, s_1, s_2) + §5.6 (odd case zero-padding)")
        print("reduce the mult count to Islam Prop 1's (n³+12n²+11n)/3 / (n³+15n²+14n-6)/3.")
        print("Implementing those count-saving optimizations is Phase B.")
    else:
        print("[Phase A step 4] FAIL — Eq 5.3 doesn't match T(Ã, B̃, C̃).")
