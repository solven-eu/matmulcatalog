"""
Islam 2009 §5.3 — Phase A, step 2: verify Lemma 1.

Lemma 1 is the algebraic heart of Pan's trilinear aggregation: it shows
that the polynomial T(A, B, C, U, V, W, X, Y, Z) that captures THREE
disjoint matrix products A·B, U·V, X·Y can be expressed in terms of a
single "aggregated" sum S plus a small number of "correction" sums
(R-terms and U-terms). Concretely:

    T(A,B,C,U,V,W,X,Y,Z) = S(A,B,C,U,V,W,X,Y,Z)
                          − R(A,V,Z) − R(U,Y,C) − R(X,B,W)
                          − U(A,Y,C,W,Z) − U(U,B,W,Z,C) − U(X,V,Z,C,W)

with definitions (sums over 1 ≤ i,j,k ≤ m):

    T(A,B,C) = Σ a_{i,j} b_{j,k} c_{k,i}              (polynomial form of A·B)
    S_{i,j,k}     = (a_{i,j}+u_{j,k}+x_{k,i})
                    · (b_{j,k}+v_{k,i}+y_{i,j})
                    · (c_{k,i}+w_{i,j}+z_{j,k})
    S      = Σ_{i,j,k} S_{i,j,k}
    R(A,V,Z)      = Σ a_{i,j} v_{k,i} z_{j,k}
    U(A,Y,C,W,Z)  = Σ a_{i,j} y_{i,j} (c_{k,i} + w_{i,j} + z_{j,k})

The identity holds <strong>iff all row sums and all column sums of A,
B, U, V, X, Y are zero</strong> (proof in §5.3 — many cross-terms vanish
under that condition). This is exactly why the §5.4 padding scheme,
which forces those sums to zero, makes Lemma 1 applicable.

This script verifies the identity SYMBOLICALLY by computing both sides
in sympy for a small m (m = 3) and checking they expand to the same
polynomial in the entries of C, W, Z.

Run: python3 02_lemma1_disjoint_products.py
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


def symbolic_matrix(prefix, m):
    """m × m symbolic matrix with entries prefix_<row><col>."""
    return sp.Matrix(m, m, lambda r, c: sp.Symbol(f"{prefix}_{r+1}{c+1}"))


def padded_zero_sum(prefix, n):
    """
    Build an (n+1)×(n+1) symbolic matrix with all row + col sums zero, by
    padding a generic n × n matrix via L · X · R (size grows by 1 in both
    dimensions). The result is fully symbolic in the n² generic entries.
    """
    X = symbolic_matrix(prefix, n)
    L, R = build_padding(n)
    return L * X * R


def T_polynomial(A, B, C):
    """Polynomial representation of A·B with C symbolic: Σ a_{i,j} b_{j,k} c_{k,i}."""
    m = A.rows
    assert A.cols == m and B.rows == m and B.cols == m and C.rows == m and C.cols == m
    out = sp.Integer(0)
    for i in range(1, m + 1):
        for j in range(1, m + 1):
            for k in range(1, m + 1):
                out += A[i-1, j-1] * B[j-1, k-1] * C[k-1, i-1]
    return out


def S_polynomial(A, B, C, U, V, W, X, Y, Z):
    m = A.rows
    out = sp.Integer(0)
    for i in range(1, m + 1):
        for j in range(1, m + 1):
            for k in range(1, m + 1):
                f1 = A[i-1, j-1] + U[j-1, k-1] + X[k-1, i-1]
                f2 = B[j-1, k-1] + V[k-1, i-1] + Y[i-1, j-1]
                f3 = C[k-1, i-1] + W[i-1, j-1] + Z[j-1, k-1]
                out += f1 * f2 * f3
    return out


def R_polynomial(A, V, Z):
    m = A.rows
    out = sp.Integer(0)
    for i in range(1, m + 1):
        for j in range(1, m + 1):
            for k in range(1, m + 1):
                out += A[i-1, j-1] * V[k-1, i-1] * Z[j-1, k-1]
    return out


def U_polynomial(A, Y, C, W, Z):
    m = A.rows
    out = sp.Integer(0)
    for i in range(1, m + 1):
        for j in range(1, m + 1):
            for k in range(1, m + 1):
                out += A[i-1, j-1] * Y[i-1, j-1] * (
                    C[k-1, i-1] + W[i-1, j-1] + Z[j-1, k-1])
    return out


def test(n):
    """Test with input matrices padded from n×n (so the m×m matrices used in
    the identity have m = n + 1)."""
    m = n + 1
    print(f"\n── n = {n}   (m = {m})  generic blocks padded to zero row/col sums ──")

    # A, B, U, V, X, Y: zero-sum padded.
    A = padded_zero_sum("a", n)
    B = padded_zero_sum("b", n)
    U = padded_zero_sum("u", n)
    V = padded_zero_sum("v", n)
    X = padded_zero_sum("x", n)
    Y = padded_zero_sum("y", n)
    # C, W, Z: symbolic m×m matrices (output variables — NOT padded).
    C = symbolic_matrix("c", m)
    W = symbolic_matrix("w", m)
    Z = symbolic_matrix("z", m)

    print("   computing T(A,B,C) + T(U,V,W) + T(X,Y,Z) ...")
    T_total = T_polynomial(A, B, C) + T_polynomial(U, V, W) + T_polynomial(X, Y, Z)
    T_total = sp.expand(T_total)

    print("   computing S(...) − R(A,V,Z) − R(U,Y,C) − R(X,B,W) − U(A,Y,C,W,Z) − U(U,B,W,Z,C) − U(X,V,Z,C,W) ...")
    rhs = (S_polynomial(A, B, C, U, V, W, X, Y, Z)
           - R_polynomial(A, V, Z) - R_polynomial(U, Y, C) - R_polynomial(X, B, W)
           - U_polynomial(A, Y, C, W, Z) - U_polynomial(U, B, W, Z, C) - U_polynomial(X, V, Z, C, W))
    rhs = sp.expand(rhs)

    diff = sp.expand(T_total - rhs)
    if diff == 0:
        print("   OK — Lemma 1 holds for m =", m)
        return True
    else:
        print("   FAIL — residual is non-zero. First terms of T_total − rhs:")
        terms = diff.as_ordered_terms()[:5]
        for t in terms:
            print("     ", t)
        return False


if __name__ == "__main__":
    print("Islam 2009 §5.3: Lemma 1 verification")
    print("=" * 70)
    ok2 = test(2)   # m = 3
    print("\n" + "=" * 70)
    print(f"   m=3: {'OK' if ok2 else 'FAIL'}")
    if ok2:
        print("\n[Phase A step 2] PASS — Lemma 1 holds under zero row/col-sum input.")
    else:
        print("\n[Phase A step 2] FAIL — re-examine S / R / U formulas in §5.3.")
