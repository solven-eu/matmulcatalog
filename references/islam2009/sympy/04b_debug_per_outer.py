"""
Phase A step 4 debug — verify Lemma 1 expansion of EACH outer T-term
individually, then verify the SUM matches T(Ã, B̃, C̃).

If each per-outer-T expansion checks out but the sum doesn't, the error
is in the R-cancellation (Lemma 2). If a per-outer-T expansion is off,
the error is in argument mapping for that T.
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


def symblock(prefix, n, i, j):
    return sp.Matrix(n, n, lambda r, c: sp.Symbol(f"{prefix}_{i}{j}_{r+1}{c+1}"))


def T_disjoint(A, B, C, U, V, W, X, Y, Z, m):
    """T(A,B,C,U,V,W,X,Y,Z) = T(A,B,C) + T(U,V,W) + T(X,Y,Z)."""
    out = sp.Integer(0)
    for i in range(1, m + 1):
        for j in range(1, m + 1):
            for k in range(1, m + 1):
                out += A[i-1, j-1] * B[j-1, k-1] * C[k-1, i-1]
                out += U[j-1, k-1] * V[k-1, i-1] * W[i-1, j-1]
                out += X[k-1, i-1] * Y[i-1, j-1] * Z[j-1, k-1]
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


def R_op(A, V, Z, m):
    out = sp.Integer(0)
    for i in range(1, m + 1):
        for j in range(1, m + 1):
            for k in range(1, m + 1):
                out += A[i-1, j-1] * V[k-1, i-1] * Z[j-1, k-1]
    return out


def U_op(A, Y, C, W, Z, m):
    out = sp.Integer(0)
    for i in range(1, m + 1):
        for j in range(1, m + 1):
            for k in range(1, m + 1):
                out += A[i-1, j-1] * Y[i-1, j-1] * (
                    C[k-1, i-1] + W[i-1, j-1] + Z[j-1, k-1])
    return out


def verify_per_outer():
    """Verify Lemma 1 expansion of each of the 4 outer T-terms."""
    n = 2
    m = n + 1
    L, R = build_padding(n)
    A = {(i, j): symblock("A", n, i, j) for i in (1, 2) for j in (1, 2)}
    B = {(i, j): symblock("B", n, i, j) for i in (1, 2) for j in (1, 2)}
    At = {ij: L * A[ij] * R       for ij in A}
    Bt = {ij: L * B[ij] * L.T     for ij in B}
    Ct = {(i, j): sp.Matrix(m, m, lambda r, c, _i=i, _j=j:
                            sp.Symbol(f"ct_{_i}{_j}_{r+1}{c+1}"))
          for i in (1, 2) for j in (1, 2)}

    # 4 outer T-terms (from Eq 5.2 — page 44, after the (1/3) regrouping).
    # Each is T(A_arg, B_arg, C_arg, U_arg, V_arg, W_arg, X_arg, Y_arg, Z_arg).
    A11, A12, A21, A22 = At[(1, 1)], At[(1, 2)], At[(2, 1)], At[(2, 2)]
    B11, B12, B21, B22 = Bt[(1, 1)], Bt[(1, 2)], Bt[(2, 1)], Bt[(2, 2)]
    C11, C12, C21, C22 = Ct[(1, 1)], Ct[(1, 2)], Ct[(2, 1)], Ct[(2, 2)]

    # Each tuple: (label, (A,B,C,U,V,W,X,Y,Z), coef).
    outer = [
        ("diag1",   (A11, B11, C11,  A11, B11, C11,  A11, B11, C11),   sp.Rational(1, 3)),
        ("off1",    (A12, B21, C11, -A11, B12, -C12, A21, B11, C21),   sp.Integer(1)),
        ("off2",    (A12, B22, C12,  A21, B12, C22, -A22, B21, -C21),  sp.Integer(1)),
        ("diag2",   (A22, B22, C22,  A22, B22, C22,  A22, B22, C22),   sp.Rational(1, 3)),
    ]

    for label, args, coef in outer:
        lhs = sp.expand(coef * T_disjoint(*args, m))
        # Lemma 1: T(A,B,C,U,V,W,X,Y,Z) = S - R(A,V,Z) - R(U,Y,C) - R(X,B,W)
        #                                     - U(A,Y,C,W,Z) - U(U,B,W,Z,C) - U(X,V,Z,C,W)
        A_, B_, C_, U_, V_, W_, X_, Y_, Z_ = args
        rhs = coef * (
            S_full(*args, m)
            - R_op(A_, V_, Z_, m) - R_op(U_, Y_, C_, m) - R_op(X_, B_, W_, m)
            - U_op(A_, Y_, C_, W_, Z_, m)
            - U_op(U_, B_, W_, Z_, C_, m)
            - U_op(X_, V_, Z_, C_, W_, m)
        )
        rhs = sp.expand(rhs)
        diff = sp.expand(lhs - rhs)
        print(f"  {label:6s}: {'OK' if diff == 0 else f'FAIL ({len(diff.as_ordered_terms())} terms)'}")

    # Now: sum of 4 outer T-terms (= T(Ã, B̃, C̃) by Eq 5.2).
    print()
    full_sum_outers = sp.Integer(0)
    for _, args, coef in outer:
        full_sum_outers += coef * T_disjoint(*args, m)
    full_sum_outers = sp.expand(full_sum_outers)
    # Direct T(Ã, B̃, C̃).
    At_full = sp.Matrix.vstack(
        sp.Matrix.hstack(A11, A12), sp.Matrix.hstack(A21, A22))
    Bt_full = sp.Matrix.vstack(
        sp.Matrix.hstack(B11, B12), sp.Matrix.hstack(B21, B22))
    Ct_full = sp.Matrix.vstack(
        sp.Matrix.hstack(C11, C12), sp.Matrix.hstack(C21, C22))
    direct = sp.Integer(0)
    for i in range(2 * m):
        for j in range(2 * m):
            for k in range(2 * m):
                direct += At_full[i, j] * Bt_full[j, k] * Ct_full[k, i]
    direct = sp.expand(direct)
    diff = sp.expand(direct - full_sum_outers)
    print(f"  Eq 5.2 (sum of 4 outer T = T(Ã,B̃,C̃)): "
          f"{'OK' if diff == 0 else f'FAIL ({len(diff.as_ordered_terms())} terms)'}")


if __name__ == "__main__":
    print("Phase A step 4 debug — per-outer Lemma 1 + Eq 5.2 check")
    print("=" * 70)
    verify_per_outer()
