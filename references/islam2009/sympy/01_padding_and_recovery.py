"""
Islam 2009 §5.4 — Phase A, step 1: verify zero-sum padding + recovery.

The first non-trivial claim of the TA construction (Section 5.4) is:

    Given A, B of size 2n with arbitrary entries, subdivide into 2×2 blocks
    of size n, pad each block to size (n+1) via L = [I; -u] (left) and
    R = [I - (1/(n+1))u^T u, -(1/(n+1))u^T] (right), so that

        Ã^{i,j} = L · A^{i,j} · R          (size (n+1)×(n+1))
        B̃^{i,j} = L · B^{i,j} · L^T       (size (n+1)×(n+1))

    Then ALL rows and ALL columns of Ã^{i,j} and B̃^{i,j} sum to zero, and
    the matrix product C = AB can be recovered from C̃ = ÃB̃ by RESIZING
    each C̃^{i,j} block from size (n+1) back to size n.

Section 5.4's exact recovery rule is:

    "C^{i,j} is obtained by resizing C̃^{i,j} in size (n,n)"

This script verifies BOTH claims symbolically for small n (n = 2, n = 3),
which is what step 1 of Phase A asks for. Successful run prints "OK".

Run: python3 01_padding_and_recovery.py
"""
import sympy as sp


def build_padding(n):
    """Build L (size (n+1)×n) and R (size n×(n+1)) for zero-sum padding."""
    I_n = sp.eye(n)
    u = sp.ones(n, 1)         # column vector of n ones
    # L = [I_n ; -u^T]          (n+1) × n
    L = sp.Matrix.vstack(I_n, -u.T)
    # R = [I_n - (1/(n+1)) u u^T  |  -(1/(n+1)) u]      n × (n+1)
    Rblock = I_n - sp.Rational(1, n + 1) * u * u.T
    Rcol = -sp.Rational(1, n + 1) * u
    R = sp.Matrix.hstack(Rblock, Rcol)
    return L, R


def has_zero_row_col_sums(M, name):
    """Return None if M has all row sums and all column sums == 0; else error description."""
    row_sums = [sp.simplify(sum(M.row(i))) for i in range(M.rows)]
    col_sums = [sp.simplify(sum(M.col(j))) for j in range(M.cols)]
    bad_rows = [i for i, s in enumerate(row_sums) if s != 0]
    bad_cols = [j for j, s in enumerate(col_sums) if s != 0]
    if bad_rows or bad_cols:
        return f"{name}: bad rows {bad_rows}, bad cols {bad_cols}"
    return None


def symbolic_block(prefix, n, i, j):
    """Build a symbolic n×n matrix called e.g. 'A_12' with entries prefix_ij_<row><col>."""
    return sp.Matrix(n, n, lambda r, c:
                     sp.Symbol(f"{prefix}_{i}{j}_{r+1}{c+1}"))


def test(n):
    print(f"\n── n = {n}   (block size, full matrices are {2*n} × {2*n}) ──")
    L, R = build_padding(n)
    print(f"   L shape = {L.shape},  R shape = {R.shape}")

    # The 4 blocks of A and B (symbolic, fully generic).
    A_blocks = {(i, j): symbolic_block("A", n, i, j) for i in (1, 2) for j in (1, 2)}
    B_blocks = {(i, j): symbolic_block("B", n, i, j) for i in (1, 2) for j in (1, 2)}

    # Tilde blocks.
    At = {(i, j): L * A_blocks[(i, j)] * R       for (i, j) in A_blocks}
    Bt = {(i, j): L * B_blocks[(i, j)] * L.T     for (i, j) in B_blocks}

    # CLAIM 1: zero row / column sums.
    print("\n   [claim 1]  Ã and B̃ have zero row + column sums:")
    named = {}
    for (i, j), M in At.items():
        named[f"Ã^{i}{j}"] = M
    for (i, j), M in Bt.items():
        named[f"B̃^{i}{j}"] = M
    for name, M in named.items():
        err = has_zero_row_col_sums(M, name)
        print(f"      {name}: {'OK' if err is None else err}")
        if err is not None:
            return False

    # CLAIM 2: C = AB is recoverable from C̃ = ÃB̃ by taking the top-left
    # n × n sub-block of each C̃^{i,j}.
    print("\n   [claim 2]  Recovery from C̃ by resizing to size n:")

    # Direct C = AB via 2×2 block product (the truth).
    A_full = sp.Matrix.vstack(
        sp.Matrix.hstack(A_blocks[(1, 1)], A_blocks[(1, 2)]),
        sp.Matrix.hstack(A_blocks[(2, 1)], A_blocks[(2, 2)]),
    )
    B_full = sp.Matrix.vstack(
        sp.Matrix.hstack(B_blocks[(1, 1)], B_blocks[(1, 2)]),
        sp.Matrix.hstack(B_blocks[(2, 1)], B_blocks[(2, 2)]),
    )
    C_direct = sp.expand(A_full * B_full)

    # Ã × B̃ via block product, then resize each C̃^{i,j} (top-left n×n).
    # C̃^{1,1} = Ã^{1,1} B̃^{1,1} + Ã^{1,2} B̃^{2,1}
    # C̃^{1,2} = Ã^{1,1} B̃^{1,2} + Ã^{1,2} B̃^{2,2}
    # C̃^{2,1} = Ã^{2,1} B̃^{1,1} + Ã^{2,2} B̃^{2,1}
    # C̃^{2,2} = Ã^{2,1} B̃^{1,2} + Ã^{2,2} B̃^{2,2}
    Ct = {
        (1, 1): At[(1, 1)] * Bt[(1, 1)] + At[(1, 2)] * Bt[(2, 1)],
        (1, 2): At[(1, 1)] * Bt[(1, 2)] + At[(1, 2)] * Bt[(2, 2)],
        (2, 1): At[(2, 1)] * Bt[(1, 1)] + At[(2, 2)] * Bt[(2, 1)],
        (2, 2): At[(2, 1)] * Bt[(1, 2)] + At[(2, 2)] * Bt[(2, 2)],
    }
    # Resize each C̃^{i,j} to size n (top-left).
    C_recovered_blocks = {(i, j): Ct[(i, j)][:n, :n] for (i, j) in Ct}
    C_recovered = sp.Matrix.vstack(
        sp.Matrix.hstack(C_recovered_blocks[(1, 1)], C_recovered_blocks[(1, 2)]),
        sp.Matrix.hstack(C_recovered_blocks[(2, 1)], C_recovered_blocks[(2, 2)]),
    )
    C_recovered = sp.expand(C_recovered)

    diff = sp.simplify(C_direct - C_recovered)
    nonzero = sum(1 for c in diff if c != 0)
    if nonzero == 0:
        print("      OK — C_recovered == AB (all 4n² entries match)")
        return True
    print(f"      FAIL — {nonzero} entries differ between C_recovered and AB")
    # Show the first few discrepancies for debugging.
    for r in range(diff.rows):
        for c in range(diff.cols):
            if diff[r, c] != 0:
                print(f"         C[{r},{c}]_direct − C[{r},{c}]_recovered = {diff[r, c]}")
                break
    return False


if __name__ == "__main__":
    print("Islam 2009 §5.4: padding + recovery verification")
    print("=" * 70)
    ok2 = test(2)
    ok3 = test(3)
    print("\n" + "=" * 70)
    print(f"   n=2: {'OK' if ok2 else 'FAIL'}")
    print(f"   n=3: {'OK' if ok3 else 'FAIL'}")
    if ok2 and ok3:
        print("\n[Phase A step 1] PASS — padding scheme is correct.")
    else:
        print("\n[Phase A step 1] FAIL — re-examine padding formulas in §5.4.")
