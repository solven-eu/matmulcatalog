"""
Phase A step 4 (re-implementation) — port of Islam's Magma `TA.mgm`
`ProductEven` function from
[Wayback snapshot of csd.uwo.ca/~mislam63/TA.mgm](https://web.archive.org/web/20120223044300/http://www.csd.uwo.ca:80/~mislam63/TA.mgm)
(see `references/typos.md` for provenance).

This re-derives DIS09 Appendix Lemma 4 directly from Islam's reference
implementation rather than from a hand-transcription of the paper, so it
avoids the published typo we documented in `references/typos.md`
(DIS09's 7th `u`-correction has `−C̃^{1,2}` in slot 3 where the Magma
source — and Islam's thesis — have `−C̃^{2,1}`).

What this script verifies:

1. **Correctness** — apply the construction to two symbolic input
   matrices `A, B` of size `2n × 2n` (full input), compute the
   reconstructed `C` block-by-block, and check that it equals the
   direct matrix product `A·B`. Symbolic equality, no numerics.

2. **Multiplication count** — count the rank-1 atoms (each
   `(linear_form_A) × (linear_form_B) × (output_c_variable)` triple
   that the construction emits) and check that the total matches
   Islam Prop 1: `(n³ + 12n² + 11n) / 3` for `n = 2n_b` even.

Test cases: `n_b = 2` (full input 4×4, target cost 100) and
`n_b = 3` (full input 6×6, target cost 238).

Run: `python3 04c_lemma4_from_magma.py`
"""
import sympy as sp


# ──────────────────────────────────────────────────────────────────────
# Multiplication counter (the Magma `Mul` function)
# ──────────────────────────────────────────────────────────────────────

class MulCounter:
    """Counts the rank-1 multiplications emitted by the construction.

    A "multiplication" is one `(linear_form_in_A) × (linear_form_in_B)`
    product, accumulated into an arbitrary linear combination of output
    `c` variables. The Magma's `Mul(a, b, c)` increments a counter
    whenever `a·b·c` is non-zero.
    """
    def __init__(self):
        self.count = 0

    def mul(self, a, b, c):
        if a == 0 or b == 0 or c == 0:
            return sp.Integer(0)
        self.count += 1
        return a * b * c


# ──────────────────────────────────────────────────────────────────────
# Padding matrices L, R  (zero-row/col sum padding from §5.4)
# ──────────────────────────────────────────────────────────────────────

def build_L(m):
    """L = [I_m ; -1^T] of size (m+1) × m."""
    I_m = sp.eye(m)
    u_row = sp.ones(1, m)   # 1 × m row of ones
    return sp.Matrix.vstack(I_m, -u_row)


def build_R(m):
    """R = [I - (1/(m+1)) 1·1^T  | -(1/(m+1)) 1]  of size m × (m+1)."""
    I_m = sp.eye(m)
    u_col = sp.ones(m, 1)
    u_row = sp.ones(1, m)
    block = I_m - sp.Rational(1, m + 1) * u_col * u_row
    last_col = -sp.Rational(1, m + 1) * u_col
    return sp.Matrix.hstack(block, last_col)


def enlarge_left(M):
    """Ã^{i,j} = L · A^{i,j} · R."""
    m = M.rows
    n = M.cols
    return build_L(m) * M * build_R(n)


def enlarge_right(M):
    """B̃^{i,j} = L · B^{i,j} · L^T."""
    n = M.rows
    p = M.cols
    return build_L(n) * M * build_L(p).T


# ──────────────────────────────────────────────────────────────────────
# Atomic operators s_0, s_1, s_2, u_1, u_2, u_2', u_3, u_4
# (faithful ports of `TA.mgm` lines 116-192)
# ──────────────────────────────────────────────────────────────────────

def prop_s1(i, j, k, m):
    """(i, j, k) ∈ S_1 (uniform-coincidence index set)."""
    b1 = (0 <= i) and (i <= j) and (j < k)
    b2 = (0 <= k) and (k < j) and (j <= i)
    return b1 or b2


def prop_size_m(i, j, k, m):
    """At most ONE of {i, j, k} equals m — the §5.5 improvement."""
    return sum(1 for x in (i, j, k) if x == m) < 2


def prop_S1_improved(i, j, k, m):
    return prop_s1(i, j, k, m) and prop_size_m(i, j, k, m)


def prop_S2_improved(i, j, k, m):
    return prop_size_m(i, j, k, m) and not (i == j == k)


def s0(A, B, C, mul):
    """s_0(U,V,W) = Σ_{1 ≤ i < m} u_{i,i} v_{i,i} w_{i,i}."""
    m = A.rows
    out = sp.Integer(0)
    for i in range(1, m):   # 1 ≤ i < m
        out += mul.mul(A[i-1, i-1], B[i-1, i-1], C[i-1, i-1])
    return out


def s1(A, B, C, mul):
    """s_1 — diagonal-uniform: A=B=C in S/S_0/S_1 sense, restricted to S_1∩(≤1 index = m)."""
    m = A.rows
    out = sp.Integer(0)
    for i in range(1, m + 1):
        for j in range(1, m + 1):
            for k in range(1, m + 1):
                if not prop_S1_improved(i, j, k, m):
                    continue
                f1 = A[i-1, j-1] + A[j-1, k-1] + A[k-1, i-1]
                f2 = B[j-1, k-1] + B[k-1, i-1] + B[i-1, j-1]
                f3 = C[k-1, i-1] + C[i-1, j-1] + C[j-1, k-1]
                out += mul.mul(f1, f2, f3)
    return out


def s2(A, B, C, U, V, W, X, Y, Z, mul):
    """s_2 — the full S-style operator restricted to S_2 (≠ diagonal) ∩ (≤1 index = m)."""
    m = A.rows
    out = sp.Integer(0)
    for i in range(1, m + 1):
        for j in range(1, m + 1):
            for k in range(1, m + 1):
                if not prop_S2_improved(i, j, k, m):
                    continue
                f1 = A[i-1, j-1] + U[j-1, k-1] + X[k-1, i-1]
                f2 = B[j-1, k-1] + V[k-1, i-1] + Y[i-1, j-1]
                f3 = C[k-1, i-1] + W[i-1, j-1] + Z[j-1, k-1]
                out += mul.mul(f1, f2, f3)
    return out


def u1(A, B, C, W, Z, mul):
    """u_1 — i ≠ j case of the U-correction (off-diagonal in i,j)."""
    m = A.rows
    out = sp.Integer(0)
    for j in range(1, m):       # 1 ≤ j < m
        for i in range(1, m):   # 1 ≤ i < m, i ≠ j
            if i == j:
                continue
            outer_a = A[i-1, j-1]
            outer_b = B[i-1, j-1]
            inner = m * W[i-1, j-1]
            for k in range(1, m + 1):
                inner += C[k-1, i-1] + Z[j-1, k-1]
            out += mul.mul(outer_a, outer_b, inner)
    return out


def u2(A, B, D, W, Z, mul):
    """u_2 — i = j case (diagonal) of the U-correction."""
    m = A.rows
    out = sp.Integer(0)
    for i in range(1, m):
        outer_a = A[i-1, i-1]
        outer_b = B[i-1, i-1]
        inner = sp.Integer(0)
        for k in range(1, m + 1):
            inner += D[k-1, i-1]
        inner += m * W[i-1, i-1]
        for k in range(1, m + 1):
            inner += Z[i-1, k-1]
        out += mul.mul(outer_a, outer_b, inner)
    return out


def u2_prime(A, B, C, mul):
    """u'_2 — modified u_2 that absorbs the s_0(diag) contribution.

    Coefficient is (m − 9) on the central c_{i,i} (this is where DIS09
    page 28's apparent ‘9’ in s_0 actually lives — it’s tucked into
    `u2_prime`, not `s_0`)."""
    m = A.rows
    out = sp.Integer(0)
    for i in range(1, m):
        outer_a = A[i-1, i-1]
        outer_b = B[i-1, i-1]
        inner = (m - 9) * C[i-1, i-1]
        for k in range(1, m + 1):
            inner += C[k-1, i-1] + C[i-1, k-1]
        out += mul.mul(outer_a, outer_b, inner)
    return out


def u3(A, Y, C, mul):
    """u_3 — i index pinned to m (last row of A-block)."""
    m = A.rows
    out = sp.Integer(0)
    for i in range(1, m):
        outer_a = A[i-1, m-1]
        outer_b = Y[i-1, m-1]
        inner = sp.Integer(0)
        for k in range(1, m + 1):
            inner += C[k-1, i-1]
        out += mul.mul(outer_a, outer_b, inner)
    return out


def u4(A, Y, Z, mul):
    """u_4 — j index pinned to m (last column of A-block)."""
    m = A.rows
    out = sp.Integer(0)
    for j in range(1, m):
        outer_a = A[m-1, j-1]
        outer_b = Y[m-1, j-1]
        inner = sp.Integer(0)
        for k in range(1, m + 1):
            inner += Z[j-1, k-1]
        out += mul.mul(outer_a, outer_b, inner)
    return out


# ──────────────────────────────────────────────────────────────────────
# ProductEven — the full Lemma 4 assembly (port of `TA.mgm` lines 200-296)
# ──────────────────────────────────────────────────────────────────────

def product_even(a, b):
    """Compute C = a·b using Islam's Pan-TA Lemma 4 construction.

    Inputs `a, b` are symbolic Q-coefficient matrices of size 2n × 2n.
    Returns (C, mul_count) where C is the n_recovered × n_recovered
    reconstructed product and mul_count is the count of rank-1 atoms.
    """
    nn = a.rows
    assert nn % 2 == 0, "even case requires even input size"
    n = nn // 2
    m = nn + 2  # padded size

    # C is a symbolic m × m matrix of c_{ij} indeterminates. Per the
    # Magma, the last-but-one row+col of each padded block are zeroed
    # before the sum is computed; the §5.5 improvement.
    C = sp.Matrix(m, m, lambda r, c: sp.Symbol(f"c_{r+1}_{c+1}"))
    # Zero the "drop" positions exactly as the Magma does (lines 211-216).
    for i in range(1, nn + 1):
        C[i-1, n + 1 - 1] = sp.Integer(0)   # col n+1
        C[i-1, m - 1] = sp.Integer(0)       # col m
        C[n + 1 - 1, i-1] = sp.Integer(0)   # row n+1
        C[m - 1, i-1] = sp.Integer(0)       # row m

    # Block decomposition of a, b.
    a11 = a[0:n, 0:n];        a12 = a[0:n, n:2*n]
    a21 = a[n:2*n, 0:n];      a22 = a[n:2*n, n:2*n]
    b11 = b[0:n, 0:n];        b12 = b[0:n, n:2*n]
    b21 = b[n:2*n, 0:n];      b22 = b[n:2*n, n:2*n]

    # C is split into four (n+1) × (n+1) blocks.
    C11 = C[0:n+1, 0:n+1];        C12 = C[0:n+1, n+1:m]
    C21 = C[n+1:m, 0:n+1];        C22 = C[n+1:m, n+1:m]

    # Zero-sum padding.
    A11 = enlarge_left(a11);   A12 = enlarge_left(a12)
    A21 = enlarge_left(a21);   A22 = enlarge_left(a22)
    B11 = enlarge_right(b11);  B12 = enlarge_right(b12)
    B21 = enlarge_right(b21);  B22 = enlarge_right(b22)

    mul = MulCounter()
    product = sp.Integer(0)

    # 4 × s_0
    product += s0(A12 - A11 + A21, B21 + B12 + B11, C11 - C12 + C21, mul)
    product += s0(A12 + A21 - A22, B22 + B12 + B21, C12 + C22 - C21, mul)
    # 2 × s_1 (diagonal)
    product += s1(A11, B11, C11, mul)
    product += s1(A22, B22, C22, mul)
    # 2 × s_2 (off-diagonal triples)
    product += s2(A12, B21, C11,  -A11, B12, -C12,  A21, B11, C21, mul)
    product += s2(A12, B22, C12,  A21, B12, C22,  -A22, B21, -C21, mul)
    # 8 × u_1
    product -= u1( A11, B11,  C11, C11, C11, mul)
    product -= u1(-A11, B21, -C12, C21, C11, mul)
    product -= u1( A12, B21,  C12, C22, -C21, mul)
    product -= u1(-A22, B12, -C21, C12, C22, mul)
    product -= u1( A12, B11,  C11, -C12, C21, mul)
    product -= u1( A21, B12,  C21, C11, -C12, mul)
    product -= u1( A21, B22,  C22, -C21, C12, mul)
    product -= u1( A22, B22,  C22, C22, C22, mul)
    # 2 × u'_2 + 6 × u_2
    product -= u2_prime(A11, B11, C11, mul)
    product -= u2(A11, B21, C12, -C21, -C11, mul)
    product -= u2(A12, B11, C11, -C12, C21, mul)
    product -= u2(A21, B12, C21, C11, -C12, mul)
    product -= u2(A22, B12, C21, -C12, -C22, mul)
    product -= u2(A21, B22, C22, -C21, C12, mul)
    product -= u2(A12, B21, C12, C22, -C21, mul)
    product -= u2_prime(A22, B22, C22, mul)
    # 4 × u_3
    product -= u3(A11 + A12, B11, C11, mul)
    product -= u3(A11 + A12, B21, C12, mul)
    product -= u3(A21 + A22, B12, C21, mul)
    product -= u3(A21 + A22, B22, C22, mul)
    # 4 × u_4
    product -= u4(A11, B11 - B21, C11, mul)
    product -= u4(A12, B11 - B21, C21, mul)
    product -= u4(A21, -B12 + B22, C12, mul)
    product -= u4(A22, -B12 + B22, C22, mul)

    # Extract the m × m coefficient matrix.
    product = sp.expand(product)
    out_C = sp.Matrix(m, m, lambda r, c: sp.Integer(0))
    for i in range(m):
        for j in range(m):
            sym = C[i, j]
            if sym == 0:
                continue
            out_C[i, j] = product.coeff(sym)

    # Top-left n×n of each padded (n+1)×(n+1) block is the recovered
    # original-size product, transposed per the Magma's `Transpose(extr(...))`.
    extr_11 = out_C[0:n, 0:n].T
    extr_12 = out_C[0:n, n+1:n+1+n].T
    extr_21 = out_C[n+1:n+1+n, 0:n].T
    extr_22 = out_C[n+1:n+1+n, n+1:n+1+n].T

    top = sp.Matrix.hstack(extr_11, extr_12)
    bot = sp.Matrix.hstack(extr_21, extr_22)
    return sp.Matrix.vstack(top, bot), mul.count


def expected_mul_count(n_full):
    """Islam Prop 1 — even case."""
    return (n_full**3 + 12 * n_full**2 + 11 * n_full) // 3


def test(n_block):
    n_full = 2 * n_block
    print(f"\n── n_block = {n_block}  (full input {n_full}×{n_full},  "
          f"m = {n_full + 2}) ──")
    A = sp.Matrix(n_full, n_full,
                  lambda r, c: sp.Symbol(f"a_{r+1}_{c+1}"))
    B = sp.Matrix(n_full, n_full,
                  lambda r, c: sp.Symbol(f"b_{r+1}_{c+1}"))
    print("   running ProductEven ...")
    C_TA, mul_count = product_even(A, B)
    print("   computing direct A·B ...")
    C_direct = sp.expand(A * B)

    diff_C = sp.expand(C_TA - C_direct)
    nonzero_entries = sum(1 for c in diff_C if c != 0)
    correctness_ok = (nonzero_entries == 0)

    expected = expected_mul_count(n_full)
    count_ok = (mul_count == expected)

    print(f"   mult count:  got {mul_count},  expected (Islam Prop 1) "
          f"= ({n_full}³ + 12·{n_full}² + 11·{n_full}) / 3 = {expected}   "
          f"→ {'OK' if count_ok else 'FAIL'}")
    print(f"   correctness: {'OK — C_TA == A·B' if correctness_ok else f'FAIL ({nonzero_entries} differing entries)'}")
    if not correctness_ok:
        # First few diffs.
        shown = 0
        for r in range(diff_C.rows):
            for c in range(diff_C.cols):
                if diff_C[r, c] != 0 and shown < 3:
                    print(f"     C[{r+1},{c+1}]: TA − direct = "
                          f"{str(sp.simplify(diff_C[r, c]))[:120]} ...")
                    shown += 1

    return correctness_ok and count_ok


if __name__ == "__main__":
    print("Islam 2009 / DIS09 Appendix Lemma 4 — port from Magma")
    print("=" * 70)
    results = [test(n_block) for n_block in (2, 3)]
    print("\n" + "=" * 70)
    if all(results):
        print("[Phase A step 4] PASS — Lemma 4 construction is correct.")
        print("   sympy verifies BOTH (a) C_TA = A·B symbolically,")
        print("   AND  (b) multiplication count = Islam Prop 1 formula.")
    else:
        print("[Phase A step 4] FAIL — see diagnostics above.")
