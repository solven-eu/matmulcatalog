"""
Phase A step 5 — port of Islam's Magma `ProductOdd` (TA.mgm lines 304-440).

The odd case (input size `nn = 2n+1`) reuses the same operator zoo as
the even case (`s_0, s_1, s_2, u_1, u_2, u'_2, u_3, u_4`) but with two
differences:

1. **Rectangular block decomposition**: a/b are split as
   `a^{1,1}: n×n`, `a^{1,2}: n×(n+1)`, `a^{2,1}: (n+1)×n`, `a^{2,2}: (n+1)×(n+1)`
   (and similarly for b). Each block is enlarged to size (n+1)×(n+1) by
   zero-row/col insertion **at the second-last position** (not the last),
   because the last position is already cheap from earlier savings.
2. **Argument permutations** in `s_2` (compared to the even case, the 9
   args are cyclically rotated — equivalent under S₂'s permutation
   symmetry).

The expected multiplication count for input size `nn` (odd) is
**`(nn³ + 15·nn² + 14·nn − 6) / 3`** — Islam Prop 1 odd case.

Run: `python3 05_lemma4_odd_from_magma.py`
"""
import sympy as sp

# Re-use the operator definitions from the even-case file.
import importlib.util, sys, os
_spec = importlib.util.spec_from_file_location(
    "step4", os.path.join(os.path.dirname(__file__), "04c_lemma4_from_magma.py"))
_step4 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_step4)

MulCounter = _step4.MulCounter
enlarge_left = _step4.enlarge_left
enlarge_right = _step4.enlarge_right
s0 = _step4.s0
s1 = _step4.s1
s2 = _step4.s2
u1 = _step4.u1
u2 = _step4.u2
u2_prime = _step4.u2_prime
u3 = _step4.u3
u4 = _step4.u4


def swap_cols(M, c1, c2):
    """Return M with columns c1 and c2 swapped (1-based indices)."""
    out = M.copy()
    for r in range(M.rows):
        out[r, c1 - 1], out[r, c2 - 1] = M[r, c2 - 1], M[r, c1 - 1]
    return out


def swap_rows(M, r1, r2):
    out = M.copy()
    for c in range(M.cols):
        out[r1 - 1, c], out[r2 - 1, c] = M[r2 - 1, c], M[r1 - 1, c]
    return out


def insert_block(target_size_rc, block, top_left=(1, 1)):
    """Build a target_size_rc zero matrix with `block` inserted at top_left (1-based)."""
    rows, cols = target_size_rc
    out = sp.zeros(rows, cols)
    tr, tc = top_left
    for r in range(block.rows):
        for c in range(block.cols):
            out[tr - 1 + r, tc - 1 + c] = block[r, c]
    return out


def product_odd(a, b):
    """Compute C = a·b for ODD input size `nn = 2n+1` via Islam's Lemma 4 odd-case."""
    nn = a.rows
    assert nn % 2 == 1, "odd case requires odd input size"
    n = (nn + 1) // 2     # block size (the "+1" handles the odd input via padding)
    m = nn + 2

    # Symbolic C of size m × m.
    C = sp.Matrix(m, m, lambda r, c: sp.Symbol(f"c_{r+1}_{c+1}"))

    # Cut a into rectangular blocks per the Magma:
    #   a11: n × n,  a12: n × (nn-n),  a21: (nn-n) × n,  a22: (nn-n) × (nn-n).
    a11 = a[0:n, 0:n];           a12 = a[0:n, n:nn]
    a21 = a[n:nn, 0:n];          a22 = a[n:nn, n:nn]

    b11 = b[0:n, 0:n];           b12 = b[0:n, n:nn]
    b21 = b[n:nn, 0:n];          b22 = b[n:nn, n:nn]

    # C blocks (rectangular slices of C).
    c11 = C[0:n+1,  0:n+1]
    c12 = C[0:n+1,  n+1:nn-n+1+n+1]
    c21 = C[n+1:nn-n+1+n+1, 0:n+1]
    c22 = C[n+1:nn-n+1+n+1, n+1:nn-n+1+n+1]

    # Enlarge a blocks via L · a · R (sizes grow by +1 in each dimension).
    a11_en = enlarge_left(a11)            # (n+1) × (n+1)
    a12_en = enlarge_left(a12)            # (n+1) × (nn-n+1)
    a21_en = enlarge_left(a21)            # (nn-n+1) × (n+1)
    a22_en = enlarge_left(a22)            # (nn-n+1) × (nn-n+1)

    # Pad each enlarged block to (n+1)×(n+1) with zeros, then push the
    # zero rows/cols from the END to the LAST-BUT-ONE position via swaps.
    target = (n + 1, n + 1)
    A11 = insert_block(target, a11_en)                # no swap (already square right size)
    A12 = swap_cols(insert_block(target, a12_en), n, n + 1)
    A21 = swap_rows(insert_block(target, a21_en), n, n + 1)
    A22 = swap_cols(swap_rows(insert_block(target, a22_en), n, n + 1), n, n + 1)

    b11_en = enlarge_right(b11)
    b12_en = enlarge_right(b12)
    b21_en = enlarge_right(b21)
    b22_en = enlarge_right(b22)

    B11 = insert_block(target, b11_en)
    B12 = swap_cols(insert_block(target, b12_en), n, n + 1)
    B21 = swap_rows(insert_block(target, b21_en), n, n + 1)
    B22 = swap_cols(swap_rows(insert_block(target, b22_en), n, n + 1), n, n + 1)

    # C blocks padded to (n+1)×(n+1) with zeros at last rows/cols
    # (Magma inserts at (1,1), so trailing rows/cols are zero by default).
    C11 = insert_block(target, c11)
    C12 = insert_block(target, c12)
    C21 = insert_block(target, c21)
    C22 = insert_block(target, c22)

    mul = MulCounter()
    product = sp.Integer(0)

    product += s0(A12 - A11 + A21, B21 + B12 + B11, C11 - C12 + C21, mul)
    product += s0(A12 + A21 - A22, B22 + B12 + B21, C12 + C22 - C21, mul)
    product += s1(A11, B11, C11, mul)
    product += s1(A22, B22, C22, mul)
    # NOTE: arg order in s_2 differs from even case (Magma lines 401-402).
    product += s2(-A11, B12, -C12,  A21, B11, C21,  A12, B21, C11, mul)
    product += s2(-A22, B21, -C21,  A12, B22, C12,  A21, B12, C22, mul)

    # In Magma odd case u_2 block appears BEFORE u_1.
    product -= u2_prime(A11, B11, C11, mul)
    product -= u2(-A11, B21, -C12, C21, C11, mul)
    product -= u2( A12, B11,  C11, -C12, C21, mul)
    product -= u2( A21, B12,  C21, C11, -C12, mul)
    product -= u2(-A22, B12, -C21, C12, C22, mul)
    product -= u2( A21, B22,  C22, -C21, C12, mul)
    product -= u2( A12, B21,  C12, C22, -C21, mul)
    product -= u2_prime(A22, B22, C22, mul)

    product -= u1( A11, B11,  C11, C11, C11, mul)
    product -= u1(-A11, B21, -C12, C21, C11, mul)
    product -= u1( A12, B21,  C12, C22, -C21, mul)
    product -= u1(-A22, B12, -C21, C12, C22, mul)
    product -= u1( A12, B11,  C11, -C12, C21, mul)
    product -= u1( A21, B12,  C21, C11, -C12, mul)
    product -= u1( A21, B22,  C22, -C21, C12, mul)
    product -= u1( A22, B22,  C22, C22, C22, mul)

    product -= u3(A11 + A12, B11, C11, mul)
    product -= u3(A11 + A12, B21, C12, mul)
    product -= u3(A21 + A22, B12, C21, mul)
    product -= u3(A21 + A22, B22, C22, mul)

    product -= u4(A11, B11 - B21, C11, mul)
    product -= u4(A12, B11 - B21, C21, mul)
    product -= u4(A21, -B12 + B22, C12, mul)
    product -= u4(A22, -B12 + B22, C22, mul)

    product = sp.expand(product)

    # Extract coefficient matrix entry-by-entry.
    out_C = sp.Matrix(m, m, lambda r, c: sp.Integer(0))
    for i in range(m):
        for j in range(m):
            sym = C[i, j]
            if sym == 0:
                continue
            out_C[i, j] = product.coeff(sym)

    # Each block extracted as n×n (Magma uses the same `extr=Submatrix(.,1,1,n,n)`
    # for the odd case too) — uniform size makes the assembly trivial, and the
    # "extra" rows/cols beyond the original rectangular block sizes hold zero
    # coefficients (padding c-variables that the algorithm never writes to).
    # After assembling to 2n × 2n we trim to nn × nn.
    extr_11 = out_C[0:n, 0:n].T
    extr_12 = out_C[0:n, n+1:n+1+n].T
    extr_21 = out_C[n+1:n+1+n, 0:n].T
    extr_22 = out_C[n+1:n+1+n, n+1:n+1+n].T

    top = sp.Matrix.hstack(extr_11, extr_12)
    bot = sp.Matrix.hstack(extr_21, extr_22)
    assembled = sp.Matrix.vstack(top, bot)
    return assembled[0:nn, 0:nn], mul.count


def expected_mul_count_odd(nn):
    return (nn**3 + 15 * nn**2 + 14 * nn - 6) // 3


def test(nn):
    print(f"\n── nn = {nn}  (odd, n = (nn+1)/2 = {(nn+1)//2},  m = {nn + 2}) ──")
    A = sp.Matrix(nn, nn, lambda r, c: sp.Symbol(f"a_{r+1}_{c+1}"))
    B = sp.Matrix(nn, nn, lambda r, c: sp.Symbol(f"b_{r+1}_{c+1}"))
    print("   running ProductOdd ...")
    C_TA, mul_count = product_odd(A, B)
    print("   computing direct A·B ...")
    C_direct = sp.expand(A * B)
    diff_C = sp.expand(C_TA - C_direct)
    nonzero = sum(1 for e in diff_C if e != 0)
    correctness_ok = (nonzero == 0)
    expected = expected_mul_count_odd(nn)
    count_ok = (mul_count == expected)
    print(f"   mult count:  got {mul_count}, expected (Islam Prop 1 odd) "
          f"= ({nn}³ + 15·{nn}² + 14·{nn} − 6) / 3 = {expected}   "
          f"→ {'OK' if count_ok else 'FAIL'}")
    print(f"   correctness: {'OK — C_TA == A·B' if correctness_ok else f'FAIL ({nonzero} differing entries)'}")
    if not correctness_ok:
        shown = 0
        for r in range(diff_C.rows):
            for c in range(diff_C.cols):
                if diff_C[r, c] != 0 and shown < 3:
                    print(f"     C[{r+1},{c+1}]: TA − direct = "
                          f"{str(sp.simplify(diff_C[r, c]))[:120]} ...")
                    shown += 1
    return correctness_ok and count_ok


if __name__ == "__main__":
    print("Islam 2009 / DIS09 Appendix Lemma 4 ODD CASE — port from Magma")
    print("=" * 70)
    results = [test(nn) for nn in (3, 5)]
    print("\n" + "=" * 70)
    if all(results):
        print("[Phase A step 5] PASS — odd-case Lemma 4 construction is correct.")
    else:
        print("[Phase A step 5] FAIL — see diagnostics above.")
