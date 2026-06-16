"""
Islam 2009 §5.4 — Phase A, step 3: verify Lemmas 3 and 4.

Lemma 3 (used for diagonal products `T(Ã^{i,i}, B̃^{i,i}, C̃^{i,i})`):
    S(U, V, W, U, V, W, U, V, W)  =  27 · S_0(U, V, W)  +  3 · S_1(U, V, W)

Lemma 4 (used for off-diagonal disjoint triples in Proposition 2):
    S(A, B, C, U, V, W, X, Y, Z)
        =  S_2(A, B, C, U, V, W, X, Y, Z)  +  S_0(A + U + X, B + V + Y, C + W + Z)

with the index-set decomposition (m × m matrices):

    S        = {(i,j,k) : 1 ≤ i,j,k ≤ m}
    S_0      = {(i,i,i) : 1 ≤ i ≤ m}                  (diagonal)
    S_1      = {(i,j,k) : 1 ≤ i ≤ j < k ≤ m  or
                                1 ≤ k < j ≤ i ≤ m}    (ordered triples)
    S_2      = {(j,k,i) : (i,j,k) ∈ S_1}
    S_3      = {(k,i,j) : (i,j,k) ∈ S_1}
    S'_2     = S \\ S_0                                (off-diagonal)

V_α(U,V,W) is the restricted sum over S_α with all three triples coinciding.

Run: python3 03_lemmas_3_and_4.py
"""
import sympy as sp


def symbolic_matrix(prefix, m):
    return sp.Matrix(m, m, lambda r, c: sp.Symbol(f"{prefix}_{r+1}{c+1}"))


def S_ijk(A, B, C, U, V, W, X, Y, Z, i, j, k):
    """Single S_{i,j,k} product. Indices are 1-based."""
    f1 = A[i-1, j-1] + U[j-1, k-1] + X[k-1, i-1]
    f2 = B[j-1, k-1] + V[k-1, i-1] + Y[i-1, j-1]
    f3 = C[k-1, i-1] + W[i-1, j-1] + Z[j-1, k-1]
    return f1 * f2 * f3


def S_full(A, B, C, U, V, W, X, Y, Z, m):
    out = sp.Integer(0)
    for i in range(1, m + 1):
        for j in range(1, m + 1):
            for k in range(1, m + 1):
                out += S_ijk(A, B, C, U, V, W, X, Y, Z, i, j, k)
    return out


def S0_full(U, V, W, m):
    """Σ_i u_{i,i} v_{i,i} w_{i,i}."""
    out = sp.Integer(0)
    for i in range(1, m + 1):
        out += U[i-1, i-1] * V[i-1, i-1] * W[i-1, i-1]
    return out


def s1_indices(m):
    """The index set S_1 = {(i,j,k): 1 ≤ i ≤ j < k ≤ m  OR  1 ≤ k < j ≤ i ≤ m}."""
    out = []
    for i in range(1, m + 1):
        for j in range(1, m + 1):
            for k in range(1, m + 1):
                if (1 <= i <= j < k <= m) or (1 <= k < j <= i <= m):
                    out.append((i, j, k))
    return out


def S1_uniform(U, V, W, m):
    """S_1 restricted to the 'all-three-same' case (A=U=X, B=V=Y, C=W=Z)."""
    out = sp.Integer(0)
    for (i, j, k) in s1_indices(m):
        out += S_ijk(U, V, W, U, V, W, U, V, W, i, j, k)
    return out


def Sprime2_full(A, B, C, U, V, W, X, Y, Z, m):
    """S restricted to S' = S \\ S_0 = all (i,j,k) except (i,i,i)."""
    out = sp.Integer(0)
    for i in range(1, m + 1):
        for j in range(1, m + 1):
            for k in range(1, m + 1):
                if i == j == k:
                    continue
                out += S_ijk(A, B, C, U, V, W, X, Y, Z, i, j, k)
    return out


def test_lemma3(m):
    print(f"\n── Lemma 3, m = {m} ──")
    U = symbolic_matrix("u", m)
    V = symbolic_matrix("v", m)
    W = symbolic_matrix("w", m)
    lhs = sp.expand(S_full(U, V, W, U, V, W, U, V, W, m))
    rhs = sp.expand(27 * S0_full(U, V, W, m) + 3 * S1_uniform(U, V, W, m))
    diff = sp.expand(lhs - rhs)
    if diff == 0:
        print(f"   OK — S(U,V,W,U,V,W,U,V,W) = 27·S_0 + 3·S_1   for m = {m}")
        return True
    print(f"   FAIL — residual = {diff.as_ordered_terms()[:3]}")
    return False


def test_lemma4(m):
    print(f"\n── Lemma 4, m = {m} ──")
    A = symbolic_matrix("a", m); B = symbolic_matrix("b", m); C = symbolic_matrix("c", m)
    U = symbolic_matrix("u", m); V = symbolic_matrix("v", m); W = symbolic_matrix("w", m)
    X = symbolic_matrix("x", m); Y = symbolic_matrix("y", m); Z = symbolic_matrix("z", m)
    lhs = sp.expand(S_full(A, B, C, U, V, W, X, Y, Z, m))
    # S_0(A+U+X, B+V+Y, C+W+Z) — diagonal at the summed matrices
    APU = sp.expand(A + U + X)
    BPV = sp.expand(B + V + Y)
    CPW = sp.expand(C + W + Z)
    rhs = sp.expand(
        Sprime2_full(A, B, C, U, V, W, X, Y, Z, m)
        + S0_full(APU, BPV, CPW, m)
    )
    diff = sp.expand(lhs - rhs)
    if diff == 0:
        print(f"   OK — S(A,B,C,U,V,W,X,Y,Z) = S'_2 + S_0(A+U+X,B+V+Y,C+W+Z)   for m = {m}")
        return True
    print(f"   FAIL — residual = {diff.as_ordered_terms()[:3]}")
    return False


if __name__ == "__main__":
    print("Islam 2009 §5.4: Lemmas 3 and 4")
    print("=" * 70)
    results = []
    for m in (3, 4):
        results.append(test_lemma3(m))
        results.append(test_lemma4(m))
    print("\n" + "=" * 70)
    if all(results):
        print("[Phase A step 3] PASS — Lemmas 3 + 4 hold for m ∈ {3, 4}.")
    else:
        print("[Phase A step 3] FAIL — re-examine S_0 / S_1 / S'_2 definitions.")
