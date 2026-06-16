"""Symbolic verification of Hopcroft-Kerr ⟨3,2,3⟩=15 example from page 4 of the paper."""
import sympy as sp

# Symbolic A (3x2) and X (2x3)
a = sp.symbols('a1_1 a1_2 a2_1 a2_2 a3_1 a3_2', commutative=True)
A11, A12, A21, A22, A31, A32 = a
x = sp.symbols('x1_1 x2_1 x1_2 x2_2 x1_3 x2_3', commutative=True)
X11, X21, X12, X22, X13, X23 = x

# A(a_i, x_l) = a_{i,2} * (x_{1,l} + x_{2,l})
def A_(ai2, x1l, x2l): return ai2 * (x1l + x2l)
def B_(ai1, ai2, x1l): return (ai1 - ai2) * x1l
def C_(aj1, aj2, x2l): return (aj1 - aj2) * x2l
def D_(aj1, x1l, x2l): return aj1 * (x1l + x2l)
def E_(ai2, aj2, x2k, x2l): return (ai2 + aj2) * (x2k + x2l)
def F_(ai1, aj1, x1k, x1l): return (ai1 + aj1) * (x1k + x1l)
def G_(aj1, ai2, x1k, x2l): return (aj1 + ai2) * (x1k - x2l)

# Page-4 formulas for ⟨3,2,3⟩
# Method (1) diagonal: y_{ii} = A(a_i, x_i) + B(a_i, x_i) — y_{11}
# Method (2) diagonal: y_{ii} = -C(a_i, x_i) + D(a_i, x_i) — y_{22}
# Method (3) diagonal: y_{ii} = E(a_i, x_i) + F(a_i, x_i) — y_{33}
# (where E/F applied to a single (a_i, x_i) means the "sum" is just a_i with a_j=0)

# Single-arg E and F: E(a, x) = a_2 * x_2, F(a, x) = a_1 * x_1
def E1_(ai2, x2l): return ai2 * x2l
def F1_(ai1, x1l): return ai1 * x1l

# Single-arg A and B (when the "subscript" is a substitution):
# A((-a_i+a_j), (-x_k+x_l)) = (-a_{i,2}+a_{j,2}) * ((-x_{1,k}+x_{1,l}) + (-x_{2,k}+x_{2,l}))
# This is just the standard A definition applied to the substituted indeterminate.

# Diagonals
y11 = A_(A12, X11, X21) + B_(A11, A12, X11)             # method (1)
y22 = -C_(A21, A22, X22) + D_(A21, X12, X22)             # method (2)
y33 = E1_(A32, X23) + F1_(A31, X13)                       # method (3)

# Page-4 off-diagonal formulas:
# y_{12} = -B(a_1, x_1) - D(a_2, x_2) + F(a_1+a_2, x_1+x_2) - G(a_1, a_2, x_1, x_2)
y12 = -B_(A11, A12, X11) - D_(A21, X12, X22) + F_(A11, A21, X11, X12) - G_(A21, A12, X11, X22)
# y_{21} = -A(a_1, x_1) + C(a_2, x_2) + E(a_1+a_2, x_1+x_2) + G(a_1, a_2, x_1, x_2)
y21 = -A_(A12, X11, X21) + C_(A21, A22, X22) + E_(A12, A22, X21, X22) + G_(A21, A12, X11, X22)

# y_{13} = A(a_1, x_1) - D(-a_1+a_3, -x_1+x_3) + F(a_3, x_3) - G(a_1, -a_1+a_3, x_1, -x_1+x_3)
y13 = (A_(A12, X11, X21)
       - D_(-A11+A31, -X11+X13, -X21+X23)
       + F1_(A31, X13)
       - G_(-A11+A31, A12, X11, -X21+X23))
# y_{31} = B(a_1, x_1) + D(-a_1+a_3, -x_1+x_3) + E(a_3, x_3) + G(a_1, -a_1+a_3, x_1, -x_1+x_3)
y31 = (B_(A11, A12, X11)
       + D_(-A11+A31, -X11+X13, -X21+X23)
       + E1_(A32, X23)
       + G_(-A11+A31, A12, X11, -X21+X23))

# y_{23} = -A(-a_2+a_3, -x_2+x_3) + D(a_2, x_2) + E(a_3, x_3) + G(-a_2+a_3, a_2, -x_2+x_3, x_2)
y23 = (-A_(-A22+A32, -X12+X13, -X22+X23)
       + D_(A21, X12, X22)
       + E1_(A32, X23)
       + G_(-A21+A31, A22, -X11+X12, X22))  # Hmm let me re-check args
# Wait the paper formula: G(-a_2+a_3, a_2, -x_2+x_3, x_2). G takes (a_j, a_i, x_k, x_l) so:
# a_j = -a_2+a_3 (subscript 1 = -a_{2,1}+a_{3,1}), a_i = a_2 (subscript 2 = a_{2,2})
# x_k = -x_2+x_3, x_l = x_2
# G = (-a_{2,1}+a_{3,1} + a_{2,2}) * ((-x_{1,2}+x_{1,3}) - x_{2,2})
y23 = (-A_(-A22+A32, -X12+X13, -X22+X23)
       + D_(A21, X12, X22)
       + E1_(A32, X23)
       + G_(-A21+A31, A22, -X12+X13, X22))

# y_{32} = -B(-a_2+a_3, -x_2+x_3) - C(a_2, x_2) + F(-a_2+a_3, x_2) - G(-a_2+a_3, a_2, -x_2+x_3, x_2)
# Hmm — F(-a_2+a_3, x_2) is single-arg form. F(a, x) = a_1 * x_1.
# So F(-a_2+a_3, x_2) = (-a_{2,1}+a_{3,1}) * x_{1,2}
y32 = (-B_(-A21+A31, -A22+A32, -X12+X13)
       - C_(A21, A22, X22)
       + F1_(-A21+A31, X12)
       - G_(-A21+A31, A22, -X12+X13, X22))

# Expected outputs
e11 = A11*X11 + A12*X21
e12 = A11*X12 + A12*X22
e13 = A11*X13 + A12*X23
e21 = A21*X11 + A22*X21
e22 = A21*X12 + A22*X22
e23 = A21*X13 + A22*X23
e31 = A31*X11 + A32*X21
e32 = A31*X12 + A32*X22
e33 = A31*X13 + A32*X23

# Verify (y_{33} is computed as E(a_3,x_3)+F(a_3,x_3) = a_{3,2}x_{2,3} + a_{3,1}x_{1,3} = e33)
for name, computed, expected in [
    ('y11', y11, e11), ('y22', y22, e22), ('y33', y33, e33),
    ('y12', y12, e12), ('y21', y21, e21),
    ('y13', y13, e13), ('y31', y31, e31),
    ('y23', y23, e23), ('y32', y32, e32),
]:
    diff = sp.expand(computed - expected)
    status = 'PASS' if diff == 0 else f'FAIL diff={diff}'
    print(f'  {name}: {status}')
