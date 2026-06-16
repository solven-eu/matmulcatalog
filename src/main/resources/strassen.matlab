# https://fr.mathworks.com/help/symbolic/solve-a-system-of-linear-equations.html

syms A11B11 A11B12 A11B21 A11B22 A12B11 A12B12 A12B21 A12B22 A21B11 A21B12 A21B21 A21B22 A22B11 A22B12 A22B21 A22B22 a b c d e f g h I11 J11 K11 L11 I12 J12 K12 L12 I21 J21 K21 L21 I22 J22 K22 L22

syms A11 A12 A21 A22 B11 B12 B21 B22 a b c d e f g h I11 J11 K11 L11 I12 J12 K12 L12 I21 J21 K21 L21 I22 J22 K22 L22
eqn1 = A11*B11 + A21*B12 == I11 * (a * A11 + b * A12 + c * A21 + d * A22) * (e * A11 + f * A12 + g * A21 + h * A22);
eqn2 = A11*B21 + A21*B22 == J11 * (a * A11 + b * A12 + c * A21 + d * A22) * (e * A11 + f * A12 + g * A21 + h * A22);
eqn3 = A12*B11 + A22*B12 == K11 * (a * A11 + b * A12 + c * A21 + d * A22) * (e * A11 + f * A12 + g * A21 + h * A22);
eqn4 = A12*B21 + A22*B22 == L11 * (a * A11 + b * A12 + c * A21 + d * A22) * (e * A11 + f * A12 + g * A21 + h * A22);

[A,B] = equationsToMatrix([eqn1, eqn2, eqn3, eqn4], [a b c d e f g h I11 J11 K11 L11 I12 J12 K12 L12 I21 J21 K21 L21 I22 J22 K22 L22])
# -> Unable to convert to matrix form because the system does not seem to be linear


# https://fr.mathworks.com/help/optim/ug/fsolve.html

function F = myfun(x)
F = [A11*B11 + A21*B12 - I11 * (a * A11 + b * A12 + c * A21 + d * A22) * (e * A11 + f * A12 + g * A21 + h * A22);
     A11*B21 + A21*B22 - J11 * (a * A11 + b * A12 + c * A21 + d * A22) * (e * A11 + f * A12 + g * A21 + h * A22);
     A12*B11 + A22*B12 - K11 * (a * A11 + b * A12 + c * A21 + d * A22) * (e * A11 + f * A12 + g * A21 + h * A22);
     A12*B21 + A22*B22 == L11 * (a * A11 + b * A12 + c * A21 + d * A22) * (e * A11 + f * A12 + g * A21 + h * A22);
    ]
    
    
    x0 = [-5;-5];
options = optimoptions('fsolve','Display','iter');

[x,fval] = fsolve(@myfun,x0,options)



syms a b c x
eqn = a*x^2 + b*x + c == 0;
sol = solve(eqn)


# https://fr.mathworks.com/help/symbolic/solve.html
syms A11 A12 A21 A22 B11 B12 B21 B22 
syms a1 b1 c1 d1 e1 f1 g1 h1
syms a2 b2 c2 d2 e2 f2 g2 h2
syms a3 b3 c3 d3 e3 f3 g3 h3
syms a4 b4 c4 d4 e4 f4 g4 h4
syms I11 J11 K11 L11 I12 J12 K12 L12 I21 J21 K21 L21 I22 J22 K22 L22
eqns = [A11*B11 + A21*B12 == I11 * (a1 * A11 + b1 * A12 + c1 * A21 + d1 * A22) * (e1 * B11 + f1 * B12 + g1 * B21 + h1 * B22),
        A11*B21 + A21*B22 == J11 * (a2 * A11 + b2 * A12 + c2 * A21 + d2 * A22) * (e2 * B11 + f2 * B12 + g2 * B21 + h2 * B22),
        A12*B11 + A22*B12 == K11 * (a3 * A11 + b3 * A12 + c3 * A21 + d3 * A22) * (e3 * B11 + f3 * B12 + g3 * B21 + h3 * B22),
        A12*B21 + A22*B22 == L11 * (a4 * A11 + b4 * A12 + c4 * A21 + d4 * A22) * (e4 * B11 + f4 * B12 + g4 * B21 + h4 * B22)];
vars = [a b c d e f g h I11 J11 K11 L11 I12 J12 K12 L12 I21 J21 K21 L21 I22 J22 K22 L22];
[sola, solb, solc, sold, sole, solf, solg, solh, solI11, solJ11, solK11, solL11, solI12, solJ12, solK12, solL12, solI21, solJ21, solK21, solL21, solI22, solJ22, solK22, solL22] = solve(eqns, vars)