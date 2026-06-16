% execute with: 
% - swipl (to start prolog)
% - ['/Users/benoit/workspace/strassen/src/main/resources/strassen7.prolog'].
% - strassen(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,A1,A2,A3,A4,B1,B2,B3,B4,C1,C2,C3,C4,D1,D2,D3,D4,E1,E2,E3,E4,F1,F2,F3,F4,G1,G2,G3,G4,H1,H2,H3,H4).
% - time(strassen(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,A1,A2,A3,A4,B1,B2,B3,B4,C1,C2,C3,C4,D1,D2,D3,D4,E1,E2,E3,E4,F1,F2,F3,F4,G1,G2,G3,G4,H1,H2,H3,H4)).
% - time(strassen(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,A1,A2,A3,A4,B1,B2,B3,B4,C1,C2,C3,C4,D1,D2,D3,D4,E1,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 )).
% - time(strassen(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,A1,A2,A3,A4,B1,B2,B3,B4,C1,C2,C3,C4,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 )).
% - time(strassen(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,A1,A2,A3,A4,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 )).
% - time(strassen(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 )).
% trace. notrace.

% http://www.swi-prolog.org/pldoc/man?section=tabling
% Tabling is prolog term for predicate-result caching
% use_module(library(tabling)).

:- table is_one/3.
% :- table is_zero/3.

% solve(I, J, K, L, A, B, C, D, E, F, G, H).
% I = A, A = [0, 0, 0, 0, -1],
% J = D, D = [0, 0, 0, -1, 0],
% K = [_2290, _2296, _2302, _2308, _2314],
% L = [_2320, _2326, _2332, _2338, _2344],
% B = [_2380, _2386, _2392, _2398, _2404],
% C = [_2410, _2416, _2422, _2428, _2434],
% E = [0, 0, 0, 0, 1],
% F = [0, 0, 0, 1, 0],
% G = [_2530, _2536, _2542, _2548, _2554],
% H = [_2560, _2566, _2572, _2578, _2584] .

restrict(A):-
	nth1(1,A,U),
	member(U, [0, -1, 1]),
	nth1(2,A,V),
	member(V, [0, -1, 1]),
	nth1(3,A,W),
	member(W, [0, -1, 1]),
	nth1(4,A,X),
	member(X, [0, -1, 1]),
	nth1(5, A, Y),
	member(Y, [0, -1, 1]).
	

restrict(I, A, E, O):-
	nth1(O,I,I1),
	member(I1, [0, -1, 1]),
	nth1(O,A,A1),
	member(A1, [0, -1, 1]),
	nth1(O,E,E1),
	member(E1, [0, -1, 1]).

mulAll(I, A, E, MT):-
	member(O1, [0, -1, 1]),
	member(O2, [0, -1, 1]),
	member(O3, [0, -1, 1]),
	member(O4, [0, -1, 1]),
	member(O5, [0, -1, 1]),
	
	sumAll(O1, O2, O3, O4, O5, MT),
	
	restrict(I, A, E, 1),
	nth1(1,I,I1),
	nth1(1,A,A1),
	nth1(1,E,E1),
	O1 is I1*A1*E1,
	
	restrict(I, A, E, 2),
	nth1(2,I,I2),
	nth1(2,A,A2),
	nth1(2,E,E2),
	O2 is I2*A2*E2,

	restrict(I, A, E, 3),	
	nth1(3,I,I3),
	nth1(3,A,A3),
	nth1(3,E,E3),
	O3 is I3*A3*E3,
	
	restrict(I, A, E, 4),
	nth1(4,I,I4),
	nth1(4,A,A4),
	nth1(4,E,E4),
	O4 is I4*A4*E4,
	
	restrict(I, A, E, 5),
	nth1(5,I,I5),
	nth1(5,A,A5),
	nth1(5,E,E5),
	O5 is I5*A5*E5.
	
	
mulAll2(I, A, E, MT):-
	nth1(1,I,I1),
	nth1(1,A,A1),
	nth1(1,E,E1),
	O1 is I1*A1*E1,
	
	nth1(2,I,I2),
	nth1(2,A,A2),
	nth1(2,E,E2),
	O2 is I2*A2*E2,

	nth1(3,I,I3),
	nth1(3,A,A3),
	nth1(3,E,E3),
	O3 is I3*A3*E3,
	
	nth1(4,I,I4),
	nth1(4,A,A4),
	nth1(4,E,E4),
	O4 is I4*A4*E4,
	
	nth1(5,I,I5),
	nth1(5,A,A5),
	nth1(5,E,E5),
	O5 is I5*A5*E5,
	
	sumAll(O1, O2, O3, O4, O5, MT).
	
sumAll(O1, O2, O3, O4, O5, ST):-
	ST is O1+O2+O3+O4+O5.
	
is_one(I, A, E):-
	mulAll(I, A, E, 1).
	
is_zero(I, A, E):-
	mulAll2(I,A,E, 0).
	
solve(I, J, K, L, A, B, C, D, E, F, G, H):-
	P_SIZE is 5,
	length(I,P_SIZE),
	length(J,P_SIZE),
	length(K,P_SIZE),
	length(L,P_SIZE),
	length(A,P_SIZE),
	length(B,P_SIZE),
	length(C,P_SIZE),
	length(D,P_SIZE),
	length(E,P_SIZE),
	length(F,P_SIZE),
	length(G,P_SIZE),
	length(H,P_SIZE),
	is_one(I,A,E),
	is_one(J,D,F).
	%is_zero(I,D,F).
	%is_zero(J,A,E),
	%is_zero(I,A,F),
	%is_zero(J,A,F),
	%is_zero(I,D,E),
	%is_zero(J,D,E).
