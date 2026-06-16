% execute with: 
% - swipl (to start prolog)
% - ['/Users/benoit/workspace/strassen/src/main/resources/strassen.prolog'].
% - strassen(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,A1,A2,A3,A4,B1,B2,B3,B4,C1,C2,C3,C4,D1,D2,D3,D4,E1,E2,E3,E4,F1,F2,F3,F4,G1,G2,G3,G4,H1,H2,H3,H4).
% - time(strassen(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,A1,A2,A3,A4,B1,B2,B3,B4,C1,C2,C3,C4,D1,D2,D3,D4,E1,E2,E3,E4,F1,F2,F3,F4,G1,G2,G3,G4,H1,H2,H3,H4)).
% - time(strassen(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,A1,A2,A3,A4,B1,B2,B3,B4,C1,C2,C3,C4,D1,D2,D3,D4,E1,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 )).
% - time(strassen(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,A1,A2,A3,A4,B1,B2,B3,B4,C1,C2,C3,C4,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 )).
% - time(strassen(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,A1,A2,A3,A4,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 )).
% - time(strassen(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 )).

strassen(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,A1,A2,A3,A4,B1,B2,B3,B4,C1,C2,C3,C4,D1,D2,D3,D4,E1,E2,E3,E4,F1,F2,F3,F4,G1,G2,G3,G4,H1,H2,H3,H4):-

%	is_one(I1,A1,E1,J1,A2,E2,K1,A3,E3,L1,A4,E4),
	ijklIsInversible(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,A1,A2,A3,A4,B1,B2,B3,B4,C1,C2,C3,C4,D1,D2,D3,D4,E1,E2,E3,E4,F1,F2,F3,F4,G1,G2,G3,G4,H1,H2,H3,H4),

I1 is I1.

restrict(A):-
	member(A, [0, -1, 1]).
restrict(I1,J1,K1,L1):-	
	restrict(I1),
	restrict(J1),
	restrict(K1),
	restrict(L1).

	
is_one(I1,A1,E1,J1,A2,E2,K1,A3,E3,L1,A4,E4):-
	restrict(I1),
	restrict(J1),
	restrict(K1),
	restrict(L1),
	
	restrict(A1),
	restrict(A2),
	restrict(A3),
	restrict(A4),
	
	restrict(E1),
	restrict(E2),
	restrict(E3),
	restrict(E4),
	
	1 is I1*A1*E1+J1*A2*E2+K1*A3*E3+L1*A4*E4.
	
is_zero(I1,A1,E1,J1,A2,E2,K1,A3,E3,L1,A4,E4):-
	restrict(I1),
	restrict(J1),
	restrict(K1),
	restrict(L1),
	
	restrict(A1),
	restrict(A2),
	restrict(A3),
	restrict(A4),
	
	restrict(E1),
	restrict(E2),
	restrict(E3),
	restrict(E4),
	
	0 is I1*A1*E1+J1*A2*E2+K1*A3*E3+L1*A4*E4.

ijklIsInversible(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,A1,A2,A3,A4,B1,B2,B3,B4,C1,C2,C3,C4,D1,D2,D3,D4,E1,E2,E3,E4,F1,F2,F3,F4,G1,G2,G3,G4,H1,H2,H3,H4):-
	restrict(I1,J1,K1,L1),
	restrict(I2,J2,K2,L2),
	restrict(I3,J3,K3,L3),
	restrict(I4,J4,K4,L4),
	
	
	% https://stackoverflow.com/questions/1148309/inverting-a-4x4-matrix
	
	% https://stackoverflow.com/questions/30464504/how-to-find-the-nth-element-of-a-list-in-prolog
	%inv([Z0,Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,Z10,Z11,Z12,Z13,Z14,Z15]),
	
	
	Z0 is J2  * K3 * L4 - 
             J2  * L3 * K4 - 
             J3  * K2  * L4 + 
             J3  * L2  * K4 +
             J4 * K2  * L3 - 
             J4 * L2  * K3,

    Z4 is -I2  * K3 * L4 + 
              I2  * L3 * K4 + 
              I3  * K2  * L4 - 
              I3  * L2  * K4 - 
              I4 * K2  * L3 + 
              I4 * L2  * K3,

    Z8 is I2  * J3 * L4 - 
             I2  * L3 * J4 - 
             I3  * J2 * L4 + 
             I3  * L2 * J4 + 
             I4 * J2 * L3 - 
             I4 * L2 * J3,

    Z12 is -I2  * J3 * K4 + 
               I2  * K3 * J4 +
               I3  * J2 * K4 - 
               I3  * K2 * J4 - 
               I4 * J2 * K3 + 
               I4 * K2 * J3,

    Z1 is -J1  * K3 * L4 + 
              J1  * L3 * K4 + 
              J3  * K1 * L4 - 
              J3  * L1 * K4 - 
              J4 * K1 * L3 + 
              J4 * L1 * K3,

    Z5 is I1  * K3 * L4 - 
             I1  * L3 * K4 - 
             I3  * K1 * L4 + 
             I3  * L1 * K4 + 
             I4 * K1 * L3 - 
             I4 * L1 * K3,

    Z9 is -I1  * J3 * L4 + 
              I1  * L3 * J4 + 
              I3  * J1 * L4 - 
              I3  * L1 * J4 - 
              I4 * J1 * L3 + 
              I4 * L1 * J3,

    Z13 is I1  * J3 * K4 - 
              I1  * K3 * J4 - 
              I3  * J1 * K4 + 
              I3  * K1 * J4 + 
              I4 * J1 * K3 - 
              I4 * K1 * J3,

    Z2 is J1  * K2 * L4 - 
             J1  * L2 * K4 - 
             J2  * K1 * L4 + 
             J2  * L1 * K4 + 
             J4 * K1 * L2 - 
             J4 * L1 * K2,

    Z6 is -I1  * K2 * L4 + 
              I1  * L2 * K4 + 
              I2  * K1 * L4 - 
              I2  * L1 * K4 - 
              I4 * K1 * L2 + 
              I4 * L1 * K2,

    Z10 is I1  * J2 * L4 - 
              I1  * L2 * J4 - 
              I2  * J1 * L4 + 
              I2  * L1 * J4 + 
              I4 * J1 * L2 - 
              I4 * L1 * J2,

    Z14 is -I1  * J2 * K4 + 
               I1  * K2 * J4 + 
               I2  * J1 * K4 - 
               I2  * K1 * J4 - 
               I4 * J1 * K2 + 
               I4 * K1 * J2,

    Z3 is -J1 * K2 * L3 + 
              J1 * L2 * K3 + 
              J2 * K1 * L3 - 
              J2 * L1 * K3 - 
              J3 * K1 * L2 + 
              J3 * L1 * K2,

    Z7 is I1 * K2 * L3 - 
             I1 * L2 * K3 - 
             I2 * K1 * L3 + 
             I2 * L1 * K3 + 
             I3 * K1 * L2 - 
             I3 * L1 * K2,

    Z11 is -I1 * J2 * L3 + 
               I1 * L2 * J3 + 
               I2 * J1 * L3 - 
               I2 * L1 * J3 - 
               I3 * J1 * L2 + 
               I3 * L1 * J2,

    Z15 is I1 * J2 * K3 - 
              I1 * K2 * J3 - 
              I2 * J1 * K3 + 
              I2 * K1 * J3 + 
              I3 * J1 * K2 - 
              I3 * K1 * J2,

	0 =\= I1 * Z0 + J1 * Z4 + K1 * Z8 + L1 * Z12,
	
	between(-100,100,det),
    det is I1 * Z0 + J1 * Z4 + K1 * Z8 + L1 * Z12
    .
