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
	is_one(I1,A1,E1,J1,A2,E2,K1,A3,E3,L1,A4,E4),
	
	is_zero(I2,A1,E1,J2,A2,E2,K2,A3,E3,L2,A4,E4),
	is_zero(I3,A1,E1,J3,A2,E2,K3,A3,E3,L3,A4,E4),
	is_zero(I4,A1,E1,J4,A2,E2,K4,A3,E3,L4,A4,E4),
	
	is_zero(I1,A1,F1,J1,A2,F2,K1,A3,F3,L1,A4,F4),
	is_zero(I2,A1,F1,J2,A2,F2,K2,A3,F3,L2,A4,F4),
	is_zero(I3,A1,F1,J3,A2,F2,K3,A3,F3,L3,A4,F4),
	is_zero(I4,A1,F1,J4,A2,F2,K4,A3,F3,L4,A4,F4),
	
	is_zero(I1,A1,G1,J1,A2,G2,K1,A3,G3,L1,A4,G4),
	is_zero(I2,A1,G1,J2,A2,G2,K2,A3,G3,L2,A4,G4),
	is_one(I3,A1,G1,J3,A2,G2,K3,A3,G3,L3,A4,G4),
	is_zero(I4,A1,G1,J4,A2,G2,K4,A3,G3,L4,A4,G4),
		
	is_zero(I1,A1,H1,J1,A2,H2,K1,A3,H3,L1,A4,H4),
	is_zero(I2,A1,H1,J2,A2,H2,K2,A3,H3,L2,A4,H4),
	is_zero(I3,A1,H1,J3,A2,H2,K3,A3,H3,L3,A4,H4),
	is_zero(I4,A1,H1,J4,A2,H2,K4,A3,H3,L4,A4,H4),
		
	is_zero(I1,B1,E1,J1,B2,E2,K1,B3,E3,L1,B4,E4),
	is_one(I2,B1,E1,J2,B2,E2,K2,B3,E3,L2,B4,E4),
	is_zero(I3,B1,E1,J3,B2,E2,K3,B3,E3,L3,B4,E4),
	is_zero(I4,B1,E1,J4,B2,E2,K4,B3,E3,L4,B4,E4),
	
	is_zero(I1,B1,F1,J1,B2,F2,K1,B3,F3,L1,B4,F4),
	is_zero(I2,B1,F1,J2,B2,F2,K2,B3,F3,L2,B4,F4),
	is_zero(I3,B1,F1,J3,B2,F2,K3,B3,F3,L3,B4,F4),
	is_zero(I4,B1,F1,J4,B2,F2,K4,B3,F3,L4,B4,F4),
	
	is_zero(I1,B1,G1,J1,B2,G2,K1,B3,G3,L1,B4,G4),
	is_zero(I2,B1,G1,J2,B2,G2,K2,B3,G3,L2,B4,G4),
	is_zero(I3,B1,G1,J3,B2,G2,K3,B3,G3,L3,B4,G4),
	is_one(I4,B1,G1,J4,B2,G2,K4,B3,G3,L4,B4,G4),
		
	is_zero(I1,B1,H1,J1,B2,H2,K1,B3,H3,L1,B4,H4),
	is_zero(I2,B1,H1,J2,B2,H2,K2,B3,H3,L2,B4,H4),
	is_zero(I3,B1,H1,J3,B2,H2,K3,B3,H3,L3,B4,H4),
	is_zero(I4,B1,H1,J4,B2,H2,K4,B3,H3,L4,B4,H4),
	
	
	is_zero(I1,C1,E1,J1,C2,E2,K1,C3,E3,L1,C4,E4),
	is_zero(I2,C1,E1,J2,C2,E2,K2,C3,E3,L2,C4,E4),
	is_zero(I3,C1,E1,J3,C2,E2,K3,C3,E3,L3,C4,E4),
	is_zero(I4,C1,E1,J4,C2,E2,K4,C3,E3,L4,C4,E4),
	
	is_one(I1,C1,F1,J1,C2,F2,K1,C3,F3,L1,C4,F4),
	is_zero(I2,C1,F1,J2,C2,F2,K2,C3,F3,L2,C4,F4),
	is_zero(I3,C1,F1,J3,C2,F2,K3,C3,F3,L3,C4,F4),
	is_zero(I4,C1,F1,J4,C2,F2,K4,C3,F3,L4,C4,F4),
		
	is_zero(I1,C1,G1,J1,C2,G2,K1,C3,G3,L1,C4,G4),
	is_zero(I2,C1,G1,J2,C2,G2,K2,C3,G3,L2,C4,G4),
	is_zero(I3,C1,G1,J3,C2,G2,K3,C3,G3,L3,C4,G4),
	is_zero(I4,C1,G1,J4,C2,G2,K4,C3,G3,L4,C4,G4),
	
	is_zero(I1,C1,H1,J1,C2,H2,K1,C3,H3,L1,C4,H4),
	is_zero(I2,C1,H1,J2,C2,H2,K2,C3,H3,L2,C4,H4),
	is_one(I3,C1,H1,J3,C2,H2,K3,C3,H3,L3,C4,H4),
	is_zero(I4,C1,H1,J4,C2,H2,K4,C3,H3,L4,C4,H4),
		
		
	is_zero(I1,D1,E1,J1,D2,E2,K1,D3,E3,L1,D4,E4),
	is_zero(I2,D1,E1,J2,D2,E2,K2,D3,E3,L2,D4,E4),
	is_zero(I3,D1,E1,J3,D2,E2,K3,D3,E3,L3,D4,E4),
	is_zero(I4,D1,E1,J4,D2,E2,K4,D3,E3,L4,D4,E4),
	
	is_zero(I1,D1,F1,J1,D2,F2,K1,D3,F3,L1,D4,F4),
	is_one(I2,D1,F1,J2,D2,F2,K2,D3,F3,L2,D4,F4),
	is_zero(I3,D1,F1,J3,D2,F2,K3,D3,F3,L3,D4,F4),
	is_zero(I4,D1,F1,J4,D2,F2,K4,D3,F3,L4,D4,F4),
		
	is_zero(I1,D1,G1,J1,D2,G2,K1,D3,G3,L1,D4,G4),
	is_zero(I2,D1,G1,J2,D2,G2,K2,D3,G3,L2,D4,G4),
	is_zero(I3,D1,G1,J3,D2,G2,K3,D3,G3,L3,D4,G4),
	is_zero(I4,D1,G1,J4,D2,G2,K4,D3,G3,L4,D4,G4),
	
	is_zero(I1,D1,H1,J1,D2,H2,K1,D3,H3,L1,D4,H4),
	is_zero(I2,D1,H1,J2,D2,H2,K2,D3,H3,L2,D4,H4),
	is_zero(I3,D1,H1,J3,D2,H2,K3,D3,H3,L3,D4,H4),
	is_one(I4,D1,H1,J4,D2,H2,K4,D3,H3,L4,D4,H4),
	
I1 is I1.

	
is_one(I1,A1,E1,J1,A2,E2,K1,A3,E3,L1,A4,E4):-
	member(I1, [0, -1, 1]),
	member(J1, [0, -1, 1]),
	member(K1, [0, -1, 1]),
	member(L1, [0, -1, 1]),
	
	member(A1, [0, -1, 1]),
	member(A2, [0, -1, 1]),
	member(A3, [0, -1, 1]),
	member(A4, [0, -1, 1]),
	
	member(E1, [0, -1, 1]),
	member(E2, [0, -1, 1]),
	member(E3, [0, -1, 1]),
	member(E4, [0, -1, 1]),
	
	1 is I1*A1*E1+J1*A2*E2+K1*A3*E3+L1*A4*E4.
	
is_zero(I1,A1,E1,J1,A2,E2,K1,A3,E3,L1,A4,E4):-
	member(I1, [0, -1, 1]),
	member(J1, [0, -1, 1]),
	member(K1, [0, -1, 1]),
	member(L1, [0, -1, 1]),
	
	member(A1, [0, -1, 1]),
	member(A2, [0, -1, 1]),
	member(A3, [0, -1, 1]),
	member(A4, [0, -1, 1]),
	
	member(E1, [0, -1, 1]),
	member(E2, [0, -1, 1]),
	member(E3, [0, -1, 1]),
	member(E4, [0, -1, 1]),
	
	0 is I1*A1*E1+J1*A2*E2+K1*A3*E3+L1*A4*E4.
