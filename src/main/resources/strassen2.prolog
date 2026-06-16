% execute with: 
% - swipl (to start prolog)
% - ['/Users/benoit/workspace/strassen/src/main/resources/strassen.prolog'].
% - strassen(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,A1,A2,A3,A4,B1,B2,B3,B4,C1,C2,C3,C4,D1,D2,D3,D4,E1,E2,E3,E4,F1,F2,F3,F4,G1,G2,G3,G4,H1,H2,H3,H4)
% - time(strassen(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 )).


strassen(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,A1,A2,A3,A4,B1,B2,B3,B4,C1,C2,C3,C4,D1,D2,D3,D4,E1,E2,E3,E4,F1,F2,F3,F4,G1,G2,G3,G4,H1,H2,H3,H4):-
	strassen_AE(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,A1,A2,A3,A4,E1,E2,E3,E4),
	strassen_AF(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,A1,A2,A3,A4,F1,F2,F3,F4),
	strassen_AG(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,A1,A2,A3,A4,G1,G2,G3,G4),
	strassen_AH(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,A1,A2,A3,A4,H1,H2,H3,H4),
	
	strassen_BE(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,B1,B2,B3,B4,E1,E2,E3,E4),
	strassen_BF(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,B1,B2,B3,B4,F1,F2,F3,F4),
	strassen_BG(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,B1,B2,B3,B4,G1,G2,G3,G4),
	strassen_BH(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,B1,B2,B3,B4,H1,H2,H3,H4),
	
	strassen_CE(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,C1,C2,C3,C4,E1,E2,E3,E4),
	strassen_CF(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,C1,C2,C3,C4,F1,F2,F3,F4),
	strassen_CG(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,C1,C2,C3,C4,G1,G2,G3,G4),
	strassen_CH(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,C1,C2,C3,C4,H1,H2,H3,H4),
	
	strassen_DE(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,D1,D2,D3,D4,E1,E2,E3,E4),
	strassen_DF(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,D1,D2,D3,D4,F1,F2,F3,F4),
	strassen_DG(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,D1,D2,D3,D4,G1,G2,G3,G4),
	strassen_DH(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,D1,D2,D3,D4,H1,H2,H3,H4),
	
	I1 is I1.

strassen_AE(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,A1,E1,A2,E2,A3,E3,A4,E4):-
	member(I1, [0, -1, 1]),
	member(I2, [0, -1, 1]),
	member(I3, [0, -1, 1]),
	member(I4, [0, -1, 1]),
	
	member(J1, [0, -1, 1]),
	member(J2, [0, -1, 1]),
	member(J3, [0, -1, 1]),
	member(J4, [0, -1, 1]),
	
	member(K1, [0, -1, 1]),
	member(K2, [0, -1, 1]),
	member(K3, [0, -1, 1]),
	member(K4, [0, -1, 1]),
	
	member(L1, [0, -1, 1]),
	member(L2, [0, -1, 1]),
	member(L3, [0, -1, 1]),
	member(L4, [0, -1, 1]),
	
	member(A1, [0, -1, 1]),
	member(A2, [0, -1, 1]),
	member(A3, [0, -1, 1]),
	member(A4, [0, -1, 1]),
	
	member(E1, [0, -1, 1]),
	member(E2, [0, -1, 1]),
	member(E3, [0, -1, 1]),
	member(E4, [0, -1, 1]),
	
	1 is I1*A1*E1+J1*A2*E2+K1*A3*E3+L1*A4*E4,
	0 is I2*A1*E1+J2*A2*E2+K2*A3*E3+L2*A4*E4,
	0 is I3*A1*E1+J3*A2*E2+K3*A3*E3+L3*A4*E4,
	0 is I4*A1*E1+J4*A2*E2+K4*A3*E3+L4*A4*E4,
	
	I1 is I1.

	
strassen_AF(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,A1,A2,A3,A4,F1,F2,F3,F4):-
	member(I1, [0, -1, 1]),
	member(I2, [0, -1, 1]),
	member(I3, [0, -1, 1]),
	member(I4, [0, -1, 1]),
	
	member(J1, [0, -1, 1]),
	member(J2, [0, -1, 1]),
	member(J3, [0, -1, 1]),
	member(J4, [0, -1, 1]),
	
	member(K1, [0, -1, 1]),
	member(K2, [0, -1, 1]),
	member(K3, [0, -1, 1]),
	member(K4, [0, -1, 1]),
	
	member(L1, [0, -1, 1]),
	member(L2, [0, -1, 1]),
	member(L3, [0, -1, 1]),
	member(L4, [0, -1, 1]),

	member(A1, [0, -1, 1]),
	member(A2, [0, -1, 1]),
	member(A3, [0, -1, 1]),
	member(A4, [0, -1, 1]),
	member(F1, [0, -1, 1]),
	member(F2, [0, -1, 1]),
	member(F3, [0, -1, 1]),
	member(F4, [0, -1, 1]),
	
	0 is I1*A1*F1+J1*A2*F2+K1*A3*F3+L1*A4*F4,
	0 is I2*A1*F1+J2*A2*F2+K2*A3*F3+L2*A4*F4,
	0 is I3*A1*F1+J3*A2*F2+K3*A3*F3+L3*A4*F4,
	0 is I4*A1*F1+J4*A2*F2+K4*A3*F3+L4*A4*F4,
	
	I1 is I1.
	

strassen_AF(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,A1,A2,A3,A4,F1,F2,F3,F4):-
I1 is I1.
strassen_AG(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,A1,A2,A3,A4,G1,G2,G3,G4):-I1 is I1.
strassen_AH(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,A1,A2,A3,A4,H1,H2,H3,H4):-I1 is I1.

strassen_BE(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,B1,B2,B3,B4,E1,E2,E3,E4):-I1 is I1.
strassen_BF(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,B1,B2,B3,B4,F1,F2,F3,F4):-I1 is I1.
strassen_BG(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,B1,B2,B3,B4,G1,G2,G3,G4):-I1 is I1.
strassen_BH(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,B1,B2,B3,B4,H1,H2,H3,H4):-I1 is I1.

strassen_CE(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,C1,C2,C3,C4,E1,E2,E3,E4):-I1 is I1.
strassen_CF(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,C1,C2,C3,C4,F1,F2,F3,F4):-I1 is I1.
strassen_CG(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,C1,C2,C3,C4,G1,G2,G3,G4):-I1 is I1.
strassen_CH(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,C1,C2,C3,C4,H1,H2,H3,H4):-
I1 is I1.

strassen_DE(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,D1,D2,D3,D4,E1,E2,E3,E4):-
I1 is I1.
strassen_DF(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,D1,D2,D3,D4,F1,F2,F3,F4):-
I1 is I1.
strassen_DG(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,D1,D2,D3,D4,G1,G2,G3,G4):-
I1 is I1.
strassen_DH(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,D1,D2,D3,D4,H1,H2,H3,H4):-
	member(I1, [0, -1, 1]),
	member(I2, [0, -1, 1]),
	member(I3, [0, -1, 1]),
	member(I4, [0, -1, 1]),
	
	member(J1, [0, -1, 1]),
	member(J2, [0, -1, 1]),
	member(J3, [0, -1, 1]),
	member(J4, [0, -1, 1]),
	
	member(K1, [0, -1, 1]),
	member(K2, [0, -1, 1]),
	member(K3, [0, -1, 1]),
	member(K4, [0, -1, 1]),
	
	member(L1, [0, -1, 1]),
	member(L2, [0, -1, 1]),
	member(L3, [0, -1, 1]),
	member(L4, [0, -1, 1]),
	
	member(D1, [0, -1, 1]),
	member(D2, [0, -1, 1]),
	member(D3, [0, -1, 1]),
	member(D4, [0, -1, 1]),
	
	member(H1, [0, -1, 1]),
	member(H2, [0, -1, 1]),
	member(H3, [0, -1, 1]),
	member(H4, [0, -1, 1]),

	0 is I1*D1*H1+J1*D2*H2+K1*D3*H3+L1*D4*H4,
	0 is I2*D1*H1+J2*D2*H2+K2*D3*H3+L2*D4*H4,
	0 is I3*D1*H1+J3*D2*H2+K3*D3*H3+L3*D4*H4,
	1 is I4*D1*H1+J4*D2*H2+K4*D3*H3+L4*D4*H4,
	
	!,
I1 is I1.
