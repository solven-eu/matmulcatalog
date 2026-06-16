% execute with: 
% - swipl (to start prolog)
% - ['/Users/benoit/workspace/strassen/src/main/resources/strassen.fake.prolog'].
% -      strassen_fake(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,A1,A2,A3,A4,B1,B2,B3,B4,C1,C2,C3,C4,D1,D2,D3,D4,E1,E2,E3,E4,F1,F2,F3,F4,G1,G2,G3,G4,H1,H2,H3,H4).
% - time(strassen_fake(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,A1,A2,A3,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 )).

strassen_fake(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,A1,A2,A3,A4,B1,B2,B3,B4,C1,C2,C3,C4,D1,D2,D3,D4,E1,E2,E3,E4,F1,F2,F3,F4,G1,G2,G3,G4,H1,H2,H3,H4):-
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
	member(B1, [0, -1, 1]),
	member(C1, [0, -1, 1]),
	member(D1, [0, -1, 1]),
	member(E1, [0, -1, 1]),
	member(F1, [0, -1, 1]),
	member(G1, [0, -1, 1]),
	member(H1, [0, -1, 1]),
	
	member(A2, [0, -1, 1]),
	member(B2, [0, -1, 1]),
	member(C2, [0, -1, 1]),
	member(D2, [0, -1, 1]),
	member(E2, [0, -1, 1]),
	member(F2, [0, -1, 1]),
	member(G2, [0, -1, 1]),
	member(H2, [0, -1, 1]),
	
	member(A3, [0, -1, 1]),
	member(B3, [0, -1, 1]),
	member(C3, [0, -1, 1]),
	member(D3, [0, -1, 1]),
	member(E3, [0, -1, 1]),
	member(F3, [0, -1, 1]),
	member(G3, [0, -1, 1]),
	member(H3, [0, -1, 1]),
	
	member(A4, [0, -1, 1]),
	member(B4, [0, -1, 1]),
	member(C4, [0, -1, 1]),
	member(D4, [0, -1, 1]),
	member(E4, [0, -1, 1]),
	member(F4, [0, -1, 1]),
	member(G4, [0, -1, 1]),
	member(H4, [0, -1, 1]),
	
	E111_1 is I1*A1*E1,
	E111_2 is J1*A2*E2,
	E111_3 is K1*A3*E3,
	E111_4 is L1*A4*E4,
	E111 is E111_1*E111_2*E111_3*E111_4,
	E111 is 1,
	
	E114 is E114.
