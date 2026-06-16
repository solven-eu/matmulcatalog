% execute with: 
% - swipl (to start prolog)
% - ['/Users/benoit/workspace/strassen/src/main/resources/strassen.prolog'].
% - time(strassen(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,K3,K4,L1,L2,L3,L4,A1,A2,A3,A4,B1,B2,B3,B4,C1,C2,C3,C4,D1,D2,D3,D4,E1,E2,E3,E4,F1,F2,F3,F4,G1,G2,G3,G4,H1,H2,H3,H4)).
% - time(strassen(I1,I2,I3,I4,J1,J2,J3,J4,K1,K2,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 )).


strassen(I1,I2):-
	between(-1,1,I2),
	between(-1,1,I1),
	I1 is I1.
