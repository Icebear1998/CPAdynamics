% % Here are the equation codes for different Es:
% % Equations for 1 E:
% eq1 = Kp*R1 - R0 == 0;
% eq2 = R0 + Koff*R1E - (Kon*Ef + Kp)*R1 == 0;
% eq3 = Kon*Ef*R1 - Koff*R1E == 0;
% 
% % Equations for two Es:
% eq1 = Kp*R1 - R0 == 0;
% eq2 = R0 + Koff*R1E + Kp*R2 - (Kon*Ef + Kp + 1)*R1 == 0;
% eq3 = Kon*Ef*R1 + Kp*R2E - (Koff + 1)*R1E == 0;
% eq4 = R1 + Koff*R2E - (2*Kon*Ef + Kp)*R2 == 0;
% eq5 = 2*Kon*Ef*R2 + R1E + 2*Koff*R2E2 - (Kon*Ef + Kp + Koff)*R2E == 0;
% eq6 = Kon*Ef*R2E - 2*Koff*R2E2 == 0;
% 
% % Equations for three Es:
% eq1 = Kp*R1 - R0 == 0;
% eq2 = R0 + Koff*R1E + Kp*R2 - (Kon*Ef + Kp + 1)*R1 == 0;
% eq3 = Kon*Ef*R1 + Kp*R2E - (Koff + 1)*R1E == 0;
% eq4 = R1 + Koff*R2E + Kp*R3 - (2*Kon*Ef + Kp + 1)*R2 == 0;
% eq5 = 2*Kon*Ef*R2 + R1E + Koff*R2E2 + Kp*R3E - (Kon*Ef + Kp + Koff + 1)*R2E == 0;
% eq6 = Kon*Ef*R2E + Kp*R3E2 - (Koff + 1)*R2E2 == 0;
% eq7 = R2 + Koff*R3E - (3*Kon*Ef + Kp)*R3 == 0;
% eq8 = 3*Kon*Ef*R3 + R2E + Koff*R3E2 - (2*Kon*Ef + Kp + Koff)*R3E == 0;
% eq9 = 2*Kon*Ef*R3E + R2E2 + Koff*R3E3 - (Kon*Ef + Kp + Koff)*R3E2 == 0;
% eq10 = Kon*Ef*R3E2 - Koff*R3E3 == 0;
% 
% % Equations for four Es:
% eq1 = Kp*R1 - R0 == 0;
% eq2 = R0 + Koff*R1E + Kp*R2 - (Kon*Ef + Kp + 1)*R1 == 0;
% eq3 = Kon*Ef*R1 + Kp*R2E - (Koff + 1)*R1E == 0;
% eq4 = R1 + Koff*R2E + Kp*R3 - (2*Kon*Ef + Kp + 1)*R2 == 0;
% eq5 = 2*Kon*Ef*R2 + R1E + Koff*R2E2 + Kp*R3E - (Kon*Ef + Kp + Koff + 1)*R2E == 0;
% eq6 = Kon*Ef*R2E + Kp*R3E2 - (Koff + 1)*R2E2 == 0;
% eq7 = R2 + Koff*R3E + Kp*R4 - (3*Kon*Ef + Kp + 1)*R3 == 0;
% eq8 = 3*Kon*Ef*R3 + R2E + Koff*R3E2 + Kp*R4E - (2*Kon*Ef + Kp + Koff + 1)*R3E == 0;
% eq9 = 2*Kon*Ef*R3E + R2E2 + Koff*R3E3 + Kp*R4E2 - (Kon*Ef + Kp + Koff + 1)*R3E2 == 0;
% eq10 = Kon*Ef*R3E2 + Kp*R4E3 - (Koff + 1)*R3E3 == 0;
% eq11 = R3 + Koff*R4E - (4*Kon*Ef + Kp)*R4 == 0;
% eq12 = 4*Kon*Ef*R4 + R3E + Koff*R4E2 - (3*Kon*Ef + Kp + Koff)*R4E == 0;
% eq13 = 3*Kon*Ef*R4E + R3E2 + Koff*R4E3 - (2*Kon*Ef + Kp + Koff)*R4E2 == 0;
% eq14 = 2*Kon*Ef*R4E2 + R3E3 + Koff*R4E4 - (Kon*Ef + Kp + Koff)*R4E3 == 0;
% eq15 = Kon*Ef*R4E3 - Koff*R4E4 == 0;
