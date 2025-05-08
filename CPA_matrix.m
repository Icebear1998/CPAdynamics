% Define symbolic variables
syms Kp %Kon Koff Ef

Kon = 1;
Koff = 2;
Ef = 1;

% r = [R0 R1 R1E];
K1 = [-1, Kp, 0;
      1, -(Kon*Ef + Kp), Koff;
      0, Kon*Ef, -Koff];
  



% r = [R0 R1 R1E R2 R2E R2E2]
K2 = [
    -1, Kp, 0, 0, 0, 0;                                % eq1
    1, -(Kon*Ef + Kp + 1), Koff, Kp, 0, 0;             % eq2
    0, Kon*Ef, -(Koff + 1), 0, Kp, 0;                  % eq3
    0, 1, 0, -(2*Kon*Ef + Kp), Koff, 0;                % eq4
    0, 0, 1, 2*Kon*Ef, -(Kon*Ef + Kp + Koff), 2*Koff;  % eq5
    0, 0, 0, 0, Kon*Ef, -2*Koff                        % eq6
];

% x= [R0 R1 R1E R2 R2E R2E2 R3 R3E R3E2 R3E3]
K3 = [
    -1, Kp, 0, 0, 0, 0, 0, 0, 0, 0;                                % eq1
    1, -(Kon*Ef + Kp + 1), Koff, Kp, 0, 0, 0, 0, 0, 0;             % eq2
    0, Kon*Ef, -(Koff + 1), 0, Kp, 0, 0, 0, 0, 0;                  % eq3
    0, 1, 0, -(2*Kon*Ef + Kp + 1), Koff, 0, Kp, 0, 0, 0;           % eq4
    0, 0, 1, 2*Kon*Ef, -(Kon*Ef + Kp + Koff + 1), Koff, 0, Kp, 0, 0; % eq5
    0, 0, 0, 0, Kon*Ef, -(Koff + 1), 0, 0, Kp, 0;                  % eq6
    0, 0, 0, 1, 0, 0, -(3*Kon*Ef + Kp), Koff, 0, 0;                % eq7
    0, 0, 0, 0, 1, 0, 3*Kon*Ef, -(2*Kon*Ef + Kp + Koff), Koff, 0;  % eq8
    0, 0, 0, 0, 0, 1, 0, 2*Kon*Ef, -(Kon*Ef + Kp + Koff), Koff;    % eq9
    0, 0, 0, 0, 0, 0, 0, 0, Kon*Ef, -Koff                         % eq10
];


% r = [R0,R1,R1E,R2,R2E,R2E2,R3,R3E,R3E2,R3E3,R4,R4E,R4E2,R4E3,R4E4]
K4 = [
    -1, Kp, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;                                % eq1
     1, -(Kon*Ef + Kp + 1), Koff, Kp, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;           % eq2
     0, Kon*Ef, -(Koff + 1), 0, Kp, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;                % eq3
     0, 1, 0, -(2*Kon*Ef + Kp + 1), Koff, 0, Kp, 0, 0, 0, 0, 0, 0, 0, 0;         % eq4
     0, 0, 1, 2*Kon*Ef, -(Kon*Ef + Kp + Koff + 1), Koff, 0, Kp, 0, 0, 0, 0, 0, 0, 0; % eq5
     0, 0, 0, 0, Kon*Ef, -(Koff + 1), 0, 0, Kp, 0, 0, 0, 0, 0, 0;                % eq6
     0, 0, 0, 1, 0, 0, -(3*Kon*Ef + Kp + 1), Koff, 0, 0, Kp, 0, 0, 0, 0;         % eq7
     0, 0, 0, 0, 1, 0, 3*Kon*Ef, -(2*Kon*Ef + Kp + Koff + 1), Koff, 0, 0, Kp, 0, 0, 0; % eq8
     0, 0, 0, 0, 0, 1, 0, 2*Kon*Ef, -(Kon*Ef + Kp + Koff + 1), Koff, 0, 0, Kp, 0, 0; % eq9
     0, 0, 0, 0, 0, 0, 0, 0, Kon*Ef, -(Koff + 1), 0, 0, 0, Kp, 0;                % eq10
     0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -(4*Kon*Ef + Kp), Koff, 0, 0, 0;              % eq11
     0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 4*Kon*Ef, -(3*Kon*Ef + Kp + Koff), Koff, 0, 0; % eq12
     0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 3*Kon*Ef, -(2*Kon*Ef + Kp + Koff), Koff, 0; % eq13
     0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2*Kon*Ef, -(Kon*Ef + Kp + Koff), Koff;   % eq14
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Kon*Ef, -Koff                         % eq15
];


% K*r = 0 => r = null(K)
% Compute the null space
null_space = null(K2);
disp(null_space);


