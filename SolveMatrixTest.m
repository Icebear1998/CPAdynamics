P.kEon = 1;
P.kEoff = 2;
P.Ef = 1;
P.kHon = 0.2;
P.kHoff = 0.1;
P.kc = 0.1;
P.kPmax = 4;
n = 2;

[r_E_BeforePas, r_P, r_k_AfterPas] = compute_steady_states(P, n);

disp(vpa(r_E_BeforePas, 4));  % Display results with symbolic kP
disp(vpa(r_P, 4));