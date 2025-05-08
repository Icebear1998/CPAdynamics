%% Define Symbolic Variables
syms R0 R1 RE0 RE1 REL1 real

% Symbolic expressions for free enzyme (E_f) and free ligand (L_f)
syms E_f L_f real

%% 1. Define Model Parameters
timescale = 65;
L_a = 100;  % segment length in bp

% Rate parameters scaled by the timescale
P.k_in    = 2/timescale;           % Initiation rate
P.k_c     = 0.8/timescale;         % Cleavage rate
kE_on_min = 0.00001/timescale;     % Minimum E binding rate
kE_on_max = 0.00025/timescale;     % Maximum E binding rate
P.kE_on = 0.00025/timescale;
P.kE_off  = 10/timescale;          % E unbinding rate
P.kL_on   = 0.00025/timescale;     % L binding rate
P.kL_off  = 0.001/timescale;       % L unbinding rate
P.k_e     = 65/timescale/L_a;      % Elongation rate before PAS
P.k_e2    = 30/timescale/L_a;      % Elongation rate after PAS
P.E_total = 70000;                 % Total E molecules
P.L_total = 100000;                % Total L molecules

% Gene and PAS specifications
geneLength_bp = 25000;             % Total gene length in bp
PASposition   = 20000;             % PAS position in bp
N = floor(geneLength_bp / L_a);    % Number of nodes along the gene
PAS = floor(PASposition / L_a);    % Node index corresponding to PAS
Pol_total = 70000;                 % Total Pol II number

%% Express Free E_f and L_f in Terms of Other Variables
E_f_expr = P.E_total - (RE0 + RE1 + REL1);
L_f_expr = P.L_total - REL1;

%% Define Steady-State Equations
eq1 = P.k_in - P.k_e*R0 - P.kE_on*E_f_expr*R0 + P.kE_off*RE0 == 0;
eq2 = - P.k_e*RE0 + P.kE_on*E_f_expr*R0 - P.kE_off*RE0 == 0;
eq3 = P.k_e*R0 - P.k_e*R1 - P.kE_on*E_f_expr*R1 + P.kE_off*RE1 == 0;
eq4 = P.k_e*RE0 - P.k_e*RE1 + P.kE_on*E_f_expr*R1 - P.kE_off*RE1 ...
    - P.kL_on*L_f_expr*RE1 + P.kL_off*REL1 == 0;
eq5 = P.kL_on*L_f_expr*RE1 - P.kL_off*REL1 - P.k_c*REL1 == 0;

%% Solve for Steady State Using vpasolve
eqs = [eq1, eq2, eq3, eq4, eq5];
vars = [R0, R1, RE0, RE1, REL1];
initial_guess = [10, 10, 10, 10, 10];  % Adjusted as necessary

sol = solve(eqs, vars);

%% Display Solutions
if isempty(sol)
    disp('No solution found. Try adjusting initial guesses or checking parameter values.');
else
    disp('Steady-state solutions:');
    disp(sol);
end
