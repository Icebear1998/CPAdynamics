function P = default_parameters()
% DEFAULT_PARAMETERS  Standard parameter set for CPAdynamics model
%
% Returns a struct P containing all base simulation parameters.
% All analysis scripts should call this function to get the base
% parameter set, then override only the parameters specific to
% their analysis.
%
% Usage:
%   P = default_parameters();          % get standard set
%   P.kc = 0.2;                        % override one parameter
%
% Standard parameter set (see also Docs/claude.md):

    % --- Geometry ---
    P.L_a            = 100;        % bp per node
    P.geneLength_bp  = 25000;      % Total gene length (bp)
    P.PASposition    = 20000;      % PAS position from TSS (bp)

    % --- Pol II kinetics ---
    P.k_in   = 2;                  % Pol II initiation rate
    P.k_e    = 65 / P.L_a;        % Elongation rate (before PAS)
    P.k_e2   = 30 / P.L_a;        % Elongation rate (after PAS, in REH)

    % --- Pool sizes ---
    P.E_total   = 100000;          % Total E factor pool
    P.Pol_total = 70000;           % Total Pol II pool

    % --- E factor binding ---
    P.kEon   = 0.000001;          % E factor on-rate
    P.kEoff  = 0.2;                % E factor off-rate

    % --- PAS recognition (hexamer) ---
    P.kHon   = 1;                  % Hexamer on-rate
    P.kHoff  = 0.5;                  % Hexamer off-rate

    % --- Cleavage ---
    P.kc     = 0.13;                % Cleavage rate

    % --- Ser2P phosphorylation ---
    P.kPon_min   = 0.01;           % Min Ser2P phosphorylation rate (at TSS)
    P.kPon_slope = 0.005;          % Linear slope of kPon along gene
    P.kPoff      = 1;              % Ser2P dephosphorylation rate

end
