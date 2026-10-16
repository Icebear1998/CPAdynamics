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
%   P.ranges.kHon                     % appendix sensitivity interval [min max]
% Baseline: the MATLAB/Python comparison parameter set supplied 2026-09-15.
% Sensitivity intervals below retain the revised appendix ranges.

    % --- Geometry ---
    P.L_a            = 100;        % bp per node
    P.geneLength_bp  = 25000;      % Total gene length (bp)
    P.PASposition    = 20000;      % PAS position from TSS (bp)

    % --- Pol II kinetics ---
    P.k_in   = 2;                  % Pol II initiation rate
    P.k_e    = 65 / P.L_a;        % 65 bp/s for unrecognized Pol II
    P.k_e2   = 30 / P.L_a;        % Elongation rate (after PAS, in REH)

    % --- Pool sizes ---
    P.E_total   = 100000;          % Total E factor pool
    P.Pol_total = 100000;           % Total Pol II pool

    % --- E factor binding ---
    P.kEon   = 2.5e-6;             % molecule^-1 s^-1
    P.kEoff  = 0.5;                % E factor off-rate
    % Engaged-E CTD detachment. Conservative estimate between a firm hold
    % by H and none (= kEoff); one value for all PAS strengths (appendix C).
    P.kEoff_engaged = 0.1;         % s^-1
 
    % --- PAS recognition (hexamer) ---
    P.kHon   = 7;                  % s^-1 per bound E
    P.kHoff  = 1;                  % s^-1; weak PAS (Kd ~2 uM, appendix B)
 
    % --- Cleavage ---
    P.kc     = 0.14;               % s^-1

    % --- Ser2P phosphorylation ---
    % One site p = two adjacent Ser2P repeats (appendix F). kPoff = 0.01
    % (mark half-life ~2.3 min) sits just above the inhibitor-based estimate
    % (0.0016-0.006), whose upper end is soft (appendix F).
    P.kPon_min   = 0;              % s^-1; no pairs at the TSS (pairs need two marks)
    P.kPon_slope = 5e-4;           % s^-1 per node; 2*k_dep^2/(52*k_e), k_dep ~0.09 s^-1
    P.kPoff      = 0.01;           % s^-1 per free site (appendix F)

    % --- Sensitivity ranges [minimum, maximum] from the revised appendix ---
    % These are working sensitivity intervals, not confidence intervals.
    % No ranges are supplied for k_in or geometry; no bounds are invented
    % for those parameters here.
    P.ranges.kHon       = [4, 17];      % s^-1 (B.2)
    P.ranges.kHoff      = [0.0015, 1];     % s^-1 (B); AAUAAA Kd 3 nM to ~2 uM (Hamilton 2019; Huang 2025)
    P.ranges.kEon       = [1.66e-6, 4.15e-6]; % molecule^-1 s^-1 (C)
    P.ranges.kEoff      = [0.1, 5];           % s^-1 (C); avidity / in vivo residence to single CID
    P.ranges.kEoff_engaged = [0.01, 0.5]; % s^-1 (C); firm hold by H to no hold (kEoff)
    P.ranges.kc         = [0.015, 0.6];    % s^-1 (D); in vivo cleavage timing to holo-CstF affinity
    P.ranges.k_e        = [3000, 5000] / 65 / P.L_a; % ~46-77 bp/s (~2.8-4.6 kb/min) (E)
    P.ranges.k_e2       = [10, 50] / P.L_a;  % bp/s -> node/s (E)
    P.ranges.kPon_slope = [1e-4, 2e-3];      % s^-1 per node (F); E-site 50% at ~20 to ~4 kb
    P.ranges.kPoff      = [0.0016, 0.02];  % s^-1 (F); inhibitor lower bound to 2x the working value
    P.ranges.E_total    = [50000, 200000];   % molecules (G)
    P.ranges.Pol_total  = [25000, 140000];   % nuclear molecules (G)

    % Derived dissociation/association ratio envelopes for sensitivity maps.
    % These combine opposite endpoints of the individual on/off intervals.
    P.ranges.kHd = [P.ranges.kHoff(1)/P.ranges.kHon(2), ...
                    P.ranges.kHoff(2)/P.ranges.kHon(1)]; % dimensionless
    P.ranges.kEd = [P.ranges.kEoff(1)/P.ranges.kEon(2), ...
                    P.ranges.kEoff(2)/P.ranges.kEon(1)]; % molecules

end
