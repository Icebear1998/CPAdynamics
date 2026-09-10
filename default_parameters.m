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
% Source: Appendix_Parameter_Estimation_revised.docx, sections B-G.
% Rates use the rounded working values reported in the appendix. Geometry
% and initiation are retained where that document supplies no new baseline.


    % --- Geometry ---
    P.L_a            = 100;        % bp per node
    P.geneLength_bp  = 25000;      % Total gene length (bp)
    P.PASposition    = 20000;      % PAS position from TSS (bp)

    % --- Pol II kinetics ---
    P.k_in   = 2;                  % Pol II initiation rate
    P.k_e    = (4000 / 60) / P.L_a; % 4 kb/min before PAS (section E)
    P.k_e2   = 30 / P.L_a;        % Elongation rate (after PAS, in REH)

    % --- Pool sizes ---
    P.E_total   = 100000;          % Total E factor pool
    P.Pol_total = 70000;           % Total Pol II pool

    % --- E factor binding ---
    P.kEon   = 2.22e-6;            % molecule^-1 s^-1; 374 um^3 nucleus (C)
    P.kEoff  = 0.5;                % E factor off-rate
 
    % --- PAS recognition (hexamer) ---
    P.kHon   = 7.04;               % s^-1; effective PAS encounter rate (B.2)
    P.kHoff  = 0.5;                % s^-1; provisional working default
    % Section B.3 gives [0.05, 1.0] s^-1 but no baseline; 0.5 is a
    % working choice within that interval, not a quoted appendix estimate.
 
    % --- Cleavage ---
    P.kc     = 0.439;              % s^-1; sequential commitment estimate (D)

    % --- Ser2P phosphorylation ---
    P.kPon_min   = 0.01;           % Min Ser2P phosphorylation rate (at TSS)
    P.kPon_slope = 0.005;          % Linear slope of kPon along gene
    P.kPoff      = 1;              % Ser2P dephosphorylation rate

    % --- Sensitivity ranges [minimum, maximum] from the revised appendix ---
    % These are working sensitivity intervals, not confidence intervals.
    % No ranges are supplied for k_in, geometry, kPon_min, kPoff or
    % kEoff_engaged; no bounds are invented for those parameters here.
    P.ranges.kHon       = [4.02, 17.25];      % s^-1 (B.2)
    P.ranges.kHoff      = [0.05, 1.0];        % s^-1 (B.3)
    P.ranges.kEon       = [1.66e-6, 4.15e-6]; % molecule^-1 s^-1 (C)
    P.ranges.kEoff      = [0.5, 5];           % s^-1 (C)
    P.ranges.kc         = [0.0077, 1.06];     % s^-1 (D)
    P.ranges.k_e        = [3000, 5000] / 60 / P.L_a; % 3-5 kb/min (E)
    P.ranges.k_e2       = [10, 50] / P.L_a;  % bp/s -> node/s (E)
    P.ranges.kPon_slope = [0.001, 0.01];     % s^-1 per node (F)
    P.ranges.E_total    = [50000, 200000];   % molecules (G)
    P.ranges.Pol_total  = [25000, 140000];   % nuclear molecules (G)

    % Derived dissociation/association ratio envelopes for sensitivity maps.
    % These combine opposite endpoints of the individual on/off intervals.
    P.ranges.kHd = [P.ranges.kHoff(1)/P.ranges.kHon(2), ...
                    P.ranges.kHoff(2)/P.ranges.kHon(1)]; % dimensionless
    P.ranges.kEd = [P.ranges.kEoff(1)/P.ranges.kEon(2), ...
                    P.ranges.kEoff(2)/P.ranges.kEon(1)]; % molecules

end
