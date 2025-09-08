function save_analysis_results(analysis_type, data, parameters, varargin)
% SAVE_ANALYSIS_RESULTS - Utility function to save analysis results
%
% Usage:
%   save_analysis_results(analysis_type, data, parameters, 'Name', Value, ...)
%
% Inputs:
%   analysis_type - String identifier for the analysis (e.g., 'PASUsage', 'ParameterSweep')
%   data          - Structure containing the results data
%   parameters    - Structure containing the simulation parameters
%   
% Optional Name-Value pairs:
%   'SavePlot'    - Boolean, whether to save the current figure (default: true)
%   'PlotFormat'  - String, format for plot ('png', 'fig', 'eps') (default: 'png')
%   'CustomName'  - String, custom filename prefix (default: auto-generated)
%   'ExtraInfo'   - String, additional information for filename
%
% Examples:
%   % For PASUsageAnalysis:
%   data.results_matrix = results_matrix;
%   data.x_values = kHoff_values;
%   data.y_values = E_total_values;
%   data.x_label = 'kHoff';
%   data.y_label = 'E_total';
%   data.metric_distance = fixed_distance_bp;
%   save_analysis_results('PASUsage', data, P);
%
%   % For SweepParameterPASusage:
%   data.results_matrix = proximal_usage_results_cdf;
%   data.x_values = inter_pas_distances_bp;
%   data.sweep_values = sweep_param_values;
%   data.sweep_param = sweep_param_name;
%   save_analysis_results('ParameterSweep', data, P, 'ExtraInfo', 'EBinding5');

% Parse optional inputs
p = inputParser;
addParameter(p, 'SavePlot', true, @islogical);
addParameter(p, 'PlotFormat', 'png', @ischar);
addParameter(p, 'CustomName', '', @ischar);
addParameter(p, 'ExtraInfo', '', @ischar);
parse(p, varargin{:});

save_plot = p.Results.SavePlot;
plot_format = p.Results.PlotFormat;
custom_name = p.Results.CustomName;
extra_info = p.Results.ExtraInfo;

% Create output directory
output_dir = sprintf('SecondVersionResults/%s/', analysis_type);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Generate filename
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

if ~isempty(custom_name)
    base_filename = custom_name;
else
    switch analysis_type
        case 'PASUsage'
            base_filename = sprintf('PASUsage_kHoff%.0e-%.0e_Etotal%d-%d_%s', ...
                min(data.x_values), max(data.x_values), ...
                min(data.y_values), max(data.y_values), timestamp);
        case 'ParameterSweep'
            base_filename = sprintf('Sweep_%s_%.3g-%.3g_%s', ...
                data.sweep_param, min(data.sweep_values), max(data.sweep_values), timestamp);
        case 'CPA_multipleE_main'
            base_filename = sprintf('CPA_main_EBinding%d_%s', ...
                data.EBindingNumber, timestamp);
        case 'parameter_sweep_1D'
            base_filename = sprintf('Sweep1D_%s_EBinding%d_%s', ...
                data.sweep_param, data.EBindingNumber, timestamp);
        case 'parameter_sweep_2D'
            base_filename = sprintf('Sweep2D_%s_%s_EBinding%d_%s', ...
                data.param1, data.param2, data.EBindingNumber, timestamp);
        case 'PASUsagevsInterPASDistance'
            base_filename = sprintf('PASvsDistance_EBinding%d_%s', ...
                data.EBindingNumber, timestamp);
        otherwise
            base_filename = sprintf('%s_%s', analysis_type, timestamp);
    end
end

% Add extra info if provided
if ~isempty(extra_info)
    base_filename = sprintf('%s_%s', base_filename, extra_info);
end

% Save the plot if requested
if save_plot && ~isempty(get(0, 'CurrentFigure'))
    plot_filename = fullfile(output_dir, [base_filename '.' plot_format]);
    switch plot_format
        case 'png'
            saveas(gcf, plot_filename);
        case 'fig'
            savefig(gcf, plot_filename);
        case 'eps'
            print(gcf, plot_filename, '-depsc');
        otherwise
            saveas(gcf, plot_filename);
    end
    fprintf('Plot saved to: %s\n', plot_filename);
end

% Save raw data as text file
data_filename = fullfile(output_dir, [base_filename '_data.txt']);
fid = fopen(data_filename, 'w');

% Write header with metadata
fprintf(fid, '%% %s Analysis Results\n', analysis_type);
fprintf(fid, '%% Generated on: %s\n', datestr(now));

% Write analysis-specific information
switch analysis_type
    case 'PASUsage'
        write_pas_usage_data(fid, data, parameters);
    case 'ParameterSweep'
        write_parameter_sweep_data(fid, data, parameters);
    case 'CPA_multipleE_main'
        write_cpa_main_data(fid, data, parameters);
    case 'parameter_sweep_1D'
        write_sweep_1d_data(fid, data, parameters);
    case 'parameter_sweep_2D'
        write_sweep_2d_data(fid, data, parameters);
    case 'PASUsagevsInterPASDistance'
        write_pas_distance_data(fid, data, parameters);
    otherwise
        write_generic_data(fid, data, parameters);
end

fclose(fid);
fprintf('Data saved to: %s\n', data_filename);

end

%% Helper functions for different analysis types

function write_pas_usage_data(fid, data, P)
    fprintf(fid, '%% Fixed distance for metric: %d bp\n', data.metric_distance);
    fprintf(fid, '%% \n');
    write_base_parameters(fid, P);
    fprintf(fid, '%% \n');
    fprintf(fid, '%% %s_values (%d values):\n', data.x_label, length(data.x_values));
    fprintf(fid, '%g ', data.x_values);
    fprintf(fid, '\n');
    fprintf(fid, '%% \n');
    fprintf(fid, '%% %s_values (%d values):\n', data.y_label, length(data.y_values));
    fprintf(fid, '%g ', data.y_values);
    fprintf(fid, '\n');
    fprintf(fid, '%% \n');
    fprintf(fid, '%% Results Matrix (rows = %s, columns = %s):\n', data.y_label, data.x_label);
    fprintf(fid, '%% Cumulative Exit Percentage at %d bp\n', data.metric_distance);
    
    % Write the results matrix
    for i = 1:size(data.results_matrix, 1)
        fprintf(fid, '%.6f ', data.results_matrix(i, :));
        fprintf(fid, '\n');
    end
end

function write_parameter_sweep_data(fid, data, P)
    fprintf(fid, '%% Swept parameter: %s\n', data.sweep_param);
    fprintf(fid, '%% \n');
    write_base_parameters(fid, P);
    fprintf(fid, '%% \n');
    fprintf(fid, '%% Distance values (%d values):\n', length(data.x_values));
    fprintf(fid, '%g ', data.x_values);
    fprintf(fid, '\n');
    fprintf(fid, '%% \n');
    fprintf(fid, '%% %s values (%d values):\n', data.sweep_param, length(data.sweep_values));
    fprintf(fid, '%g ', data.sweep_values);
    fprintf(fid, '\n');
    fprintf(fid, '%% \n');
    fprintf(fid, '%% Results Matrix (rows = distance, columns = %s):\n', data.sweep_param);
    fprintf(fid, '%% Cumulative Exit Percentage\n');
    
    % Write the results matrix
    for i = 1:size(data.results_matrix, 1)
        fprintf(fid, '%.6f ', data.results_matrix(i, :));
        fprintf(fid, '\n');
    end
end

function write_generic_data(fid, data, P)
    write_base_parameters(fid, P);
    fprintf(fid, '%% \n');
    fprintf(fid, '%% Results Data:\n');
    
    % Write any numeric matrices in the data structure
    fields = fieldnames(data);
    for i = 1:length(fields)
        field = fields{i};
        value = data.(field);
        if isnumeric(value) && ~isscalar(value)
            fprintf(fid, '%% %s:\n', field);
            if isvector(value)
                fprintf(fid, '%g ', value);
                fprintf(fid, '\n');
            else
                for j = 1:size(value, 1)
                    fprintf(fid, '%.6f ', value(j, :));
                    fprintf(fid, '\n');
                end
            end
            fprintf(fid, '\n');
        end
    end
end

function write_base_parameters(fid, P)
    fprintf(fid, '%% Base Parameters:\n');
    if isfield(P, 'L_a'), fprintf(fid, '%% L_a = %g\n', P.L_a); end
    if isfield(P, 'k_in'), fprintf(fid, '%% k_in = %g\n', P.k_in); end
    if isfield(P, 'kEon'), fprintf(fid, '%% kEon = %g\n', P.kEon); end
    if isfield(P, 'kEoff'), fprintf(fid, '%% kEoff = %g\n', P.kEoff); end
    if isfield(P, 'k_e'), fprintf(fid, '%% k_e = %g\n', P.k_e); end
    if isfield(P, 'k_e2'), fprintf(fid, '%% k_e2 = %g\n', P.k_e2); end
    if isfield(P, 'E_total'), fprintf(fid, '%% E_total = %g\n', P.E_total); end
    if isfield(P, 'L_total'), fprintf(fid, '%% L_total = %g\n', P.L_total); end
    if isfield(P, 'Pol_total'), fprintf(fid, '%% Pol_total = %g\n', P.Pol_total); end
    if isfield(P, 'kHon'), fprintf(fid, '%% kHon = %g\n', P.kHon); end
    if isfield(P, 'kHoff'), fprintf(fid, '%% kHoff = %g\n', P.kHoff); end
    if isfield(P, 'kc'), fprintf(fid, '%% kc = %g\n', P.kc); end
    if isfield(P, 'kPon_min'), fprintf(fid, '%% kPon_min = %g\n', P.kPon_min); end
    if isfield(P, 'kPon_max'), fprintf(fid, '%% kPon_max = %g\n', P.kPon_max); end
    if isfield(P, 'kPoff_min'), fprintf(fid, '%% kPoff_min = %g\n', P.kPoff_min); end
    if isfield(P, 'kPoff_max'), fprintf(fid, '%% kPoff_max = %g\n', P.kPoff_max); end
    if isfield(P, 'kPoff_const'), fprintf(fid, '%% kPoff_const = %g\n', P.kPoff_const); end
    if isfield(P, 'geneLength_bp'), fprintf(fid, '%% geneLength_bp = %g\n', P.geneLength_bp); end
    if isfield(P, 'PASposition'), fprintf(fid, '%% PASposition = %g\n', P.PASposition); end
    if isfield(P, 'EBindingNumber'), fprintf(fid, '%% EBindingNumber = %g\n', P.EBindingNumber); end
end

function write_cpa_main_data(fid, data, P)
    fprintf(fid, '%% EBindingNumber: %d\n', data.EBindingNumber);
    fprintf(fid, '%% \n');
    write_base_parameters(fid, P);
    fprintf(fid, '%% \n');
    fprintf(fid, '%% Solution Data:\n');
    fprintf(fid, '%% R_sol (Polymerase concentrations at each node):\n');
    fprintf(fid, '%.6f ', data.R_sol);
    fprintf(fid, '\n');
    fprintf(fid, '%% REH_sol (Paused polymerase concentrations):\n');
    fprintf(fid, '%.6f ', data.REH_sol);
    fprintf(fid, '\n');
    if isfield(data, 'Ser2P')
        fprintf(fid, '%% Ser2P (Average phosphorylation levels):\n');
        fprintf(fid, '%.6f ', data.Ser2P);
        fprintf(fid, '\n');
    end
    if isfield(data, 'avg_E_bound')
        fprintf(fid, '%% Average E bound:\n');
        fprintf(fid, '%.6f ', data.avg_E_bound);
        fprintf(fid, '\n');
    end
    fprintf(fid, '%% Final free E: %.6f\n', data.Ef_ss);
    fprintf(fid, '%% Final free Pol II: %.6f\n', data.Pol_f_final);
end

function write_sweep_1d_data(fid, data, P)
    fprintf(fid, '%% 1D Parameter Sweep Results\n');
    fprintf(fid, '%% Swept parameter: %s\n', data.sweep_param);
    fprintf(fid, '%% EBindingNumber: %d\n', data.EBindingNumber);
    fprintf(fid, '%% \n');
    write_base_parameters(fid, P);
    fprintf(fid, '%% \n');
    fprintf(fid, '%% Parameter values (%d values):\n', length(data.param_values));
    fprintf(fid, '%g ', data.param_values);
    fprintf(fid, '\n');
    fprintf(fid, '%% \n');
    fprintf(fid, '%% Cutoff values (bp at 25%% cleavage):\n');
    fprintf(fid, '%.6f ', data.cutoff_values);
    fprintf(fid, '\n');
end

function write_sweep_2d_data(fid, data, P)
    fprintf(fid, '%% 2D Parameter Sweep Results\n');
    fprintf(fid, '%% Parameter 1: %s\n', data.param1);
    fprintf(fid, '%% Parameter 2: %s\n', data.param2);
    fprintf(fid, '%% EBindingNumber: %d\n', data.EBindingNumber);
    fprintf(fid, '%% \n');
    write_base_parameters(fid, P);
    fprintf(fid, '%% \n');
    fprintf(fid, '%% %s values (%d values):\n', data.param1, length(data.param1_values));
    fprintf(fid, '%g ', data.param1_values);
    fprintf(fid, '\n');
    fprintf(fid, '%% \n');
    fprintf(fid, '%% %s values (%d values):\n', data.param2, length(data.param2_values));
    fprintf(fid, '%g ', data.param2_values);
    fprintf(fid, '\n');
    fprintf(fid, '%% \n');
    fprintf(fid, '%% Cutoff Matrix (rows = %s, columns = %s):\n', data.param2, data.param1);
    fprintf(fid, '%% Values represent cutoff positions in bp\n');
    
    % Write the cutoff matrix
    for i = 1:size(data.cutoff_matrix, 1)
        fprintf(fid, '%.6f ', data.cutoff_matrix(i, :));
        fprintf(fid, '\n');
    end
end

function write_pas_distance_data(fid, data, P)
    fprintf(fid, '%% PAS Usage vs Inter-PAS Distance Analysis\n');
    fprintf(fid, '%% EBindingNumber: %d\n', data.EBindingNumber);
    fprintf(fid, '%% \n');
    write_base_parameters(fid, P);
    fprintf(fid, '%% \n');
    fprintf(fid, '%% Inter-PAS distances (bp) (%d values):\n', length(data.inter_pas_distances_bp));
    fprintf(fid, '%g ', data.inter_pas_distances_bp);
    fprintf(fid, '\n');
    fprintf(fid, '%% \n');
    fprintf(fid, '%% Proximal usage probability (0-1):\n');
    fprintf(fid, '%.6f ', data.proximal_usage_prob);
    fprintf(fid, '\n');
    if isfield(data, 'R_sol')
        fprintf(fid, '%% \n');
        fprintf(fid, '%% R_sol (Polymerase concentrations):\n');
        fprintf(fid, '%.6f ', data.R_sol);
        fprintf(fid, '\n');
    end
    if isfield(data, 'REH_sol')
        fprintf(fid, '%% REH_sol (Paused polymerase concentrations):\n');
        fprintf(fid, '%.6f ', data.REH_sol);
        fprintf(fid, '\n');
    end
end
