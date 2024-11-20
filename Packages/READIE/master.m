addpath('.');
addpath('./utils/');
addpath('./triallevel/');
addpath('./overall/');

%{
General Parameters
%}

% Define directory paths for data and results
SAVE_ROOT = '/Users/macbook/Desktop/Reliability_Project_Scripts/READIE_Toolbox_0806/results';
DATA_FOLDER = '/Users/macbook/Desktop/Reliability_Project_Scripts/READIE_Toolbox_0806/Preprocessed_Individual_trial_data_readieERP/readieERPs/ERP_calculatedVals';

% Extract participant name from file names using the specified divider
% e.g., for "2_191_49685484_3_20220901_110819_generatedERPvals_27-02-2024.csv",
% participant name is "2_191_49685484_3_20220901_110819"
FILENAME_DIVIDER = "_readieERPvals";

% Specify exclusion files 
IGNORE_CONTAINS = [
    "AllSubsAve", ...
    ];

% Define ALL conditions present in the dataset
% If the dataset only contains one condition, leave CONDITIONS blank
CONDITIONS = [
    % "_HUp+_", "_FUp+_", ...
    ];

% Specify the time windows of interest for data analysis
VALUE_COLUMNS = [
    "Mean Amplitude for Window 75-130",
    "Mean Amplitude for Window 100-230", 
    ];

%{
SME, reliability, and effect size iteration parameters
%}
NUM_ITERATIONS = 1000;

%{
Reliability parameters(Trial-level) 
%}
N_FROM = 10;
N_TO = 100;
N_BY = 10;


%{
other parameters,
no need to modify
%}
COLORS = [0, 0.4470, 0.7410; ... % MATLAB default blue
          0.8500, 0.3250, 0.0980; ... % MATLAB default orange
          0.9290, 0.6940, 0.1250; ... % MATLAB default yellow
          0.4940, 0.1840, 0.5560; ... % MATLAB default purple
          0.4660, 0.6740, 0.1880; ... % MATLAB default green
          0.3010, 0.7450, 0.9330; ... % MATLAB default light blue
          0.6350, 0.0780, 0.1840]; % MATLAB default red
rng(123)

warning('off', 'MATLAB:MKDIR:DirectoryExists');

% setting the number of CPUs and using parpool function 
numCores = feature('numcores');
numWorkers = numCores - 1;
if numWorkers <= 0
    numWorkers = 1; % Ensure at least one worker is used
end

% get data
fprintf(1, 'Getting Data ... \n');
% Read in all data from the specified folder for all participants 
all_data = read_data(DATA_FOLDER, FILENAME_DIVIDER, IGNORE_CONTAINS, ...
    CONDITIONS, VALUE_COLUMNS);
% Define the value sequence for trial-level reliability and effect size 
values_sequence = N_FROM:N_BY:N_TO;

%Loop for each time window 
for i = 1:length(VALUE_COLUMNS)
    vc = VALUE_COLUMNS(i);
    fprintf(1, '\nCalculating for variable %s:\n', vc);
    % Calculate SME 
    result_SME = calc_SME(all_data, vc, NUM_ITERATIONS, SAVE_ROOT);
    % Calculate trial level reliability 
    [res_rel_trial, summ_rel_trial] = ...
        calc_reliability_triallevel(all_data, vc, values_sequence, NUM_ITERATIONS, SAVE_ROOT);
    % Calculate overall reliability 
    [res_rel_overall, summ_rel_overall] = ...
        calc_reliability_overall(all_data, vc, NUM_ITERATIONS, SAVE_ROOT);
    % Plot trial level reliability 
    plot_from_summary(summ_rel_trial, "reliability", values_sequence, vc, unique(all_data.cond), COLORS, SAVE_ROOT);
    % Calculate trial level effect size 
    [res_eff_trial, summ_eff_trial] = calc_effectsize_triallevel(all_data, vc, values_sequence, NUM_ITERATIONS, SAVE_ROOT);
    % Calculate overall effect size 
    [res_eff_overall, summ_eff_overall] = ...
        calc_effectsize_overall(all_data, vc, NUM_ITERATIONS, SAVE_ROOT);
    % Plot trial level effect size
    plot_from_summary(summ_eff_trial, "effect size", values_sequence, vc, unique(all_data.cond), COLORS, SAVE_ROOT);
end
