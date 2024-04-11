%% SET FOLDERS AND PATH
% Use input from the command line to set the path to the data. If an 
% invalid path is entered, repeat until a valid path is entered.
clear ;
fprintf('Preparing to generate reliability...\n') ;
happeDir = strrep(fileparts(which(mfilename('fullpath'))), [filesep '3. ' ...
    'generate'], '') ;
addpath([happeDir filesep '3. generate'], ...
    [happeDir filesep 'files'], ...
    [happeDir filesep 'scripts'], ...
    [happeDir filesep 'scripts' filesep 'ui'], ...
    [happeDir filesep 'scripts' filesep 'support'], ...
    [happeDir filesep 'packages'], ...
    [happeDir filesep 'packages' filesep 'READIE'], ...
    [happeDir filesep 'packages' filesep 'READIE' filesep 'utils'], ...
    [happeDir filesep 'packages' filesep 'READIE' filesep 'triallevel'], ...
    [happeDir filesep 'packages' filesep 'READIE' filesep 'overall']) ;

%% SELECT PREVIOUSLY RUN SCRIPT
while true
    fprintf(['Select previous script used to create your .csv files:\n' ...
        '  1 - generateERPs\n  2 - generatePower\n  3 - Other\n'])
    prevScript = input('> ', 's') ;
    if str2double(prevScript) == 1 || str2double(prevScript) == 2; break;
    elseif str2double(prevScript) == 3
        error('Other scripts are not currently supported.') ;
    else; fprintf('Invalid input: please enter an integer between 1 and 3.\n') ;
    end
end
if str2double(prevScript) == 1; params.divider = '_generatedERPvals' ;
%  elseif str2double(prevScript) == 2; params.divider = ''
end

%% DETERMINE AND SET PATH TO THE DATA
% Use input from the command line to set the path to the data. If an 
% invalid path is entered, repeat until a valid path is entered.
while true
    srcDir = input('Enter the path to the folder containing your files:\n> ','s') ;
    if exist(srcDir, 'dir') == 7; break ;
    else; fprintf(['Invalid input: please enter the complete path to the ' ...
            'folder containing your files.\n']) ;
    end
end
cd(srcDir) ;
saveDir = [srcDir filesep 'READIE'] ;
if ~exist(saveDir, 'dir') == 7
    mkdir(saveDir) ;
end

%% DETERMINE THE INPUT PARAMETERS USED TO CREATE CSVs
while true
    % Use command window user input to collect the file name of the
    % pre-existing parameters. If the entered file is not an existing
    % file, repeats until a valid input is entered.
    fprintf(['Enter the parameter file you used to create these files, ' ...
        'including the full path and file extension:\n']) ;
    paramFile = input('> ', 's') ;
    if isfile(paramFile); break;
    else; fprintf('Invalid input: please enter the correct file\n') ;
    end
end
% LOAD THE PARAMETER FILE
fprintf('Loading parameters...\n') ;
load(paramFile) ;
paramsOrig = params ;
fprintf('Parameters loaded.\n') ;

%% SET PARAMETERS FOR RELIABILITY
% DETERMINE CONDITIONS OF INTEREST
fprintf(['Enter the conditions of interest, one at a time, pressing Enter/Return' ...
    ' between each entry.\nExample: stm+'])
params.conds = UI_cellArray(1, {}) ;

% SET ITERATION PARAMETERS
fprintf('Number of iterations to compute:\nExample: 1000\n') ;
params.its = input('> ') ;
fprintf('Subsample Start Value:\nExample: 10\n') ;
params.nFrom = input('> ') ;
fprintf('Subsample End Value:\nExample: 100\n') ;
params.nTo = input('> ') ;
fprintf('Subsample Increment Value:\nExample: 5\n') ;
params.nBy = input('> ') ;

%% OTHER PARAMETERS - NO NEED TO MODIFY
ignore_contains = "AllSubsAve" ;
COLORS = [0, 0.4470, 0.7410; ...                                            % MATLAB default blue
          0.8500, 0.3250, 0.0980; ...                                       % MATLAB default orange
          0.9290, 0.6940, 0.1250; ...                                       % MATLAB default yellow
          0.4940, 0.1840, 0.5560; ...                                       % MATLAB default purple
          0.4660, 0.6740, 0.1880; ...                                       % MATLAB default green
          0.3010, 0.7450, 0.9330; ...                                       % MATLAB default light blue
          0.6350, 0.0780, 0.1840];                                          % MATLAB default red
rng(123)

warning('off', 'MATLAB:MKDIR:DirectoryExists');

% Set number of CPUs for the parpool function 
numCores = feature('numcores');
numWorkers = numCores - 1;
if numWorkers <= 0; numWorkers = 1; end                                     % Ensure at least one worker

%% FIND MEAN AMP COL NAMES
if str2double(prevScript) == 1
    meanAmpNames = cell(1, paramsOrig.calcVals.meanAmpMethod(1)*size(paramsOrig.calcVals.windows,1)) ;
    if paramsOrig.calcVals.meanAmpMethod(1)
        for i=1:size(paramsOrig.calcVals.windows,1)
            meanAmpNames{i} = ['Mean Amplitude for Window ' ...
                paramsOrig.calcVals.windows{i,1} '-' paramsOrig.calcVals.windows{i,2}] ;
        end
    end
end
    

%% RUN THE MASTER SCRIPT
fprintf('Getting Data...\n') ;
all_data = read_data(srcDir, params.divider, ignore_contains, params.conds, ...
    meanAmpNames) ;
values_sequence = params.nFrom:params.nBy:params.nTo ;

for i = 1:length(meanAmpNames,2)
    vc = meanAmpNames{i} ;
    fprintf(1, 'Calculating for variable %s:\n', vc) ;
    % Calculate SME
    result_SME = calc_SME(all_data, vc, params.its, saveDir) ;

    % Calculate trial level reliability
    [res_rel_trial, summ_rel_trial] = calc_reliability_triallevel(alldata, ...
        vc, values_sequence, params.its, saveDir) ;

    % Calculate overall reliability 
    [res_rel_overall, summ_rel_overall] = calc_reliability_overall(all_data, ...
        vc, params.its, saveDir);

    % Plot trial level reliability 
    plot_from_summary(summ_rel_trial, "reliability", values_sequence, vc, ...
        unique(all_data.cond), COLORS, saveDir);

    % Calculate trial level effect size 
    [res_eff_trial, summ_eff_trial] = calc_effectsize_triallevel(all_data, ...
        vc, values_sequence, params.its, saveDir);

    % Calculate overall effect size 
    [res_eff_overall, summ_eff_overall] = calc_effectsize_overall(all_data, ...
        vc, params.its, saveDir);
    
    % Plot trial level effect size
    plot_from_summary(summ_eff_trial, "effect size", values_sequence, vc, ...
        unique(all_data.cond), COLORS, saveDir);
end