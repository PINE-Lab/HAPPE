%% generateNonlinearMeasures
% Generate Nonlinear Measures - a post-processing script to calculate a
% variety of nonlinear measures on HAPPE-processed data. Relies on python
% functions, and is adapted from a python script (see credit below).
%
% Developed at Northeastern University's PINE Lab
%
% Originally written in Python by William Bosl and Flemming Peck.
%
% Authors: A.D. Monachino, PINE Lab at Northeastern University, 2023
%          Gerardo Parra, Boston Children's Hospital, 2023
%
% This file is part of HAPPE.
% Copyright 2018-2023 Alexa Monachino, Kelsie Lopez, Laurel Gabard-Durnam
%
% HAPPE is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
% 
% HAPPE is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
% details.
% 
% You should have received a copy of the GNU General Public License along
% with HAPPE. If not, see <https://www.gnu.org/licenses/>.

%% ADD NECESSARY FOLDERS TO THE PATH
clear ;
fprintf('Preparing HAPPE generateNonlinearMeasures...\n') ;
happeDir = strrep(fileparts(which(mfilename('fullpath'))), [filesep '3. ' ...
    'generate'], '') ;
eeglabDir = [happeDir filesep 'packages' filesep 'eeglab2022.0'] ;
addpath([happeDir filesep '3. generate'], ...
    [happeDir filesep 'files'], ...
    [happeDir filesep 'scripts'], ...
    [happeDir filesep 'scripts' filesep 'ui'], ...
    [happeDir filesep 'scripts' filesep 'support'], ...
    [happeDir filesep 'scripts' filesep 'python'], ...
    eeglabDir, [eeglabDir filesep 'functions']) ;
rmpath(genpath([eeglabDir filesep 'functions' filesep 'octavefunc'])) ;

%% DETERMINE AND SET PATH TO THE DATA
% Use input from the Command Window to set the path to the data. If an 
% invalid path is entered, repeat until a valid path is entered.
while true
    srcDir = input(['Enter the path to the folder containing the ' ...
        'processed dataset(s):\n> '],'s') ;
    if exist(srcDir, 'dir') == 7; break ;
    else; fprintf(['Invalid input: please enter the complete path to the ' ...
            'folder containing the dataset(s).\n']) ;
    end
end
cd(srcDir) ;

%% COLLECT FILES
% Locates all .set files to run through the script. If no files meeting
% this criteria are found, throw an error and end the run.
fprintf('Gathering files...\n') ;
FileNames = {dir('*.set').name} ;
if isempty(FileNames); error('ERROR: NO .SET FILES FOUND') ; end

%% CREATE OUTPUTS FOLDER
% Create the folder to store outputs of this script.
fprintf('Creating output folder...\n') ;
if ~isfolder([srcDir filesep 'nonlinear_measures'])
    mkdir([srcDir filesep 'nonlinear_measures']) ;
end
addpath('nonlinear_measures') ;
cd([srcDir filesep 'nonlinear_measures']) ;

%% DETERMINE IF LOADING PARAMETERS
fprintf('Load a pre-existing set of parameters? [Y/N]\n') ;
preExist = choose2('N','Y') ;
% If indicated by the user, use Command Window input to collect and load
% the file of pre-existing parameters. If the entered file does not exist,
% repeat until an existing file is entered.
if preExist
    while true
        fprintf(['Enter your parameter file, including the full path ' ...
            'and file extension:\n']) ;
        paramFile = input('> ', 's') ;
        if isfile(paramFile); break;
        else; fprintf('Invalid input: please enter the correct file\n') ;
        end
    end
    fprintf('Loading parameters...\n') ; load(paramFile) ; fprintf('Parameters loaded.\n') ;
    % List the loaded parameters for user review and ask if they should be
    % changed. If the parameter set cannot be loaded, inform the user that
    % their parameter set is invalid.
    try genNlM_listParams(params) ;
    catch
        error('ERROR: This is not a valid set of parameters for this script\n') ;
    end
    fprintf('Change an existing parameter? [Y/N]\n') ;
    changedParams = choose2('n', 'y') ;
% If no existing parameters are loaded, create a set of default parameters
% to be edited/filled out during the Set Parameters step below.
else
    changedParams = 0 ;
    params = struct() ;
end

%% SET PARAMETERS
while true
    userChoice = '' ;
    %% BREAK IF NOT CHANGING EXISTING PARAMETERS
    if preExist && ~changedParams; break ; end
    
    %% IF CHANGING EXISTING PARAMETERS, LIST OPTIONS AND COLLECT USER INPUT
    if changedParams
        fprintf(['Parameter to change: features, frequency bands, ' ...
            'python version\n']) ;
        userChoice = input('> ', 's') ;
    end

    %% NONLINEAR MEASURE FEATURES:
    if ~preExist || strcmpi(userChoice, 'features')
        params.featureList = {'katz', 'higuchi', 'lzw_rho_median', ...
                'SampE', 'hurstrs', 'dfa'} ;
        % ASK TO COMPUTE ALL FEATURES OR A SUBSET
        fprintf('Compute all features or a subset of features?\n') ;
        % IF SELECTING A SUBSET, COLLECT USER INPUT ABOUT WHICH FEATURES TO
        % INCLUDE
        if choose2('all', 'subset')
            fprintf(['Enter your features of interest from the list ' ...
                    'below, exactly as displayed.\nPress Enter/Return' ...
                    'between each entry.\nWhen you have finished entering ' ...
                    'all features, enter "done" (without quotations).\n' ...
                    'Possible Features: ' sprintf('%s, ', params.featureList{1:end-1}) ...
                    params.featureList{end} '\n']) ;
            tempFeatures = {} ;
            while true
                ui = input('> ', 's') ;
                if ismember(ui, params.featureList)
                    tempFeatures = [tempFeatures, ui] ;                      %#ok<AGROW> 
                elseif strcmpi(ui, 'done')
                    if isempty(tempFeatures); fprintf('No features entered.\n') ;
                    else; break ;
                    end
                else; fprintf([ui ' is not a valid feature.\n']) ;
                end
            end
            params.featureList = intersect(params.featureList, unique(tempFeatures)) ;
            clear tempFeatures ;
        end
    end

    %% FREQUENCY BANDS:
    if ~preExist || strcmpi(userChoice, 'frequency bands')
        [params.fb.customBands, params.fb.bands] = setFreqBands() ;
    end
   
    %% PYTHON VERSION:
    if ~preExist || strcmpi(userChoice, 'python version')
        fprintf('Python version:\n') ;
        while true
            params.pyVers = input('> ','s') ;
            if str2double(params.pyVers) < 3
                error('Python 3 or above is required.\n') ;
            else; break;
            end
        end
    end

    %% CONFIRM PARAMETERS:
    if ~preExist || strcmpi('done', userChoice)
        fprintf('Please check your input parameters before continuing.\n') ;
        genNlM_listParams(params) ;
        fprintf('Are the above parameters correct? [Y/N]\n') ;
        if choose2('n','y'); break ;
        elseif ~preExist
            changedParams = 1; preExist = 1;
        end
    end
end

%% SAVE INPUT PARAMETERS
% If created a new or changed a parameter set, save as a new .mat file
if ~preExist || changedParams
    % Prompt to use a default or custom name for parameter the file. 
    % If the file exists, ask to create new file with a different name
    % or overwrite existing file.
    fprintf(['Parameter file save name:\n  default = Default name (genNlM' ...
        '_parameters_dd-mm-yyyy.mat)\n  custom = Create a custom file name' ...
        '\n']) ;
    if choose2('custom', 'default')
        paramFile = paramFile_validateExist(['genNlM_parameters_' ...
            datestr(now, 'dd-mm-yyyy') '.mat'], 'genNlM_parameters_', 2) ;
    else
        fprintf('File name (Do not include .mat):\n') ;
        paramFile = paramFile_validateExist([input('> ', 's') '.mat'], ...
            'genERP_parameters_', 0) ;
    end

    % Save the params variable to a .mat file using the name created above.
    fprintf('Saving parameters...\n') ;
    save(paramFile, 'params') ;
    fprintf('Parameters saved.\n') ;
end
clear('preExist', 'changedParams', 'paramFile') ;

%% LOAD AND CONFIGURE PYTHON
pyenv('Version', params.pyVers) ;
% IMPORT MODULES:
nolds = py.importlib.import_module('nolds') ;
antropy = py.importlib.import_module('antropy') ;
lzw_ruffini = py.importlib.import_module('LZW_Ruffini') ;

%% SET VARIABLES
decompBySub = cell(1, size(FileNames,2)) ; % Array to hold all the features

%% RUN NONLINEAR MEASURE CALCULATIONS ON FILES
cd(srcDir) ;
for currFile=1:size(FileNames,2)
    % LOAD THE FILE
    EEG = load('-mat', FileNames{currFile}) ;
     
    % RUN OVER EACH FREQUENCY BAND
    for currBand=1:size(params.fb.bands,1)
        % FILTER TO THE BAND
        EEGfilt = pop_eegfiltnew(EEG, str2double(params.fb.bands{currBand,2}), ...
            str2double(params.fb.bands{currBand,3}), [], 0, [], 0) ;
        % ITERATE THROUGH CHANNELS, skipping NaN channels
        for currChan=1:EEG.nbchan
            % ITERATE THROUGH EPOCHS
            for currEpoch=1:EEG.trials
                currData = EEGfilt.data(currChan, :, currEpoch) ;
                % 1A. SAMPLE ENTROPY
                if contains('SampE', params.featureList)
                    try D.("SampE").(params.fb.bands{currBand,1})(currChan, ...
                            currEpoch) = nolds.sampen(currData) ;
                    catch; D.("SampE").(params.fb.bands{currBand,1})(currChan, ...
                            currEpoch) = NaN ;
                    end
                end
                % 1B. HURST PARAMETER
                if contains('hurstrs', params.featureList)
                    try D.("hurstrs").(params.fb.bands{currBand,1})(currChan, ...
                            currEpoch) = nolds.hurst_rs(currData) ;
                    catch; D.("hurstrs").(params.fb.bands{currBand,1})(currChan, ...
                            currEpoch) = NaN ;
                    end
                end
                % 1C. DFA
                if contains('dfa', params.featureList)
                    try D.("dfa").(params.fb.bands{currBand,1})(currChan, ...
                            currEpoch) = nolds.dfa(currData) ;
                    catch; D.("dfa").(params.fb.bands{currBand,1})(currChan, ...
                            currEpoch) = NaN ;
                    end
                end
                % 2A. FRACTAL DIMENSION - KATZ
                if contains('katz', params.featureList)
                    try D.("katz").(params.fb.bands{currBand,1})(currChan, ...
                            currEpoch) = antropy.fractal.higuchi_fd(currData) ;
                    catch; D.("katz").(params.fb.bands{currBand,1})(currChan, ...
                            currEpoch) = NaN ;
                    end
                end
                % 2B. FRACTAL DIMENSION - HIGUCHI
                if contains('higuchi', params.featureList)
                    try D.("higuchi").(params.fb.bands{currBand,1})(currChan, ...
                            currEpoch) = antropy.fractal.higuchi_fd(currData) ;
                    catch; D.("higuchi").(params.fb.bands{currBand,1})(currChan, ...
                            currEpoch) = NaN ;
                    end
                end
                % 3. LEMPEL ZIV COMPLEXITY
                if contains('lzw_rho_median', params.featureList)
                    try D.("lzw_rho_median").(params.fb.bands{currBand,1})(currChan, ...
                            currEpoch) = lzw_ruffini.Compute_rho0(binarize(currData)) ;
                    catch; D.("lzw_rho_median").(params.fb.bands{currBand,1})(currChan, ...
                            currEpoch) = NaN ;
                    end
                end
            end
        end
    end
    decompBySub{currFile} = D ;
end

%% RECONFIGURE MATRIX TO OUTPUT FORMAT
% Average over epochs. *** Will add a means to keep all epochs
output = cell(1, size(params.featureList,2)) ;
for currFile=1:size(decompBySub,2)
    for currFeat=1:size(params.featureList,2)
        for currBand=1:size(params.fb.bands,1)
            for currChan=1:EEG.nbchan
                output{currFeat}{currBand}{currFile,currChan} = ...
                    mean(decompBySub{currFile}.(params.featureList{...
                    currFeat}).(params.fb.bands{currBand,1})(currChan, ...
                    :), 2, 'omitnan') ;
            end
        end
    end
end


%% OUTPUT TABLES AS .XLSX FILES
for currFeat=1:size(output,2) % Output a new file for each feature
    % Create Save Name
    saveName = [params.featureList{currFeat} '_generatedNonlinearMeasures' ...
        '_' datestr(now, 'dd-mm-yyyy') '.xlsx'] ;
    for currBand=1:size(output{currFeat},2)
        % for each band, append a new sheet to the file made of a table
        % from the current output
        writetable(array2table(output{currFeat}{currBand}, 'VariableNames', ...
            {EEG.chanlocs.labels}, 'RowNames', FileNames'), saveName, ...
            'WriteRowNames', true, 'Sheet', params.fb.bands{currBand,1}) ;
    end
end

%% SUPPORT FUNCTIONS
function genNlM_listParams(params)
    fprintf(['Nonlinear Measure Features: ' sprintf('%s, ', ...
        params.featureList{1:end-1}) params.featureList{end} '\n' ...
        'Frequency Bands:\n']) ;
    for i=1:size(params.fb.bands,1)
        fprintf([' - ' params.fb.bands{i,1} ': ' params.fb.bands{i,2} ...
            ' Hz - ' params.fb.bands{i,3} ' Hz\n']) ;
    end
    fprintf(['Python Version: ' params.pyVers '\n']) ;
end

% function to binarize vector
function v_b = binarize(v,method)
    % verify inptus
    if ~isvector(v), error('Input must be a vector'); end
    if ~exist('method','var'), method = 'median'; end
    % intialize vector of zeros
    v_b = zeros(1,length(v));
    % set threshold based on method
    switch method
        case 'median'
            thr = median(v);
        case 'mean'
            thr = mean(v);
    end
    % set indices above thr to one
    v_b(v>thr) = 1;
end