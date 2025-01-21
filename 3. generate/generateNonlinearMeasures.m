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
eeglabDir = [happeDir filesep 'packages' filesep 'eeglab2024.0'] ;
addpath([happeDir filesep '3. generate'], ...
    [happeDir filesep 'files'], ...
    [happeDir filesep 'scripts'], ...
    [happeDir filesep 'scripts' filesep 'ui'], ...
    [happeDir filesep 'scripts' filesep 'support'], ...
    [happeDir filesep 'scripts' filesep 'python'], ...
    eeglabDir,  genpath([eeglabDir filesep 'functions'])) ;
rmpath(genpath([eeglabDir filesep 'functions' filesep 'octavefunc'])) ;
pluginDir = dir([eeglabDir filesep 'plugins']) ;
pluginDir = strcat(eeglabDir, filesep, 'plugins', filesep, {pluginDir.name}, ';') ;
addpath([pluginDir{:}]) ;

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
    % BREAK IF NOT CHANGING EXISTING PARAMETERS
    if preExist && ~changedParams; break ; end
    
    % IF CHANGING EXISTING PARAMETERS, LIST OPTIONS AND COLLECT USER INPUT
    if changedParams
        fprintf(['Parameter to change: features, frequency bands, ' ...
            'python version\n']) ;
        userChoice = input('> ', 's') ;
    end

    % NONLINEAR MEASURE FEATURES:
    if ~preExist || strcmpi(userChoice, 'features')
        params.featureList = {'katz', 'higuchi', 'lzw_rho_median', ...
                'SampE', 'hurstrs', 'dfa'} ;
        % ASK TO COMPUTE ALL FEATURES OR A SUBSET
        fprintf('Compute all features or a subset of features?\n') ;
        % IF SELECTING A SUBSET, COLLECT USER INPUT ABOUT WHICH FEATURES TO
        % INCLUDE
        if choose2('all', 'subset')
            fprintf(['Enter your features of interest from the list ' ...
                    'below, exactly as displayed.\nPress Enter/Return ' ...
                    'between each entry.\nWhen you have finished entering ' ...
                    'all features, enter "done" (without quotations).\n' ...
                    'Possible Features: ' sprintf('%s, ', params.featureList{1:end-1}) ...
                    params.featureList{end} '\n']) ;
            tempFeatures = {} ;
            while true
                ui = input('> ', 's') ;
                if ismember(ui, params.featureList)
                    tempFeatures = [tempFeatures, ui] ;                       
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

    % FREQUENCY BANDS:
    if ~preExist || strcmpi(userChoice, 'frequency bands')
        [params.fb.customBands, params.fb.bands] = setFreqBands() ;
    end

    % CHANNELS OF INTEREST/REGIONS OF INTEREST:
    if ~preExist || strcmpi(userChoice, 'channels of interest')
        params.rois = {} ;
        indx = 1;
        while true
            [params.rois{1,indx}, params.rois{2,indx}] = determ_chanIDs() ;
            fprintf(['Enter a name for this subset:\nNOTE: We recommend ' ...
                'short names of approximately 6 characters or less\n']) ;
            params.rois{3,indx} = input('> ', 's') ;
            fprintf('Add another subset? [Y/N]\n') ;
            if ~choose2('N','Y'); break; end
            indx = indx+1;
        end
        % Remove ROIs with duplicate names.
        for i=1:size(params.rois,2)
            if sum(ismember(params.rois(3,:), params.rois{3,i})) > 1
                dupeNames = find(strcmpi(params.rois(3,:), params.rois{3,i})) ;
                params.rois = params.rois(:, setdiff(1:size(params.rois,2), ...
                    dupeNames(2:end))) ;
            end
        end
        % Remove duplicate requests for "all" channels
        if sum(ismember(params.rois(1,:), 'all')) > 1
            dupeNames = find(strcmpi(params.rois(1,:), 'all')) ;
            params.rois = params.rois(:, setdiff(1:size(params.rois,2), ...
                dupeNames(2:end))) ;
        end
    end

    % INDIVIDUAL, AVERAGE, OR BOTH TRIALS/SEGS
    fprintf(['Calculate nonlinear measures:\n  individual = By individual ' ...
        'trials/segments\n  average = For the average over trials/segments' ...
        '\n  both = Both individual trials/segments and average over ' ...
        'trials/segments\n']) ;
    while true
        ui = input('> ','s') ;
        if strcmpi(ui, 'individual'); params.method = [1,0]; break;
        elseif strcmpi(ui, 'average'); params.method = [0,1]; break;
        elseif strcmpi(ui, 'both'); params.method = [1,1]; break;
        else; fprintf(['Invalid input: please enter "individual", "average"' ...
                ', or "both" (without quotations)\n']) ;
        end
    end
   
    % PYTHON VERSION:
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

    % CONFIRM PARAMETERS:
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
            'genNlM_parameters_', 0) ;
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
cd([happeDir filesep 'scripts' filesep 'python'])
lzw_ruffini = py.importlib.import_module('LZW_Ruffini') ;

%% CREATE OUTPUTS FOLDERS
% Create the folder to store outputs of this script.
fprintf('Creating output folder...\n') ;
if ~isfolder([srcDir filesep 'nonlinear_measures'])
    mkdir([srcDir filesep 'nonlinear_measures']) ;
end
addpath([srcDir filesep 'nonlinear_measures']) ;
cd([srcDir filesep 'nonlinear_measures']) ;
for currFeat=1:size(params.featureList,2)
    if ~isfolder([srcDir filesep 'nonlinear_measures' filesep ...
            params.featureList{currFeat}])
        mkdir([srcDir filesep 'nonlinear_measures' filesep ...
            params.featureList{currFeat}]) ;
    end
end

%% SET VARIABLES
decompBySub = cell(size(FileNames,2), size(params.rois,2)) ; % Array to hold all the features

%% RUN NONLINEAR MEASURE CALCULATIONS ON FILES
cd(srcDir) ;
for currFile=1:size(FileNames,2)
    % LOAD THE FILE
    fprintf(sprintf('<strong>Loading %s...</strong>\n', FileNames{currFile})) ;
    EEG = load('-mat', FileNames{currFile}) ;
     
    % RUN OVER EACH FREQUENCY BAND
    for currBand=1:size(params.fb.bands,1)
        fprintf(sprintf('Processing %s band...\n', params.fb.bands{currBand,1})) ;
        % FILTER TO THE BAND
        EEGfilt = pop_eegfiltnew(EEG, str2double(params.fb.bands{currBand,2}), ...
            str2double(params.fb.bands{currBand,3}), [], 0, [], 0) ;
        
        % FIND CHANNELS OF INTEREST
        for currROI=1:size(params.rois,2)
            chanIndxs = [] ;
            if ~isempty(EEG.chanlocs)
                % If the ROI requests all channels, simply use all possibly
                % indicies.
                if strcmpi(params.rois{1,currROI}, 'all')
                    chanIndxs = 1:size(EEG.data,1) ;
                % If selecting to include or exclude channels, find use the
                % channel locations to get the index associated with the
                % relevant channel in the data matrix.
                else
                    for i=1:size(params.rois{2,currROI}, 2)
                        chanIndxs = [chanIndxs find(ismember({EEG.chanlocs.labels}, ...
                            params.rois{2, currROI}{i}))] ;                 %#ok<*AGROW> 
                    end
                    if strcmpi(params.rois{1,currROI}, 'coi_exclude')
                        chanIndxs = setdiff(1:size(EEG.data,1), ...
                            chanIndxs) ;
                    elseif ~strcmpi(params.rois{1,currROI}, 'coi_include')
                        error('Invalid channel COI selection.') ;
                    end
                end
            else
                fprintf(['Cannot select ROIs for data without channel names. ' ...
                    'Will use all existing channels.']) ;
                chanIndxs = 1:size(EEG.data,1) ;
            end
            saveChanIndxs{currROI} = chanIndxs ;

            % ITERATE THROUGH CHANNELS, skipping NaN channels
            for currChan=1:size(chanIndxs,2)
                % ITERATE THROUGH EPOCHS
                for currEpoch=1:EEG.trials
                    currData = EEGfilt.data(chanIndxs(currChan), :, currEpoch) ;
                    % 1A. SAMPLE ENTROPY
                    if contains('SampE', params.featureList)
                        try chanCalc.("SampE")(currEpoch) = nolds.sampen(currData) ;
                            
                            allSubs.("SampE"){currROI}{currBand}{currFile, ...
                                currChan}(currEpoch) = nolds.sampen(currData) ;
                        catch; D.("SampE").(params.fb.bands{currBand,1})(currChan, ...
                                currEpoch) = NaN ;
                        end
                    end
                    % 1B. HURST PARAMETER
                    if contains('hurstrs', params.featureList)
                        try allSubs.("SampE"){currROI}{currBand}{currFile, ...
                                currChan}{currEpoch} = nolds.hurst_rs(currData) ;
                        catch; D.("hurstrs").(params.fb.bands{currBand,1})(currChan, ...
                                currEpoch) = NaN ;
                        end
                    end
                    % 1C. DFA
                    if contains('dfa', params.featureList)
                        try allSubs.("dfa"){currROI}{currBand}{currFile, ...
                                currChan}{currEpoch} = nolds.dfa(currData) ;
                        catch; D.("dfa").(params.fb.bands{currBand,1})(currChan, ...
                                currEpoch) = NaN ;
                        end
                    end
                    % 2A. FRACTAL DIMENSION - KATZ
                    if contains('katz', params.featureList)
                        try allSubs.("katz"){currROI}{currBand}{currFile, ...
                                currChan}{currEpoch} = ...
                                antropy.fractal.katz_fd(currData) ;
                        catch; D.("katz").(params.fb.bands{currBand,1})(currChan, ...
                                currEpoch) = NaN ;
                        end
                    end
                    % 2B. FRACTAL DIMENSION - HIGUCHI
                    if contains('higuchi', params.featureList)
                        try allSubs.("higuchi"){currROI}{currBand}{currFile, ...
                                currChan}{currEpoch} = ...
                                antropy.fractal.higuchi_fd(currData) ;
                        catch; D.("higuchi").(params.fb.bands{currBand,1})(currChan, ...
                                currEpoch) = NaN ;
                        end
                    end
                    % 3. LEMPEL ZIV COMPLEXITY
                    if contains('lzw_rho_median', params.featureList)
                        try allSubs.("lzw_rho_median"){currROI}{currBand}{currFile, ...
                                currChan}{currEpoch} = ...
                                lzw_ruffini.Compute_rho0(binarize(currData)) ;
                        catch; D.("lzw_rho_median").(params.fb.bands{currBand,1})(currChan, ...
                                currEpoch) = NaN ;
                        end
                    end
                end
%                 fields = fieldnames(chanCalc) ;
%                 for currFeat = 1:length(fields)
%                     allSubs.(fields{currFeat}){currROI}{currBand}{currFile, ...
%                         currChan} = mean(chanCalc.(fields{currFeat})) ;
%                 end
            end
        end
    end
end

%% OUTPUT TABLES
chanLabels = {EEG.chanlocs.labels} ;
for currFeat = 1:size(params.featureList,2)
    currName = ['genNlM_' params.featureList{currFeat}] ;
    for currROI = 1:size(params.rois,2)
        currName = [currName '_' params.rois{3,currROI}] ;
        for currBand = 1:size(params.fb.bands,1)
            currName = [currName '_' params.fb.bands{currBand,1}] ;
            currMat = allSubs.(params.featureList{currFeat}){currROI}{currBand} ;
            for currFile = 1:size(FileNames,2)
                if params.method(1); outMat = [] ; end
                for currChan = 1:size(currMat, 2)
                    if params.method(1)
                        outMat = [outMat; currMat{currFile, currChan}] ;
                    end
                    if params.method(2)
                        allSubs.(params.featureList{currFeat}){currROI} ...
                            {currBand}{currFile, currChan} = mean(currMat{currFile, ...
                            currChan}) ;
                    end
                end
                if params.method(1)
                    colNames = {} ;
                    for i=1:size(outMat,2)
                        colNames = [colNames, ['Epoch_' num2str(i)]] ;
                    end
                    writetable(array2table(outMat, 'VariableNames', colNames, ...
                        'RowNames',  chanLabels(saveChanIndxs{currROI})'), ...
                        helpName([filesep 'indivTrials' filesep currName '-' ...
                        strrep(FileNames{currFile}, '.set', '.csv')]), ...
                        'WriteRowNames', true, 'QuoteStrings', true) ;
                end
            end
            if params.method(2)
                writetable(array2table(allSubs.(params.featureList{currFeat}){currROI}{currBand}, ...
                    'VariableNames', chanLabels(saveChanIndxs{currROI}), ...
                    'RowNames', FileNames'), helpName([currName '.csv'], '.csv'), ...
                    'WriteRowNames', true, 'QuoteStrings', true) ;
            end
        end
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