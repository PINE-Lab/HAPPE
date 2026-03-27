%% generatePower
% Generate Power - A post-processing script to calculate PSDs, bandpower,
% and specparam (as from Voytek & colleagues). Relies on functions from the
% EEGLAB specparam wrapper plugin by Donoghue et al., 2020.
% (https://github.com/bfbarry/EEGLAB-specparam), and python functionality.
%
% If using this script, please ensure you have python installed and set up
% to work with MATLAB.
%
% If using this script for analyses in publications, cite both HAPPE and
% the specparam/FOOOF EEGLAB wrapper script.
%       Donoghue T, Haller M, Peterson EJ, Varma P, Sebastian P, Gao R, 
%   Noto T, Lara AH, Wallis JD, Knight RT, Shestyuk A, & Voytek B (2020). 
%   Parameterizing neural power spectra into periodic and aperiodic
%   components. Nature Neuroscience, 23, 1655-1665. 
%   DOI: 10.1038/s41593-020-00744-x
%
% Developed at Northeastern University's PINE Lab
%
% Authors: A.D. Monachino, PINE Lab at Northeastern University, 2023
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

% NOTE: Code is not currently optimized for processing speed. Author
% intends to make these updates. Otherwise, should functionally work.


%% SET FOLDERS AND PATH
% Use input from the command line to set the path to the data. If an 
% invalid path is entered, repeat until a valid path is entered.
clear ;
fprintf('Preparing HAPPE generatePower...\n') ;
happeDir = strrep(fileparts(which(mfilename('fullpath'))), [filesep '3. ' ...
    'generate'], '') ;
addpath([happeDir filesep '3. generate'], ...
    [happeDir filesep 'files'], ...
    [happeDir filesep 'scripts'], ...
    [happeDir filesep 'scripts' filesep 'ui'], ...
    [happeDir filesep 'scripts' filesep 'support'], ...
    [happeDir filesep 'packages' filesep 'eeglab2024.0'], ...
    genpath([happeDir filesep 'packages' filesep 'eeglab2024.0' filesep ...
    'plugins' filesep 'EEGLAB-specparam-master'])) ;
pluginDir = dir([happeDir filesep 'packages' filesep 'eeglab2024.0' filesep 'plugins']) ;
pluginDir = strcat([happeDir filesep 'packages' filesep 'eeglab2024.0'], ...
    filesep, 'plugins', filesep, {pluginDir.name}, ';') ;
addpath([pluginDir{:}]) ;

%% DETERMINE AND SET PATH TO THE DATA
% Use input from the command line to set the path to the data. If an 
% invalid path is entered, repeat until a valid path is entered.
while true
    srcDir = input('Enter the path to the folder containing the processed dataset(s):\n> ','s') ;
    if exist(srcDir, 'dir') == 7; break ;
    else; disp("Invalid input: please enter the complete path to the folder containing the dataset(s).") ;
    end
end
cd(srcDir) ;

%% CREATE OUTPUT FOLDERS P.1
% Create the folder in which to store outputs.
fprintf('Creating output folders...\n') ;
if ~isfolder([srcDir filesep 'generatePower']); mkdir([srcDir filesep 'generatePower']); end
addpath('generatePower') ;

%% DETERMINE IF USING PRE-EXISTING SET OF PARAMETERS
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
    try genPower_listParams(params) ;
    catch; error('ERROR: This is not a valid set of parameters for this script\n') ;
    end
    fprintf('Change an existing parameter? [Y/N]\n') ;
    changedParams = choose2('N', 'Y') ;

% If no existing parameters are loaded, create a set of default parameters
% to be edited/filled out during the Set Parameters step below.
else; changedParams = 0; params = struct(); 
end

%% SET PARAMETERS THROUGH USER INPUT
% Set the necessary parameters to run through this script using user input
% via the Command Window.
while true
    userChoice = '' ;
    if preExist && ~changedParams; break; end

    % CHANGING PARAMETERS:
    % List out the options of parameters that can be changed by the user.
    if changedParams
        fprintf(['Parameter to change: psd indiv/ave/both, ' ...
            'limit output frequencies,\nchannels of interest, power in frequency ' ...
            'bands, specparam,\nspecparam outputs, specparam python, ' ...
            'specparam spectrum bounds,\nspecparam peak width, ' ...
            'specparam peak number, specparam peak height,\nspecparam ' ...
            'peak threshold, specparam mode, specparam peaks by freq bands\n' ...
            'Enter "done" (without quotations) when finished changing ' ...
            'parameters.\n']) ;
        userChoice = input('> ', 's') ;
    end

    % INDIVIDUAL OUTPUTS, AVERAGE OUTPUTS, OR BOTH:
    % Determine whether to write output files for each individual epoch,
    % for the average of epochs, or both. If outputting for individual
    % epochs, determine whether to use .csv or .xlsx for the output file
    % type.
    if ~preExist || strcmpi(userChoice, 'psd indiv/ave/both')
        [params.method, params.csvFormat] = determ_IndivAveBoth() ;
    end

    % LIMIT FREQUENCIES OF PSD OUTPUT:
    % Determine whether to select a subset of frequencies to include in the
    % PSD output file(s). For example, a user may choose to only examine
    % frequencies below 100 Hz.
    if ~preExist || strcmpi(userChoice, 'limit output frequencies')
        fprintf(['Limit the frequencies included in the PSD output? [Y/N]\n' ...
            'Example: include frequencies from 1 Hz to 40 Hz\n']) ;
        params.freqs.limit = choose2('N','Y') ;
        if params.freqs.limit
            params.freqs.min = input('Lower frequency bound:\nIf none, enter 0\n> ') ;
            params.freqs.max = input('Upper frequency bound:\nIf none, enter 0\n> ') ;
        end
    end

    % CHANNELS OF INTEREST/REGIONS OF INTEREST:
    % Determine the collections of channels of interest in which to examine
    % power measures. Enable the selection of all channels, or a subset,
    % and for multiple regions of interest.
    if ~preExist || strcmpi('channels of interest', userChoice)
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

    % CALCULATE POWER IN FREQUENCY BANDS:
    % Determine if calculating power by frequency bands. If so, select the
    % frequency bands - either using the default set or creating a custom
    % set. Additionally probes whether to export individual, average, or
    % both outputs. If individual is selected, determine whether to use
    % .csv or .xlsx for the file type.
    if ~preExist || strcmpi(userChoice, 'power in frequency bands')
        fprintf('Calculate power in frequency bands? [Y/N]\n') ;
        params.bands.on = choose2('N','Y') ;
        [~, params.bands.vals] = setFreqBands() ;
        [params.bands.method, params.bands.csvFormat] = determ_IndivAveBoth() ;
    end

    % CALCULATE SPECPARAM (FORMERLY KNOWN AS FOOOF):
    % Determine if calculating specparam. If so, set the relevant settings
    % including whether to enable visualizations, the bounds of the
    % spectrum to fit, and peak width limits, among others.
    if ~preExist || strcmpi(userChoice, 'specparam')
        fprintf('Calculate specparam (as from Voytek & colleagues)? [Y/N]\n') ;
        params.specparam.on = choose2('N','Y') ;
    end
    % Calculating Individual, Average, or Both Trials/Segments: 
    % If multiple ROIs, determine whether to split the output files by ROI
    % or keep all information from all ROIs in a single file. When
    % splitting, determine whether outputs should be multiple .csv files or
    % a single Excel file with multiple sheets.
    if params.specparam.on && (~preExist || strcmpi(userChoice, 'specparam outputs'))
        [params.specparam.method, params.specparam.csvFormat] = determ_IndivAveBoth() ;
        if size(params.rois,2) > 1
            fprintf(['Split specparam output files by ROIs? [Y/N]\nNOTE: ' ...
                'Choosing not to split = more columns; choosing to split ' ...
                '= more files/tabs.\n']) ;
            params.specparam.splitROI.on = choose2('N','Y') ;
            if params.specparam.splitROI.on
                fprintf(['Print the ROIs as seperate .csv files or as sheets\n' ...
                    'in a single Excel file?\n  csv = Output seperate ' ...
                    '.csv files\n  sheets = Output a single Excel file with ' ...
                    'multiple sheets\n']) ;
                while true
                    ui = input('> ','s') ;
                    if strcmpi(ui, 'csv'); params.specparam.splitROI.csvFormat = 1; break;
                    elseif strcmpi(ui, 'sheets'); params.specparam.splitROI.csvFormat = 0; break;
                    else; fprintf(['Invalid input: please enter "csv" or "sheets"' ...
                            ' (without quotations).']) ;
                    end
                end
            else; params.specparam.splitROI.csvFormat = 0 ;
            end
        else; params.specparam.splitROI.on = 0;
        end
        
        % Determine if enabling visualizations during the run. Only
        % recommended for a single file with a small ROI for testing
        % reasons as it can dramatically slow down processing and result in
        % many windows. Visuals will save regardless of choice.
        if params.specparam.method(2)
            fprintf('Enable visualizations? [Y/N]\n') ;
            params.specparam.vis = choose2('N', 'Y') ;
        else; params.specparam.vis = 0 ;
        end
    end
    % Determine Python Version:
    % Get the python version installed on the user's computer. For Windows,
    % the version number is sufficient. For Mac and Linux computers, you
    % need to enter the full path to the python executable.
    if params.specparam.on && (~preExist || strcmpi(userChoice, 'specparam python'))
        if ispc; params.specparam.pyVers = input('Python version:\n','s') ;
        else; params.specparam.pyVers = input(['Full path to your python ' ...
                'executable:\n'], 's') ;
        end
    end
    % Determine Bounds of Specparam Spectrum:
    % Determine the upper and lower limits, in Hz, for the specparam
    % spectrum calculations.
    if params.specparam.on && (~preExist || strcmpi(userChoice, 'specparam spectrum bounds'))
        fprintf('Enter specparam lower bound of spectrum to fit, in Hz:\nExample: 2\n') ;
        params.specparam.min = input('> ') ;
        fprintf('Enter specparam upper bound of spectrum to fit, in Hz:\nExample: 40\n') ;
        params.specparam.max = input('> ');
    end
    % Determine Peak Width - DO NOT RENAME VARIABLE
    if params.specparam.on && (~preExist || strcmpi(userChoice, 'specparam peak width'))
        fprintf('Enter lower peak width limit:\n') ;
        params.specparam.settings.peak_width_limits(1) = input('> ') ;
        fprintf('Enter upper peak width limit:\n') ;
        params.specparam.settings.peak_width_limits(2) = input('> ') ;
    end
    % Determine Number of Peaks to Find - DO NOT RENAME VARIABLE
    if params.specparam.on && (~preExist || strcmpi(userChoice, 'specparam peak number'))
        fprintf('Enter max number of peaks:\n') ;
        params.specparam.settings.max_n_peaks = input('> ') ;
    end
    % Determine Minimum Peak Height - DO NOT RENAME VARIABLE
    if params.specparam.on && (~preExist || strcmpi(userChoice, 'specparam peak height'))
        fprintf('Minimum peak height:\n') ;
        params.specparam.settings.min_peak_height = input('> ') ;
    end
    % Determine Peak Threshold - DO NOT RENAME VARIABLE
    if params.specparam.on && (~preExist || strcmpi(userChoice, 'specparam peak threshold'))
            fprintf('Enter peak threshold:\n') ;
            params.specparam.settings.peak_threshold = input('> ') ;
    end
    % Determine Peak Mode - DO NOT RENAME VARIABLE
    if params.specparam.on && (~preExist || strcmpi(userChoice, 'specparam mode'))
        fprintf(['Aperiodic mode?\n  fixed = fixed aperiodic mode\n  ' ...
            'knee = knee aperiodic mode\n']) ;
        while true
            ui = input('> ','s') ;
            if strcmpi(ui, 'knee')
                params.specparam.settings.aperiodic_mode = 'knee' ;
                break ;
            elseif strcmpi(ui, 'fixed')
                params.specparam.settings.aperiodic_mode = 'fixed' ;
                break ;
            else; fprintf(['Invalid Input: Please enter "fixed" or "knee"' ...
                    '(without quotations).']) ;
            end
        end
    end
    % Look for Peaks within Frequency Bands:
    % Determine whether or not to look for peaks within specific frequency
    % bands. If so, define the bands in which to look for peaks.
    if params.specparam.on && (~preExist || strcmpi(userChoice, 'specparam peaks by freq bands'))
        fprintf('Select peaks within frequency ranges? [Y/N]\n') ;
        params.specparam.bands.on = choose2('N','Y') ;
        if params.specparam.bands.on
            [~, params.specparam.bands.vals] = setFreqBands() ;
        end
    end

    % DONE:
    % When finished entering parameters, ask the user to review their
    % parameters and make any necessary corrections.
    if ~preExist || strcmpi('done', userChoice)
        fprintf('Please check your input parameters before continuing.\n') ;
        genPower_listParams(params) ;
        fprintf('Are the above parameters correct? [Y/N]\n') ;
        if choose2('N', 'Y'); break;
        elseif ~preExist; changedParams = 1; preExist = 1;
        end
    end
end

%% SAVE INPUT PARAMETERS
% If created a new or changed a parameter set, save as a new .mat file
if ~preExist || changedParams
    % Prompt to use a default or custom name for parameter the file. 
    % If the file exists, ask to create new file with a different name
    % or overwrite existing file.
    fprintf(['Parameter file save name:\n  default = Default name (genPower' ...
        '_parameters_dd-mm-yyyy.mat)\n  custom = Create a custom file name' ...
        '\n']) ;
    if choose2('custom', 'default')
        paramFile = paramFile_validateExist(['genPower_parameters_' ...
            datestr(now, 'dd-mm-yyyy') '.mat'], 'genPower_parameters_', 2) ;
    else
        fprintf('File name (Do not include .mat):\n') ;
        paramFile = paramFile_validateExist([input('> ', 's') '.mat'], ...
            'genPower_parameters_', 0) ;
    end

    % Save the params variable to a .mat file using the name created above.
    fprintf('Saving parameters...\n') ; 
    save([srcDir filesep 'generatePower' filesep paramFile], 'params') ;
    fprintf('Parameters saved.\n') ;
end
clear('preExist', 'changedParams', 'paramFile') ;

%% CREATE OUTPUT FOLDERS P.2
if params.method(1) && ~isfolder([srcDir filesep 'generatePower' filesep ...
        'PSD_IndivTrials'])
    mkdir([srcDir filesep 'generatePower' filesep 'PSD_IndivTrials']) ;
end
if params.method(2) && ~isfolder([srcDir filesep 'generatePower' filesep ...
        'PSD_AveOverTrials'])
    mkdir([srcDir filesep 'generatePower' filesep 'PSD_AveOverTrials']) ;
end
if params.bands.on && params.bands.method(1) && ~isfolder([srcDir filesep ...
        'generatePower' filesep 'bandpower'])
    mkdir([srcDir filesep 'generatePower' filesep 'bandpower_IndivTrials']) ;
end
if params.bands.on && params.bands.method(2) && ~isfolder([srcDir filesep ...
        'generatePower' filesep 'bandpower'])
    mkdir([srcDir filesep 'generatePower' filesep 'bandpower_AveOverTrials']) ;
end
if params.specparam.on && params.specparam.method(1) && ~isfolder([srcDir ...
        filesep 'generatePower' filesep 'specparam_IndivTrials'])
    mkdir([srcDir filesep 'generatePower' filesep 'specparam_IndivTrials']) ;
end
if params.specparam.on && params.specparam.method(2) && ~isfolder([srcDir ...
        filesep 'generatePower' filesep 'specparam_AveOverTrials'])
    mkdir([srcDir filesep 'generatePower' filesep 'specparam_AveOverTrials']) ;
end
fprintf('Output folders created.\n') ;

%% COLLECT FILES
% Locates all .set files to run through the script. If no files meeting
% this criteria are found, throw an error and end the run.
fprintf('Gathering files...\n') ;
FileNames = {dir('*.set').name} ;
if isempty(FileNames); error('ERROR: NO .SET FILES FOUND.'); end

%% LOAD AND CONFIGURE PYTHON (IF RUNNING SPECPARAM)
if params.specparam.on
    pyenv('Version', params.specparam.pyVers) ;
end
%% SUBJECT LEVEL MATRIX
% Create a matrix to hold the calculated measures across subjects.
allSubs = cell(size(FileNames,2),4) ;
subLvl_PSD = cell(1, size(params.rois,2)) ;
subLvl_bandpower = cell(2, size(params.rois,2)) ;
subLvl_specparam = cell(2, size(params.rois,2)) ;
for i=1:size(subLvl_specparam,2)
    subLvl_specparam{1,i} = {} ;
end

errorLog = {} ;

%% RUN CALCULATIONS ON EACH FILE
for currFile = 1:size(FileNames,2)
    try
        %% LOAD EEG
        % Try to load the processed file. If unable to load, alert the user
        % and rethrow the error.
        try 
            fprintf(['Loading ' FileNames{currFile} '...\n']) ;
            EEG = load('-mat', FileNames{currFile}) ;
        catch ME; fprintf(2, ['ERROR: Unable to load ' FileNames{currFile} ...
                '.\n ']) ;
            rethrow(ME) ;
        end

        %% FIND CHANNELS OF INTEREST
        % If the EEG has channel locations/names, use the user-provided
        % channel names in each region to determine the associated index.
        % If no channel locations are included, 
        chanIndxs = cell(1, size(params.rois,2)) ;
        for currROI=1:size(chanIndxs,2)
            if ~isempty(EEG.chanlocs)
                % If the ROI requests all channels, simply use all possibly
                % indicies.
                if strcmpi(params.rois{1,currROI}, 'all')
                    chanIndxs{1, currROI} = 1:size(EEG.data,1) ;
                % If selecting to include or exclude channels, find use the
                % channel locations to get the index associated with the
                % relevant channel in the data matrix.
                else
                    temp = [] ;
                    for i=1:size(params.rois{2,currROI}, 2)
                        temp = [temp find(ismember({EEG.chanlocs.labels}, ...
                            params.rois{2, currROI}{i}))] ;                                 
                    end
                    if strcmpi(params.rois{1,currROI}, 'coi_exclude')
                        chanIndxs{1, currROI} = setdiff(1:size(EEG.data,1), ...
                            temp) ;
                    elseif strcmpi(params.rois{1,currROI}, 'coi_include')
                        chanIndxs{1, currROI} = temp ;
                    else
                        error('Invalid channel COI selection.') ;
                    end
                end
            else
                fprintf(['Cannot select ROIs for data without channel names. ' ...
                    'Will use all existing channels.']) ;
                chanIndxs{1, currROI} = 1:size(EEG.data,1) ;
            end
        end
    
        %% FIND CHANNELS WITH ALL ZEROS
        zeroChanIndxs = [] ;
        for i=1:size(EEG.data,1)
            if all(all(EEG.data(i,:,:) == 0)); zeroChanIndxs = [zeroChanIndxs i]; end %#ok<*AGROW> 
        end

        %% CALCULATE MULTITAPER PSD AND 95% CONFIDENCE INTERVALS:
        % Use 4 multitaper windows, the default nfft, EEG's sampling rate, and
        % Slepian tapers for each segment within each channel.
        try
            clear('psd', 'psdLow', 'psdHigh') ;
            fprintf(['Calculating multitaper PSD and 95%% confidence ' ...
                'intervals...\n']) ;
            for currChan=1:EEG.nbchan
                for currSeg=1:EEG.trials
                    [pxx, f, pxxc] = pmtm(EEG.data(currChan,:,currSeg), 4, [], ...
                        EEG.srate, 'Tapers', 'slepian', 'ConfidenceLevel', 0.95) ;
                    psd(currChan,:,currSeg) = pxx' ;                        %#ok<*SAGROW> 
                    psdLow(currChan,:,currSeg) = pxxc(:,1)' ;
                    psdHigh(currChan,:,currSeg) = pxxc(:,2)' ;
                end
            end
            allSubs{currFile, 1} = psd ;
        catch ME; fprintf(2, ['ERROR: Unable to calculate multitaper PSD ' ...
                'for ' FileNames{currFile} '.\n ']) ;
            rethrow(ME) ;
        end

        %% CALCULATE AVERAGE PSD
        try
            fprintf('Calculating the average PSD...\n') ;
            psd_ave = mean(psd,3)' ;
            psdLow_ave = mean(psdLow,3)' ;
            psdHigh_ave = mean(psdHigh,3)' ;
            allSubs{currFile,2} = psd_ave ;
        catch ME; fprintf(2, ['ERROR: Unable to calculate average PSD ' ...
                'for ' FileNames{currFile} '.\n ']) ;
            rethrow(ME) ;
        end
        
        %% CREATE OUTPUT TABLE INFORMATION
        % Determine indicies for the frequencies to output, using the value
        % closest to the limit and within the specified range. If not limiting
        % the frequencies according to user-specificiation, uses them all.
        if params.freqs.limit
            [~, freqMin] = min(abs(f-params.freqs.min)) ;
            while f(freqMin) < params.freqs.min && freqMin < size(f,1)
                freqMin = freqMin+1 ;
            end
            [~, freqMax] = min(abs(f-params.freqs.max)) ;
            while f(freqMax) > params.freqs.max && freqMin < size(f,1)
                freqMax = freqMax-1 ;
            end
        else; freqMin = 1; freqMax = size(f,1) ;
        end
    
        % Channel Names: If no channel names exist, simply use 'Channel #' for
        % each channel. Otherwise, use the existing channel labels.
        if isempty(EEG.chanlocs)
            varNames = cell(1,size(psd,1)) ;
            for i=1:size(varNames,2); varNames{i} = ['Channel ' num2str(i)] ; end
        else; varNames = {EEG.chanlocs.labels} ;
        end
    
        % Confidence Interval Channel Names: Has the lower CI bound first,
        % followed by the upper CI bound for each channel. Use existing channel
        % names when applicable. Otherwise, label channels using their index.
        CIvarNames = cell(1,size(psd,1)*2) ;
        if isempty(EEG.chanlocs)
            for i=1:size(psd,1)
                CIvarNames{i*2-1} = ['Channel ' num2str(i) ' Lower CI'] ;
                CIvarNames{i*2} = ['Channel ' num2str(i) ' Upper CI'] ;
            end
        else
            for i=1:size(psd,1)
                CIvarNames{i*2-1} = [varNames{i} ' Lower CI'] ;
                CIvarNames{i*2} = [varNames{i} ' Upper CI'] ;
            end
        end
        
        %% PRINT OUT INDIVIDUAL PSD FILES (IF ENABLED)
        if params.method(1)
            fprintf('Saving Individual PSDs...\n') ;
            if params.csvFormat
                if ~isfolder([srcDir filesep 'generatePower' filesep ...
                        'PSD_IndivTrials' filesep strrep(FileNames{currFile}, ...
                        '.set', '')])
                    mkdir([srcDir filesep 'generatePower' filesep 'PSD_IndivTrials' ...
                        filesep strrep(FileNames{currFile}, '.set', '')]) ;
                end
                saveName = [srcDir filesep 'generatePower' filesep 'PSD_IndivTrials' ...
                    filesep strrep(FileNames{currFile}, '.set', '') filesep ...
                    strrep(FileNames{currFile}, '.set', '_indivPSD_epoch#.csv')] ;
                CIsaveName = [srcDir filesep 'generatePower' filesep 'PSD_IndivTrials' ...
                    filesep strrep(FileNames{currFile}, '.set', '') filesep ...
                    strrep(FileNames{currFile}, '.set', '_indivPSD_CI_epoch#.csv')] ;
            else
                saveName = helpName([srcDir filesep 'generatePower' filesep 'PSD_IndivTrials' ...
                    filesep strrep(FileNames{currFile}, '.set', '_indivPSD.xlsx')], ...
                    '.xlsx') ;
                CIsaveName = helpName([srcDir filesep 'generatePower' filesep 'PSD_IndivTrials' ...
                    filesep strrep(FileNames{currFile}, '.set', '_indivPSD_CI.xlsx')], ...
                    '.xlsx') ;
            end
            for currSeg=1:EEG.trials
                for i=1:size(psd,1)
                    CIpsd(:, i*2-1) = psdLow(i,:,currSeg) ;
                    CIpsd(:, i*2) = psdHigh(i,:,currSeg) ;
                end
                if params.csvFormat
                    saveName2 = strrep(saveName, '#', num2str(currSeg)) ;
                    CIsaveName2 = strrep(CIsaveName, '#', num2str(currSeg)) ;
                    writetable(array2table(psd(:,freqMin:freqMax,currSeg)', ...
                        'VariableNames', varNames, 'RowNames', ...
                        cellstr(num2str(f(freqMin:freqMax)))), helpName(saveName2, ...
                        '.csv'), 'WriteRowNames', true, 'QuoteStrings', true) ;
                    writetable(array2table(CIpsd(freqMin:freqMax,:), ...
                        'VariableNames', CIvarNames, 'RowNames', ...
                        cellstr(num2str(f(freqMin:freqMax)))), helpName(CIsaveName2, ...
                        '.csv'), 'WriteRowNames', true, 'QuoteStrings', true) ;
                else
                    writetable(array2table(psd(:,freqMin:freqMax,currSeg)', ...
                        'VariableNames', varNames, 'RowNames', ...
                        cellstr(num2str(f(freqMin:freqMax)))), saveName, 'Sheet', ...
                        ['epoch_' num2str(currSeg)], 'WriteRowNames', true) ;
                    writetable(array2table(CIpsd(freqMin:freqMax,:), ...
                        'VariableNames', CIvarNames, 'RowNames', ...
                        cellstr(num2str(f(freqMin:freqMax)))), CIsaveName, 'Sheet', ...
                        ['epoch_' num2str(currSeg)], 'WriteRowNames', true) ;
                end
            end
        end
    
        %% PRINT OUT THE AVERAGE PSD FILE
        if params.method(2)
            fprintf('Saving average PSD...\n') ;
            for i=1:size(psd,1)
                CIpsd_ave(:, i*2-1) = psdLow_ave(:,i) ;
                CIpsd_ave(:, i*2) = psdHigh_ave(:,i) ;
            end
    
            saveName = helpName([srcDir filesep 'generatePower' filesep ...
                'PSD_AveOverTrials' filesep strrep(FileNames{currFile}, '.set', ...
                '_avePSD.csv')], '.csv') ;
            CIsaveName = helpName([srcDir filesep 'generatePower' filesep ...
                'PSD_AveOverTrials' filesep strrep(FileNames{currFile}, '.set', ...
                '_avePSD_CI.csv')], '.csv') ;
            
            writetable(array2table(psd_ave(freqMin:freqMax,:), 'VariableNames', ...
                varNames, 'RowNames', cellstr(num2str(f(freqMin:freqMax)))), ...
                saveName, 'WriteRowNames', true, 'QuoteStrings', true) ;
            writetable(array2table(CIpsd_ave(freqMin:freqMax,:), 'VariableNames', ...
                CIvarNames, 'RowNames', cellstr(num2str(f(freqMin:freqMax)))), ...
                CIsaveName, 'WriteRowNames', true, 'QuoteStrings', true) ;
        end

        % SAVE INFO FOR ALLSUBS_PSD
        for i=1:size(params.rois,2)
            subLvl_PSD{1, i} = [subLvl_PSD{i} mean(allSubs{currFile,2}(:, ...
                chanIndxs{1,i}),2)] ;
        end
        
        %% CALCULATE PSD BY FREQUENCY BAND
        if params.bands.on
            allSubs{currFile, 3} = cell(size(params.bands.vals,1),size(chanIndxs,2)) ;
            for currBand = 1:size(params.bands.vals,1)
                fprintf(['Calculating PSD in the ' params.bands.vals{currBand,1} ...
                    ' band...\n']) ;
                % Determine indicies for frequency band:
                [~, bandMin] = min(abs(f-str2double(params.bands.vals{currBand,2}))) ;
                while f(bandMin) < str2double(params.bands.vals{currBand,2})
                    bandMin = bandMin+1 ;
                end
                [~, bandMax] = min(abs(f-str2double(params.bands.vals{currBand,3}))) ;
                while f(bandMax) > str2double(params.bands.vals{currBand,3})
                    bandMax = bandMax-1 ;
                end
                    
                for currROI=1:size(chanIndxs,2)
                    if params.bands.method(1)
                        % For each segment/trial, calculate the bandpower.
                        for currSeg=1:EEG.trials
                            powerBand = bandpower(psd(setdiff(chanIndxs{currROI}, ...
                                zeroChanIndxs), :, currSeg)', f, [f(bandMin) ...
                                f(bandMax)], 'psd') ;
                            allSubs{currFile,3}{currBand,currROI} = ...
                                [allSubs{currFile,3}{currBand, currROI}; ...
                                [powerBand mean(powerBand)]] ;
                        end
                    end
                    
                    powerBand = bandpower(psd_ave(:, setdiff(chanIndxs{currROI}, ...
                        zeroChanIndxs)), f, [f(bandMin) f(bandMax)], 'psd') ;
                    allSubs{currFile,3}{currBand,currROI} = [allSubs{currFile, ...
                        3}{currBand,currROI}; [powerBand mean(powerBand)]] ;
                    clear('powerBand') ;
                end
            end
    
            % PRINT OUT BANDS
            rowNames = {} ;
            for i=1:EEG.trials; rowNames = [rowNames; ['epoch_' num2str(i)]]; end  
            rowNames = [rowNames; 'average_epoch'] ;
            if params.bands.method(1)
                if params.bands.csvFormat
                    if ~isfolder([srcDir filesep 'generatePower' filesep ...
                            'bandpower_IndivTrials' filesep strrep(FileNames{currFile}, ...
                            '.set', '')])
                        mkdir([srcDir filesep 'generatePower' filesep ...
                            'bandpower_IndivTrials' filesep strrep(FileNames{currFile}, ...
                            '.set', '')]) ;
                    end
                    saveName = [srcDir filesep 'generatePower' filesep ...
                        'bandpower_IndivTrials' filesep strrep(FileNames{currFile}, ...
                        '.set', '') filesep strrep(FileNames{currFile}, ...
                        '.set', '_indivBandpower_FREQBAND_ROI.csv')] ;
                else
                    saveName = helpName([srcDir filesep 'generatePower' filesep ...
                        'bandpower_IndivTrials' filesep strrep(FileNames{currFile}, ...
                        '.set', '_indivBandpower.xlsx')], '.xlsx') ;
                end
                for currROI=1:size(chanIndxs,2)
                    for currBand=1:size(params.bands.vals,1)
                        if params.bands.csvFormat
                            saveName2 = strrep(strrep(saveName, 'FREQBAND',  ...
                                params.bands.vals{currBand,1}), 'ROI', ...
                                params.rois{3, currROI});
                            writetable(array2table(allSubs{currFile,3}{currBand, ...
                                currROI}, 'VariableNames', [varNames(setdiff(chanIndxs{currROI}, ...
                                zeroChanIndxs)), 'AverageChan'], 'RowNames', ...
                                rowNames), helpName(saveName2, '.csv'), ...
                                'WriteRowNames', true, 'QuoteStrings', true) ;
                        else
                            writetable(array2table(allSubs{currFile,3}{currBand, ...
                                currROI}, 'VariableNames', [varNames(setdiff(chanIndxs{currROI}, ...
                                zeroChanIndxs)), 'AverageChan'], 'RowNames', ...
                                rowNames), saveName, 'Sheet', ...
                                [params.bands.vals{currBand,1} '_' ...
                                params.rois{3,currROI}], 'WriteRowNames', true) ;
                        end
                    end
                end
            end
    
            if params.method(2)
                for currROI=1:size(params.rois,2)
                    saveName = helpName([srcDir filesep 'generatePower' filesep ...
                        'bandpower_AveOverTrials' filesep strrep(FileNames{currFile}, ...
                        '.set', ['_aveBandpower_' params.rois{3, currROI} ...
                        '.csv'])], '.csv') ;
        
                    temp = [] ;
                    for currBand=1:size(params.bands.vals,1)
                        temp = [temp; allSubs{1,3}{currBand,currROI}(end,:)] ;
                    end
        
                    writetable(array2table(temp, 'VariableNames', ...
                        [varNames(setdiff(chanIndxs{currROI}, zeroChanIndxs)), ...
                        'AverageChan'], 'RowNames', params.bands.vals(:,1)), saveName, ...
                        'WriteRowNames', true, 'QuoteStrings', true) ;
                end
            end

            % ALLSUBS_BANDPOWER
            for currROI=1:size(params.rois,2)
                tempBand = [] ;
                tempName = {};
                for i=1:size(params.bands.vals,1)
                    tempBand = [tempBand allSubs{currFile,3}{i, currROI}(end,:)] ;
                    if currFile == 1
                        tempName_chans = {} ;
                        for j=1:size(allSubs{currFile,3}{1, currROI},2)-1 % NEEDS TO BE BASED ON CHANNELS
                            tempName_chans = [tempName_chans {[params.bands.vals{i,1} ...
                                '_' varNames{chanIndxs{currROI}(j)}]}] ;
                        end
                        tempName = [tempName [tempName_chans, {[params.bands.vals{i,1} ...
                            '_ROIave']}]] ;
                    end
                end
                subLvl_bandpower{1,currROI} = [subLvl_bandpower{1,currROI}; tempBand] ;
                if currFile == 1
                    subLvl_bandpower{2,currROI} = [subLvl_bandpower{2,currROI}; tempName] ;
                end
            end
        end
    
        %% CALCULATE SPECPARAM
        if params.specparam.on
            fprintf('Calculating specparam...\n') ;
            allSubs{currFile,4} = cell(params.specparam.method(1) * ...
                EEG.trials+params.specparam.method(2), size(chanIndxs,2)) ;

            % Specparam by individual epochs
            for currROI=1:size(chanIndxs,2)
                if params.specparam.method(1)
                    for currSeg=1:EEG.trials
                        currPSD = psd(setdiff(chanIndxs{currROI}, zeroChanIndxs), ...
                            :, currSeg) ;
                        allSubs{currFile,4}{currSeg,currROI} = fooof_group(f', ...
                            [currPSD; mean(currPSD)]', [params.specparam.min ...
                            params.specparam.max], params.specparam.settings, ...
                            true) ;
                    end
                end
        
                % Specparam on average psd
                currPSD = psd_ave(:, setdiff(chanIndxs{currROI}, zeroChanIndxs)) ;
                allSubs{currFile,4}{end,currROI} = fooof_group(f', [currPSD, ...
                    mean(currPSD,2)], [params.specparam.min params.specparam.max], ...
                    params.specparam.settings, true) ;
            end
    
            % output results
            if ~isfolder([srcDir filesep 'generatePower' filesep 'specparam' ...
                    '_AveOverTrials' filesep strrep(FileNames{currFile}, '.set', ...
                    '_figures')]) && params.specparam.method(2)
                mkdir([srcDir filesep 'generatePower' filesep 'specparam_AveOverTrials' ...
                    filesep strrep(FileNames{currFile}, '.set', '_figures')]) ;
            end
            % Seperate "all" roi
%             if size(params.rois,2) > 1 && any(ismember(params.rois(1,:), 'all'))
%                 filteredROIindxs = setdiff(1:size(params.rois,2), ...
%                     find(strcmpi(params.rois(1,:), 'all'))) ;
%             else; filteredROIindxs = 1:size(params.rois,2) ;
%             end
            filteredROIindxs = 1:size(params.rois,2) ;

            % CREATE A BASE VARIABLE NAME SET THAT CAN BE USED OR ADAPTED DEPENDING ON
            % EXPORT TABLE FORMAT: [aper_offset, aper_exponent, peak#_centerFreq,
            % peak#_powerAboveAper, peak#_bandwidth, ..., error, r_squared,
            % BAND#_peaks]
            sp_outCols_orig = cell(1, params.specparam.settings.max_n_peaks*3) ;
            for i=1:params.specparam.settings.max_n_peaks
                sp_outCols_orig{1,i*3-2} = ['peak' num2str(i) '_centerFreq'] ;
                sp_outCols_orig{1,i*3-1} = ['peak' num2str(i) '_powerAboveAper'] ;
                sp_outCols_orig{1,i*3} = ['peak' num2str(i) '_bandwidth'] ;
            end
            if strcmp(params.specparam.settings.aperiodic_mode, 'knee')
                sp_outCols_orig = ['aper_offset', 'aper_knee', 'aper_exponent', sp_outCols_orig, ...
                'error', 'r_squared'] ;
            else
                sp_outCols_orig = ['aper_offset', 'aper_exponent', sp_outCols_orig, ...
                    'error', 'r_squared'] ;
            end
            if params.specparam.bands.on
                for i=1:size(params.specparam.bands.vals,1)
                    sp_outCols_orig = [sp_outCols_orig, [params.specparam.bands.vals{i,1} ...
                            '_peaks']] ;
                end
            end

            % Print out individual trials
            if params.specparam.method(1)
                sp_outRows = cell(EEG.trials+1,1) ;
                for i=1:EEG.trials
                    sp_outRows{i} = ['epoch_' num2str(i)] ;
                end
                sp_outRows{EEG.trials+1} = 'average_epoch' ;
            
                % Create base savename for the export file
                if size(params.rois,2) > 1 && params.specparam.splitROI.on
                    if params.specparam.splitROI.csvFormat
                        saveName = [srcDir filesep 'generatePower' filesep ...
                            'specparam_IndivTrials' filesep strrep(FileNames{currFile}, ...
                            '.set', '_indivTrialsSpecparam_ROI.csv')] ;
                    else
                        saveName = helpName([srcDir filesep 'generatePower' filesep ...
                            'specparam_IndivTrials' filesep strrep(FileNames{currFile}, ...
                            '.set', '_indivTrialsSpecparam.xlsx')], '.xlsx') ;
                    end
                else
                    saveName = helpName([srcDir filesep 'generatePower' filesep ...
                            'specparam_IndivTrials' filesep strrep(FileNames{currFile}, ...
                            '.set', '_indivTrialsSpecparam.csv')], '.csv') ;
                end
                if ~params.specparam.splitROI.on; sp_outCols_all = [] ; sp_allMets_all = []; end
            
                for currROI=1:size(filteredROIindxs,2)
                    roiSize = size(allSubs{currFile,4}{1, filteredROIindxs(currROI)},2) ;
                    sp_allMets = [] ;
                    sp_outCols = [] ;
                    for currChan=1:roiSize
                        sp_aper = [] ;
                        sp_peaks = [] ;
                        sp_error = cell(EEG.trials+1,1) ;
                        sp_r = cell(EEG.trials+1,1) ;
                        if params.specparam.bands.on
                            sp_bandpeaks = cell(EEG.trials+1,size(params.specparam.bands.vals,1)) ;
                        end
                        for currSeg=1:EEG.trials+1
                            sp_aper = [sp_aper; allSubs{currFile,4}{currSeg, ...
                                filteredROIindxs(currROI)}(currChan).aperiodic_params] ;
            
                            currCell = allSubs{currFile,4}{currSeg, ...
                                filteredROIindxs(currROI)}(currChan).peak_params ;
                            sp_peaks_temp = [] ;
                            for i=1:size(currCell,1)
                                sp_peaks_temp = [sp_peaks_temp currCell(i,:)] ;
                            end
                            if size(sp_peaks_temp,2) < 3*params.specparam.settings.max_n_peaks
                                sp_peaks_temp(1, size(sp_peaks_temp,2)+1:3*params.specparam.settings.max_n_peaks) = NaN ;
                            end
                            sp_peaks = [sp_peaks; sp_peaks_temp] ;
            
                            sp_error{currSeg} = allSubs{currFile,4}{currSeg, ...
                                filteredROIindxs(currROI)}(currChan).error ;
                            sp_r{currSeg} = allSubs{currFile,4}{currSeg, ...
                                filteredROIindxs(currROI)}(currChan).r_squared ;
            
                            if params.specparam.bands.on
                                for i=1:size(params.specparam.bands.vals,1)
                                    for j=1:size(currCell,1)
                                        if (str2double(params.specparam.bands.vals{i,2}) < ...
                                                currCell(j,1)) && (currCell(j,1) < ...
                                                str2double(params.specparam.bands.vals{i,3}))
                                            if isempty(sp_bandpeaks{currSeg, i})
                                                sp_bandpeaks{currSeg,i} = num2str(currCell(j,1)) ;
                                            else
                                                sp_bandpeaks{currSeg, i} = [sp_bandpeaks{currSeg, ...
                                                    i} ', ' num2str(currCell(j,1))] ;
                                            end
                                        end
                                    end
                                    if isempty(sp_bandpeaks{currSeg, i})
                                        sp_bandpeaks{currSeg, i} = 'NaN' ;
                                    end
                                end
                            end
                        end
                        sp_outCols_chan = cell(1, size(sp_outCols_orig,2)) ;
                        if currChan > roiSize-1; addOn = 'Average' ;
                        else; addOn = varNames{chanIndxs{filteredROIindxs(currROI)}(currChan)} ;
                        end
                        for i=1:size(sp_outCols_orig,2)
                            sp_outCols_chan{1,i} = [addOn '_' sp_outCols_orig{1,i}] ;
                        end
                        sp_allMets_temp = [num2cell(sp_aper) num2cell(sp_peaks) sp_error sp_r] ;
                        if params.specparam.bands.on
                            sp_allMets_temp = [sp_allMets_temp sp_bandpeaks] ;
                        end
                        sp_allMets = [sp_allMets sp_allMets_temp] ;
                        sp_outCols = [sp_outCols sp_outCols_chan] ;
                    end
                    if params.specparam.splitROI.on
                        % append tag to savename and print
                        if params.specparam.splitROI.csvFormat
                            saveName2  = strrep(saveName, 'ROI', params.rois{3, ...
                                filteredROIindxs(currROI)}) ;
                            writetable(array2table(sp_allMets, 'VariableNames', ...
                                sp_outCols, 'RowNames', sp_outRows), ...
                                helpName(saveName2, '.csv'), 'WriteRowNames', true, ...
                                'QuoteStrings', true) ;
                        else
                            writetable(array2table(sp_allMets, 'VariableNames', ...
                                sp_outCols, 'RowNames', sp_outRows), saveName, 'Sheet', ...
                                params.rois{3, filteredROIindxs(currROI)}, ...
                                'WriteRowNames', true) ;
                        end
                    else
                        % affix ROI to VarNames and concat output matrix
                        for i=1:size(sp_outCols, 2)
                            sp_outCols{1,i} = [params.rois{3, filteredROIindxs(currROI)} ...
                                '_' sp_outCols{1,i}] ;
                        end
                        sp_outCols_all = [sp_outCols_all sp_outCols] ;
                        sp_allMets_all = [sp_allMets_all sp_allMets] ;
                    end
                    if currFile == 1
                        subLvl_specparam{2,currROI} = sp_outCols ;
                    end
                    subLvl_specparam{1,currROI} = [subLvl_specparam{1,currROI}; sp_allMets(end,:)] ;
                end
                if ~params.specparam.splitROI.on
                    writetable(array2table(sp_allMets_all, 'VariableNames', ...
                        sp_outCols_all, 'RowNames', sp_outRows), saveName, ...
                        'WriteRowNames', true, 'QuoteStrings', true) ;
                end
            else
                if currFile==1
                    for currROI=1:size(filteredROIindxs,2)
                        roiSize = size(allSubs{currFile,4}{1, filteredROIindxs(currROI)},2) ;
                        for currChan=1:roiSize
                            sp_outCols_chan = cell(1, size(sp_outCols_orig,2)) ;
                            if currChan > roiSize-1; addOn = 'aveROI' ;
                            else; addOn = varNames{chanIndxs{filteredROIindxs(currROI)}(currChan)} ;
                            end
                            for i=1:size(sp_outCols_orig,2)
                                sp_outCols_chan{1,i} = [addOn '_' sp_outCols_orig{1,i}] ;
                            end
                            subLvl_specparam{2,currROI} = [subLvl_specparam{2,currROI} sp_outCols_chan] ;
                        end
                    end
                end
            end

            % Print out average over trials
            if params.specparam.method(2)
                % Create base savename for the export file
                if size(params.rois,2) > 1 && params.specparam.splitROI.on
                    if params.specparam.splitROI.csvFormat
                        saveName = [srcDir filesep 'generatePower' filesep ...
                            'specparam_AveOverTrials' filesep strrep(FileNames{currFile}, ...
                            '.set', '_aveSpecparam_ROI.csv')] ;
                    else
                        saveName = helpName([srcDir filesep 'generatePower' filesep ...
                            'specparam_AveOverTrials' filesep strrep(FileNames{currFile}, ...
                            '.set', '_aveSpecparam.xlsx')], '.xlsx') ;
                    end
                else
                    saveName = helpName([srcDir filesep 'generatePower' filesep ...
                            'specparam_AveOverTrials' filesep strrep(FileNames{currFile}, ...
                            '.set', '_aveSpecparam.csv')], '.csv') ;
                end
            
%                 if params.specparam.splitROI.on
                    sp_outMat = [] ;
%                 end
                sp_outRows = [] ;
                for currROI=1:size(filteredROIindxs,2)
                    roiSize = size(allSubs{currFile,4}{1, filteredROIindxs(currROI)},2) ;
                    % DO MATRIX CALCULATION STUFF
                    % Aperiodic Info
                    currCell = {allSubs{currFile,4}{end, filteredROIindxs(currROI)}.aperiodic_params} ;
                    sp_aper = [] ;
                    for i=1:size(currCell,2)
                        sp_aper = [sp_aper; currCell{i}] ;                      
                    end
                    
                    % Peak Info
                    currCell = {allSubs{currFile,4}{end, filteredROIindxs(currROI)}.peak_params} ;
                    sp_peaks = zeros(size(currCell,2),params.specparam.settings.max_n_peaks*3) ;
                    for i=1:size(currCell,2)
                        for j=1:params.specparam.settings.max_n_peaks
                            if j<=size(currCell{i},1)
                                sp_peaks(i,j*3-2:j*3) = currCell{i}(j,:) ;
                            else; sp_peaks(i,j*3-2:j*3) = NaN(1,3) ;
                            end
                        end
                    end
            
                    sp_outMat_temp = [sp_aper sp_peaks cell2mat({allSubs{currFile,4}{end, ...
                        filteredROIindxs(currROI)}.error}') ...
                        cell2mat({allSubs{currFile,4}{end, ...
                        filteredROIindxs(currROI)}.r_squared}')] ;
            
                    if params.specparam.bands.on
                        bandPeaks = cell(size(sp_peaks,1), size(params.specparam.bands.vals,1)) ;
                        for currBand=1:size(params.specparam.bands.vals,1)
                            for i=1:size(sp_peaks,1)
                                for j=1:params.specparam.settings.max_n_peaks
                                    if sp_peaks(i,j*3-2) > str2double(params.specparam.bands.vals{currBand,2}) && ...
                                        sp_peaks(i,j*3-2) < str2double(params.specparam.bands.vals{currBand,3})
                                        if isempty(bandPeaks{i, currBand})
                                            bandPeaks{i, currBand} = num2str(sp_peaks(i,j*3-2)) ;
                                        else; bandPeaks{i, currBand} = [bandPeaks{i, ...
                                                currBand} ', ' num2str(sp_peaks(i,j*3-2))] ;
                                        end
                                    end
                                end
                                if isempty(bandPeaks{i, currBand})
                                    bandPeaks{i, currBand} = 'NaN' ;
                                end
                            end 
                        end
                        sp_outMat_temp = [num2cell(sp_outMat_temp), bandPeaks] ;
                    end
                    % IF SPLITTING, PRINT WITH TAG APPENDED TO SAVENAME
                    if params.specparam.splitROI.on
                        sp_outRows = [varNames(setdiff(chanIndxs{filteredROIindxs(currROI)}, ...
                            zeroChanIndxs)) 'Average'] ;
                        % update savename
                        if params.specparam.splitROI.csvFormat
                            saveName2 = strrep(saveName, 'ROI', ...
                                params.rois{3, filteredROIindxs(currROI)}) ;
                            writetable(array2table(sp_outMat_temp, 'VariableNames', ...
                                sp_outCols_orig, 'RowNames', sp_outRows'), ...
                                helpName(saveName2, '.csv'), 'WriteRowNames', true, ...
                                'QuoteStrings', true) ;
                        else
                            writetable(array2table(sp_outMat_temp, 'VariableNames', ...
                                sp_outCols_orig, 'RowNames', sp_outRows'), ...
                                saveName, 'Sheet', params.rois{3, filteredROIindxs(currROI)}, ...
                                'WriteRowNames', true) ;
                        end
                    % IF NOT SPLITTING, APPEND ROI TAG TO ROW NAMES
                    else
                        % Update Row Names
                        sp_outRows_temp = [varNames(setdiff(chanIndxs{filteredROIindxs(currROI)}, ...
                            zeroChanIndxs)) 'Average'] ;
                        for i=1:size(sp_outRows_temp,2)
                            sp_outRows_temp{i} = [params.rois{3,filteredROIindxs(currROI)} ...
                                '_' sp_outRows_temp{i}] ;
                        end
                        sp_outRows = [sp_outRows sp_outRows_temp] ;
            
                        % Concat and update output matrix
                        sp_outMat = [sp_outMat; sp_outMat_temp] ;
                    end
                    if ~params.specparam.method(1)
                        tempBand = [] ;
                        for i=1:size(sp_outMat_temp,1)
                            tempBand = [tempBand sp_outMat_temp(i,:)] ;
                        end
                        subLvl_specparam{1,currROI} = [subLvl_specparam{1,...
                            currROI}; tempBand] ; % TODO: check cases where num2cell(tempBand) is valid 
                    end
                    % Visualizations for chans
                    figName = [srcDir filesep 'generatePower' filesep ...
                            'specparam_AveOverTrials' filesep ...
                            strrep(FileNames{currFile}, '.set', '_figures') ...
                            filesep 'ROI_CHAN'] ;
                    for i=1:size(allSubs{currFile,4}{end, filteredROIindxs(currROI)},2)
                        plotSpecparam(allSubs{currFile,4}{end,filteredROIindxs(currROI)}(i), ...
                            params.specparam.vis) ;
                        if i < size(allSubs{currFile,4}{end, filteredROIindxs(currROI)},2)
                            figName2 = strrep(figName, 'CHAN', ...
                                varNames{chanIndxs{filteredROIindxs(currROI)}(i)}) ;
                        else; figName2 = strrep(figName, 'CHAN', 'average') ;
                        end
                        saveas(gcf, [strrep(figName2, 'ROI', params.rois{3, ...
                            filteredROIindxs(currROI)}) '.png']) ;
                    end
                end
                if ~params.specparam.splitROI.on
                    writetable(array2table(sp_outMat, 'VariableNames', sp_outCols_orig, ...
                        'RowNames', sp_outRows'), saveName, 'WriteRowNames', true, ...
                        'QuoteStrings', true) ;
                end
            end

        end
    catch ME
        name = {ME.stack.name} ;
        line = {ME.stack.line} ;
        if size(line,2) > 1
            errList = [sprintf('Line %d in %s; ', line{1:end-1}, ...
                name{1:end-1}) 'Line ' num2str(line{end}) ' in ' name{end}] ;
        else; errList = ['Line ' num2str(line{1}) ' in ' name] ;
        end
        errorLog = [errorLog; {FileNames{currFile}, ME.message, errList}] ;
    end
end

%% SAVE ERROR LOG
% If there were any errors while running HAPPE, save an error log so the
% user can troubleshoot.
if ~isempty(errorLog)
    fprintf('Saving error log...\n') ;
    errTabName = helpName(['HAPPE_errorLog_' datestr(now, 'dd-mm-yyyy') ...
        '.csv'], '.csv') ;
    writetable(cell2table(errorLog, 'VariableNames', {'File', ...
        'Error Message' 'Stack Trace'}), errTabName) ;
end

%% PRINT OUT ALLSUBS FILE
if ~isfolder([srcDir filesep 'generatePower' filesep 'allSubs'])
    mkdir([srcDir filesep 'generatePower' filesep 'allSubs']) ;
end
fprintf('Saving information from all files...\n') ;
for currROI=1:size(params.rois,2)
    writetable(array2table(subLvl_PSD{currROI}(freqMin:freqMax,:), ...
        'VariableNames', FileNames, 'RowNames', cellstr(num2str(f(freqMin:freqMax)))), ...
        helpName([srcDir filesep 'generatePower' filesep 'allSubs' filesep ...
        'allSubs_PSD_' params.rois{3, currROI} '.csv'], '.csv'), ...
        'WriteRowNames', true, 'QuoteStrings', true) ;

    if params.bands.on
        writetable(array2table(subLvl_bandpower{1,currROI}, 'VariableNames', ...
            subLvl_bandpower{2,currROI}, 'RowNames', FileNames'), helpName([srcDir ...
            filesep 'generatePower' filesep 'allSubs' filesep 'allSubs_bandpower_' ...
            params.rois{3, currROI} '.csv'], '.csv'), 'WriteRowNames', true, ...
            'QuoteStrings', true) ;
    end

    if params.specparam.on && ismember(currROI, filteredROIindxs)
        writetable(array2table(subLvl_specparam{1,currROI}, 'VariableNames', ...
            subLvl_specparam{2,currROI}, 'RowNames', FileNames'), ...
            helpName([srcDir filesep 'generatePower' filesep 'allSubs' ...
            filesep 'allSubs_specparam_' params.rois{3,currROI} '.csv'], ...
            '.csv'), 'WriteRowNames', true, 'QuoteStrings', true) ;

    end
end

%% SUPPORT FUNCTIONS
% DETERMINE IF CALCULATING FOR INDIVIDUAL, AVERAGE, OR BOTH TRIALS
function [method, csvFormat] = determ_IndivAveBoth()
    fprintf(['Calculate power:\n  individual = By individual trials/segments\n' ...
        '  average = For the average over trials/segments\n  both = Both ' ...
        'individual trials/segments and average over trials/segments\n']) ;
    while true
        ui = input('> ','s') ;
        if strcmpi(ui, 'individual'); method = [1,0]; break;
        elseif strcmpi(ui, 'average'); method = [0,1]; break;
        elseif strcmpi(ui, 'both'); method = [1,1]; break;
        else; fprintf(['Invalid input: please enter "individual", "average"' ...
                ', or "both" (without quotations)\n']) ;
        end
    end
    if method(1)
        fprintf(['Print the individual trials/segments as seperate .csv ' ...
            'files or as\nsheets in a single Excel file?\n  csv = Output ' ...
            'seperate .csv files\n  sheets = Output a single Excel file ' ...
            'with multiple sheets\n']) ;
        while true
            ui = input('> ','s') ;
            if strcmpi(ui, 'csv'); csvFormat = 1; break;
            elseif strcmpi(ui, 'sheets'); csvFormat = 0; break;
            else; fprintf(['Invalid input: please enter "csv" or "sheets"' ...
                    ' (without quotations).']) ;
            end
        end
    else; csvFormat = 1;
    end
end

% LIST OUT PARAMETER SET
function genPower_listParams(params)
    fprintf('Calculate PSD for Individual Trials/Segments: ') ;
    if params.method(1); fprintf('On\n - Output Format: ') ;
        if params.csvFormat; fprintf('Seperate .csv files\n') ;
        else; fprintf('Sheets in an Excel file\n') ;
        end
    else; fprintf('Off\n') ;
    end
    
    fprintf('Calculate PSD for Average of Trials/Segments: ') ;
    if params.method(2); fprintf('On\n') ;
    else; fprintf('Off\n') ;
    end
    
    fprintf('Limit PSD Frequencies: ') ;
    if params.freqs.limit
        fprintf(['On\n - Minimum: ' num2str(params.freqs.min) ' Hz\n - ' ...
            'Maximum: ' num2str(params.freqs.max) ' Hz\n']) ;
    else; fprintf('Off\n') ;
    end
    
    fprintf('Channels of Interest:\n') ;
    for i=1:size(params.rois,2)
        fprintf(' - ') ;
        if strcmpi(params.rois{1,i}, 'all')
            fprintf([params.rois{3,i} ': all channels\n']) ;
        elseif strcmpi(params.rois{1,i}, 'coi_include')
            fprintf([params.rois{3,i} ': ' sprintf('%s, ', params.rois{2,...
                i}{1:end-1}) params.rois{2,i}{end} '\n']) ;
        elseif strcmpi(params.rois{1,i}, 'coi_exclude')
            fprintf([params.rois{3,i} ': all except' sprintf('%s, ', ...
                params.rois{2,i}{1:end-1}) params.rois{2,i}{end} '\n']) ;
        end
    end

    fprintf('Calculate PSD in Frequency Bands: ') ;
    if params.bands.on
        fprintf('On\n - Bands:\n') ;
        for i=1:size(params.bands.vals,1)
            fprintf(['    - ' params.bands.vals{i,1} ': ' params.bands.vals{i,2} ...
                ' Hz - ' params.bands.vals{i,3} ' Hz\n']) ;
        end
        fprintf(' - Calculate for Individual Trials: ') ;
        if params.bands.method(1); fprintf('On\n    - Output Format: ') ;
            if params.csvFormat; fprintf('Seperate .csv files\n') ;
            else; fprintf('Sheets in an Excel file\n') ;
            end
        else; fprintf('Off\n') ;
        end
        fprintf(' - Calculate for Average of Trials: ') ;
        if params.bands.method(2); fprintf('On\n') ;
        else; fprintf('Off\n') ;
        end
    else; fprintf('Off\n') ;
    end
    
    fprintf('Calculate Specparam: ') ;
    if params.specparam.on
        fprintf(['On\n - Python Version: ' num2str(params.specparam.pyVers) ...
            '\n - Calculate for Individual Trials/Segments: ']) ;
        if params.specparam.method(1); fprintf('On\n'); else; fprintf('Off\n') ; end
        fprintf(' - Calculate for Average of Trials/Segments: ');
        if params.specparam.method(2); fprintf('On\n'); else; fprintf('Off\n') ; end
        fprintf(' - Visualizations: ')
        if params.specparam.vis; fprintf('On\n'); else; fprintf('Off\n') ; end
        if ispc; fprintf([' - Python Version: ' params.specparam.pyVers '\n']) ;
        else; fprintf([' - Path to Python Executable: ' params.specparam.pyVers ...
                '/n']) ;
        end
        fprintf([' - Spectrum Limits: ' num2str(params.specparam.min) ' Hz - ' ...
            num2str(params.specparam.max) ' Hz\n - Peak Width Limits: ' ...
            num2str(params.specparam.settings.peak_width_limits(1)) ' - ' ...
            num2str(params.specparam.settings.peak_width_limits(2)) '\n - ' ...
            'Max Number of Peaks: ' num2str(params.specparam.settings.max_n_peaks) ...
            '\n - Minimum Peak Height: ' num2str(params.specparam.settings.min_peak_height) ...
            '\n - Peak Threshold: ' num2str(params.specparam.settings.peak_threshold) ...
            '\n - Aperiodic Mode: ' params.specparam.settings.aperiodic_mode '\n']) ;
    else; fprintf('Off\n') ;
    end
end