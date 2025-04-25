%% generatePAC
% Generate Phase-Amplitude Coupling (PAC) - A post-processing script to
% calculate Phase-Amplitude Coupling values for within-electrode and
% between-electrode analyses.
% Relies on functions from the EEGLAB PACTools wrapper plugin by Ramon Martinez-Cancino.
% (Martinez-Cancino et al. 2020 PMC6342492).
%
%
% If using this script for analyses in publications, cite both HAPPE and
% the PACTools plugin.
%
%       Martinez-Cancino, Heng, Delorme, Kreutz-Delgado, Sotero, Makeig (2020). 
%   Meauring transient phase-amplitude coupling using local mutual information. NeuroImage, 361-378. 
%   DOI: 10.1016/j.neuroimage.2018.10.034
%
% Developed at Northeastern University's PINE Lab
%
% Authors: Joshua Paul Rodriguez, Tess Forest, PINE Lab at Northeastern
% University, 2025
%
% This file is part of HAPPE.
% Copyright 2018-2025 
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
fprintf('Preparing HAPPE generatePAC...\n') ;
happeDir = strrep(fileparts(which(mfilename('fullpath'))), [filesep '3. ' ...
    'generate'], '') ;
eeglabDir = [happeDir filesep 'packages' filesep 'eeglab2024.0'] ;
addpath([happeDir filesep '3. generate'], ...
    [happeDir filesep 'files'], ...
    [happeDir filesep 'scripts'], ...
    [happeDir filesep 'scripts' filesep 'ui'], ...
    [happeDir filesep 'scripts' filesep 'support'], ...
    eeglabDir, ...
    genpath([eeglabDir filesep 'functions']), ...
    genpath([eeglabDir filesep ...
    'plugins' filesep 'PACTools'])) ;
pluginDir = dir([eeglabDir filesep 'plugins']) ;
pluginDir = strcat(eeglabDir, filesep, 'plugins', filesep, {pluginDir.name}, ';') ;
addpath([pluginDir{:}]) ;
rmpath(genpath([eeglabDir filesep 'functions' filesep 'octavefunc'])) ;

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
if ~isfolder([srcDir filesep 'generatePAC']); mkdir([srcDir filesep 'generatePAC']); end
addpath('generatePAC') ;

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
    try genPAC_listParams(params) ;
    catch; error('ERROR: This is not a valid set of parameters for this script\n') ;
    end
    fprintf('Change an existing parameter? [Y/N]\n') ;
    changedParams = choose2('N', 'Y') ;

% If no existing parameters are loaded, create a set of default parameters
% to be edited/filled out during the Set Parameters step below.
else; changedParams = 0;     
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
        fprintf(['Parameter to change: channels of interest, pac channels within/between, ' ...
            'low frequency phase range, high freqency amplitude range,\n' ...
            'number of frequency bins for phase range, ' ...
            'number of frequency bins for amplitude range,\n' ...
            'frequency scale mode, PAC method, number of bootstraps, alpha,\n' ...
            'bonforoni correction, time limits,\n' ...
            'Enter "done" (without quotations) when finished changing ' ...
            'parameters.\n']) ;
        userChoice = input('> ', 's') ;
    end

    % INDIVIDUAL OUTPUTS, AVERAGE OUTPUTS, OR BOTH:
    % Determine whether to calculate PAC between each channel of interest
    % or between all channels of interest.
    if ~preExist || strcmpi(userChoice, 'pac channels within/between')
        [params.method, params.csvFormat] = determ_WithinBetween() ;
    end

    % CHANNELS OF INTEREST/REGIONS OF INTEREST:
    % Determine the collections of channels of interest in which to examine
    % power measures. Enable the selection of all channels, or a subset,
    % and for multiple regions of interest.
    if ~preExist || strcmpi('channels of interest', userChoice)
        [params.chansAll, params.chanIDs] = determ_chanIDs() ;
    end

    % Determine Bounds of Phase Amplitude Spectrum:
    % Determine the upper and lower limits, in Hz, for PAC
    % spectrum calculations.
    if ~preExist || strcmpi(userChoice, 'low frequency phase range') || ...
            strcmpi(userChoice, 'phase range')           
        fprintf('Enter phase lower bound of spectrum to fit, in Hz:\nExample: 2\n') ;
        params.pac.phase_range(1)  = input('> ') ;
        fprintf('Enter phase upper bound of spectrum to fit, in Hz:\nExample: 20\n') ;
        params.pac.phase_range(2) = input('> ');
    end
    % Determine Bounds of Amplitude:
    % Determine the upper and lower limits, in Hz, for PAC
    % spectrum calculations.
    if ~preExist || strcmpi(userChoice, 'high freqency amplitude range') || ...
            strcmpi(userChoice, 'amplitude range')
        fprintf('Enter amplitude lower bound of spectrum to fit:\nExample: 30\n') ;
        params.pac.amplitude_range(1)  = input('> ') ;
        fprintf('Enter amplitude upper bound of spectrum to fit:\nExample: 90\n') ;
        params.pac.amplitude_range(2) = input('> ');
    end

    % Determine Frequency Bins
    if ~preExist || strcmpi(userChoice, 'number of frequency bins for phase range') ||...
            strcmpi(userChoice, 'phase bins')
        fprintf('Enter number of frequency bins for phase range:\n') ;
        params.pac.phase_bins = input('> ') ;
    end
    % Determine Amplitude Bins
    if ~preExist || strcmpi(userChoice, 'number of frequency bins for amplitude range') ||...
            strcmpi(userChoice, 'amplitude bins')
        fprintf('Enter number of frequency bins for amplitude range:\n') ;
        params.pac.amplitude_bins = input('> ') ;
    end

    % Determine Frequency Scale Mode
    if ~preExist || strcmpi(userChoice, 'frequency scale mode')
        fprintf(['Frequency scale mode?\n  linear = linear mode\n  ' ...
            'log = logarithmic mode\n']) ;
        while true
            ui = input('> ','s') ;
            if strcmpi(ui, 'linear')
                params.pac.scale_mode = 'linear' ;
                break ;
            elseif strcmpi(ui, 'log')
                params.pac.scale_mode = 'log' ;
                break ;
            else; fprintf(['Invalid Input: Please enter "linear" or "log"' ...
                    '(without quotations).']) ;
            end
        end
    end

    % Determine Frequency Scale Mode
    if ~preExist || strcmpi(userChoice, 'PAC method')
        fprintf(['PAC method?\n  mvlmi = Mean Vector Length Modulation Index\n  ' ...
            'klmi = Kullback-Leibler Modulation Index\n  ' ...
            'glm = Generalized Linear Model\n  ' ...
            'ermipac = Event related MIPAC\n  ' ...
            'instmipac = Instantaneous MIPAC\n']) ;
        while true
            ui = input('> ','s') ;
            if strcmpi(ui, 'mvlmi')
                params.pac.method = 'mvlmi' ;
                break ;
            elseif strcmpi(ui, 'klmi')
                params.pac.method = 'klmi' ;
                break ;
            elseif strcmpi(ui, 'glm')
                params.pac.method = 'glm' ;
                break ;
            elseif strcmpi(ui, 'ermipac')
                params.pac.method = 'ermipac' ;
                break ;
            elseif strcmpi(ui, 'instmipac')
                params.pac.method = 'instmipac' ;
                break ;
            else; fprintf(['Invalid Input: Please enter "mvlmi", "klmi", "glm", "ermipac", or "instmipac"' ...
                    '(without quotations).']) ;
            end
        end
    end

    % Determine Number of Bootstraps
    if ~preExist || strcmpi(userChoice, 'number of bootstraps')
            fprintf('Enter number of bootstraps:\n') ;
            params.pac.n_boots = input('> ') ;
    end

    % Determine Alpha
    if ~preExist || strcmpi(userChoice, 'alpha')
            fprintf('Enter alpha value:\n') ;
            params.pac.alpha = input('> ') ;
    end

    % Determine whether to use bonforoni correction
    if ~preExist || strcmpi(userChoice, 'bonforoni correction')
        fprintf(['Apply Bonforoni Correction? [Y/N]\n']) ;
        params.pac.bonfcorr = choose2('N','Y') ;
    end

    % Determine Time Limits
    if ~preExist || strcmpi(userChoice, 'time limits')
        fprintf('Time Limits\nUse defaults? ([0 number of time points / sampling rate]) [Y/N]\n') ;
        if choose2('N','Y')
            params.pac.time_limits = {} ;
        else
            clear params.pac.time_limits
            fprintf('Enter lower time limit:\n') ;
            params.pac.time_limits(1) = input('> ') ;
            fprintf('Enter upper time limit:\n') ;
            params.pac.time_limits(2) = input('> ') ;
        end
    end

    % DONE:
    % When finished entering parameters, ask the user to review their
    % parameters and make any necessary corrections.
    if ~preExist || strcmpi('done', userChoice)
        fprintf('Please check your input parameters before continuing.\n') ;
        genPAC_listParams(params) ;
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
        paramFile = paramFile_validateExist(['genPAC_parameters_' ...
            datestr(now, 'dd-mm-yyyy') '.mat'], 'genPAC_parameters_', 2) ;
    else
        fprintf('File name (Do not include .mat):\n') ;
        paramFile = paramFile_validateExist([input('> ', 's') '.mat'], ...
            'genPAC_parameters_', 0) ;
    end

    % Save the params variable to a .mat file using the name created above.
    fprintf('Saving parameters...\n') ; 
    save([srcDir filesep 'generatePAC' filesep paramFile], 'params') ;
    fprintf('Parameters saved.\n') ;
end
clear('preExist', 'changedParams', 'paramFile') ;

% %% CREATE OUTPUT FOLDERS P.2
% if params.method(1) && ~isfolder([srcDir filesep 'generatePAC' filesep ...
%         'IndivTrials'])
%     mkdir([srcDir filesep 'generatePAC' filesep 'IndivTrials']) ;
% end
% if params.method(2) && ~isfolder([srcDir filesep 'generatePAC' filesep ...
%         'AveOverTrials'])
%     mkdir([srcDir filesep 'generatePAC' filesep 'AveOverTrials']) ;
% end
% fprintf('Output folders created.\n') ;

%% COLLECT FILES
% Locates all .set files to run through the script. If no files meeting
% this criteria are found, throw an error and end the run.
fprintf('Gathering files...\n') ;
FileNames = {dir('*.set').name} ;
if isempty(FileNames); error('ERROR: NO .SET FILES FOUND.'); end

%% SUBJECT LEVEL MATRIX
% Create a matrix to hold the calculated measures across subjects.
allSubs = cell(size(FileNames,2),1) ;
% subLvl_PAC = cell(1, 1) ;
allSub_freqs_amp = cell(1) ;
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
        if ~isempty(EEG.chanlocs)
            % If the ROI requests all channels, simply use all possibly
            % indicies.
            if strcmpi(params.chansAll, 'all')
                chanIndxs = 1:size(EEG.data,1) ;
            % If selecting to include or exclude channels, find use the
            % channel locations to get the index associated with the
            % relevant channel in the data matrix.
            else
                temp = [] ;
                for i=1:size(params.chanIDs, 2)
                    temp = [temp find(ismember({EEG.chanlocs.labels}, ...
                        params.chanIDs{i}))] ;                                 
                end
                if strcmpi(params.chansAll, 'coi_exclude')
                    chanIndxs = setdiff(1:size(EEG.data,1), ...
                        temp) ;
                elseif strcmpi(params.chansAll, 'coi_include')
                    chanIndxs = temp ;
                else
                    error('Invalid channel COI selection.') ;
                end
            end
        else
            fprintf(['Cannot select ROIs for data without channel names. ' ...
                'Will use all existing channels. TODO:CHANGE']) ;
            chanIndxs = 1:size(EEG.data,1) ;
        end

        % n_epochs = EEG.trials ;
        n_epochs = 2 ;


        %% RUN POP PAC

        % WITHIN CHANNELS
        if params.method(1) && ~params.method(2)
            within_pac_values = cell(n_epochs, 1) ;
            % Per epoch
            for ep = 1:n_epochs

                data = pop_select(EEG, 'trial', ep) ;

                fprintf(['Calculating within channels PAC, Epoch:' num2str(ep) ...
                    '/' num2str(n_epochs) ' ...\n']) ;
                if size(params.pac.time_limits, 2) <= 1
                    within_pac_values_struct = pop_pac(data,'Channels', ...
                        params.pac.phase_range, params.pac.amplitude_range, chanIndxs, chanIndxs, ...
                        'nfreqs1', params.pac.phase_bins, ...
                        'nfreqs2', params.pac.amplitude_bins, ...
                        'freqscale', params.pac.scale_mode, ...
                        'method', params.pac.method, ...
                        'nboot', params.pac.n_boots, ...
                        'alpha', params.pac.alpha, ...
                        'bonfcorr', params.pac.bonfcorr).etc.eegpac ;
                else
                    within_pac_values_struct = pop_pac(data,'Channels', ...
                        params.pac.phase_range, params.pac.amplitude_range, chanIndxs, chanIndxs, ...
                        'nfreqs1', params.pac.phase_bins, ...
                        'nfreqs2', params.pac.amplitude_bins, ...
                        'freqscale', params.pac.scale_mode, ...
                        'method', params.pac.method, ...
                        'nboot', params.pac.n_boots, ...
                        'alpha', params.pac.alpha, ...
                        'bonfcorr', params.pac.bonfcorr, ...
                        'tlimits', params.pac.time_limits).etc.eegpac ;
                end

                % Convert to common form
                within_pac_values_at_ep = cell(1, length(within_pac_values_struct)) ;
                for i = 1:length(within_pac_values_struct)
                    within_pac_values_at_ep{i} = within_pac_values_struct(i);
                end

                within_pac_values{ep} = within_pac_values_at_ep ;
            end
            pac_values = within_pac_values ;
            fprintf('Within channels PAC finished.\n') ;
        end


        % BETWEEN CHANNELS
        if params.method(2)
            between_pac_values = cell(n_epochs, 1) ;
            % Per epoch
            for ep = 1:n_epochs

                fprintf(['Calculating between channels PAC, epoch:' num2str(ep) ...
                    '/' num2str(n_epochs) '...\n']) ;

                pac_value_combo = cell((length(chanIndxs)*length(chanIndxs)), 1) ;
                acc_indx = 1 ;
                data = pop_select(EEG, 'trial', ep) ;
                for chanIndx_phase = chanIndxs
                    for chanIndx_amp = chanIndxs
                        fprintf(['Set-' num2str(currFile) ' Epoch ' num2str(ep) ...
                            '/' num2str(n_epochs) ' | Phase: ' ...
                            num2str(chanIndx_phase) ', Amp: ' num2str(chanIndx_amp) '\n']) ;
                        if size(params.pac.time_limits, 2) <= 1
                            pac_value_combo{acc_indx} = pop_pac(data,'Channels', ...
                                params.pac.phase_range, params.pac.amplitude_range, ...
                                chanIndx_phase, chanIndx_amp, ...
                                'nfreqs1', params.pac.phase_bins, ...
                                'nfreqs2', params.pac.amplitude_bins, ...
                                'freqscale', params.pac.scale_mode, ...
                                'method', params.pac.method, ...
                                'nboot', params.pac.n_boots, ...
                                'alpha', params.pac.alpha, ...
                                'bonfcorr', params.pac.bonfcorr).etc.eegpac ;
                        else
                            pac_value_combo{acc_indx} = pop_pac(data,'Channels', ...
                                params.pac.phase_range, params.pac.amplitude_range, ...
                                chanIndx_phase, chanIndx_amp, ...
                                'nfreqs1', params.pac.phase_bins, ...
                                'nfreqs2', params.pac.amplitude_bins, ...
                                'freqscale', params.pac.scale_mode, ...
                                'method', params.pac.method, ...
                                'nboot', params.pac.n_boots, ...
                                'alpha', params.pac.alpha, ...
                                'bonfcorr', params.pac.bonfcorr, ...
                                'tlimits', params.pac.time_limits).etc.eegpac ;
                        end

                        acc_indx = acc_indx + 1 ;
                    end
                end
                between_pac_values{ep} = pac_value_combo ;
            end
            pac_values = between_pac_values ;
            fprintf('Between channels PAC finished.\n') ;
        end


        % FORMAT MATRIX
        num_epairs = size(pac_values{1, 1}, 1) ;
        num_trs = length(pac_values) ;
        pac_file_mat = cell((n_epochs * num_epairs ...
            * params.pac.phase_bins), (5 + 2 * params.pac.amplitude_bins)) ;
        for tr = 1:n_epochs
            pac_value = pac_values{tr, 1} ;
            for e = 1:num_epairs
                epair = pac_value{e} ; % this gets the current line of the pac_values matrix
                % Extract phase and amplitude electrode indices
                phase_electrode = epair.dataindx(1,1) ; % first value from data index of current line
                amp_electrode = epair.dataindx(1,2) ;
                % Extract PAC matrix for this electrode pairing
                pac_matrix = epair.klmi.pacval ; 
                nbins_matrix = epair.klmi.nbinskl ;
                % Extract frequency labels for this specific pairing
                freqs_phase = epair.params.freqs_phase ; % Row names
                freqs_amp = epair.params.freqs_amp ; % Column names

                % Initialize subject level freq cols if haven't
                if size(allSub_freqs_amp, 2) <= 2
                    allSub_freqs_amp = freqs_amp ;
                end

                if size(pac_matrix, 1) ~= size(nbins_matrix, 1)
                    error('Incongruent lengh of frenquency phases');
                end
                pac_pair_mat =  [repmat(tr, size(pac_matrix, 1), 1) ...
                    repmat(phase_electrode, size(pac_matrix, 1), 1) ...
                    repmat(amp_electrode, size(pac_matrix, 1), 1) ...
                    freqs_phase(:) ...
                    pac_matrix nbins_matrix] ;
                pac_file_mat((tr-1)*num_epairs*params.pac.phase_bins ...
                    + (e-1)*params.pac.phase_bins + 1:(tr-1)*num_epairs*params.pac.phase_bins ...
                    + e*params.pac.phase_bins, 2:end) = ...
                    num2cell(pac_pair_mat) ;

            end
        end
        [pac_file_mat{:,1}] = deal(FileNames{currFile}) ;

        % PRINT OUT INDIVIDUAL FILE
        fprintf('Saving information from current files...\n') ;
        saveName = helpName([srcDir filesep 'generatePAC' filesep ...
            strrep(FileNames{currFile}, ...
            '.set', '_pac.csv')], '.csv') ;
        writetable(array2table(pac_file_mat(:,2:end), ...
                'VariableNames', ['epoch' 'phase_electrode' 'freq_phase' ...
                    'amp_electrode' compose('PAC_Amp_%gHz', allSub_freqs_amp) ...
                    compose('nbins_Amp_%gHz', allSub_freqs_amp)]), ...
            saveName, 'WriteRowNames', true, 'QuoteStrings', true) ;

        allSubs{currFile, 1} =  pac_file_mat ;

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
if ~isfolder([srcDir filesep 'generatePAC' filesep 'allSubs'])
    mkdir([srcDir filesep 'generatePAC' filesep 'allSubs']) ;
end
fprintf('Saving information from all files...\n') ;
writetable(array2table(vertcat(allSubs{:}), ...
        'VariableNames', ['filename' 'epoch' 'phase_electrode' 'freq_phase' ...
            'amp_electrode' compose('PAC_Amp_%gHz', allSub_freqs_amp) ...
            compose('nbins_Amp_%gHz', allSub_freqs_amp)]), ...
    helpName([srcDir filesep 'generatePAC' filesep 'allSubs' filesep ...
    'allSubs_pac.csv'], '.csv'), ...
    'WriteRowNames', true, 'QuoteStrings', true) ;

%% SUPPORT FUNCTIONS
% DETERMINE IF CALCULATING FOR INDIVIDUAL, AVERAGE, OR BOTH TRIALS
function [method, csvFormat] = determ_WithinBetween()
    fprintf(['Calculate PAC:\n  within = within each channel\n' ...
        '  between = between each channel\n']) ;
    while true
        ui = input('> ','s') ;
        if strcmpi(ui, 'within'); method = [1,0]; break;
        elseif strcmpi(ui, 'between'); method = [0,1]; break;
        else; fprintf(['Invalid input: please enter "within" or "between"' ...
                ' (without quotations)\n']) ;
        end
    end
    % if method(1)
    %     fprintf(['Print the individual trials/segments as seperate .csv ' ...
    %         'files or as\nsheets in a single Excel file?\n  csv = Output ' ...
    %         'seperate .csv files\n  sheets = Output a single Excel file ' ...
    %         'with multiple sheets\n']) ;
    %     while true
    %         ui = input('> ','s') ;
    %         if strcmpi(ui, 'csv'); csvFormat = 1; break;
    %         elseif strcmpi(ui, 'sheets'); csvFormat = 0; break;
    %         else; fprintf(['Invalid input: please enter "csv" or "sheets"' ...
    %                 ' (without quotations).']) ;
    %         end
    %     end
    % else; 
    csvFormat = 1;
    % end
end

% LIST OUT PARAMETER SET
function genPAC_listParams(params)
    formatOffOnArray = {'Off\n' 'On\n'} ;
    if size(params.pac.time_limits, 2) <= 1
        timeLimitStr = 'Default' ;
    else
        timeLimitStr = [num2str(params.pac.time_limits(1)) ... 
        ' - ' num2str(params.pac.time_limits(2))] ;
    end

    % fprintf('Calculate PAC for Individual Trials/Segments: ') ;
    % fprintf(formatOffOnArray{params.method(1)+1})
    % if params.method(1); fprintf(' - Output Format: ') ;
    %     if params.csvFormat; fprintf('Seperate .csv files\n') ;
    %     else; fprintf('Sheets in an Excel file\n') ;
    %     end
    % end
    
    % fprintf('Calculate PAC for Average of Trials/Segments: ') ;
    % fprintf(formatOffOnArray{params.method(2)+1})
    
    fprintf('Channels of Interest:\n') ;
    fprintf(' - ') ;
    if strcmpi(params.chansAll, 'all')
        fprintf('all channels\n') ;
    elseif strcmpi(params.chansAll, 'coi_include')
        fprintf([sprintf('%s, ', ...
            params.chanIDs{1:end-1}) params.chanIDs{end} '\n']) ;
    elseif strcmpi(params.chansAll, 'coi_exclude')
        fprintf(['all except ' sprintf('%s, ', ...
            params.chanIDs{1:end-1}) params.chanIDs{end} '\n']) ;
    end

    fprintf('Calculate PAC: On\n') ;
    fprintf([' - Phase Range: ' ...
          num2str(params.pac.phase_range(1)) ' Hz - ' num2str(params.pac.phase_range(2)) ' Hz\n' ...
        ' - Amplitude Range: ' num2str(params.pac.amplitude_range(1)) ' - ' ...
          num2str(params.pac.amplitude_range(2)) '\n' ...
        ' - Number of Phase Bins: ' num2str(params.pac.phase_bins) '\n' ...
        ' - Number of Amplitude Bins: ' num2str(params.pac.amplitude_bins) '\n' ...
        ' - Scale Mode: ' params.pac.scale_mode '\n' ...
        ' - PAC method: ' params.pac.method '\n' ...
        ' - Number of bootstraps: ' num2str(params.pac.n_boots) '\n' ...
        ' - Alpha: ' num2str(params.pac.alpha) '\n' ...
        ' - Bonforoni Correction: ' formatOffOnArray{params.pac.bonfcorr+1} ...
        ' - Time Limits: ' timeLimitStr '\n']) ;
end