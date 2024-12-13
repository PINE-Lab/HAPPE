%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% HAPPE Version 3.3
%
% Developed at Northeastern University's PINE Lab
%
% For a detailed description of the pipeline and user options, please see 
% the following manuscripts:
%   Gabard-Durnam, et al., 2018 - HAPPE v1.
%   Lopez, et al. (2022) - HAPPILEE
%   Monachino, et al. (2022) - HAPPE+ER
%   ...
%
% Contributors to HAPPE:
%   Laurel Joy Gabard-Durnam (laurel.gabarddurnam@gmail.com)
%   Adriana S. Mendez Leal (asmendezleal@gmail.com)
%   Carol L. Wilkinson (carol.wilkinson@childrens.harvard.edu)
%   April R. Levin (april.levin@childrens.harvard.edu)
%   ----
%   Alexa D. Monachino (alexamonachino@gmail.com)
%   Kelsie L. Lopez (kelsie_lopez@alumni.brown.edu)
%
% HAPPE includes code that is dependent on the following third-party 
% software. Please reference this third-party software in any manuscripts 
% making use of HAPPE as below:
%   EEGLAB: A Delorme & S Makeig (2004) EEGLAB: an open source toolbox for
%       analysis of single-trial EEG dynamics. Journal of Neuroscience 
%       Methods 134:9-21
%   Cleanline by Tim Mullen (as an EEGLAB plug-in): Mullen, T. (2012). 
%       NITRC: CleanLine: Tool/Resource Info. Available online at: 
%       http://www.nitrc.org/projects/cleanline.
%   ECG Reduction (if using ECGone): Isler JR, Pini N, Lucchini M, 
%       Shuffrey LC, Mitsuyama M, Welch MG, Fifer WP, Stark RI, Myers MM. 
%       A Novel Method for ECG Artifact Removal from EEG without 
%       Simultaneous ECG. Annu Int Conf IEEE Eng Med Biol Soc. 2022 
%       %Jul;2022:3582-3585. doi: 10.1109/EMBC48229.2022.9871252. 
%       PMID: 36086135.
%   FASTER segment-level channel interpolation code (if interpolating 
%       within segments): Nolan, H., Whelan, R., & Reilly, R.B. (2010). 
%       FASTER: Fully Automated Statistical Thresholding for EEG artifact
%       Rejection. Journal of Neuroscience Methods, 192, 152-162.
%   ERPLAB (for Butterworth filter only): Lopez-Calderon J, Luck SJ. 
%       ERPLAB: an open-source toolbox for the analysis of event-related 
%       potentials. Front Hum Neurosci. 2014 Apr 14;8:213.
%   ICLabel EEGLAB plugin (for muscIL option only) - Pion-Tonachini et al. 
%       (2019). ICLabel: An automated electroencephalographic independent
%       component classifier, dataset, and website. Neuroimage, 198, 181–197.
%   REST EEGLAB plugin (for REST re-reference only): Li Dong*, Fali Li, 
%       Qiang Liu, Xin Wen, Yongxiu Lai, Peng Xu and Dezhong Yao*. MATLAB 
%       Toolboxes for Reference Electrode Standardization Technique (REST) 
%       of Scalp EEG. Frontiers in Neuroscience, 2017:11(601).
%
% HAPPE includes code adapted from third-party software/scripts. 
% See the relevant function for further documentation:
%   butterFilt.m adapted from ERPLAB (Lopez-Calderon & Luck, 2014)
%   adaptedREST.m adapted from REST (Yao, 2001), (Dong et al., 2019)
%   ECGone.m (Isler et al., 2022)
%
% Any code that is not part of the third-party dependencies is released
% under the GNU General Public License version 3.
%
% Copyright 2018-2023 Alexa Monachino, Kelsie Lopez, Laurel Gabard-Durnam
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License (version 3) as
% published by the Free Software Foundation.
%
% This software is being distributed with the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See GNU General
% Public License for more details.
%
% In no event shall Boston Children’s Hospital (BCH), the BCH Division of 
% Developmental Medicine, the Laboratories of Cognitive Neuroscience (LCN),
% Northeastern University, the Center for Cognitive and Brain Health (CCBH),
% the Department of Psychology at Northeastern University, the Plasticity 
% in Neurodevelopment (PINE) Lab, or software contributors be liable to any
% party for direct, indirect, special, incidental, or consequential damages,
% including lost profits, arising out of the use of this software and its 
% documentation, even if any of the above listed entities and/or contributors
% have been advised of the possibility of such damage. Software and 
% documentation is provided “as is.” The listed entities and/or software 
% contributors are under no obligation to provide maintenance, support, 
% updates, enhancements, or modifications.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ATTENTION: WE STRONGLY ADVISE AGAINST EDITING THE CODE UNLESS YOU HAVE 
% ADVANCED KNOWLEDGE OF BOTH MATLAB CODING AND SIGNAL PROCESSING.
%
% If you do not have knowledge of both MATLAB coding and signal processing,
% instead interface with this script by running it and entering your 
% parameters in the command window as detailed in the HAPPE User Guide.
%
% If you edit the code we CANNOT GUARANTEE the code will function as
% intended and may negatively affect HAPPE's performance.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% CLEAR THE WORKSPACE
clear ;

%% SET PATH TO INCLUDE NECESSARY FOLDERS
fprintf('Preparing HAPPE...\n') ;
% SET HAPPE AND EEGLAB PATHS USING THE RUNNING SCRIPT
happeDir = strrep(fileparts(which(mfilename('fullpath'))), [filesep '1. ' ...
    'pre-process'], '') ;
eeglabDir = [happeDir filesep 'packages' filesep 'eeglab2024.0'] ;

% ADD HAPPE AND REQUIRED FOLDERS TO PATH
addpath([happeDir filesep '1. pre-process'], ...
    [happeDir filesep 'files'], ...
    [happeDir filesep 'files' filesep 'acquisition_layout_information'], ...
    [happeDir filesep 'scripts'], ...
    [happeDir filesep 'scripts' filesep 'pipeline_steps'], ...
    [happeDir filesep 'scripts' filesep 'ui'], ...
    [happeDir filesep 'scripts' filesep 'support'], ...
    [happeDir filesep 'packages'], ...
    eeglabDir, genpath([eeglabDir filesep 'functions'])) ;
rmpath(genpath([eeglabDir filesep 'functions' filesep 'octavefunc'])) ;

% ADD EEGLAB PLUGIN FOLDERS TO PATH
pluginDir = dir([eeglabDir filesep 'plugins']) ;
pluginDir = strcat(eeglabDir, filesep, 'plugins', filesep, {pluginDir.name}, ';') ;
addpath([pluginDir{:}]) ;

% ADD CLEANLINE FOLDERS TO PATH
if exist('cleanline', 'file')
    cleanlineDir = which('eegplugin_cleanline.m') ;
    cleanlineDir = cleanlineDir(1:strfind(cleanlineDir, 'eegplugin_cleanline.m')-1) ;
    addpath(genpath(cleanlineDir)) ;
else; error('Please make sure cleanline is on your path') ;
end

%% ENSURE TOOLBOXES ARE INSTALLED
toolboxes = matlab.addons.installedAddons ;
if ~ismember('Wavelet Toolbox', toolboxes.Name)
    error('ERROR: WAVELET TOOLBOX NOT INSTALLED. INSTALL TOOLBOX TO RUN HAPPE!') ;
end
if ~ismember('Signal Processing Toolbox', toolboxes.Name)
    error(['ERROR: SIGNAL PROCESSING TOOLBOX NOT INSTALLED. INSTALL ' ...
        'TOOLBOX TO RUN HAPPE!']) ;
end
if ~ismember('Statistics and Machine Learning Toolbox', toolboxes.Name)
    error(['ERROR: STATISTICS AND MACHINE LEARNING TOOLBOX NOT INSTALLED. ' ...
        'INSTALL TOOLBOX TO RUN HAPPE!']) ;
end
if ~ismember('Optimization Toolbox', toolboxes.Name)
    error('ERROR: OPTIMIZATION TOOLBOX NOT INSTALLED. INSTALL TOOLBOX TO RUN HAPPE!') ;
end
clear('toolboxes') ;

%% DETERMINE AND SET PATH TO DATA
% Use input from the Command Window to set the path to the data. If an 
% invalid path is entered, repeat until a valid path is entered.
while true
    srcDir = input('Enter the path to the folder containing the dataset(s):\n> ','s') ;
    if exist(srcDir, 'dir') == 7; break ;
    else; fprintf(['Invalid input: please enter the complete path to the ' ...
            'folder containing the dataset(s).\n']) ;
    end
end
cd (srcDir) ;

%% DETERMINE IF REPROCESSING DATA
[reprocess, ranMuscIL, rerunExt] = isReprocessed() ;

%% DETERMINE IF USING PRESET PARAMETERS
ver = '3_3_0' ;
[preExist, params, changeParams] = isPreExist(reprocess, ver) ;

%% SET PARAMETERS
params = setParams(params, preExist, reprocess, changeParams, happeDir) ;

%% SAVE INPUT PARAMETERS
% If the parameter set was newly created or if a loaded parameter set was
% changed, save the new parameter set as a .mat file. Prompt to use the
% default or a custom name for this saved file. If the file already exists,
% ask to create a new file with a different name or overwrite the existing
% file.
if ~preExist || changeParams
    % CREATE "input_parameters" FOLDER AND ADD IT TO PATH, unless it
    % already exists.
    if ~isfolder([srcDir filesep 'input_parameters'])
        mkdir([srcDir filesep 'input_parameters']) ;
    end
    addpath('input_parameters') ;
    cd input_parameters ;
    
    % DETERMINE PARAMETER FILE NAME: Prompt to use a default or custom name
    % for parameter file. If file exists, ask to create new file with a 
    % different name or overwrite existing file.
    fprintf(['Parameter file save name:\n  default = Default name (input' ...
        'Parameters_dd-mm-yyyy.mat).\n  custom = Create a custom file name' ...
        '\n']) ;
    if choose2('custom', 'default')
        paramFile = paramFile_validateExist(['inputParameters_' ...
            datestr(now, 'dd-mm-yyyy') '.mat'], 'inputParameters_', 2) ;
    else
        fprintf('File name (Do not include .mat):\n') ;
        paramFile = paramFile_validateExist([input('> ', 's') '.mat'], ...
            'inputParameters_', 0) ;
    end

    % SAVE PARAMETERS: Save the params variable to a .mat file using the
    % name created above.
    params.HAPPEver = ver ;
    fprintf('Saving parameters...') ; save(paramFile, 'params') ;
    fprintf('Parameters saved.') ;
end

%% CREATE OUTPUT FOLDERS
% Create the folders in which to store intermediate and final outputs.
% If you aren't running ERP analyses, don't add the ERP_filtered folder.
cd(srcDir) ;
fprintf('Creating output folders...\n') ;
allDirNames = {'intermediate_processing', 'wavelet_cleaned_continuous', ...
    'muscIL', 'ERP_filtered', 'segmenting', 'processed', ...
    'quality_assessment_outputs'} ;
if ~params.paradigm.ERP.on; allDirNames(ismember(allDirNames, 'ERP_filtered')) = []; end
if ~params.muscIL; allDirNames(ismember(allDirNames, 'muscIL')) = []; end
dirNames = cell(1,size(allDirNames,2)) ;

% If reprocessing with a set of steps that has a different number of
% folders than the original run, adapt the folder structure to reflect the
% new set of folders. For example, if reprocessing with muscIL, add a
% 'muscIL' folder.
if reprocess
    currDirNames = dir(srcDir) ;
    currDirNames = {currDirNames([dir(srcDir).isdir]).name} ;
    if any(contains(allDirNames, 'muscIL')) && ~any(contains(currDirNames, 'muscIL'))
        for i=1:length(allDirNames)
            dirNames{i} = [num2str(i) ' - ' allDirNames{i}] ;
            if ~isfolder(dirNames{i}) && ~any(contains(currDirNames, allDirNames{i}))
                mkdir([srcDir filesep dirNames{i}]) ;
            elseif ~isfolder(dirNames{i}) && any(contains(currDirNames, allDirNames{i}))
                movefile(currDirNames{contains(currDirNames, allDirNames{i})}, dirNames{i}) ;
            end
        end
    elseif ~any(contains(allDirNames, 'muscIL')) && any(contains(currDirNames, 'muscIL'))
        wavInd = find(contains(setdiff(currDirNames, {'.', '..'}),'wavelet')) ;
        allDirNames = [allDirNames(1:wavInd) 'muscIL' allDirNames(wavInd+1:end)] ;
        for i=1:length(allDirNames)
            dirNames{i} = [num2str(i) ' - ' allDirNames{i}] ;
        end
        clear('wavInd') ;
    else
        for i=1:length(allDirNames)
            dirNames{i} = [num2str(i) ' - ' allDirNames{i}] ;
        end
    end
    clear('currDirNames') ;
else
    for i=1:length(allDirNames)
        dirNames{i} = [num2str(i) ' - ' allDirNames{i}] ;
        if ~isfolder([srcDir filesep dirNames{i}])
            mkdir([srcDir filesep dirNames{i}]) ;
        end
    end
end
clear('allDirNames') ;
fprintf('Output folders created.\n') ;

%% COLLECT DATA TO RUN
fprintf('Gathering files...\n') ;
cd(srcDir) ;
% DETERMINE FILE EXTENSION USING USER LOAD INFO
switch params.loadInfo.inputFormat
    case 1; inputExt = '.mat' ;
    case 2; inputExt = '.raw' ;
    case 3; inputExt = '.set' ;
    case 4; inputExt = '.cdt' ;
    case 5; inputExt = '.mff' ;
    case 6; inputExt = '.EDF' ;
    case 7; inputExt = '.set' ;
end
% COLLECT FILE NAMES: gather the names of all the files with the relevant
% file extension in the directory indicated by the user. If no files are
% found, end the run via error.
FileNames = {dir(['*' inputExt]).name} ;
if isempty(FileNames); error(['ERROR: No ' inputExt ' files detected!']) ; end
% LOCATE STIM FILE AND NAMES
if params.loadInfo.inputFormat == 1 && params.paradigm.task
    % ***
end

%% INITIALIZE QUALITY REPORT METRICS
fprintf('Initializing report metrics...\n') ;
% DATA QUALITY METRICS: create a variable holding all the names of each
% metric and a variable to hold the metrics for each file.
dataQCnames = {'File_Length_in_Seconds', 'Number_User-Selected_Chans', ...
    'Number_Good_Chans_Selected', 'Percent_Good_Chans_Selected', 'Bad_Chan_IDs', ...
    'Percent_Var_Retained_Post-Wav', 'Number_ICs_Rej', 'Percent_ICs_Rej', ...
    'Chans_Interpolated_per_Seg', 'Number_Segs_Pre-Seg_Rej', ...
    'Number_Segs_Post-Seg_Rej', 'Percent_Segs_Post-Seg_Rej', 'Kept_Segs_Indxs'} ;
dataQC = cell(length(FileNames), length(dataQCnames)) ;
% If processing for tasks, create an additional variable to hold specific
% data metrics for each onset tag.
if params.paradigm.task
    dataQC_task = cell(length(FileNames), length(params.paradigm.onsetTags)*5) ;
    dataQCnames_task = cell(1, length(params.paradigm.onsetTags)*5) ;
    for i=1:size(params.paradigm.onsetTags,2)
        dataQCnames_task{i*5-4} = ['Number_' params.paradigm.onsetTags{i} ...
            '_Segs_Pre-Seg_Rej'] ;
        dataQCnames_task{i*5-3} = ['Number_' params.paradigm.onsetTags{i} ...
            '_Segs_Post-Seg_Rej'] ;
        dataQCnames_task{i*5-2} = ['Percent_' params.paradigm.onsetTags{i} ...
            '_Segs_Post-Seg_Rej'] ;
        dataQCnames_task{i*5-1} = ['Kept_' params.paradigm.onsetTags{i} ...
            '_Segs_Indxs_allSegs'] ;
        dataQCnames_task{i*5} = ['Kept_' params.paradigm.onsetTags{i} ...
            '_Segs_Indxs_byTag'] ;
    end
    
    % If grouping any tags by condition, create additional variable to hold
    % specific data metrics for each condition.
    if params.paradigm.conds.on
        dataQC_conds = cell(length(FileNames), ...
            size(params.paradigm.conds.groups,1)*3) ;
        dataQCnames_conds = cell(1, size(params.paradigm.conds.groups,1)*3) ;
        for i = 1:size(params.paradigm.conds.groups,1)
            dataQCnames_conds{i*3-2} = ['Number_' params.paradigm.conds.groups{i, ...
                1} '_Segs_Pre-Seg_Rej'] ;
            dataQCnames_conds{i*3-1} = ['Number_' params.paradigm.conds.groups{i, ...
                1} '_Segs_Post-Seg_Rej'] ;
            dataQCnames_conds{i*3} = ['Percent_' params.paradigm.conds.groups{i, ...
                1} '_Segs_Post-Seg_Rej'] ;
        end
    end
end
errorLog = {} ;

% PIPELINE QUALITY METRICS: initialize the variables to hold the pipeline
% quality metrics for assessing HAPPE performance on line noise reduction
% and waveleting.
if ~reprocess; lnMeans = [] ; wavMeans = [] ; end
if params.muscIL; icaMeans = [] ; end

%% RUN THE PROCESSING PIPELINE OVER EACH FILE
for currFile = 1:length(FileNames)
    cd(srcDir) ;
    try
        if ~reprocess
            %% LOAD AND VALIDATE RAW FILE
            fprintf(['Loading ' FileNames{currFile} '...\n']) ;
            try
                switch params.loadInfo.inputFormat
                    % .MAT
                    case 1
                        % NetStation Format:
                        if params.loadInfo.NSformat
                            load(FileNames{currFile}) ;
                            % Use known information about Net Station sampling 
                            % rate variable names to determine the sampling 
                            % rate of the file.
                            srate = intersect(who, {'samplingRate', 'EEGSamplingRate'}) ;
                            srate = double(eval(srate{1})) ;
                            % Use the potential variable names given by the user to
                            % determine the data variable.
                            eegVName = intersect(who, params.loadInfo.NSvarNames) ;
                            % Try to load the EEG using an EGI Net Station EEGLAB
                            % plugin. If unable to load the file, throws an error.
                            % Specifically looks for MATLAB:badsubscript.
                            try EEGraw = pop_importegimat(FileNames{currFile}, ...
                                    srate, 0.00, eegVName{1}) ;
                            catch ME
                                if strcmp(ME.identifier, 'MATLAB:badsubscript')
                                    fprintf(['Sorry, could not read the variable' ...
                                        ' name of the EEG data. Please check your' ...
                                        ' file.\n']) ;
                                end
                                rethrow(ME) ;
                            end
                        % Matrix:
                        else
                            % If all the files have the same sampling rate, use
                            % the single user specified value. Otherwise, load 
                            % the table of sampling rates specified by the user
                            % and locate the corresponding value from the table
                            % using the file name.
                            if params.loadInfo.srate.same
                                srate = params.loadInfo.srate.val ;
                            else
                                srateTable = table2cell(readtable(params.loadInfo.srate.file)) ;
                                srate = srateTable{contains(list, ...
                                    FileNames{currFile}), 2} ;
                            end
                            % Use the built-in EEGLAB function to load the
                            % data. Uses the sampling rate determined above and
                            % the user-specified channel locations. If the user
                            % did not specify any channel locations, is still
                            % able to load successfully.
                            EEGraw = pop_importdata('dataformat', 'matlab', ...
                                'data', FileNames{currFile}, 'srate', srate, ...
                                'chanlocs', params.loadInfo.chanlocs.file) ;
                        end
                        % If task data, import event information using the file
                        % indicated by the user.
                        if params.paradigm.task
                            cd(params.loadInfo.eventLoc) ;
                            EEGraw = pop_importevent(EEGloaded, 'append', 'no', ...
                                'event', StimNames{currFile}, 'fields', {'type', ...
                                'latency', 'status'}, 'skipline', 1, 'timeunit', ...
                                1E-3, 'align', NaN) ;
                            origEvents = EEGraw.urevent ;
                            events = EEGraw.event ;
                        end
                    % .RAW FILES
                    case 2; EEGraw = pop_readegi(FileNames{currFile}, [], ...
                            [], 'auto') ;
                    % .SET FILES
                    case 3
                        EEGraw = load('-mat', FileNames{currFile}) ;
                        if isfield(EEGraw, 'EEG'); EEGraw = EEGraw.EEG; end
                    % .CDT FILES
                    case 4; EEGraw = loadcurry([srcDir filesep ...
                            FileNames{currFile}], 'CurryLocation', 'false') ;
                    % .MFF FILES
                    case 5; EEGraw = pop_mffimport(FileNames{currFile} , ...
                            params.loadInfo.typeFields, 0, 0) ;
                    % .EDF FILES
                    case 6; EEGraw = pop_biosig(FileNames{currFile}, ...
                            'blockepoch', 'off') ;
                    % .BDF->.SET
                    case 7; EEGraw = pop_loadset('filename', FileNames{currFile}) ;
                end

                % Determine the sampling rate from the loaded EEG. If task
                % related, save the events from the loaded file in its own
                % variable.
                srate = double(EEGraw.srate) ;
                if params.paradigm.task; events = EEGraw.event ; end

                % Validate the loaded file using EEGLAB's checkset function
                fprintf('Validating file...\n') ;
                EEGraw.setname = 'rawEEG' ;
                EEG = eeg_checkset(EEGraw) ;
                dataQC{currFile, 1} = EEG.xmax ;

            % If unable to load the file, indicate this to the user via the
            % command line and throw an error.
            catch ME
                error('HAPPE:LoadFail', ...
                    ['ERROR: Unable to load ' FileNames{currFile} '\n']) ;
            end
            
            %% CORRECT CHANNELS, IF NEEDED
            % Channels need to be corrected for EGI files, regardless of
            % input format. Please note that for .raw files, if exported
            % incorrectly from netstation, this will not work!
            if params.loadInfo.correctEGI
                fprintf('Correcting EGI Channel Locations...\n') ;
                switch params.loadInfo.layout(1)
                    % LOAD THE CHANNEL LOCATIONS FILE USING THE LAYOUT
                    % INFORMATION PARAMETERS
                    case 1
                        if params.loadInfo.layout(2) == 64
                        params.loadInfo.chanlocs.file = [happeDir filesep ...
                            'files' filesep 'acquisition_layout_information' ...
                            filesep 'GSN65v2_0.sfp'] ;
                        else; error('Invalid sensor layout selection.') ;
                        end
                    case 2
                        if ismember(params.loadInfo.layout(2), [32, 64, 128, 256])
                            params.loadInfo.chanlocs.file = [happeDir filesep ...
                                'files' filesep 'acquisition_layout_information' ...
                                filesep 'GSN-HydroCel-' ...
                                num2str(params.loadInfo.layout(2)+1) '.sfp'] ;
                        else; error('Invalid sensor layout selection.') ;     
                        end
                end
                % CORRECT THE CHANNELS USING EEGLAB FUNCTIONALITY
                EEG = pop_chanedit(EEG, EEG.chanlocs, 'load', ...
                    {params.loadInfo.chanlocs.file 'filetype' 'autodetect'}) ;
                % VALIDATE THAT THE EEG WAS CORRECTED CORRECTLY
                EEG = eeg_checkset(EEG) ;
            end
            
            %% SET THE FLATLINE REFERENCE CHANNEL
            % If the reference channel is included in the data, remove it
            % from the data channels. This will prevent it from being 
            % falsely flagged as a bad channel.
            if params.loadInfo.chanlocs.inc && params.reref.on && ...
                    params.reref.flat
                if ismember(params.reref.chan, {EEG.chanlocs.labels})
                    indx = find(strcmpi({EEG.chanlocs.labels}, ...
                        params.reref.chan)) ;
                    EEG.data(indx,:) = [] ;
                    EEG.nbchan = EEG.nbchan - 1;
                    EEG.chaninfo.refChan = EEG.chanlocs(indx) ;
                    EEG.chanlocs(indx) = [] ;
                elseif ~isempty(EEG.chaninfo.nodatchans) && ...
                        ismember(params.reref.chan, {EEG.chaninfo.nodatchans.labels})
                    indx = find(strcmpi({EEG.chaninfo.nodatchans.labels}, ...
                        params.reref.chan)) ;
                    EEG.chaninfo.refChan = EEG.chaninfo.nodatchans(indx) ;
                    EEG.chaninfo.refChan = rmfield(EEG.chaninfo.refChan, ...
                        'datachan') ;
                    warning(['Located reference channel in nodatchans.' ...
                        ' Channel locations for this channel may ' ...
                        'not be correct!']) ;
                else
                    warning(['Unable to locate your defined reference ' ...
                        'channel! Will run as if no flatline reference ' ...
                        'is included.']) ;
                    params.reref.flat = 0 ;
                end 
            end

            %% UPDATE CHANNEL IDS AND FILTER TO THE CHANNELS OF INTEREST
            % IF CHANNEL LOCATIONS ARE INCLUDED, update the channel IDs in
            % the EEG structure to match the channels of interest.
            if params.loadInfo.chanlocs.inc
                % Updates chanIDs (a variable holding the channels of
                % interest for the current subject) using information
                % provided by the user. If for some reason one or more
                % channels are missing in a file, this will protect against
                % the channel IDs being impacted in other files.
                fprintf('Selecting channels of interest...\n') ;
                if strcmpi(params.chans.subset, 'coi_exclude')
                    chanIDs = setdiff({EEG.chanlocs.labels}, params.chans.IDs) ;
                elseif strcmpi(params.chans.subset, 'coi_include')
                    chanIDs = intersect({EEG.chanlocs.labels}, params.chans.IDs) ;
                elseif strcmpi(params.chans.subset, 'all')
                    chanIDs = {EEG.chanlocs.labels} ;
                end
                % Filter to the channels of interest
                fprintf('Filtering to channels of interest...\n') ;
                EEG = pop_select(EEG, 'channel', chanIDs) ;
                EEG.setname = [EEG.setname '_cs'] ;
                % Report how many channels are present in the dataset as a
                % data metric.
                dataQC{currFile,2} = size(chanIDs,2) ;
            
            % IF NO CHANNEL LOCATIONS ARE INCLUDED, report how many
            % channels are present in the dataset as a data metric.
            else
                dataQC{currFile,2} = size(EEG.data,1) ;
                chanIDs = 1:size(EEG.data,1) ;
            end

            % Enter the intermediate_processing folder
            cd([srcDir filesep dirNames{contains(dirNames, ...
                'intermediate_processing')}]) ;

            %% REDUCE LINE NOISE
            % Attempt to reduce line noise using the happe_reduceLN
            % function (see function for further documentation). If HAPPE
            % fails during this step, alert the user that processing failed
            % during the line noise step in the command window and rethrow
            % the error to proceed to the next file.
            try [EEG, lnMeans] = happe_reduceLN(EEG, params.lineNoise, ...
                    lnMeans) ;
                EEG.setname = [EEG.setname '_ln'] ;
            catch ME
                fprintf(2, 'ERROR: Line noise reduction failed.\n') ;
                rethrow(ME) ;
            end

            pop_saveset(EEG, 'filename', strrep(FileNames{currFile}, ...
                inputExt, '_lnreduced.set')) ;
            
            %% RESAMPLE DATA
            % If the user enabled downsampling, resample the data to the
            % user-specified frequency using EEGLAB functionality.
            if params.downsample
                fprintf(['Resampling the data to ' num2str(params.downsample) ...
                    ' Hz...\n']) ;
                EEG = pop_resample(EEG, params.downsample) ;
            end

            %% FILTER AND SAVE INTERMEDIATE FILE
            % If not an ERP paradigm, filter the data using user-specified
            % high- and low-passes. Save the filtered dataset to the
            % intermediate_processing folder.
            if ~params.paradigm.ERP.on && params.filt.on
                EEG = pop_eegfiltnew(EEG, params.filt.highpass, ...
                    params.filt.lowpass, [], 0, [], 0) ;
                EEG.setname = [EEG.setname '_filt'] ;
                pop_saveset(EEG, 'filename', strrep(FileNames{currFile}, ...
                    inputExt, '_filtered.set')) ;
            end
            
            %% DETECT BAD CHANNELS
            % If bad channel detection is on, detect and reject bad
            % channels using the happe_detectBadChans function (see
            % function for further documentation). Regardless if detecting
            % bad channels, update the number of channels included in
            % processing (good channels), the percent of good channels
            % selected (100% if bad channel detection is disabled), and the
            % list of channels, if any, marked as bad.
            origChans = EEG.chanlocs ;
            if params.badChans.rej
                [EEG, dataQC] = happe_detectBadChans(EEG, params, dataQC, ...
                    chanIDs, currFile) ;
                EEG.setname = [EEG.setname 'badc'] ;
                pop_saveset(EEG, 'filename', strrep(FileNames{currFile}, ...
                    inputExt, '_badchanrej.set')) ;
            else
                dataQC{currFile,3} = size(chanIDs,2) ;
                dataQC{currFile,4} = dataQC{currFile,3}/dataQC{currFile,2}*100;
                dataQC{currFile,5} = 'NA' ;
            end

            %% ECGONE
            % If ECGone is enabled, remove significant ECG artifact through
            % the ECGone function (see function for further documentation),
            % as adapted from the method created by Isler et al., 2022.
            try EEG = ECGone(EEG, params.ecgone) ;
            catch ME
                fprintf(2, ['ERROR: ECGone failed. Proceeding to wavelet ' ...
                    'thresholding...\n']) ;
%                 rethrow(ME) ;
            end
            
            %% WAVELET THRESHOLDING
            % Attempt to wavelet threshold using the happe_wavThresh
            % function (see function for further documentation). If HAPPE
            % fails during this step, alert the user that processing failed
            % during wavelet thresholding in the command window and rethrow
            % the error to proceed to the next file.
            cd([srcDir filesep dirNames{contains(dirNames, ...
                'wavelet_cleaned_continuous')}]) ;
            pop_saveset(EEG, 'filename', strrep(FileNames{currFile}, ...
                inputExt, '_prewav.set')) ;

            try [EEG, wavMeans, dataQC] = happe_wavThresh(EEG, params, ...
                    wavMeans, dataQC, currFile) ;
            catch ME
                fprintf(2, 'ERROR: Wavelet thresholding failed.\n') ;
                rethrow(ME) ;
            end
            
            % Save the wavelet-thresholded EEG as an intermediate output.
            pop_saveset(EEG, 'filename', strrep(FileNames{currFile}, ...
                inputExt, '_wavclean.set')) ;
            
        else
            %% LOAD AND VALIDATE PROCESSED FILE
            try
                fprintf(['Loading ' FileNames{currFile} '...\n']) ;
                % LOAD FILTERED & LINE NOISE REDUCED FILE:
                % Use the raw file name. Collect the original number of
                % user-selected channels from this file, then clear the
                % dataset.
                EEGloaded = pop_loadset('filename', strrep(FileNames{currFile}, ...
                    inputExt, '_lnreduced.set'), 'filepath', [srcDir filesep ...
                    dirNames{contains(dirNames, 'intermediate_processing')}]) ;
                dataQC{currFile, 2} = size(EEGloaded.data,1) ;
                origChans = EEGloaded.chanlocs ;
                clear('EEGloaded') ;
 
                % LOAD PRE-SEGMENTED PROCESSED DATA:
                % Load the wavelet-thresholded data. Using the raw file to 
                % locate and load the correct file.
                if ranMuscIL; EEGloaded = pop_loadset('filename', ...
                        strrep(FileNames{currFile}, inputExt, ...
                        '_muscILcleaned.set'), 'filepath', [srcDir filesep ...
                        dirNames{contains(dirNames, 'muscIL')}]) ;
                else; EEGloaded = pop_loadset('filename', ...
                        strrep(FileNames{currFile}, inputExt, ...
                        '_wavclean.set'), 'filepath', [srcDir filesep ...
                        dirNames{contains(dirNames, ...
                        'wavelet_cleaned_continuous')}]) ;
                end
                        
                % VALIDATE LOADED FILE
                EEG = eeg_checkset(EEGloaded) ;             
                
                % FILL IN DATA QUALITY METRICS:
                % Using the loaded data, fill in the following data quality
                % metrics about the current file: length of the EEG
                % recording, the number of good channels, and the percent 
                % of good channels.
                dataQC{currFile,1} = EEG.xmax ;
                dataQC{currFile,3} = size(EEG.data,1) ;
                dataQC{currFile,4} = dataQC{currFile,3}/dataQC{currFile,2}*100 ;
                dataQC{currFile,5} = 'NA' ;
                dataQC{currFile,6} = 'NA' ;
                
            % If unable to load the file, indicate this to the user via the
            % command line and throw an error.
            catch ME
                error('HAPPE:LoadFail', ['ERROR: Unable to load ' ...
                    FileNames{currFile} '\n']) ;
            end
        end
        
        %% MUSCIL
        if ~params.paradigm.ERP.on && params.muscIL && ~ranMuscIL
            cd([srcDir filesep dirNames{contains(dirNames, 'muscIL')}]) ;
            try [EEG, dataQC] = happe_muscIL(EEG, dataQC, FileNames, ...
                    currFile, inputExt) ;
            catch ME; rethrow(ME) ;
            end
        else
            dataQC{currFile, 7} = 'NA' ;
            dataQC{currFile, 8} = 'NA' ;
        end
        
        %% FILTER USING ERP CUTOFFS AND SAVE
        % If performing preprocessing for ERPs, attempts to filter 
        % using the user-defined ERP cutoffs. If HAPPE fails during 
        % this step, alert the user that processing failed during 
        % filtering to ERP cutoffs in the command window and rethrow
        % the error to proceed to the next file.
        if params.paradigm.ERP.on
            fprintf('Filtering using ERP cutoffs...\n') ;
            cd([srcDir filesep dirNames{contains(dirNames, ...
                'ERP_filtered')}]) ;
            try
                % If the user indicated to use ERPLAB's butterworth
                % filter, filter using the butterFilt function (adapted
                % from ERPLAB's functions - see butterFilt
                % documentation for more information). Currently only
                % available for ERP paradigms.
                if params.filt.butter
                    EEG = butterFilt(EEG, params.filt.lowpass, ...
                        params.filt.highpass) ;
                % If the user indicated to use EEGLAB's FIR filter,
                % filter using the EEGLAB function.
                else
                    EEG = pop_eegfiltnew(EEG, params.filt.highpass, ...
                        params.filt.lowpass, [], 0, [], 0) ;
                end
                EEG.setname = [EEG.setname '_forERP'] ;
            catch ME
                fprintf(2, 'ERROR: Filtering to ERP cutoffs failed.\n') ;
                rethrow(ME) ;
            end

            % Save the filtered EEG as an intermediate output.
            pop_saveset(EEG, 'filename', strrep(FileNames{currFile}, ...
                inputExt, ['_wavcleaned_filteredERP' rerunExt '.set'])) ;
        end
        
        %% RE-ADD FLATLINE RE-REFERENCE
        % If there was a flatline channel, return it to the data. If HAPPE
        % is unable to find the designated channel, throw an error to
        % proceed to the next file and alert the user via the command
        % window.
        if params.loadInfo.chanlocs.inc && params.reref.on && params.reref.flat
            try
                EEG.data(EEG.nbchan+1,:) = 0 ;
                EEG.nbchan = EEG.nbchan + 1 ;
                EEG.chanlocs(EEG.nbchan) = EEG.chaninfo.refChan ;
                origChans(length(origChans)+1) = EEG.chaninfo.refChan ;
            catch ME
                fprintf(2, 'ERROR: No flatline channel detected.') ;
                rethrow(ME) ;
            end
        end
        
        %% SEGMENT DATA
        % If the data is not already segmented, use the user-specified
        % parameters to segment the data. If it is segmented, alert the
        % user via the command window that the data could not be segmented.
        cd([srcDir filesep dirNames{contains(dirNames, 'segmenting')}]) ;
        if params.segment.on
            try EEG = happe_segment(EEG, params) ;

                % SAVE THE SEGMENTED DATASET AS AN INTERMEDIATE FILE
                pop_saveset(EEG, 'filename', strrep(FileNames{currFile}, inputExt, ...
                    ['_segmented' rerunExt '.set'])) ;
            catch ME
                fprintf(2, 'ERROR: Segmenting failed.\n') ;
                rethrow(ME) ;
            end
        end
        
        % Add the number of segments to the dataQC matrix. In the case of
        % ERP paradigms, add the number of segments for each event tag and
        % each condition to the dataQC_task and dataQC_conds matrices.
        dataQC{currFile, 10} = EEG.trials ;
        if params.paradigm.task
            % Segments per Event Tag:
            for i=1:length(params.paradigm.onsetTags)
                try dataQC_task{currFile, i*5-4} = length(pop_selectevent(EEG, ...
                       'type', params.paradigm.onsetTags{i}).epoch) ;
                catch
                    dataQC_task{currFile, i*5-4} = 0 ;
                    fprintf('No instances of "%s" in this file.\n', ...
                        params.paradigm.onsetTags{i}) ;
                end
            end
            
            % Segments per Condition:
            if params.paradigm.conds.on
                for i=1:size(params.paradigm.conds.groups,1)
                    try 
                        currGroup = params.paradigm.conds.groups(i,2:end) ;
                        currGroup = currGroup(~cellfun('isempty',currGroup)) ;
                        dataQC_conds{currFile, i*3-2} = length(pop_selectevent(EEG, ...
                            'type', currGroup).epoch) ;
                    catch
                        dataQC_conds{currFile, i*3-2} = 0 ;
                    end
                end
            end
        end
         
        %% BASELINE CORRECT AND SAVE (ERPs)
        % If pre-processing ERP data and baseline correction is enabled,
        % baseline correct the data using the baseline window specified by
        % the user.
        if params.paradigm.ERP.on && params.baseCorr.on
            fprintf('Correcting baseline...\n') ;
            try EEG = pop_rmbase(EEG, [params.baseCorr.start ...
                    params.baseCorr.end]) ;

                % SAVE BASELINE CORRECTED DATA AS AN INTERMEDIATE FILE
                pop_saveset(EEG, 'filename', strrep(FileNames{currFile}, inputExt, ...
                    ['_segmented_BLcorr' rerunExt '.set'])) ;
            catch ME
                fprintf(2, 'ERROR: Baseline correction failed.\n') ;
                rethrow(ME) ;
            end
        end
        
        
        %% INTERPOLATE BAD DATA USING CHANNELS
        % If enabled, interpolate data within segments by channels. Uses
        % code adapted from the FASTER program (Nolan et al., 2010),
        % interpolating channels scoring above/below z threshold of 3 for a
        % segment. Not available if channels are not included or if bad
        % channel detection was not performed. Adds the segments and 
        % interpolated channels, if any, to the dataquality metrics.
        if params.segment.interp
            try [EEG, dataQC] = happe_interpChanBySeg(EEG, dataQC, ...
                    currFile, params) ;

                % SAVE INTERPOLATED DATA AS INTERMEDIATE FILE
                pop_saveset(EEG, 'filename', strrep(FileNames{currFile}, inputExt, ...
                    ['_segmented_interp' rerunExt '.set'])) ;
            catch ME
                fprintf(2, 'ERROR: Interpolation within segments failed.\n') ;
                rethrow(ME) ;
            end
        else; dataQC{currFile, 9} = 'NA' ;
        end
        
        %% REJECT SEGMENTS
        % If the data is segmented and segment rejection is enabled, reject
        % segments containing artifact according to the criteria as set by
        % the user. Can reject segments based on a region of interest or by
        % the entire set of channels present in the data. If the data is
        % continuous or is a single trial, cannot reject segments.
        if params.segRej.on && EEG.trials > 1
            try 
                if params.paradigm.task && size(params.paradigm.onsetTags,2) > 1
                    origTrials = EEG.epoch ;
                end
                [EEG, keptTrials] = happe_segRej(EEG, params) ;
                
                % SAVE FILE WITH SEGMENTS REJECTED
                pop_saveset(EEG, 'filename', strrep(FileNames{currFile}, ...
                    inputExt, ['_segments_postrej' rerunExt '.set'])) ;
            catch ME; rethrow(ME) ;
            end
            if size(keptTrials,2)==dataQC{currFile,11}
                dataQC{currFile,13} = 'All' ;
            else; dataQC{currFile, 13} = [sprintf('%i, ', keptTrials(1:end-1)), ...
                    num2str(keptTrials(end))] ;
            end
            if params.paradigm.task && size(params.paradigm.onsetTags,2) > 1
                tempTrials = cell(1, size(params.paradigm.onsetTags,2)) ;
                for i=1:size(params.paradigm.onsetTags,2)
                    for j=1:size(origTrials,2)
                        if ismember(params.paradigm.onsetTags{i}, ...
                                origTrials(j).eventtype)
                            tempTrials{i} = [tempTrials{i} j] ;
                        end
                    end
                end

                for i=1:size(params.paradigm.onsetTags,2)
                    for j=1:size(keptTrials,2)
                        if ismember(params.paradigm.onsetTags{i}, ...
                                origTrials(keptTrials(j)).eventtype)
                            if isempty(dataQC_task{currFile, i*5})
                                dataQC_task{currFile, i*5-1} = num2str(keptTrials(j)) ;
                                dataQC_task{currFile, i*5} = num2str(find(tempTrials{i}==keptTrials(j))) ;
                            else
                                dataQC_task{currFile, i*5-1} = [dataQC_task{currFile, ...
                                    i*5-1} ', ' num2str(keptTrials(j))] ;
                                dataQC_task{currFile, i*5} = [dataQC_task{currFile, ...
                                    i*5} ', ' num2str(find(tempTrials{i}==keptTrials(j)))] ;
                            end
                        end
                    end
                end
                
                clear('origTrials') ;
            end
            clear('keptTrials') ;
        elseif params.segRej.on && EEG.trials == 1
            error('HAPPE:rejOneSeg', ['ERROR: Cannot reject segments from' ...
                ' data with a single epoch.']) ;
        else; dataQC{currFile,13} = 'All' ;
        end
        
        % UPDATE DATA QUALITY METRICS: Add the number and percent of trials
        % retained post-segment rejection.
        dataQC{currFile, 11} = EEG.trials ;
        dataQC{currFile, 12} = dataQC{currFile, 11}/dataQC{currFile, 10}*100 ;
 
        %% INTERPOLATE BAD CHANNELS
        % If the file has channel locations, interpolate the bad channels
        % to enable proper re-referencing.
        try
            if params.loadInfo.chanlocs.inc
                fprintf('Interpolating bad channels...\n') ;
                EEG = pop_interp(EEG, origChans, 'spherical') ;
                EEG.setname = [EEG.setname '_interp'] ;
            end
        catch ME
            fprintf(2, 'ERROR: Bad channel interpolation failed.\n') ;
            rethrow(ME) ;
        end
        
        %% RE-REFERENCE DATA
        % If enabled, re-reference the data. Will perform an average
        % rereference or re-reference to a subset based on the parameters
        % set by the user. If a reference electrode is present in the data,
        % will keep it.
        if params.reref.on
            fprintf('Re-referencing...\n') ;
            if strcmpi(params.reref.method, 'average')
                EEG = pop_reref(EEG, [], 'keepref', 'on') ;
            elseif strcmpi(params.reref.method, 'subset')
                EEG = pop_reref(EEG, intersect({EEG.chanlocs.labels}, ...
                    params.reref.subset, 'stable'), 'keepref', 'on') ;
                EEG.setname = [EEG.setname '_reref'] ;
            elseif strcmpi(params.reref.method, 'rest')
                if strcmpi(params.reref.rest.ref, 'mastoid')
                    params.reref.rest.chanL = find(ismember({EEG.chanlocs.labels}, ...
                        params.reref.rest.chanL)) ;
                    params.reref.rest.chanR = find(ismember({EEG.chanlocs.labels}, ...
                        params.reref.rest.chanR)) ;
                else
                    params.reref.rest.chanlist = 1:size({EEG.chanlocs.labels},2) ;
                end
                EEG = adaptedREST(EEG, params.reref.rest) ;
            end
        end
        
        %% SPLIT EEG BY EVENT TAGS AND CONDITIONS
        % In the case of multiple event tags in a dataset, create a unique
        % EEG structure containing every instance of a single event tag
        % for each event tag. If a tag does not appear in the original EEG
        % dataset, alerts the user via the command window. Does the same 
        % per condition. Additionally updates the data quality metrics 
        % specific to tags and conditions with the number and percent 
        % trials retained.
        if params.paradigm.task && size(params.paradigm.onsetTags,2) > 1
            % Split by Event Tags:
            fprintf('Creating EEGs by event tags...\n') ;
            eegByTags = cell(1,size(params.paradigm.onsetTags,2)) ;
            for i=1:size(params.paradigm.onsetTags,2)
                try eegByTags{i} = pop_selectevent(EEG, 'type', ...
                        params.paradigm.onsetTags{i}) ;                   
                    dataQC_task{currFile, i*5-3} = length(eegByTags{i}.epoch) ;
                    dataQC_task{currFile, i*5-2} = dataQC_task{currFile, i*5-3} ...
                        / dataQC_task{currFile, i*5-4}*100 ;
                catch
                    fprintf('No instances of tag %s appear in this file.\n', ...
                        params.paradigm.onsetTags{i}) ;
                    dataQC_task{currFile, i*5-3} = 0 ;
                    dataQC_task{currFile, i*5-2} = 'NA' ;
                end
            end
            
            % Split by Conditions:
            if params.paradigm.conds.on
                eegByConds = cell(1, size(params.paradigm.conds.groups,1)) ;
                for i=1:size(params.paradigm.conds.groups,1)
                    eegByConds{i} = pop_selectevent(EEG, 'type', ...
                        params.paradigm.conds.groups(i,2:end)) ;
                    dataQC_conds{currFile, i*3-1} = length(eegByConds{i}.epoch) ;
                    dataQC_conds{currFile, i*3} = dataQC_conds{currFile, i*3-1} ...
                        / dataQC_conds{currFile, i*3-2}*100 ;
                end
            end
        end
        
        %% GENERATE VISUALIZATIONS
        % If visualizations are enabled, generate topoplots for the current
        % file. If preprocessing for ERPs, this will be the processed ERP
        % spectrum.
        cd([srcDir filesep dirNames{contains(dirNames, 'processed')}]) ;
        if params.vis.enabled
           if params.paradigm.ERP.on
               figure; pop_timtopo(EEG, [params.vis.min params.vis.max], ...
                   params.vis.toPlot) ;
               saveas(gcf, strrep(FileNames{currFile}, inputExt, ...
                   ['_processedERPspectrum' rerunExt '.jpg'])) ;
           else
               figure; pop_spectopo(EEG, 1, [], 'EEG', 'freq', ...
                   [params.vis.toPlot], 'freqrange', [[params.vis.min] ...
                   [params.vis.max]], 'electrodes', 'off') ;
               saveas(gcf, strrep(FileNames{currFile}, inputExt, ...
                   ['_processedSpectrum' rerunExt '.jpg'])) ;
           end
        end
        
        %% SAVE PRE-PROCESSED DATASET(S)
        fprintf('Saving pre-processed dataset(s)...\n') ;
        if params.outputFormat == 1 || params.outputFormat == 3
            pop_saveset(EEG, 'filename', strrep(FileNames{currFile}, ...
                inputExt, ['_processed' rerunExt '.set'])) ;
            if params.paradigm.task && size(params.paradigm.onsetTags,2) > 1
                for i=1:size(eegByTags, 2)
                    if ~isempty(eegByTags{i})
                        pop_saveset(eegByTags{i}, 'filename', ...
                            strrep(FileNames{currFile}, inputExt, ...
                            ['_processed_' params.paradigm.onsetTags{i} ...
                            rerunExt '.set'])) ;
                    end
                end
                if params.paradigm.conds.on
                    for i=1:size(eegByConds,2)
                        if ~isempty(eegByConds{i})
                            pop_saveset(eegByConds{i}, 'filename', ...
                                strrep(FileNames{currFile}, inputExt, ...
                                ['_processed_' params.paradigm.conds.groups{i,1} ...
                                rerunExt '.set'])) ;
                        end
                    end
                end
            end
        end
        if params.outputFormat == 2
            save(strrep(FileNames{currFile}, inputExt, ['_processed' ...
                rerunExt '.mat']), 'EEG') ;
            if params.paradigm.task && size(params.paradigm.onsetTags, ...
                    2) > 1
                for i=1:size(eegByTags, 2)
                    currEEG = eegByTags{i} ;
                    if ~isempty(currEEG)
                        save(strrep(FileNames{currFile}, inputExt, ...
                            ['_processed_' params.paradigm.onsetTags{i} ...
                            rerunExt '.mat']), 'currEEG') ;
                    end
                end
                if params.paradigm.conds.on
                    for i=1:size(eegByConds,2)
                        currEEG = eegByConds{i} ;
                        if ~isempty(currEEG)
                            save(strrep(FileNames{currFile}, inputExt, ...
                                ['_processed_' params.paradigm.conds.groups{i,1} ...
                                rerunExt '.mat']), 'currEEG') ;
                        end
                    end
                end
            end
        end
        if params.outputFormat == 1
            pop_export(EEG, strrep(FileNames{currFile}, inputExt, ...
                ['_processed_IndivTrial' rerunExt '.txt']), ...
                'transpose', 'on', 'precision', 8) ;
            pop_export(EEG, strrep(FileNames{currFile}, inputExt, ...
                ['_processed_AveOverTrials' rerunExt '.txt']), ...
                'transpose', 'on', 'erp', 'on', 'precision', 8) ;
            if size(params.paradigm.onsetTags, 2) > 1
                for i=1:size(eegByTags, 2)
                    if ~isempty(eegByTags{i})
                        pop_export(eegByTags{i}, strrep(FileNames{currFile}, ...
                            inputExt, ['_processed_IndivTrial_' ...
                            params.paradigm.onsetTags{i} rerunExt '.txt']), ...
                            'transpose', 'on', 'precision', 8) ;
                        pop_export(eegByTags{i}, strrep(FileNames{currFile}, ...
                            inputExt, ['_processed_AveOverTrials_' ...
                            params.paradigm.onsetTags{i} rerunExt '.txt']), ...
                            'transpose', 'on', 'erp', 'on', 'precision', 8) ;
                    end
                end
                if params.paradigm.conds.on
                    for i=1:size(eegByConds,2)
                        if ~isempty(eegByConds{i})
                            pop_export(eegByConds{i}, strrep(FileNames{currFile}, ...
                                inputExt, ['_processed_IndivTrial_' ...
                                params.paradigm.conds.groups{i,1} rerunExt, ...
                                '.txt']), 'transpose', 'on', 'precision', 8) ;
                            pop_export(eegByConds{i}, strrep(FileNames{currFile}, ...
                                inputExt, ['_processed_AveOverTrials_' ...
                                params.paradigm.conds.groups{i,1} rerunExt, ...
                                '.txt']), 'transpose', 'on', 'erp', 'on', ...
                                'precision', 8) ;
                        end
                    end
                end
            end
        end
        
    %% ERRORS
    % If HAPPE errors out, check first for common failures, and indicate
    % the cause in the data and pipeline quality assessments. If an
    % uncommon error, will simply fill with error. This allows HAPPE to
    % continue to run over the whole dataset.
    catch ME
        % ADD ERROR TO THE ERROR LOG:
        name = {ME.stack.name} ;
        line = {ME.stack.line} ;
        if size(line,2) > 1
            errList = [sprintf('Line %d in %s; ', line{1:end-1}, ...
                name{1:end-1}) 'Line ' num2str(line{end}) ' in ' name{end}] ;
        else; errList = ['Line ' num2str(line{1}) ' in ' name] ;
        end
        errorLog = [errorLog; {FileNames{currFile}, ME.message, errList}] ;  %#ok<AGROW> 
        
        % CHECK FOR COMMON ERRORS:
        if strcmp(ME.identifier, 'HAPPE:LoadFail'); fill = 'LOAD_FAIL' ;
        elseif strcmp(ME.identifier, 'HAPPE:allICsRej'); fill = 'ALL_IC_REJ' ;
        elseif strcmp(ME.identifier, 'HAPPE:noTags'); fill = 'NO_TAGS' ;
        elseif strcmp(ME.identifier, 'HAPPE:rejOneSeg'); fill = 'REJ_ONE_SEG' ;
        elseif strcmp(ME.identifier, 'HAPPE:AllTrialsRej'); fill = 'ALL_SEG_REJ' ;
        else; fill = 'ERROR' ;
        end
        
        % FILL DATA QC MATRICES WITH ERROR INDICATOR
        for i=1:size(dataQC, 2)
            if isempty(dataQC{currFile,i}); dataQC{currFile,i} = fill ; end
        end
        if params.paradigm.ERP.on
            for i=1:size(dataQC_task, 2)
                if isempty(dataQC_task{currFile, i})
                    dataQC_task{currFile, i} = fill ;
                end
            end
        end
        
        % FILL PIPELINE QC MATRICES WITH ERROR INDICATOR
        if ~reprocess
            lnMeans = qm_ERROR(lnMeans, 1+length(params.lineNoise.neighbors) + ...
                length(params.lineNoise.harms.freqs), currFile) ;
            wavMeans = qm_ERROR(wavMeans, 5+length(params.QCfreqs), ...
                currFile) ;
        end
    end
end

%% GENERATE OUTPUT TABLES
fprintf('Generating quality assessment outputs...\n') ;
cd([srcDir filesep dirNames{contains(dirNames, ...
    'quality_assessment_outputs')}]) ;
rmpath(genpath(cleanlineDir)) ;
try
    % CREATE AND SAVE PIPELINE QUALITY ASSESSMENT
    if ~reprocess
        % Create line noise reduction names.
        lnNames = {'r all freqs pre/post linenoise reduction'} ;
        for i=2:size(params.lineNoise.neighbors, 2)+1
            lnNames{i} = ['r ' num2str(params.lineNoise.neighbors(i-1)) ...
                ' hz pre/post linenoise reduction'] ;
        end
        for i=1:size(params.lineNoise.harms.freqs, 2)
            lnNames{i+size(params.lineNoise.neighbors, 2)+1} = ['r ' ...
                num2str(params.lineNoise.harms.freqs(i)) ...
                ' hz pre/post harmonic reduction'] ;
        end
        
        % Create wavelet thresholding names.
        wavNames = {'RMSE post/pre waveleting', 'MAE post/pre waveleting', ...
            'SNR post/pre waveleting', 'PeakSNR post/pre waveleting', ...
            'r alldata post/pre waveleting'} ;
        for i=1:size(params.QCfreqs,2)
            wavNames{i+5} = ['r ' num2str(params.QCfreqs(i)) ' hz post/pre ' ...
                'waveleting'] ;
        end

        % Concat the Names and QC matrices.
        pipelineQC_names = [(lnNames) (wavNames)] ;
        pipelineQC = [(lnMeans) (wavMeans)] ;
        
        % Save the pipeline QC table.
        pipelineQC_saveName = helpName(['HAPPE_pipelineQC' rerunExt '_' ...
            datestr(now, 'dd-mm-yyyy') '.csv'], '.csv') ;
        writetable(array2table(pipelineQC, 'VariableNames', pipelineQC_names, ...
            'RowNames', FileNames), pipelineQC_saveName, 'WriteRowNames', ...
            true, 'QuoteStrings', true);
    end

    % CREATE AND SAVE DATA QUALITY ASSESSMENT
    % Concat Names and QC matrices according to the presence or absence of
    % multiple onset tags and conditions.
    if params.paradigm.task && size(params.paradigm.onsetTags,2) > 1
        dataQCnames = [dataQCnames dataQCnames_task] ;
        dataQC = [dataQC dataQC_task] ;
        if params.paradigm.conds.on
            dataQCnames = [dataQCnames dataQCnames_conds] ;
            dataQC = [dataQC dataQC_conds] ;
        end
    end
    
    % Save the data QC table.
    dataQC_saveName = helpName(['HAPPE_dataQC' rerunExt '_' datestr(now, ...
        'dd-mm-yyyy') '.csv'], '.csv') ;
    writetable(cell2table(dataQC, 'VariableNames', dataQCnames, 'RowNames', ...
        FileNames), dataQC_saveName, 'WriteRowNames', true, 'QuoteStrings', ...
        true) ;

% ERRORS IN WRITING OUTPUTS:
catch ME
    name = {ME.stack.name} ;
    line = {ME.stack.line} ;
    if size(line,2) > 1
        errList = [sprintf('Line %d in %s; ', line{1:end-1}, ...
            name{1:end-1}) 'Line ' num2str(line{end}) ' in ' name{end}] ;
    else; errList = ['Line ' num2str(line{1}) ' in ' name] ;
    end
    errorLog = [errorLog; {FileNames{currFile}, ME.message, errList}] ;

    % Check for a common error, usually caused by the inability to process
    % any files at all or to completion.
    if strcmp(ME.identifier, 'MATLAB:table:IncorrectNumberOfVarNames')
        fprintf(2, ['ERROR: HAPPE was unable to process any of your files.' ...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
            '\nTo troubleshoot, check your command window and error log.\n']) ;
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