% addSimERP - a data validation script in association with HAPPE+ER to
% add known ERP waveforms to continuous baseline/resting data.
%
% Developed at Northeastern University's PINE Lab
%
% For a detailed description of this script and user options, please see 
% the following manuscript(s):
%   Monachino, et al., (----)
%
% Author: A.D. Monachino, PINE Lab at Northeastern University, 2022
%
% This file is part of HAPPE.
% Copyright 2018, 2021 Alexa Monachino, Kelsie Lopez, Laurel Gabard-Durnam
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

%% CLEAR THE WORKSPACE
clear ;

%% ADD NECESSARY FOLDERS TO THE PATH
clear ;
fprintf('Preparing HAPPE''s addSimERP script...\n') ;
happeDir = strrep(fileparts(which(mfilename('fullpath'))), [filesep '4. ' ...
    'validate'], '') ;
eeglabDir = [happeDir filesep 'packages' filesep 'eeglab2024.0'] ;
addpath([happeDir filesep '4. validate'], ...
    [happeDir filesep 'files'], ...
    [happeDir filesep 'files' filesep 'ERPs'], ...
    [happeDir filesep 'scripts'], ...
    [happeDir filesep 'scripts' filesep 'ui'], ...
    [happeDir filesep 'scripts' filesep 'support'], ...
    eeglabDir, [eeglabDir filesep 'functions']) ;
rmpath(genpath([eeglabDir filesep 'functions' filesep 'octavefunc'])) ;

%% DETERMINE AND SET PATH TO DATA
% Use input from the command line to set the path to the data. If an 
% invalid path is entered, repeat until a valid path is entered.
while true
    srcDir = input('Enter the path to the folder containing the dataset(s):\n> ','s') ;
    if exist(srcDir, 'dir') == 7; break ;
    else; fprintf(['Invalid input: please enter the complete path to the ' ...
            'folder containing the dataset(s).\n']) ;
    end
end

%% DETERMINE CHANNELS OF INTEREST
[chansAll, chanIDs] = determ_chanIDs() ;

%% DETERMINE ERP TO INSERT
fprintf(['Which ERP waveform would you like to add to your data?\n  VEP = ' ...
    'Visual Evoked Potential\n  p3-10 = P3 with Amplitude 10, Variation 4\n' ...
    '  p3-20 = P3 with Amplitude 20, Variation 7\n']) ;
suppERPs = {'VEP', 'p3-10', 'p3-20'} ;
while true
    selERP = upper(input('> ', 's')) ;
    if any(strcmpi(selERP, suppERPs)); break ;
    else; fprintf(['Invalid input: please enter one of the following: ' ...
            sprintf('%s, ', suppERPs{1:end-1}), suppERPs{end} '.\n']) ;
    end
end

%% CREATE OUTPUT FOLDERS
% Create the folder in which to store final outputs.
fprintf(['Creating output folder "+' selERP '"...\n']) ;
if ~isfolder([srcDir filesep '+' selERP])
    mkdir([srcDir filesep '+' selERP]) ;
    fprintf(['Output folder "+' selERP '" created.\n']) ;
else; fprintf(['Output folder "+' selERP '" already exists.\n']) ;
end

%% LOAD ERP TIMESERIES AND EVENT FILES
fprintf('Loading ERP timeseries and events...\n') ;
load([happeDir filesep 'files' filesep 'ERPs' filesep selERP '_ts.mat']) ;
load([happeDir filesep 'files' filesep 'ERPs' filesep selERP '_events.mat']) ;

cd(srcDir) ;
FileNames = {dir('*.set').name} ;

%% LOOP OVER EACH FILE
cd(srcDir) ;
for currFile=1:length(FileNames)
    fullERP = ERP ;
    tempevent = event ;
    
    %% LOAD AND VALIDATE THE FILE
    try
        EEGraw = load('-mat', FileNames{currFile}) ;
        if isfield(EEGraw, 'EEG'); EEGraw = EEGraw.EEG; end
        EEG = eeg_checkset(EEGraw) ;
    catch
        fprintf(2, ['ERROR: unable to load ' FileNames{currFile} ' ...\n']) ;
        continue
    end
    
    %% CHECK THAT THE DATA IS CONTINUOUS
    if ndims(EEG.data) > 2 %#ok<ISMAT>
        fprintf(2, 'ERROR: data must be continuous!\n') ;
        continue
    end
    
    %% RESAMPLE TO 250
    EEG = pop_resample(EEG, 250) ;
    
    %% DETERMINE MULTIPLIER FOR ERP
    if size(EEG.data,2) > size(ERP,2)
        lngthmult = floor(size(EEG.data,2)/size(ERP,2)) ;
    else
        fprintf(2, ['ERROR: Data cannot have fewer timepoints than the' ...
            ' ERP!\n']) ;
        continue ;
    end
    
    %% ADD EVENTS TO THE EVENT LIST
    lats = {event.latency} ;
    inits = {event.init_time} ;
    for i=1:lngthmult-1
        fullERP = [fullERP ERP] ;
        for j=1:length(event)
            tempevent(i*length(event)+j).type = event(j).type ;
            tempevent(i*length(event)+j).latency = (lats{2}-lats{1})*((i*length(event)-1)+j)+1;
            tempevent(i*length(event)+j).duration = 1 ;
            tempevent(i*length(event)+j).init_index = 1 ;
            tempevent(i*length(event)+j).init_time = ((i*length(event)+j)-1)*(inits{2}-inits{1}) ;
        end
    end
    
    %% REMOVE REMAINDER DATA POINTS
    tempdata = EEG.data ;
    tempdata(:,end-(size(EEG.data,2)-size(fullERP,2)-1):end) = [] ;
    
    %% GET COLUMN NUMBERS FOR COI:
    if strcmpi(chansAll, 'all')
        chanInds = 1:size(EEG.data, 1) ;
    else
        chanInds = [] ;
        for i=1:length(chanIDs)
            chanInds = [chanInds find(ismember({EEG.chanlocs.labels}, ...
                chanIDs(i)))] ;                                             %#ok<*AGROW> 
        end
    end
    
    %% ADD ERP TO EACH CHANNEL
    for i=1:size(chanInds,2)
        tempdata(chanInds(i),:) = tempdata(chanInds(i),:) + fullERP ;
    end
    
    %% ADD UPDATED EVENTS AND DATA TO EEG
    EEG.data = tempdata ;
    EEG.event = tempevent ;
    EEG = eeg_checkset(EEG) ;
    
    %% SAVE UPDATED EEG
    pop_saveset(EEG, 'filename', strrep(FileNames{currFile}, '.set', ...
        ['+' selERP '.set']), 'filepath', [srcDir filesep '+' selERP]) ;
end