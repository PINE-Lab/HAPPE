% Generate Microstates - a post-processing script in association with HAPPE
% to obtain various microstate features from spontaneous EEG Data.
%
% Relies on functions from:
% Poulsen, A. T., Pedroni, A., Langer, N., &  Hansen, L. K. (2018).
% Microstate EEGlab toolbox: An introductionary guide. bioRxiv.
%
% Note: Function pop_micro_selectNmicro.m in the above toolbox (MST1.0) has been 
% modified (line 308. %close(h)) to facilitate saving a figure for users.
%
% If using this script for analyses in publications, cite both HAPPE and
% Microstate EEGlab toolbox.
%
% Developed at Northeastern University's PINE Lab
%
% Author: Priyanka Ghosh, Josh Rodriguez. PINE Lab at Northeastern University, 2023
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUT: Preprocessed '.set' baseline/spontaneous EEG files.
% Before setting up the folders, please ensure that the sampling rate (srate) and the
% number of channels (nbchan) is consistent across all the input '.set' EEG Files. 
% For instance, use the 'resample' feature of HAPPE preprocessing software to
% up/downsample the files that may have a different sampling rate.

%% SET FOLDERS AND PATH
% Use input from the command line to set the path to the data. If an 
% invalid path is entered, repeat until a valid path is entered.
clear ;
fprintf('Preparing HAPPE generateMicroStates...\n') ;
happeDir = strrep(fileparts(which(mfilename('fullpath'))), [filesep '3. ' ...
    'generate'], '') ;
addpath([happeDir filesep '3. generate'], ...
    [happeDir filesep 'files'], ...
    [happeDir filesep 'scripts'], ...
    [happeDir filesep 'scripts' filesep 'ui'], ...
    [happeDir filesep 'scripts' filesep 'support'], ...
    [happeDir filesep 'packages' filesep 'eeglab2024.0']);
    genpath([happeDir filesep 'eeglab2024.0' filesep 'plugins' filesep ...
   'MST1.0']) ;
    pluginDir = dir([happeDir filesep 'packages' filesep 'eeglab2024.0' filesep 'plugins']) ;
    pluginDir = strcat([happeDir filesep 'packages' filesep 'eeglab2024.0'], ...
    filesep, 'plugins', filesep, {pluginDir.name}, ';') ;
    addpath([pluginDir{:}]) ;

%% DETERMINE AND SET PATH TO THE DATA
% Use input from the command line to set the path to the data. If an 
% invalid path is entered, repeat until a valid path is entered.
while true
    srcDir = input('Enter the path to the folder containing the processed dataset(s):\n> ','s') ;
    if isfolder(srcDir); break;
    else; disp("Invalid input: please enter the complete path to the folder containing the dataset(s).") ;
    end
end
cd(srcDir) ;

%% CREATE OUTPUT FOLDERS
% Create the folders in which to store outputs
fprintf('Creating output folders...\n') ;
if ~isfolder([srcDir filesep 'generateMicroStates']); mkdir([srcDir filesep 'generateMicrostates']) ;
end
addpath('generateMicroStates') ;

%% SELECT ACTIVE NUMBER OF MICROSTATES
prompt = 'Do you want to preset the no. of microstates(P) or select interactively based on measures of fit(I)? [P/I]: ';
choice = input(prompt, 's');
if strcmpi(choice, 'P')
Nmicro = input('Enter the number of microstates: ');  % Prompt the user to enter Nmicro
end

%% LOAD DATASETS IN EEGLAB
eeglab
EEGFiles = dir('*.set'); % retrieve a list of all EEG Files in srcDir
lastIndex = -1;  % Initialize the index of the last qualifying dataset
% exclude HAPPE processed EEG files with number of segments less than 15 
for i = 1:length(EEGFiles)
    EEG = pop_loadset('filename', EEGFiles(i).name, 'filepath', srcDir);
    if EEG.trials >= 15 
        EEG.data = reshape(EEG.data, size(EEG.data, 1), []);
        EEG.trials = 1;
        EEG.pnts = size(EEG.data,2);
        EEG.times = 0:size(EEG.data,2)-1;
        [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0);
        lastIndex = i;  % Update the index of the last qualifying dataset
    end
end

if lastIndex > -1
    EEG = pop_loadset('filename', EEGFiles(lastIndex).name, 'filepath', srcDir);
    EEG.data = reshape(EEG.data, size(EEG.data, 1), []);
    EEG.trials = 1;  EEG.pnts = size(EEG.data,2); EEG.times = 0:size(EEG.data,2)-1;                           
end
clearvars -except ALLEEG EEG CURRENTSET prompt choice Nmicro srcDir
cd([srcDir filesep 'generateMicroStates']);

%% SELECT DATA TYPE FOR MICROSTATE ANALYSIS
[EEG, ALLEEG] = pop_micro_selectdata( EEG, ALLEEG, 'datatype', 'spontaneous',...
'avgref', 0, ...
'normalise', 0, ...
'MinPeakDist', 10, ...
'Npeaks', 1000, ...
'GFPthresh', 1, ...
'dataset_idx', 1:length(ALLEEG));

% store data in a new EEG structure
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET); %updates EEGLAB datasets

%% MICROSTATE SEGMENTATION
% select the "GFPpeak" dataset and make it the active set
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, length(ALLEEG)-1, 'retrieve', length(ALLEEG),'study',0);

% Perform the microstate segmentation
EEG = pop_micro_segment( EEG, 'algorithm', 'modkmeans', ...
'sorting', 'Global explained variance', ...
'Nmicrostates', 2:8, ...
'verbose', 1, ...
'normalise', 0, ...
'Nrepetitions', 50, ...
'max_iterations', 1000, ...
'threshold', 1e-06, ...
'fitmeas', 'CV',...
'optimised',1);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

%% PLOT MICROSTATE PROTOTYPE TOPOGRAPHIES
figure; MicroPlotTopo(EEG, 'plot_range', []); savefig('Nmicrostates');

%% REVIEW AND SELECT THE NUMBER OF MICROSTATES  
if strcmpi(choice, 'P')
    EEG = pop_micro_selectNmicro_local(EEG, 'Measures', {'CV','GEV'},'Nmicro', Nmicro);
elseif strcmpi(choice, 'I') %comment line 309 in script pop_micro_selectNmicro.m to save plot
    EEG = pop_micro_selectNmicro_local(EEG, 'Measures', {'CV','GEV'}, 'do_subplots', 0); hFig = gcf; saveas(hFig, 'GEVvsCV.fig'); close(hFig);
else
    disp('Invalid choice. Please enter "P" or "I".');
end

[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

%% IMPORT MICROSTATE PROTOTYPES
for i = 1:length(ALLEEG)-1
 fprintf('Importing prototypes and backfitting for dataset %i\n',i)
 [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',i,'study',0);
 EEG = pop_micro_import_proto(EEG, ALLEEG, length(ALLEEG));

 %% BACKFIT MICROSTATES TO EEG DATASETS
 EEG = pop_micro_fit(EEG, 'polarity', 0 );

 %% TEMPORALLY SMOOTHEN MICROSTATE LABELS
 EEG = pop_micro_smooth( EEG, 'label_type', 'backfit', ...
 'smooth_type', 'reject segments', ...
 'minTime', 30, ...
 'polarity', 0 );

 %% CALCULATE MICROSTATE STATISTICS
 EEG = pop_micro_stats( EEG, 'label_type', 'backfit', ...
 'polarity', 0 );
 [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
 i
end

%% ILLUSTRATING MICROSTATE SEGMENTATION
% Plotting GFP of active microstates for the last subject
 [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'retrieve' ,1,'study',0);
 figure; MicroPlotSegments( EEG, 'label_type' , 'backfit', ...
 'plotsegnos', 'first', 'plot_time', [200 1800], 'plottopos', 1 );
 savefig('microstatesGFP');

 %save('generateMicrostate_outputs.mat');

%% FETCHING MICROSTATE VARIABLES FROM .MAT FILE INTO TABLE FORMAT
generateMicrostate_outputs=[];
for i=1:length(ALLEEG)-1;generateMicrostate_outputs(i).filename=ALLEEG(i).filename; end
for i=1:length(ALLEEG)-1;generateMicrostate_outputs(i).occurence=ALLEEG(i).microstate.stats.Occurence; end
for i=1:length(ALLEEG)-1;generateMicrostate_outputs(i).duration=ALLEEG(i).microstate.stats.Duration; end
for i=1:length(ALLEEG)-1;generateMicrostate_outputs(i).coverage=ALLEEG(i).microstate.stats.Coverage; end
for i=1:length(ALLEEG)-1;generateMicrostate_outputs(i).GEV=ALLEEG(i).microstate.stats.GEV; end
for i=1:length(ALLEEG)-1;
    for j=1:length(ALLEEG(i).microstate.stats.TP);generateMicrostate_outputs(i).TP{j}=ALLEEG(i).microstate.stats.TP(j,:);end
end
table = struct2table(generateMicrostate_outputs);
% Write the table to a CSV file
writetable(table, 'generateMicrostate_outputs.csv');