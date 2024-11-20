% Wakeman & Henson Data analysis: Extract EEG data and import events and channel location
% Data available at https://sccn.ucsd.edu/eeglab/download/ds000117_sub-01.zip
% (subject 1 of https://nemar.org/dataexplorer/detail?dataset_id=ds000117)
%
% EEGLAB plugin/extensions required
% - File-io
% - Fieldtrip
% - bids-matlab-tools
% - picard
%
% Authors: Arnaud Delorme, SCCN, 2022

clear;
clear globals;
[ALLEEG, EEG, CURRENTSET] = eeglab; % start EEGLAB

% Comment one of the two lines below to process EEG or MEG data
%chantype = { 'megmag' }; % process MEG megmag channels
%chantype = { 'megplanar' }; % process MEG megplanar channels
chantype = { 'eeg' }; % process EEG

% Paths below must be updated to the files on your environment.
dataPath = '/Users/arno/Downloads/sub-01';
if ~exist(dataPath), error('Data not found; download it from the web site and update the path'); end
filenameEEG = fullfile( dataPath, 'ses-meg','meg','sub-01_ses-meg_task-facerecognition_run-01_meg.fif');
filenameFID = fullfile( dataPath, 'ses-meg','meg','sub-01_ses-meg_coordsystem.json');
filenameMRI = fullfile( dataPath, 'ses-mri','anat','sub-01_ses-mri_acq-mprage_T1w.nii.gz');

% Step 1: Importing data with FileIO
EEG = pop_fileio(filenameEEG);
EEG = pop_select(EEG, 'chantype', chantype);
EEG = pop_select(EEG, 'rmchannel', { 'EEG061' 'EEG062' 'EEG063' 'EEG064' }); % remove EOG and EKG channels

% Preprocess and run ICA (so one may be localized)
EEG = pop_resample(EEG, 100);
EEG = pop_eegfiltnew(EEG, 1, 0);
EEG = pop_reref(EEG, []);
EEG = pop_runica( EEG , 'picard', 'maxiter', 500, 'pca', 20); % PCA not recommended if you have enough data

% Import MRI. The MRI has an associated file with the coordinates of the
% fiducials in MRI space. These are automatically imported as well. 
% The alignment should always be checked and if necessary adjusted.
EEG = pop_dipfit_headmodel( EEG, filenameMRI, 'plotmesh', 'scalp');
EEG = pop_dipfit_settings( EEG, 'coord_transform', 'alignfiducials');

% Localize first 10 component for speed
EEG = pop_multifit(EEG, 1:20,'threshold', 100, 'dipplot','off'); 
pop_dipplot(EEG, [], 'normlen', 'on');
