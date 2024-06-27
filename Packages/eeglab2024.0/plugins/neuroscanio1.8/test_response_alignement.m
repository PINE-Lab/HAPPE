% This script checks the alignemnt of reaction time
% stored in the .dat file and potential reaction time events
% stored in the CNT file. We see that in this case, there
% is a misalignemnt of about 4 samples that needs to be corrected
% The reason for this misalignment is unknown
%
% Arnaud Delorme - Dec 2023

eeglab;

filePath = fileparts(which('test_response_alignement.m'));
EEGkeystroke = pop_loadcnt(     fullfile(filePath, 'data', 'test.cnt'), 'dataformat', 'auto', 'keystroke', 'on');
EEG = pop_loadcnt(     fullfile(filePath, 'data', 'test.cnt'), 'dataformat', 'auto');
EEG = pop_loaddat(EEG, fullfile(filePath, 'data', 'test.dat'), 1600, 4); % the number 4 is a correction factor of 4 sample to align the events

% copy the keystrokes
EEGkeystroke.event(1).eegtype = [];
EEGkeystroke.event(1).response = [];
EEGkeystroke.event(1).resptype = [];
EEGkeystroke.event(1).rtlatency = [];
EEGkeystroke = pop_selectevent(EEGkeystroke, 'type', 'keypad2', 'deleteevents', 'on')
EEG.event = [ EEG.event EEGkeystroke.event];

% look for event "keypad2" overlapped with event "rt2"
pop_eegplot(EEG);
