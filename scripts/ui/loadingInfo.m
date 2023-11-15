% loadingInfo() - A helper script for HAPPE. Uses user input to determine
%                 and set the appropriate parameters to successfully load
%                 the data in EEGLAB (Delorme &  Makeig, 2004) format. Will
%                 restrict certain formats and nets depending on what HAPPE
%                 currently has the capacity to properly process.
%
% Usage: 
%   >> loadInfo = loadingInfo(params, happeDir)
%
% Inputs:
%   params   - The existing parameters to assist in ensuring the proper
%              loading information is gathered.
%   happeDir - The HAPPE directory path. Needed to access the aquisition
%              layouts.
%
% Outputs:
%   loadInfo - A struct containing all the relevant information to 
%
% Author: A.D. Monachino, PINE Lab at Northeastern University, 2021
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

function loadInfo = loadingInfo(params, happeDir)
loadInfo = params.loadInfo ;
% DETERMINE INPUT FORMAT (FILE TYPE OF RAW DATA):
loadInfo.inputFormat = setFormat() ;
% SET THE LAYOUT (AQUISITION LAYOUT OR NET TYPE):
[loadInfo.layout, loadInfo.correctEGI] = setLayout(loadInfo.inputFormat) ;

% Assume that the channel locations are included. Because this is true for
% most layouts/nets, it is easier to change in in the no-location layouts.
loadInfo.chanlocs.inc = 1;
loadInfo.sys = 0 ;

% IF INPUT FORMAT IS .MAT:
if loadInfo.inputFormat == 1
    % Choose between netstation format or a MATLAB matrix:
    fprintf(['What best describes your data:\n' ...
        '  net station = Data exported using Net Station with additional ' ...
        'information\n  matrix = MATLAB matrix of the data\n']) ;
    loadInfo.NSformat = choose2('matrix', 'net station') ;
    
    % If loading a Netstation layout, confirm that HAPPE has the
    % appropriate channel locations for the given layout. If so, collect
    % the potential EEG variable names. If not, produce an error.
    if loadInfo.NSformat
        if (loadInfo.layout(1) == 1 && loadInfo.layout(2) == 64) || ...
                (loadInfo.layout(1) == 2 && ismember(loadInfo.layout(2), ...
                [32,64,128,256]))
            fprintf(['Enter the potential EEG variable names, one at a ' ...
                'time.\nPress enter/return between each entry.\nWhen ' ...
                'you are finished, enter "done" (without quotations).\nNOTE: ' ...
                'variable names containing "segment" may cause issues.\n']) ;
            loadInfo.NSvarNames = UI_cellArray(1, {}) ;
        else
            error(['HAPPE does not support this layout as a Net Station ' ...
                '.mat file.']) ;
        end
        
    % If loading a MATLAB matrix, collect the following information:
    % whether channel locations are present, and if not, whether to run
    % without chanlocs or to import a chanlocs file; sampling rate for each
    % file (either collectively the same or import a .csv with the sampling
    % rates).
    else
        List = {[1,64], [2,32], [2,64], [2,128], [2,256]} ;
        % Collect channel location information (described above):
        if ~any(cellfun(@(m)isequal(m,loadInfo.layout),List(1,:)))
            [loadInfo.chanlocs.inc, loadInfo.chanlocs.file] = determ_chanLocs() ;
        else
            loadInfo.chanlocs.file = [happeDir filesep 'files' filesep ...
                'acquisition_layout_information' filesep] ;
            if loadInfo.layout(1) == 1
                if loadInfo.layout(2) == 64
                    loadInfo.chanlocs.file = [loadInfo.chanlocs.file 'GSN65v2_0.sfp'] ;
                else; error(['Invalid sensor layout selection.', newline, 'Users ' ...
                        'wishing to use an unsupported layout can run HAPPE through' ...
                        'BEAPP, as described in the HAPPE manuscript.']) ;
                end
            elseif loadInfo.layout(1) == 2
                if ismember(loadInfo.layout(2), [32, 64, 128, 256])
                    loadInfo.chanlocs.file = [loadInfo.chanlocs.file ...
                        'GSN-HydroCel-' num2str(loadInfo.layout(2)+1) '.sfp'] ;
                else; error(['Invalid sensor layout selection.', newline, 'Users ' ...
                        'wishing to use an unsupported layout can run HAPPE through' ...
                        'BEAPP, as described in the HAPPE manuscript.']) ;
                end
            end
        end
        clear('List') ;
        
        % Collect sampling rate information (described above):
        fprintf('Do all your files share the same sampling rate? [Y/N]\n') ;
        loadInfo.srate.same = choose2('n','y') ;
        if loadInfo.srate.same
            fprintf('Sampling rate:\n') ;
            while true
                loadInfo.srate.val = input('> ') ;
                if isnumeric(loadInfo.srate.val); break;
                else; disp('Invalid input: please enter a real number.') ;
                end
            end
        else
            fprintf(['Enter the name of the file containing the sampling ' ...
                'rates for each data file,\nincluding the path and file ' ...
                'extension:\nSee the HAPPE User Guide for how this table ' ...
                'should be formatted.\n']) ;
            while true
                loadInfo.srate.file = input('> ', 's') ;
                if isfile(loadInfo.srate.file); break;
                else; disp('Invalid input: please enter an existing file.');
                end
            end
        end 
    end
    if params.paradigm.task
        loadInfo.eventLoc = formatTask() ;
    end

% .raw
elseif loadInfo.inputFormat == 2
    if params.lowDensity == 1
        error('HAPPE does not currently support .raw low density files.') ;
    else
        if (loadInfo.layout(1) == 1 && loadInfo.layout(2) ~= 64) || ...
                (loadInfo.layout(1) == 2 && ~ismember(loadInfo.layout(2), ...
                [32, 64, 128, 256]))
            error(['The entered number of channels is not supported for this' ...
                ' net as a .raw file.']) ;
        end
    end

% .set
elseif loadInfo.inputFormat == 3

% .cdt
elseif loadInfo.inputFormat == 4

% .mff
elseif loadInfo.inputFormat == 5
    if params.lowDensity; error(['HAPPE does not currently support .mff ' ...
            'low density files.']) ;
    else
    loadInfo.typeFields = {'code'} ;
        fprintf('Do you have additional type fields besides "code"? [Y/N]\n') ;
        if choose2('n','y')
            fprintf(['Enter your type field names, one at a time. When ' ...
                'you are finished entering type fields, enter "done" (without' ...
                ' quotations).\n']) ;
            loadInfo.typeFields = UI_cellArray(2, {'code'}) ;
        end
    end
% .edf
elseif loadInfo.inputFormat == 6
    while true
        fprintf(['Select your system:\n  1 = Mentalab\n  2 = EMOTIV\n' ...
            '  3 = Other\n'])
        loadInfo.sys = input('> ') ;
        if ismember(loadInfo.sys, [1,2,3]); break ;
        else; fprintf('Please enter a valid selection.') ;
        end
    end
    if loadInfo.sys == 2
        fprintf(['Enter the channels in your aquisition set-up, one at' ...
            ' a time,\npressing Enter/Return between entries.\nWhen you ' ...
            'are finished entering channels, enter "done" (without quotations).\n']) ;
        loadInfo.chanlocs.expected = UI_cellArray(1, {}) ;
        loadInfo.chanlocs.locs = fixLocs(loadInfo.chanlocs.expected, happeDir) ;
    else
        [loadInfo.chanlocs.inc, loadInfo.chanlocs.file] = determ_chanLocs() ;
    end
elseif loadInfo.inputFormat == 7
    % Will probably need to load the chanlocs ***
end
end

%-------------------------------------------------------------------------%
%% HELPER FUNCTIONS:
% setFormat() - A helper function for loadingInfo.m from HAPPE.
%               Determines the file format of the data to be loaded in 
%               HAPPE through user input via the command line. Does not 
%               throw errors for or accept invalid input.
function inputFormat = setFormat()
    fprintf(['File format:\n  1 = .mat (MATLAB file)\n  2 = .raw' ...
        ' (Net Station simple binary)\n  3 = .set (EEGLAB format)\n' ...
        '  4 = .cdt (Neuroscan)\n  5 = .mff (EGI)\n  6 = .edf\n' ...
        '  7 = .bdf->.set (Mentalab)\n']) ;
    while true
        inputFormat = input('> ') ;
        if ismember(inputFormat, 1:7); break ;
        else; disp("Invalid input: please enter an integer between 1 and 7.") ;
        end
    end
end

% SetLayout() - A helper function for loadingInfo.m from HAPPE.
%               Determines the net layout, including company and number of
%               electrodes. Does not accept an invalid company, but will
%               accept any number of channels.
function [layout, correctEGI] = setLayout(inputFormat)
layout = [0,0] ;
fprintf(['Acquisition layout type:\n  1 = EGI Geodesic Sensor ' ...
    'Net\n  2 = EGI HydroCel Geodesic Sensor Net\n  3 = Neuroscan Quik-Cap' ...
    '\n  4 = Mentalab Explore\n  5 = Other\n']) ;
while true
    layout(1) = input('> ') ;
    if ismember(layout(1), 1:5); break;
    else; fprintf('Invalid input: please enter an integer between 1 and 5.\n') ;
    end
end

while true
    if inputFormat == 2
        if layout(1) == 1; fprintf(['For .raw EGI GSN files, HAPPE ' ...
                'supports 64 channels.\n']) ;
        elseif layout(1) == 2; fprintf(['For .raw EGI HydroCel GSN ' ...
                'files, HAPPE supports 32, 64, 128, and 256 channels.\n']) ;
        else; error('HAPPE does not support this net for .raw files.') ;
        end
    elseif inputFormat == 4
        if layout(1) ~= 3; error(['To run a .cdt file, the net must be' ...
                ' a NeuroScan layout.']) ; end
    elseif inputFormat == 5
        if layout(1) == 1; fprintf(['For .mff EGI GSN, HAPPE supports ' ...
                '64 channels.\n']) ;
        elseif layout(1) == 2; fprintf(['For .mff EGI HydroCel GSN, ' ...
                'HAPPE supports 32, 64, 128, and 256 channels.\n']) ;
        else
            error('To run a .mff file, the net must be an EGI layout.') ;
        end
    end
    
    fprintf('Number of channels: \n') ;
    layout(2) = input('> ') ;
    
    % Automatically assume that any EGI layout needs correcting
    if ismember(layout(1), [1,2]); correctEGI = 1;
    else; correctEGI = 0;
    end
    
    if inputFormat == 2 && ((layout(1) == 1 && layout(2) ~= 64) || ...
            (layout(1) == 2 && ~ismember(layout(2), [32,64,128,256])))
        fprintf('Invalid input: please enter a supported number of channels.\n') ;
    else; break;
    end  
end
end

% formatTask() - A helper function for setParams.m and HAPPE that
%                collects the path to where task information is stored
%                through user input in the command line.    
function eventLoc = formatTask()
    while true
        eventLoc = input(['Path to .txt files containing task event ' ...
            'info:\n> '], 's') ;
        if exist(eventLoc, 'dir') == 7; break ;
        else; fprintf(['Invalid input: please enter the correct path to ' ...
                'your task event info.\n']) ; 
        end
    end
end

function [inc, chanlocsFile] = determ_chanLocs()
    fprintf(['Do you have a channel locations file for your ' ...
        'data? [Y/N]\nNOTE: A list of supported file types can be ' ...
        'found in the HAPPE User Guide.\n']) ;
    inc = choose2('n', 'y') ;
    if inc
        chanlocsFile = input(['Enter the name of the ' ...
            'channel locations file, including the full path and ' ...
            'file extension.\n> '], 's') ;
    else
        fprintf(['HAPPE functionality is limited without channel ' ...
            'locations. Continue anyway? [Y/N]\n'])
        if choose2('n','y'); chanlocsFile = [] ;
        else; error(['User Termination: No channel locations for .mat ' ...
                'file.']) ;
        end
    end
end