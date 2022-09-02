% determ_chanIDs() - A helper script for setParams.m from HAPPE. Sets the
%                    channels of interest for the run through user input
%                    from the command line. Allows the user to select all
%                    channels, or to choose which subset of channels to
%                    include/exclude.
%
% Usage: 
%   >> [chansAll, chanIDs] = determ_chanIDs()
%
% Inputs:
%
% Outputs:
%   chansAll - String indicating whether the user is using all channels or
%              a subset of channels. If selecting a subset, additionally
%              indicates if the chan_IDs should be included or excluded.
%   chanIDs  - Cell array, with each cell being the name of a channel of
%              interest. The actual contents will vary based on the input
%              parameters and the chans_all selection to best function
%              with HAPPE's current code.
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

function [chansAll, chanIDs] = determ_chanIDs()
fprintf(['Select channels of interest:\n  all = Select all channels in the ' ...
    'data\n  coi = Select a user-specified subset of channels\n']) ;
while true
    % Collect and store user input
    chansAll = input('> ', 's') ;
    
    % If the user requests all channels, set chanIDs to an empty cell array
    % to be filled during the first file run.
    if strcmpi(chansAll, 'all'); chanIDs = {} ; break ;
        
    % If the user wants a subset, request user input to collect channels of interest
    elseif strcmpi(chansAll, "coi")
        fprintf(['Choose an option for entering channels:\n  include = ' ...
            'Include ONLY the entered channel names\n  exclude = Include ' ...
            'every channel EXCEPT the entered channel names\n']) ;
        while true
            chansAll = ['coi_' input('> ', 's')] ;
            if strcmpi(chansAll, 'coi_include') || strcmpi(chansAll, ...
                    'coi_exclude'); break ;
            else; fprintf(['Invalid input: please enter "include" or ' ...
                    '"exclude" (without quotations).\n']) ;
            end
        end
        % Collect channel names using user input
        fprintf(['Enter channels, including the preceding letter, one at ' ...
            'a time.\nPress enter/return between each entry.\nExamples: ' ...
            'E17\n          M1\nWhen you have entered all channels, input ' ...
            '"done" (without quotations).\n']) ;
        chanIDs = unique(UI_cellArray(1, {}), 'stable') ;
        break ;
    else; fprintf(['Invalid input: please enter "all" or "coi" (without ' ...
            'quotations)\n']) ;    
    end
end