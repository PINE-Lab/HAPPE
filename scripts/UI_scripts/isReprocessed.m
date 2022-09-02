% isReprocessed() - A helper script for HAPPE. Uses user input to
%                   determine whether the data is raw or has already been
%                   processed. Additionally sets a rerun text extension to
%                   aid in proper saving of HAPPE outputs.
%
% Usage: 
%   >> [reprocessing, rerunExt] = isReprocessed()
%
% Inputs:
%
% Outputs:
%   reprocessing - A binary number [0|1] indicating whether or not data is
%                  being reprocessed. 
%   rerunExt     - A string/character array that is added to the end of
%                  file name to distinguish it from potentially existing
%                  files. Is empty if not reprocessing.
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

function [reprocessing, ranMuscIL, rerunExt] = isReprocessed()
% ASK WHETHER DATA IS RAW OR ALREADY PROCESSED VIA COMMAND LINE INPUT
fprintf(['Select an option:\n  raw = Run on raw data from the start\n' ...
    '  reprocess = Run on HAPPE-processed data starting post-artifact ' ...
    'reduction\n']) ;
reprocessing = choose2('raw', 'reprocess') ;

% IF REPROCESSING...
if reprocessing
    % ASK TO RUN POST-WAV OR POST-MUSCIL
    fprintf(['Run from wavelet-cleaned data or muscIL-cleaned data?\n  ' ...
        'wavclean = Run starts after wavelet-cleaning\n  muscIL = Run ' ...
        'starts after muscIL (only available if muscIL-cleaned data exists)\n']) ;
    ranMuscIL = choose2('wavclean', 'muscIL') ;

    % ASK TO EITHER OVERWRITE FILES OR SAVE NEW ONES.
    fprintf(['Files, such as processed data and quality metrics may already' ...
        ' exist for this dataset.\n  overwrite = Overwrite existing files\n' ...
        '  new = Save new files\n']) ;
    % IF SAVING NEW FILES, create a rerun tag either using the default or
    % through user input.
    if choose2('overwrite', 'new')
        fprintf(['Use a default or custom label for processed data?\n' ...
            '  default = Default name (_rerun_dd-mm-yyyy)\n  custom =' ...
            ' Create your own file label, starting with an underscore\n']) ;
        if choose2('custom', 'default')
            rerunExt = ['_rerun_' datestr(now, 'dd-mm-yyyy')] ;
        else; rerunExt = input('Enter custom label:\n> ', 's') ;
        end
    % IF OVERWRITING FILES, set the rerun tag to empty.
    else; rerunExt = '' ;    
    end
    
% IF RUNNING FROM THE START, set the rerun tag to empty.
else; rerunExt = '' ; ranMuscIL = 0;
end
end