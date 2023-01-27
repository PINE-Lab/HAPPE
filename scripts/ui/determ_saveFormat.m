% determ_saveFormat() - A helper function for setParams.m from HAPPE. Sets
%                       the save format for the processed files based on
%                       user input through the command line. Options
%                       include .mat, .set, .and .txt. Does not throw
%                       and error on or accept invalid inputs.
%
% Usage: 
%   >> saveFormat = determ_saveFormat()
%
% Outputs:
%   saveFormat - An integer (1, 2, or 3) representing the file format in
%                which to save the processed data.
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

function saveFormat = determ_saveFormat()
fprintf(['Format to save processed data:\n  1 = .txt file (electrodes as ' ...
    'columns, time as rows) - Choose this for ERP timeseries\n  2 = .mat' ...
    ' file (MATLAB format)\n  3 = .set file (EEGLAB format)\n']) ;
while true
    saveFormat = input('> ') ;
    if ismember(saveFormat, [1, 2, 3]); break ;
    else; disp("Invalid input: please enter 1, 2, or 3.") ;
    end
end
end