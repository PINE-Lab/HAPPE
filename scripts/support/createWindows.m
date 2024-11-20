% createWindows() - A helper function for HAPPE's generateERPs script that
%                   creates arrays containing the data points for the 
%                   windows of interest as indicated by the user in set 
%                   parameters.
%
% Usage: 
%   >> [currSubWindows, currSubGlobal] = createWindows(lats, subAve_noBL, windows)
%
% Inputs:
%   lats           - An array containing a list of the unique latencies in
%                    the dataset.
%   subAve_noBL    - The average ERP across the channels of interest for
%                    the current subject/trial, excluding the baseline
%                    period.
%   windows        - An array containing the user-specified parameters for
%                    the windows of interest.
%
% Outputs:
%   currSubWindows - An array of arrays, with each cell being the data
%                    matrix for each user-specified window.
%   currSubGlobal  - The array of data for the current subject.
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

function [currSubWindows, currSubGlobal] = createWindows(lats, subAve_noBL, windows)
lats = round(lats) ;
currSubGlobal = [lats(find(lats==0):size(lats,1),:) subAve_noBL] ;
currSubWindows =  {};
for i=1:size(windows, 1)
    currSubWindows{i} = currSubGlobal(find(currSubGlobal==str2double(windows{i,1})):find(currSubGlobal==str2double(windows{i,2})),:) ; %#ok<*AGROW> 
end
end