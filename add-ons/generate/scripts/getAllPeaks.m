% getAllPeaks() - A helper function for HAPPE's generateERPs script that
%                 identifies all peaks in the dataset.
%
% Usage: 
%   >> [listedmaxes, listedmins] = getAllPeaks(currSubGlobal)
%
% Inputs:
%   currSubGlobal - The global timeseries for the current trial in which to
%                   locate the peaks.
%
% Outputs:
%   listedmaxes   - A character array containing a list of all the maxima
%                   in the data and the latency at which they occur.
%   listedmins    - A character array containing a list of all the minima
%                   in the data and the latency at which they occur.
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

% FIND ALL PEAKS IN THE DATA
function [listedmaxes, listedmins] = getAllPeaks(currSubGlobal)
maxes = [currSubGlobal(:,1) islocalmax(currSubGlobal(:,2))] ;
mins = [currSubGlobal(:,1) islocalmin(currSubGlobal(:,2))] ;

listedmaxes = '' ;
listedmins = '' ;
for i=1:size(maxes,1)
    if maxes(i,2)
        listedmaxes = [listedmaxes num2str(maxes(i,1)) '(' ...
            num2str(round(currSubGlobal(i,2),2)) ') '] ;
    end
    if mins(i,2)
        listedmins = [listedmins num2str(mins(i,1)) '( ' ...
            num2str(round(currSubGlobal(i,2),2)) ') '] ;
    end
end
end