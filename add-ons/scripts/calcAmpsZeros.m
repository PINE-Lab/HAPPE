% calcAmpsZeros() - A helper function for HAPPE's generateERPs script that
%                 calculates the mean amplitude for each window in a set of
%                 zero-bound windows. Lists each window's bounds in curly 
%                 brackets ({}), along with the associated mean amplitude 
%                 for that window.
%
% Usage: 
%   >> meanAmps = calcAmpsZeros(windows)
%
% Inputs:
%   windows - An array containing the zero-bound windows in which to 
%             calculate the mean amplitude.
%
% Outputs:
%   meanAmps - A character array containing the mean amplitude for each
%              zero-bound window in windows.
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

function meanAmps = calcAmpsZeros(windows)
meanAmps = '' ;
for i=1:size(windows,2)
    currWin = windows{i} ;
    if size(currWin, 1) == 1; meanAmps = currWin(1,2) ;
    else; meanAmps = mean(currWin(:,2)) ;
    end
    meanAmps = [meanAmps '{' num2str(currWin(1,1)) '-' ...
        num2str(currWin(size(currWin,1), 1)) ': ' ...
        num2str(meanAmps(length(meanAmps))) '} '] ;                         %#ok<*AGROW> 
end
end