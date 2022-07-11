% createZeroWins() - A helper function for HAPPE's generateERPs script that
%                    creates arrays containing the data points for windows
%                    defined by zero-crossings.
%
% Usage: 
%   >> zeroCrossWins = createZeroWins(currSubGlobal)
%
% Inputs:
%   currSubGlobal - The array of data for the current subject.
%
% Outputs:
%   zeroCrossWins - An array containing the data for each window, with each
%                   cell being an individual zero-bound window.
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

function zeroCrossWins = createZeroWins(currSubGlobal)
% Find the Zero Crossings
zeroCrossings = 0 ;
for i=2:size(currSubGlobal,1)
    prev = currSubGlobal(i-1,2) ;
    curr = currSubGlobal(i,2) ;
    if curr == 0
        zeroCrossings = [zeroCrossings currSubGlobal(i,1)] ;                %#ok<*AGROW> 
    elseif (prev < 0 && curr > 0) || (prev > 0 && curr < 0)
        zeroCrossings = [zeroCrossings currSubGlobal(i-1,1)] ;
    end
end
zeroCrossings = [zeroCrossings currSubGlobal(length(currSubGlobal),1)] ;

% Create Windows from the Zero Crossings
zeroCrossWins =  {} ;
for i=2:size(zeroCrossings, 2)
    zeroCrossWins{i-1} = currSubGlobal(find(currSubGlobal==zeroCrossings(i-1)):find(currSubGlobal==zeroCrossings(i)),:) ;
end
end