% calcPeaks() - A helper function for HAPPE's generateERPs script that
%               identifies maxima/minima in user-specified windows (max or 
%               min per window, as indicated by the user).
%
% Usage: 
%   >> peakVals = calcPeaks(windows, subGlobal)
%
% Inputs:
%   windows   - An array containing the data for each user-defined window,
%               with each cell being a window from the subject's data.
%   subGlobal - The global timeseries for the current trial in which to
%               locate the peaks.
%
% Outputs:
%   peakVals  - An array containing all the peaks per window in the dataset.
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

function peakVals = calcPeaks(setWindows, subWindows, subGlobal)
peakVals = cell(1,2*size(setWindows,1)+4) ;

% Find Peak (Max or Min) in Windows
for i=1:size(setWindows, 1)
    currWin = subWindows{i} ;
    if strcmpi(setWindows{i, 3}, "max")
        [~, currPeakIndx] = max(currWin(:,2)) ;
    elseif strcmpi(setWindows{i,3}, "min")
        [~, currPeakIndx] = min(currWin(:,2)) ;
    end
    peakVals{1, i*2-1} = currWin(currPeakIndx,2) ;
    peakVals{1, i*2} = currWin(currPeakIndx,1) ;
end
% Find Global Max and Min
[~, globalMaxIndx] = max(subGlobal(:,2)) ;
[~, globalMinIndx] = min(subGlobal(:,2)) ;
% Store Info
peakVals{1, size(peakVals,2)-3} = subGlobal(globalMaxIndx,2) ;
peakVals{1, size(peakVals,2)-2} = subGlobal(globalMaxIndx,1) ;
peakVals{1, size(peakVals,2)-1} = subGlobal(globalMinIndx,2) ;
peakVals{1, size(peakVals,2)} = subGlobal(globalMinIndx,1) ;
end