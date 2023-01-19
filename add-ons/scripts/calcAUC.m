% calcAUC() - A helper function for HAPPE's generateERPs script that
%             calculates the area under the curve and 50% area under the
%             curve for a series of windows using user-specified latency
%             values as the window boundaries. Additionally calculates the
%             global area under the curve and global 50% area under the
%             curve.
%
% Usage: 
%   >> [aucVals, fiftyAL] = calcAUC(windows, subGlobal)
%
% Inputs:
%   windows   - An array containing the data for each user-defined window,
%               with each cell being a window from the subject's data.
%   subGlobal - The global timeseries for the current trial on which to
%               calculate the global area under the curve and global 50%
%               area under the curve.
%
% Outputs:
%   aucVals   - The area under the curve for each window and the global
%               area under the curve.
%   fiftyAL   - The 50% area under the curve, and the latency at which it
%               was reached for each window and the global timeseries.
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

function [aucVals, fiftyAL] = calcAUC(windows, subGlobal)
aucVals = cell(1, size(windows,2)+1) ;
fiftyAL = cell(1, 2*size(windows,2)+2) ;
% Calculate AUC & 50% AUC
for i=1:size(windows,2)
    currWin = windows{i} ;
    auc = cumtrapz(currWin(:,1), abs(currWin(:,2))) ;
    aucVals{1,i} = auc(length(auc)) ;
    [~, closestIndx] = min(abs(auc-auc(length(auc))/2)) ;
    fiftyAL{1,i*2} = currWin(closestIndx,1) ;
    fiftyAL{1,i*2-1} = auc(closestIndx) ;
end

% Calculate Global AUC and 50% AUC
auc = cumtrapz(subGlobal(:,1), abs(subGlobal(:,2))) ;
aucVals{1, end} = auc(length(auc)) ;
[~, closestIndx] = min(abs(auc-auc(length(auc))/2)) ;
fiftyAL{1, end} = subGlobal(closestIndx,1) ;
fiftyAL{1, end-1} = auc(closestIndx) ;
end