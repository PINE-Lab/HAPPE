% calcAUCZeros() - A helper function for HAPPE's generateERPs script that
%                  calculates the area under the curve and 50% area under 
%                  the curve for a series of windows using zero-bound
%                  latency values as the window boundaries.
%
% Usage: 
%   >> [aucValsZeros, fiftyALZeros] = calcAUCZeros(zeroCrossWins)
%
% Inputs:
%   zeroCrossWins - An array containing the data for each zero-bound
%                   window, with each cell being a window from the 
%                   subject's data.
%
% Outputs:
%   aucVals   - The area under the curve for each window listed as a
%               character array.
%   fiftyAL   - The 50% area under the curve, and the latency at which it
%               was reached for each window, listed as a character array.
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

function [aucVals, fiftyAL] = calcAUCZeros(zeroCrossWins)
% Calculate AUC & 50% AUC for each Window
aucVals = '' ;
fiftyAL = '' ;
for i=1:size(zeroCrossWins,2)
    currWin = zeroCrossWins{i} ;
    if size(currWin, 1) == 1; auc = currWin(1,2) ;
    else; auc = cumtrapz(currWin(:,1), abs(currWin(:,2))) ;
    end
    [minVal, closestIndx] = min(abs(auc-auc(length(auc))/2)) ;
    aucVals = [aucVals '{' num2str(currWin(1,1)) '-' num2str(currWin(size(currWin,1), ...
        1)) ': ' num2str(auc(length(auc))) '} '] ;
    fiftyAL = [fiftyAL '{' num2str(currWin(1,1)) '-' ...
        num2str(currWin(size(currWin,1),1)) ': ' num2str(auc(closestIndx)) ...
        ' at ' num2str(currWin(closestIndx,1)) '} '] ;
end
end