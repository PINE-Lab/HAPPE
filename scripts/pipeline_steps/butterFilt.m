% butterFilt() - A stripped version of pop_basicfilter and its subroutines 
%                from ERPLAB (Lopez-Calderon & Luck, 2014). Uses just the 
%                bandpass butterworth filter from the original function and
%                applies it to all channels across the entire timeseries. 
%                Equivalent to pop_basicfilter(EEG, [1:length(EEG.chanlocs)],
%                'Filter', 'bandpass', 'Design', 'butter', 'Cutoff', 
%                [highpass lowpass], 'RemoveDC', 'off', 'Boundary', []) ;
%                See ERPLAB for more documentation.
%
% Usage: 
%   >> EEG = butterFilt(EEG, lowPass, highPass) ;
%
% Inputs:
%   EEG      - The EEGLAB EEG struct prior to filtering
%   lowPass  - The low-pass filter to use on the data (typically 30 - 40
%              Hz)
%   highPass - The high-pass filter to use on the data (typically 30 - 40
%              Hz)
%
% Outputs:
%   EEG      - The EEGLAB EEG struct after filtering
%
% Orignial Authors: Javier Lopez-Calderon and Steven Luck, UC Davis, 2009
% Copyright 2007 The Regents of the University of California
%
% Adapted for HAPPE by: A.D. Monachino, PINE Lab at Northeastern 
%                       University, 2022
%
% This file is part of HAPPE.
% Copyright 2018-2022 Alexa Monachino, Kelsie Lopez, Laurel Gabard-Durnam
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

function EEG = butterFilt(EEG, lowPass, highPass)
% Defines the filter order
filtOrder = 2 ;
% Creates an array of the indicies of all channels included in the data
chanArray = [1:length(EEG.chanlocs)] ;

% ERROR CHECKING: confirms that the values of the filter order, lowpass,
% and highpass are all valid to perform a butterworth bandpass filter. If
% any of the values do not meet the necissary criteria, throw an error.
if highPass >= EEG.srate/2 || lowPass >= EEG.srate/2
    error('ERROR: Cutoff frequency cannot be >= srate/2');
end
% if filtOrder*3 > EEG.pnts
%     error('ERROR: filter order too high. Samples must be 3x filter order') ;
% end
% if mod(filtOrder,2) ~= 0; error('ERROR: Filter order must be even'); end
if lowPass/(EEG.srate/2) >= 1 || highPass/(EEG.srate/2) >= 1 || ...
        lowPass < 0 || highPass < 0
    error('ERROR: Cutoff frequencies must be >0 and <srate/2') ;
end

% FILTER AS WITH ERPLAB'S filter_tf FUNCTION:
[b1, a1] = butter(filtOrder/2, lowPass/(EEG.srate/2)) ;
[b2, a2] = butter(filtOrder/2, highPass/(EEG.srate/2), 'high') ;
b = [b1; b2] ;
a = [a1; a2] ;

warning off MATLAB:singularMatrix
fprintf('Filtering ...') ;
if size(b,1) > 1
    EEG.data(chanArray, 1:EEG.pnts, 1) = filtfilt(b(1,:), a(1,:), ...
        EEG.data(chanArray, 1:EEG.pnts, 1)')';
    EEG.data(chanArray, 1:EEG.pnts, 1) = filtfilt(b(2,:), a(2,:), ...
        EEG.data(chanArray, 1:EEG.pnts, 1)')';
end
end