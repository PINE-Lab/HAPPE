% happe_detectBadChans() - A helper script for HAPPE that detects and
%                          rejects bad channels from the data.
%
% Usage: 
%   >> [EEG, dataQC] = happe_detectBadChans(EEG, params, pt)
%
% Inputs:
%   EEG    - The EEG, in EEGLAB (Delorme & Makeig, 2004) format, from which
%            to remove bad channels.
%   params - The struct of the HAPPE params struct containing all the
%            relevent information for detecting bad channels, including 
%            whether to perform the legacy or default method of rejection.
%   pt     - An indicator of whether performing bad channel rejection
%            before or after wavelet thresholding to assist in setting the
%            correct parameters for the detection functions.
%
% Outputs:
%   EEG    - The EEG, in EEGLAB (Delorme & Makeig, 2004) format, after bad
%            channel rejection.
%
% Author: A.D. Monachino, PINE Lab at Northeastern University, 2021
%
% This file is part of HAPPE.
% Copyright 2018, 2023 Alexa Monachino, Kelsie Lopez, Laurel Gabard-Durnam
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

function EEG = happe_detectBadChans(EEG, params, pt)
disp('Detecting bad channels...') ;
% If low-density data, detect bad channels using the methods optimized 
% in HAPPILEE.
if params.lowDensity
     EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion', ...
         'off', 'ChannelCriterion', .7, 'LineNoiseCriterion', ...
         2.5, 'Highpass', 'off', 'BurstCriterion', 'off', ...
         'WindowCriterion', 'off', 'BurstRejection', ...
         'off', 'Distance', 'Euclidian') ;
    EEG = pop_rejchan(EEG, 'elec', 1:EEG.nbchan, ...
        'threshold', [-2.75 2.75], 'norm', 'on', 'measure', ...
        'spec', 'freqrange', [1 100]) ;

% Otherwise, detect bad channels using methods optimized for HAPPE v2
else
    switch pt
        case 1; upThresh = 1.8935; chanCrit = 0.485; lnCrit = 7.1;
        case 2; upThresh = 2.1316; chanCrit = 'off'; lnCrit = 'off' ;
    end
    EEG = pop_rejchan(EEG, 'elec', 1:EEG.nbchan, 'threshold', [-5 upThresh], ...
        'norm', 'on', 'measure', 'spec', 'freqrange', [1 100]) ;
    EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion', 'off', 'ChannelCriterion', ...
        chanCrit, 'LineNoiseCriterion', lnCrit, 'Highpass', 'off', ...
        'BurstCriterion', 'off', 'WindowCriterion', 'off', 'BurstRejection', ...
        'off', 'Distance', 'Euclidian') ;
end
end