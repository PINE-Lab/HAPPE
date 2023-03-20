% happe_detectBadChans() - A helper script for HAPPE that detects and
%                          rejects bad channels from the data and assesses
%                          the data quality of this file.
%
% Usage: 
%   >> [EEG, dataQC] = happe_detectBadChans(inEEG, params, dataQC, chanIDs, currFile)
%
% Inputs:
%   inEEG    - The EEG, in EEGLAB (Delorme & Makeig, 2004) format, from 
%              which to remove bad channels.
%   params   - The struct of the HAPPE params struct containing all the
%              relevent information for detecting bad channels, including 
%              whether to perform the legacy or default method of 
%              rejection.
%   dataQC   - The matrix holding the data quality assessment metrics.
%   chanIDs  - The matrix holding the channels of interest for the current 
%              subject.
%   currFile - The index of the current file in the list of file names.
%              Used in updating the dataQC variable.
%
% Outputs:
%   EEG      - The EEG, in EEGLAB (Delorme & Makeig, 2004) format, after 
%              bad channel rejection.
%   dataQC   - The matrix containing the data quality assessment metrics
%              updated with the metrics for the current file.
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

function [EEG, dataQC] = happe_detectBadChans(inEEG, params, dataQC, ...
    chanIDs, currFile)
disp('Detecting bad channels...') ;
% LEGACY DETECTION (from HAPPE 1.0 - NOT RECOMMENDED):
% Conducts crude bad channel detection using spectrum criteria and
% 3 SDeviations as channel outlier threshold (twice). This option is 
% not available for low density layouts.
if params.badChans.legacy
    EEG = pop_rejchan(inEEG, 'elec', [1:length(chanIDs)], 'threshold', [-3 3], ...
        'norm', 'on', 'measure', 'spec', 'freqrange', [1 125]) ;
    EEG = pop_rejchan(EEG, 'elec', [1:inEEG.nbchan], 'threshold', ...
        [-3 3], 'norm', 'on', 'measure', 'spec', 'freqrange', [1 125]);

% DEFAULT DETECTION:
% Conducts bad channel detection optimized for HAPPE v2, HAPPE+ER, and 
% HAPPILEE.
else
    % Check for flatline channels.
    EEG = pop_clean_rawdata(inEEG, 'FlatlineCriterion', 3, ...
        'ChannelCriterion', .1, 'LineNoiseCriterion', ...
        20, 'Highpass', 'off', 'BurstCriterion', 'off', ...
        'WindowCriterion', 'off', 'BurstRejection', 'off', ...
        'Distance', 'Euclidian') ;
    
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
        EEG = pop_rejchan(EEG, 'elec', 1:EEG.nbchan, 'threshold', [-5 3.5], ...
            'norm', 'on', 'measure', 'spec', 'freqrange', [1 100]) ;
        EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion', 'off', 'ChannelCriterion', ...
            .8, 'LineNoiseCriterion', 6, 'Highpass', 'off', 'BurstCriterion', ...
            'off', 'WindowCriterion', 'off', 'BurstRejection', 'off', 'Distance', ...
            'Euclidian') ;
    end
end

% BAD CHANNEL DETECTION QC: Assesses the performance of bad channel 
% rejection. Metrics include the number of channels included in processing
% (good channels), the percent of good channels selected from the original 
% user-selected set, and the list of channels, if any, marked as bad.
disp('Evaluating bad channel detection...') ;
dataQC{currFile,3} = size(EEG.chanlocs,2) ;
dataQC{currFile, 4} = dataQC{currFile,3}/dataQC{currFile,2}*100 ;
if ~isempty(EEG.chanlocs)
    badChans = setdiff(chanIDs, {EEG.chanlocs.labels}, ...
        'stable') ;
    if isempty(badChans); dataQC{currFile, 5} = 'None' ;
    else; dataQC{currFile,5} = [sprintf('%s ', badChans{1:end-1}), badChans{end}] ;
    end
else; dataQC{currFile,5} = 'NA' ;
end
end