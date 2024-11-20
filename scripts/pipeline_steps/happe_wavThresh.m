% happe_wavThresh() - A helper function for HAPPE that applies wavelet
%                     thresholding to the data, then assesses pipeline and
%                     data quality on this file.
%
% Usage: 
%   >> [EEG, wavMeans, dataQC] = happe_wavThresh(EEG, params, wavMeans, dataQC, currFile)
%
% Inputs:
%   EEG      - The EEG, in EEGLAB (Delorme & Makeig, 2004) format, that 
%              will have wavelet thresholding applied.
%   params   - The struct containing all the relevent information for 
%              wavelet thresholding, including whether to perform the 
%              legacy or default method and whether performing ERP
%              preprocessing.
%   wavMeans - The matrix holding the pipleine quality assessment metrics.
%   dataQC   - The matrix holding the data quality assessment metrics.
%   currFile - The index of the current file in the file names list.
%
% Outputs:
%   outEEG   - The EEG, in EEGLAB format, after wavelet thresholding.
%   wavMeans - The matrix containing the pipeline quality assessment
%              metrics updated with the metrics for the current file.
%   dataQC   - The matrix containing the data quality assessment metrics
%              updated with the metrics for the current file.
%
% Author: L.J. Gabard-Durnam, PINE Lab at Northeastern University, 2021
%         A.D. Monachino, PINE Lab at Northeastern University, 2021
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

function [EEG, wavMeans, dataQC] = happe_wavThresh(EEG, params, wavMeans, ...
    dataQC, currFile) 
disp('Wavelet thresholding...') ;

%% LEGACY WAVELET - NOT RECOMMENDED
% Kept from HAPPEv1 so comparison/unfinished analyses can be run without 
% needing to switch versions. ICA for clustering data. Uses a soft, global 
% threshold for the wavelets. The wavelet family is coiflet (level 5). 
% Threshold multiplier is used to remove more high frequency noise or for 
% ERP analyses. For details, refer to wICA.m.
if params.wavelet.legacy
    if params.paradigm.ERP.on; threshmultiplier = 3;
    else; threshmultiplier = 1; end
    try
        [wIC, A, ~, ~] = wICA(EEG, 'runica', threshmultiplier, 0, ...
            EEG.srate, 5, 'coif5') ;
    catch wica_err
        if strcmp (['Output argument "wIC" (and maybe others) ' ...
                'not assigned during call to "wICA".'], wica_err.message)
            error(['Error during wICA, most likely due to memory ' ...
                'settings. Please confirm your EEGLAB memory settings ' ...
                'are set according to the description in the HAPPE ReadMe'])
        else; rethrow(wica_err)
        end
    end
    artifacts = A * wIC ;

%% DEFAULT WAVELET:
% Uses a global threshold for the wavelets. Wavelet family is coiflet 
% (level depending). Threshold the wavelet coefficients to generate 
% artifact signals, reconstructing signal as channels x samples format.
else
    % Set wavelet and decomposition level depending on the sampling rate 
    % and paradigm.
    if params.paradigm.ERP.on
        wavFam = 'coif4' ;
        if EEG.srate > 500; wavLvl = 11;
        elseif EEG.srate > 250 && EEG.srate <= 500; wavLvl = 10;
        elseif EEG.srate <=250; wavLvl = 9;
        end
    else
        wavFam = 'bior4.4' ;
        if EEG.srate > 500; wavLvl = 10;
        elseif EEG.srate > 250 && EEG.srate <= 500; wavLvl = 9;
        elseif EEG.srate <=250; wavLvl = 8;
        end
    end

    % Set the threshold rule depending on user input (only available for
    % ERP paradigms).
    if params.paradigm.ERP.on && params.wavelet.softThresh
        ThresholdRule = 'Soft' ;
    else; ThresholdRule = 'Hard' ;
    end
    
    if ~isa(EEG.data, 'double')
        EEG.data = double(EEG.data) ;
    end

    % Use wavelet thresholding to determine artifacts in the data.
    artifacts = wdenoise(reshape(EEG.data, size(EEG.data, 1), [])', wavLvl, ...
        'Wavelet', wavFam, 'DenoisingMethod', 'Bayes', 'ThresholdRule', ...
        ThresholdRule, 'NoiseEstimate', 'LevelDependent')' ;
end

% REMOVE ARTIFACT FROM DATA: Subtract out the wavelet artifact signal from 
% the EEG signal and save the wavcleaned data into an EEGLAB structure. If 
% conducting ERP analyses, filter the data to the user-specified frequency 
% range for analyses purposes only.
if params.paradigm.ERP.on
    preEEG = reshape(pop_eegfiltnew(EEG, params.filt.highpass, ...
        params.filt.lowpass, [], 0, [], 0).data, size(EEG.data, 1), ...
        []) ;
    EEG.data = reshape(EEG.data, size(EEG.data,1), []) - artifacts ;
    postEEG = reshape(pop_eegfiltnew(EEG, params.filt.highpass, ...
        params.filt.lowpass, [], 0, [], 0).data, size(EEG.data, 1), ...
        []) ;
else
    preEEG = reshape(EEG.data, size(EEG.data,1), []) ;
    postEEG = preEEG - artifacts ;
    EEG.data = postEEG ;
end
EEG.setname = 'wavcleanedEEG' ;

% WAVELETING QC METRICS: Assesses the performance of wavelet thresholding.
wavMeans = assessPipelineStep('wavelet thresholding', preEEG, postEEG, ...
    wavMeans, EEG.srate, params.QCfreqs) ;
% DATA QC METRICS: Add the variance retained to the data QC matrix.
dataQC{currFile, 6} = var(postEEG, 1, 'all')/var(preEEG, 1, 'all')*100 ;
end