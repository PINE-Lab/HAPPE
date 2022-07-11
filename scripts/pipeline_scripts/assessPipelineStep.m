% assessPipelineStep() - A helper function for HAPPE that assesses the
%                        performance of a given pipeline step on a data
%                        file by calculating values such as
%                        cross-correlation, signal to noise ratio, peak
%                        signal to noise ratio, root mean squared error,
%                        and mean average error. Note that for line noise
%                        reduction, only cross-correlation is performed.
%
% Usage: 
%   >> means = assessPipelineStep(key, preEEG, postEEG, means, srate, freqsofinterest)
%
% Inputs:
%   key             - A string indicating what pipeline step is being
%                     assessed (e.g., line noise reduction, wavelet
%                     thresholding).
%   preEEG          - The EEG data matrix prior to the pipeline step.
%   postEEG         - The EEG data matrix post the pipeline step.
%   means           - The matrix holding the values of previous files (if
%                     any) to add the newly calculated values to.
%   srate           - The sampling rate of the current EEG file.
%   freqsofinterest - The frequencies of interest on which to calculate
%                     values.
%
% Outputs:
%   means           - The updated matrix containing the newly calculated
%                     values assessing the current pipeline step.
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

function means = assessPipelineStep(key, preEEG, postEEG, means, srate, ...
    freqsofinterest)
fprintf(['Evaluating pipeline performance on ' key '...\n']) ;
% Calculate Cross Correlation:
cc = corrcoef(postEEG, preEEG) ;
ccFreqs = mean(mscohere(postEEG', preEEG', 1000, 0, freqsofinterest, ...
    srate), 2, 'omitnan')' ;

% LINE NOISE REDUCTION: If assessing line noise reduction, simply
% update the means matrix with the cross-correlation values.
if strcmpi(key, 'line noise reduction')
    means = [means; cc(1,2) ccFreqs] ;

% OTHER STEPS: For other pipeline steps, calculate additional measures
% (RMSE, MAE, SNR, PeakSNR) and add all values plus cross-correlation
% values to the means matrix.
else
    % Calculate Signal to Noise Ratio (SNR) and Peak Signal to Noise
    % Ratio (PeakSNR)
    [SNR, pSNR] = calcSNR_PSNR(preEEG, postEEG, 2) ;

    % Update means matrix with RMSE, MAE, SNR, peakSNR, and
    % cross-correlation across data and channels. RMSE and MAE are
    % calculated within this step.
    means = [means; mean(realsqrt(mean((preEEG - postEEG) .^ 2, 2))) ...
        mean(mean(abs(preEEG - postEEG), 2)) SNR pSNR cc(1,2) ccFreqs] ;
end
end