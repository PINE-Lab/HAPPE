% happe_reduceLN() - A helper script for HAPPE that reduces line noise
%                    and assesses HAPPE's performance on this file. Uses
%                    Cleanline (Mullen, 2012) as an EEGLAB (Delorme &
%                    Makeig, 2004) plugin.
%
% Usage: 
%   >> [outEEG, lnMeans] = happe_reduceLN(inEEG, lnParams, lnMeans)
%
% Inputs:
%   inEEG    - The EEG, in EEGLAB (Delorme & Makeig, 2004) format, from 
%              which line noise is reduced.
%   lnParams - The substruct of the HAPPE params struct containing all the
%              relevent information for processing line noise, including 
%              whether to perform the legacy or default method and what
%              frequency the line noise occurs at.
%   lnMeans  - The matrix holding the pipleine quality assessment metrics.
%
% Outputs:
%   outEEG   - The EEG, in EEGLAB (Delorme & Makeig, 2004) format, after 
%              line noise reduction.
%   lnMeans  - The matrix containing the pipeline quality assessment
%              metrics updated with the metrics for the current file.
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

function [outEEG, lnMeans] = happe_reduceLN(inEEG, lnParams, lnMeans)
fprintf('Reducing line noise...\n') ;
if lnParams.cl
    % If reducing harmonics, add the user-specified harmonics to the
    % frequencies to examine. Otherwise, look at just the user-specified
    % frequency and that frequency x2.
    if lnParams.harms.on
        lineFreqs = unique([lnParams.freq, lnParams.freq*2, ...
            lnParams.harms.freqs]) ;
    else; lineFreqs = [lnParams.freq, lnParams.freq*2] ;
    end
    
    % LEGACY LINE NOISE: If the user indicated to use the legacy method of line
    % noise reduction, apply the code from HAPPE v1 to the data.
    if lnParams.legacy
        outEEG = pop_cleanline(inEEG, 'bandwidth', 2, 'chanlist', ...
            1:inEEG.nbchan, 'computepower', 1, 'linefreqs', lineFreqs, ...
            'normSpectrum', 0, 'p', 0.01, 'pad', 2, 'plotfigures', 0, ...
            'scanforlines', 1, 'sigtype', 'Channels', 'tau', 100, 'verb', 0, ...
            'winsize', 4, 'winstep', 1, 'ComputeSpectralPower', 'False');
    
    % DEFAULT LINE NOISE REDUCTION: If the user indicated to use the default
    % method of line noise reduction, apply the optimized method to the data.
    % This is the updated version of CleanLine (Mullen, 2012).
    else
        [outEEG, ~] = cleanLineNoise(inEEG, struct('lineNoiseMethod', 'clean', ...
            'lineNoiseChannels', 1:inEEG.nbchan, 'Fs', inEEG.srate, ...
            'lineFrequencies', lineFreqs, 'p', 0.01, 'fScanBandWidth', 2, ...
            'taperBandWidth', 2, 'taperWindowSize', 4, 'taperWindowStep', 4, ...
            'tau', 100, 'pad', 2, 'fPassBand', [0 inEEG.srate/2], ...
            'maximumIterations', 10)) ;
    end
else
    outEEG = pop_eegfiltnew(inEEG, 'locutoff', lnParams.low, 'hicutoff', ...
        lnParams.high, 'revfilt', 1, 'plotfreqz', 0);
end

% LINE NOISE REDUCTION QM: Assesses the performance of line noise reduction.
lnMeans = assessPipelineStep('line noise reduction', reshape(inEEG.data, ...
    size(inEEG.data, 1), []), reshape(outEEG.data, size(outEEG.data,1), ...
    []), lnMeans, inEEG.srate, [lnParams.neighbors lnParams.harms.freqs]) ;
end