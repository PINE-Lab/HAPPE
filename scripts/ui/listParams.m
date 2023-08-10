% listParams() - A helper function for setParams.m and HAPPE that prints
%                out the user-entered parameters to the command window in a
%                way that is easy to read. Allows the user to check their
%                parameters prior to confirming them.
%
% Usage: 
%   >> listParams(params)
%
% Inputs:
%   params - A struct containing all of the parameters needed to run HAPPE.
%
% Outputs:
%
% NOTE: Changes to the way any parameters are encoded, or changes in the
% number of parameters (adding/subtracting) will result in this code
% needing to be updated to reflect those changes. There is no easy way to
% make it reliant on other scripts... Sorry :(
%
% Author: A.D. Monachino, PINE Lab at Northeastern University, 2022
%
% This file is part of HAPPE.
% Copyright 2018, 2021 Alexa Monachino, Kelsie Lopez, Laurel Gabard-Durnam
%
% HAPPE is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% HAPPE is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with HAPPE. If not, see <https://www.gnu.org/licenses/>.

function listParams(params)
disp("---------------------------------------------") ;
disp("PARAMETER SETTINGS:") ;
%% Density
fprintf('Density: ') ;
if params.lowDensity; fprintf('Low (<= 30 channels)\n') ;
else; fprintf('High (> 30 channels)\n') ;
end

%% Rest VS Task
fprintf('Resting State or Task: ') ;
if params.paradigm.task
    fprintf('Task\n - Task Onset Tags: ') ;
    fprintf([sprintf('%s, ', params.paradigm.onsetTags{1:end-1}), ...
        params.paradigm.onsetTags{end} '\n']) ;
    fprintf('    - Conditions: ') ;
    if params.paradigm.conds.on
        for i=1:size(params.paradigm.conds.groups,1)
            currGroup = params.paradigm.conds.groups(i,2:end) ;
            currGroup = currGroup(~cellfun('isempty',currGroup)) ;
            fprintf([params.paradigm.conds.groups{i,1} ' = ' ...
                sprintf('%s, ', currGroup{1:end-1}) currGroup{end} '\n']);
            if i~=size(params.paradigm.conds.groups,1)
                fprintf('                  ') ;
            end
        end
    else; fprintf('NA\n') ;
    end
    fprintf(' - ERP Analysis: ') ;
    if params.paradigm.ERP.on; fprintf('Yes\n') ;
    else; fprintf('No\n') ;
    end
else; fprintf('Resting State\n') ;
end

%% Data File Format
fprintf('Data File Format: ') ;
if params.loadInfo.inputFormat == 1
    fprintf('.mat (Net Station output or MATLAB matrix)\n') ;
    if params.loadInfo.NSformat
        fprintf(' - Potential EEG Variable Names: ') ;
        fprintf([sprintf('%s, ', params.loadInfo.NSvarNames{1:end-1}), ...
            params.loadInfo.NSvarNames{end} '\n']) ;
    else
        fprintf(' - Channel Locations Provided: ') ;
        if params.loadInfo.chanlocs.inc
            fprintf('Yes\n    - Channel Locations File: %s\n', ...
                params.loadInfo.chanlocs.file) ;
        else; fprintf('No\n') ;
        end
        
        fprintf(' - Same Sampling Rate Across Files: ');
        if params.loadInfo.srate.same; fprintf('Yes\n') ;
        else
            fprintf('No\n    - List of Sampling Rates: %s\n', ...
                params.loadInfo.srate.file) ;
        end
    end
    if params.paradigm.task
        fprintf(' - Task Event Info: %s\n', params.loadInfo.eventLoc) ;
    end
elseif params.loadInfo.inputFormat == 2; fprintf('.raw (Netstation simple binary)\n') ;
elseif params.loadInfo.inputFormat == 3; fprintf('.set (EEGLab)\n') ;
elseif params.loadInfo.inputFormat == 4; fprintf('.cdt (Neuroscan)\n') ;
elseif params.loadInfo.inputFormat == 5
    fprintf('.mff (EGI)\n - Type Fields: ') ;
    if length(params.loadInfo.typeFields) > 1
        fprintf([sprintf('%s, ', params.loadInfo.typeFields{1:end-1}), ...
            params.loadInfo.typeFields{end} '\n']) ;
    else; fprintf([params.loadInfo.typeFields{1} '\n']) ;
    end
elseif params.loadInfo.inputFormat == 6; fprintf('.edf\n') ;
elseif params.loadInfo.inputFormat == 7
    fprintf('.bdf->.set (Mentalab)\n') ;
end

%% Aquisition Layout
fprintf('Acquisition Layout: ') ;
if params.loadInfo.layout(1) == 1; fprintf(['%i channel EGI Geodesic Sensor ' ...
        'Net\n'], params.loadInfo.layout(2)) ;
elseif params.loadInfo.layout(1) == 2; fprintf(['%i channel EGI HydroCel ' ...
        'Geodesic Sensor Net\n'], params.loadInfo.layout(2)) ;
elseif params.loadInfo.layout(1) == 3; fprintf('%i channel Neuroscan Quik-Cap\n', ...
        params.loadInfo.layout(2)) ;
elseif params.loadInfo.layout(1) == 4; fprintf('%i channel Mentalab Explore\n', ...
        params.loadInfo.layout(2)) ;
elseif params.loadInfo.layout(1) == 5; fprintf('Unspecified\n') ;
end

%% Channels of Interest
fprintf('Channels: ') ;
if strcmpi(params.chans.subset, 'all'); fprintf('All\n') ;
elseif strcmpi(params.chans.subset, 'coi_include')
    fprintf([sprintf('%s, ', params.chans.IDs{1:end-1}), ...
        params.chans.IDs{end} '\n']) ;
elseif strcmpi(params.chans.subset, 'coi_exclude')
    fprintf('All except ') ;
    fprintf([sprintf('%s, ', params.chans.IDs{1:end-1}), ...
        params.chans.IDs{end} '\n']) ;    
end

%% Line Noise
if isfield(params, 'lineNoise')
    % FREQUENCY
    fprintf('Line Noise Frequency: %i Hz\n', params.lineNoise.freq) ;
    if params.lineNoise.harms.on
        fprintf([' - Additional Frequencies to Reduce: ' sprintf('%i, ', ...
            params.lineNoise.harms.freqs(1:end-1)) ...
            num2str(params.lineNoise.harms.freqs(end)) '\n']) ;
    end
    % LEGACY
    fprintf('Line Noise Reduction Method: ') ;
    if params.lineNoise.cl
        fprintf('CleanLine - ') ;
        if params.lineNoise.legacy; fprintf('Legacy\n') ;
        else; fprintf('Default\n') ;
        end
    else
        fprintf(['Notch Filter\n - Low Cutoff: ' num2str(params.lineNoise.low) ...
            '\n - High Cutoff: ' num2str(params.lineNoise.high) '\n']) ;
    end
end

%% Resampling
if isfield(params, 'downsample')
    fprintf('Resample: ')
    if params.downsample == 0; fprintf('Off\n') ;
    else; fprintf(sprintf('To %i Hz\n', params.downsample)) ;
    end
end

%% Filtering
fprintf('Filter: ') ;
if params.filt.on
    fprintf([' - Lowpass Cutoff: ' num2str(params.filt.lowpass) '\n' ...
        ' - Highpass Cutoff: ' num2str(params.filt.highpass) '\n' ...
        ' - Type: ']) ;
    if params.filt.butter; fprintf('ERPLAB''s bandpass Butterworth\n') ;
    else; fprintf('EEGLAB''s FIR\n') ;
    end
else; fprintf('Off\n') ;
end

%% Bad Channel Detection
if isfield(params, 'badChans')
    fprintf('Bad Channel Detection: ') ;
    if params.badChans.rej
        fprintf('On\n - Bad Channel Detection Order: ') ;
        if params.badChans.order; fprintf('After Wavelet Thresholding\n') ;
        else; fprintf('Before Wavelet Thresholding\n') ;
        end
    else; fprintf('Off\n') ;
    end  
end

%% ECGone
if isfield(params, 'ecgone')
    fprintf('ECGone: ') ;
    if params.ecgone.on
        fprintf('On\n - ECG Channel: ') ;
        if params.ecgone.ECGchan.inc
            fprintf([params.ecgone.ECGchan.ID '\n']) ;
        else
            fprintf(['Proxy made from artifact in channel(s) ' ...
                sprintf('%s, ', params.ecgone.ECGchan.ID{1:end-1}), ...
                params.ecgone.ECGchan.ID{end} '\n - Artifact Threshold:' ...
                num2str(params.ecgone.peaky) '\n']) ;
        end
        fprintf([' - Peak Detection Window Length: ' ...
            num2str(params.ecgone.peakWinSize) ' seconds\n - Template ' ...
            'Length: ' num2str(params.ecgone.procEpoch) ' seconds\n']) ;
    else; fprintf('Off\n') ;
    end
end

%% Legacy Wavelet
if isfield(params, 'wavelet')
    fprintf('Wavelet Thresholding: ') ;
    if params.wavelet.legacy; fprintf('Legacy\n') ;
    else; fprintf('Default\n') ;
        if params.paradigm.ERP.on
            fprintf(' - Threshold Rule: ') ;
            if params.wavelet.softThresh; fprintf('Soft\n') ;
            else; fprintf('Hard\n') ;
            end
        end
    end
end

%% MUSCIL
fprintf('MuscIL: ') ;
if params.muscIL; fprintf('On\n') ; else; fprintf('Off\n') ; end

%% Segmentation
fprintf('Segmentation: ') ;
if params.segment.on
    fprintf('On\n') ;
    if params.paradigm.task
        fprintf(' - Starting Parameter for Stimulus: %g seconds\n', params.segment.start) ;
        fprintf(' - Ending Parameter for Stimulus: %g seconds\n', params.segment.end) ;
        if params.paradigm.ERP.on
            fprintf(' - Task Offset: %g milliseconds\n - Baseline Correction: ', params.segment.offset) ;
            if params.baseCorr.on
                fprintf('On\n    - Baseline Correction Start: %g milliseconds\n', params.baseCorr.start) ;
                fprintf('    - Baseline Correction End: %g milliseconds\n', params.baseCorr.end) ;
            else
                fprintf('Off\n') ;
            end
        end
    else; fprintf(' - Segment Length: %g seconds\n', params.segment.length) ;
    end
end

%% Interpolation
fprintf('Interpolation: ') ;
if params.segment.interp; fprintf('On\n') ;
else; fprintf('Off\n') ;
end

%% Segment Rejection
fprintf('Segment Rejection: ') ;
if params.segRej.on
    fprintf('On\n - Segment Rejection Method: ') ;
    if strcmpi(params.segRej.method, 'both'); fprintf('Both amplitude and similarity criteria\n') ;
    elseif strcmpi(params.segRej.method, 'amplitude'); fprintf('Amplitude criteria only\n');
    elseif strcmpi(params.segRej.method, 'similarity'); fprintf('Similarity criteria only\n') ;
    end
    if strcmpi(params.segRej.method, 'both') || ...
            strcmpi(params.segRej.method, 'amplitude')
        fprintf('    - Minimum Segment Rejection Threshold: %i\n', ...
            params.segRej.minAmp) ;
        fprintf('    - Maximum Segment Rejection Threshold: %i\n', ...
            params.segRej.maxAmp) ;
        fprintf('    - Segment Rejection based on All Channels or ROI: ') ;
        if params.segRej.ROI.on
            fprintf('ROI\n     - ROI Channels: ') ;
            if ~params.segRej.ROI.include; fprintf('All except ') ; end
            fprintf([sprintf('%s, ', params.segRej.ROI.chans{1:end-1}), ...
                params.segRej.ROI.chans{end} '\n']) ;
        else
            fprintf('All\n') ;
        end
        if params.paradigm.task && params.loadInfo.inputFormat == 1
            fprintf('Restrict Analysis Using Pre-Selected Trials: ') ;
            if params.segRej.selTrials; fprintf('On\n') ;
            else; fprintf('Off\n') ;
            end
        end
    end
else; fprintf('Off\n') ;
end

%% Re-Reference
fprintf('Re-Referencing: ') ;
if params.reref.on
    fprintf(['On\n - Re-Reference Method: ' params.reref.method '\n']) ;
    if strcmpi(params.reref.method, 'subset')
        if size(params.reref.subset,2) > 1
            fprintf([' - Subset: ' sprintf('%s, ', params.reref.subset{1:end-1}) ...
                params.reref.subset{end} '\n']) ;
        else; fprintf([' - Subset: ' params.reref.subset{end} '\n']) ;
        end
    elseif strcmpi(params.reref.method, 'rest')
        fprintf([' - Reference Type: ' params.reref.rest.ref '\n']) ;
        if strcmpi(params.reref.rest.ref, 'mastoid')
            fprintf(['    - Left Channels: ' sprintf('%s, ', params.reref.rest.chanL{1:end-1}) ...
                params.reref.rest.chanL{end} '\n    - Right Channels: ' ...
                sprintf('%s, ', params.reref.rest.chanR{1:end-1}) ...
                params.reref.rest.chanR{end} '\n']) ;
        end
        fprintf(' - Load or calculate leadfield: ');
        if params.reref.rest.calc
            fprintf('calculate\n') ;
            if strcmpi(params.reref.rest.ref, 'mastoid')
                fprintf(['    - Left Reference Coordinates: ' ...
                    num2str(params.reref.rest.coordLeft(1)) ' ' ...
                    num2str(params.reref.rest.coordLeft(2)) ' ' ...
                    num2str(params.reref.rest.coordLeft(3)) '\n    - Right ' ...
                    'Reference Coordinates: ' ...
                    num2str(params.reref.rest.coordRight(1)) ' ' ...
                    num2str(params.reref.rest.coordRight(2)) ' ' ...
                    num2str(params.reref.rest.coordRight(3)) '\n']) ;
            end
        else
            fprintf(['load\n    - Leadfield File: ' params.reref.rest.file '\n']) ;
        end
    end
else; fprintf('Off\n') ;
end

%% Visualizations
fprintf('Visualizations: ') ;
if params.vis.enabled
    fprintf('On\n') ;
    if params.paradigm.ERP.on
        fprintf(' - Start Time: %i milliseconds\n', params.vis.min) ;
        fprintf(' - End Time: %i milliseconds\n', params.vis.max) ;
        fprintf(' - Times to Plot: ');
        for i=1:length(params.vis.toPlot)
            if i == length(params.vis.toPlot)
                fprintf('%i\n', params.vis.toPlot(i)) ;
            else
                fprintf('%i, ', params.vis.toPlot(i)) ;
            end
        end
    else
        fprintf(' - Power Spectrum Minimum: %i\n', params.vis.min) ;
        fprintf(' - Power Spectrum Maximum: %i\n', params.vis.max) ;
        fprintf(' - Frequencies to Plot: ');
        for i=1:length(params.vis.toPlot)
            if i == length(params.vis.toPlot)
                fprintf('%i\n', params.vis.toPlot(i)) ;
            else
                fprintf('%i, ', params.vis.toPlot(i)) ;
            end
        end
    end
else
    fprintf('Off\n') ;
end

%% Save Format
fprintf('Save Format: ') ;
if params.outputFormat == 1; fprintf('.txt file\n') ;
elseif params.outputFormat == 2; fprintf('.mat file\n') ;
elseif params.outputFormat == 3; fprintf('.set file\n') ;
end
disp("---------------------------------------------") ;    
end