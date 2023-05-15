% happe_segRej() - A helper script for HAPPE that rejects segments from the
%                  EEG data that meet certain criteria as specified by the
%                  user. Can use amplitude and/or similarity criteria.
%                  Checks for common errors such as rejecting all segments
%                  or attempting to reject data with a single epoch.
%
% Usage: 
%   >> EEG = happe_segRej(EEG, params)
%
% Inputs:
%   EEG    - The EEG, in EEGLAB (Delorme & Makeig, 2004) format, to
%            segment.
%   params - The struct containing all the relevent information for 
%            segment rejection, including rejection criteria.
%
% Outputs:
%   EEG    - The EEG, in EEGLAB (Delorme & Makeig, 2004) format, after 
%            segment rejection.
%
% Author: L.J. Gabard-Durnam, PINE Lab at Northeastern University, 2021
%         A.D. Monachino, PINE Lab at Northeastern University, 2022
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

function [EEG, keptTrials] = happe_segRej(EEG, params) 
fprintf('Rejecting segments...\n') ;
% GATHER ROI INDXS: If rejecting based on a ROI, collect the indxs for the
% channels in the ROI.
if params.segRej.ROI.on
    if params.segRej.ROI.include
        params.segRej.ROI.chans = intersect({EEG.chanlocs.labels}, ...
            params.segRej.ROI.chans) ;
    else
        params.segRej.ROI.chans = setdiff({EEG.chanlocs.labels}, ...
            params.segRej.ROI.chans) ;
    end
    ROI_indxs = [] ;
    for i=1:size(params.segRej.ROI.chans, 2)
        ROI_indxs = [ROI_indxs find(strcmpi({EEG.chanlocs.labels}, ...
            params.segRej.ROI.chans{i}))] ;                                 %#ok<*AGROW> 
    end
end

% MARK FOR REJECTION USING AMPLITUDE CRITERIA: Mark segments for rejection 
% that have an amplitude change greater than the user-defined thresholds.
if ismember(params.segRej.method, {'amplitude', 'both'})
    if params.segRej.ROI.on; EEG = pop_eegthresh(EEG, 1, ROI_indxs, ...
            [params.segRej.minAmp], [params.segRej.maxAmp], [EEG.xmin], ...
            [EEG.xmax], 2, 0) ;
    else; EEG = pop_eegthresh(EEG, 1, 1:EEG.nbchan, ...
            [params.segRej.minAmp], [params.segRej.maxAmp], [EEG.xmin], ...
            [EEG.xmax], 2, 0) ;
    end
end

% MARK FOR REJECTION USING SIMILARITY CRITERIA: Mark segments
% for rejection that don't meet the joint-probability criteria.
% For more on this criteria, consult the HAPPE manuscripts.
if ismember(params.segRej.method, {'similarity', 'both'})
    if params.lowDensity; num = 2; else; num = 3; end
    
    % If the flatline channel is included in the data, temporarily remove
    % it to properly assess joint-probability.
    if params.loadInfo.chanlocs.inc && params.reref.on && ...
        params.reref.flat
        EEGtemp = pop_select(EEG, 'channel', ...
            setdiff({EEG.chanlocs.labels}, params.reref.chan)) ;
    else
        zeroChanIndxs = [] ;
        for i=1:size(EEG.data,1)
            if all(all(EEG.data(i,:,:) == 0)); zeroChanIndxs = [zeroChanIndxs i]; end
        end
        if isempty(zeroChanIndxs); EEGtemp = EEG ;
        else
            chans = {EEG.chanlocs.labels} ;
            EEGtemp = pop_select(EEG, 'channel', chans(setdiff(1:size(chans,2), ...
                zeroChanIndxs))) ;
        end
    end
    
    % Assess the joint-probability depending on whether assessing a ROI or
    % all channels.
    if params.segRej.ROI.on; EEGtemp = pop_jointprob(EEGtemp, 1, ...
            ROI_indxs, num, num, params.vis.enabled, 0, ...
            params.vis.enabled, [], params.vis.enabled) ;
    else
        EEGtemp = pop_jointprob(EEGtemp, 1, 1:EEGtemp.nbchan, num, ...
            num, params.vis.enabled, 0, params.vis.enabled, [], ...
            params.vis.enabled) ;
    end
    
    % Apply rejection information to the EEG, adding the flatline channel
    % into the matrix if necessary.
    if params.loadInfo.chanlocs.inc && params.reref.on && ...
        params.reref.flat
        EEG.reject.rejjpE = [EEGtemp.reject.rejjpE; zeros(1, ...
            size(EEGtemp.reject.rejjpE,2))] ;
    else
        if isempty(zeroChanIndxs)
            EEG.reject.rejjpE = EEGtemp.reject.rejjpE ;
        else
            rejjpE_temp = [] ;
            indx = 1;
            for i=1:EEG.nbchan
                if ismember(i, zeroChanIndxs)
                    rejjpE_temp = [rejjpE_temp; zeros(1, ...
                        size(EEGtemp.reject.rejjpE,2))] ;
                else
                    rejjpE_temp = [rejjpE_temp; EEGtemp.reject.rejjpE(indx,:)] ;
                    indx = indx+1 ;
                end
            end
            EEGtemp.reject.rejjpE = rejjpE_temp ;
        end
    end
    EEG.reject.rejjp = EEGtemp.reject.rejjp ;
end

% SELECT SEGMENTS FOR REJECTION
if params.paradigm.task && params.loadInfo.inputFormat == 1 && ...
        params.segRej.selTrials
    EEG = pop_selectevent(EEG, 'status', 'good', 'deleteevents', ...
        'on', 'deleteepochs', 'on', 'invertepochs', 'off') ;
end
EEG = eeg_rejsuperpose(EEG, 1, 0, 1, 1, 1, 1, 1, 1) ;
if isfield(EEG, 'reject'); keptTrials = find(EEG.reject.rejglobal==0); end
if isfield(EEG, 'rej'); keptTrials = find(EEG.rej.rejglobal==0); end

% REJECT TRIALS: Prior to rejection, check to see if
% all segments have been rejected, and if so, alert the user
% that this has occured then throw an error to proceed to the
% next file. If not all segments have been rejected, remove the
% marked segments from the data.
if (isfield(EEG, 'reject') && all(EEG.reject.rejglobal)) || ...
        (isfield(EEG, 'rej') && all(EEG.rej.rejglobal))
    fprintf('All trials rejected. No further processing possible.\n') ;
    error('HAPPE:AllTrialsRej', 'All trials rejected') ;
else; EEG = pop_rejepoch(EEG, [EEG.reject.rejglobal], 0) ;
end
end