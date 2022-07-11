% std_interp() - interpolate, if needed, a list of named data channels 
%                for all datasets included in a STUDY. Currently assumes
%                that all channels have uniform locations across datasets.
%
% Usage:      >> [STUDY ALLEEG] = std_interp(STUDY, ALLEEG, chans, method);
%
% Inputs: 
%     STUDY    - EEGLAB STUDY structure
%     ALLEEG   - EEGLAB vector containing all EEG dataset structures in the STUDY.
%     chans    - [Cell array] cell array of channel names (labels) to interpolate 
%                into the data if they are missing from one of the datasets.
%     method   - [string] griddata() method to use for interpolation.
%                See  >> help eeg_interp() {default:'spherical'}
%
% Important limitation:
%     This function currently presuposes that all datasets have the same channel 
%     locations (with various channels from a standard set possibly missing). 
%     If this is not the case, the interpolation will not be performed.
%
% Output: 
%     STUDY    - study structure.
%     ALLEEG   - updated datasets.
%
% Author: Arnaud Delorme, CERCO, CNRS, August 2006-
%
% See also: eeg_interp()

% Copyright (C) Arnaud Delorme, CERCO, 2006, arno@salk.edu
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

% Deprecated:
%              - [chanlocs structure] channel location structure containing
%                a full channel structure (missing channels in the current 
%                dataset are interpolated).

function [STUDY, ALLEEG] = std_interp(STUDY, ALLEEG, chans, method);

if nargin < 2
    help std_interp;
    return;
end
if nargin < 3
    chans = [];
end
if nargin < 4
    method = 'spherical';
end

% union of all channel structures
% -------------------------------
alllocs = eeg_mergelocs(ALLEEG.chanlocs);

% check electrode names to interpolate
% ------------------------------------
if iscell(chans)
    alllabs = lower({ alllocs.labels });
    for index = 1:length(chans)
        tmpind = strmatch(lower(chans{index}), alllabs, 'exact');
        if isempty(tmpind)
            error( sprintf('Channel named ''%s'' not found in any dataset', chans{index}));
        end
    end
end

% read all STUDY datasets and interpolate electrodes
% ---------------------------------------------------
for index = 1:length(STUDY.datasetinfo)
   tmpind  = STUDY.datasetinfo(index).index;
   tmplocs = ALLEEG(tmpind).chanlocs;

   % build electrode location structure for interpolation
   % ----------------------------------------------------
   [tmp tmp2 id1] = intersect_bc({tmplocs.labels}, {alllocs.labels});
   if isempty(chans)
       interplocs = alllocs;
   elseif iscell(chans)
       [tmp tmp2 id2] = intersect_bc( chans, {alllocs.labels});
       interplocs = alllocs(union(id1, id2));
   else
       interplocs = chans;
   end
   
   if length(interplocs) ~= length(tmplocs)
       
       % search for position of electrode in backup structure
       % ----------------------------------------------
       extrachans = [];
       if isfield(ALLEEG(tmpind).chaninfo, 'nodatchans')
            if isfield(ALLEEG(tmpind).chaninfo.nodatchans, 'labels')
                extrachans = ALLEEG(tmpind).chaninfo.nodatchans;
            end
       end
       tmplabels = { tmplocs.labels };
       for i=1:length(interplocs)
           ind = strmatch( interplocs(i).labels, tmplabels, 'exact');
           if ~isempty(ind)
               interplocs(i) = tmplocs(ind); % this is necessary for polhemus
           elseif ~isempty(extrachans)
               ind = strmatch( interplocs(i).labels, { extrachans.labels }, 'exact');
               if ~isempty(ind)
                    fprintf('Found position of %s in chaninfo structure\n', interplocs(i).labels);
                    interplocs(i) = extrachans(ind);
               end
           end
       end
       
        % perform interpolation
        % ---------------------
        EEG = eeg_retrieve(ALLEEG, index);
        EEG = eeg_checkset(EEG);
        EEG = eeg_interp(EEG, interplocs, method);
        EEG.saved = 'no';
        EEG = pop_saveset(EEG, 'savemode', 'resave');
        
        % update dataset in EEGLAB
        % ------------------------
        if ischar(ALLEEG(tmpind).data)
            tmpdata = ALLEEG(tmpind).data;
            [ ALLEEG EEG ] = eeg_store(ALLEEG, EEG, tmpind);
            ALLEEG(tmpind).data  = tmpdata;
            ALLEEG(tmpind).saved = 'yes';
            clear EEG;
        else
            [ ALLEEG EEG ] = eeg_store(ALLEEG, EEG, tmpind);
            ALLEEG(tmpind).saved = 'yes';
        end
    else
        fprintf('No need for interpolation for dataset %d\n', tmpind);
    end
end
           
       
function checkchans(STUDY, ALLEEG)

    % Check to see if all the channels have the same coordinates
    % (check only the theta field)
    % ----------------------------------------------------------
    for index = 1:length(STUDY.datasetinfo)
       tmpind  = STUDY.datasetinfo(index).index;
       tmplocs = ALLEEG(tmpind).chanlocs;
       [tmp id1 id2] = intersect_bc({tmplocs.labels}, {alllocs.labels});
       for ind = 1:length(id1)
           if tmplocs(id1(ind)).theta ~= alllocs(id2(ind)).theta

               % find datasets with different coordinates
               % ----------------------------------------
               for ind2 = 1:length(STUDY.datasetinfo)
                   tmplocs2 = ALLEEG(ind2).chanlocs;
                   tmpmatch = strmatch(alllocs(id2(ind)).labels, { tmplocs2.labels }, 'exact');
                   if ~isempty(tmpmatch) 
                       if alllocs(id2(ind)).theta == tmplocs2(tmpmatch).theta
                           datind = ind2;
                           break;
                       end
                   end
               end

               error(sprintf( [ 'Dataset %d and %d do not have the same channel location\n' ...
                   'for electrode ''%s''' ], datind, tmpind, tmplocs(id1(ind)).labels));
           end
       end
    end

