% EEG_CHECKCHANLOCS - Check the consistency of the channel locations structure 
%                  of an EEGLAB dataset.
%
% Usage:
%  >> EEG = eeg_checkchanlocs(EEG); 
%  >> [chanlocs chaninfo] = eeg_checkchanlocs( chanlocs, chaninfo);
%
% Inputs:
%   EEG      - EEG dataset
%   chanlocs - EEG.chanlocs structure
%   chaninfo - EEG.chaninfo structure
%
% Outputs:
%   EEG        - new EEGLAB dataset with updated channel location structures 
%                EEG.chanlocs, EEG.urchanlocs, EEG.chaninfo
%   chanlocs   - updated channel location structure
%   chaninfo   - updated chaninfo structure
%
% Author: Arnaud Delorme, SCCN/INC/UCSD, March 2, 2011

% Copyright (C) SCCN/INC/UCSD, March 2, 2011, arno@salk.edu
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

% Hey Arno -- this is a quick fix to make an analysis work for Makoto
% I think the old version had a bug...

function [chans, chaninfo, chanedit]= eeg_checkchanlocs(chans, chaninfo);

if nargin < 1 
    help eeg_checkchanlocs;
    return;
end

if nargin < 2
    chaninfo = [];
end

processingEEGstruct = 0;
if isfield(chans, 'data')
    processingEEGstruct = 1;
    tmpEEG = chans;
    [chans, chaninfo] = insertchans(tmpEEG.chanlocs, tmpEEG.chaninfo);

    % force Nosedir to +X (done here because of DIPFIT)
    % -------------------
    if isfield(tmpEEG.chaninfo, 'nosedir')
        if ~strcmpi(tmpEEG.chaninfo.nosedir, '+x') && all(isfield(tmpEEG.chanlocs,{'X','Y','theta','sph_theta'}))
            disp(['Note for expert users: Nose direction is now set from ''' upper(tmpEEG.chaninfo.nosedir)  ''' to default +X in EEG.chanlocs']);
            [~, chaninfo, chans] = eeg_checkchanlocs(tmpEEG.chanlocs, tmpEEG.chaninfo); % Merge all channels for rotation (FID and data channels)
            if strcmpi(chaninfo.nosedir, '+y')
                rotate = 270;
            elseif strcmpi(chaninfo.nosedir, '-x')
                rotate = 180;
            else
                rotate = 90;
            end
            for index = 1:length(chans)
                rotategrad = rotate/180*pi;
                coord = (chans(index).Y + chans(index).X*sqrt(-1))*exp(sqrt(-1)*-rotategrad);
                chans(index).Y = real(coord);
                chans(index).X = imag(coord);

                if ~isempty(chans(index).theta)
                    chans(index).theta     = chans(index).theta    -rotate;
                    chans(index).sph_theta = chans(index).sph_theta+rotate;
                    if chans(index).theta    <-180, chans(index).theta    =chans(index).theta    +360; end
                    if chans(index).sph_theta>180 , chans(index).sph_theta=chans(index).sph_theta-360; end
                end
            end

            if isfield(tmpEEG, 'dipfit')
                if isfield(tmpEEG.dipfit, 'coord_transform')
                    if isempty(tmpEEG.dipfit.coord_transform)
                        tmpEEG.dipfit.coord_transform = [0 0 0 0 0 0 1 1 1];
                    end
                    tmpEEG.dipfit.coord_transform(6) = tmpEEG.dipfit.coord_transform(6)+rotategrad;
                end
            end

            chaninfo.originalnosedir = chaninfo.nosedir;
            chaninfo.nosedir = '+X';
        end
    end
    chanedit = chans;
    complicated = true;
else
    if ~isfield(chans, 'datachan')
        [chanedit,dummy,complicated] = insertchans(chans, chaninfo);
    else
        chanedit = chans;
        complicated = true;
    end
end

nosevals       = { '+X' '-X' '+Y' '-Y' };
if ~isfield(chaninfo, 'plotrad'), chaninfo.plotrad = []; end
if ~isfield(chaninfo, 'shrink'),  chaninfo.shrink = [];  end
if ~isfield(chaninfo, 'nosedir'), chaninfo.nosedir = nosevals{1}; end

% handles deprecated fields
% -------------------------
plotrad  = [];
if isfield(chanedit, 'plotrad'),
    plotrad = chanedit(1).plotrad;
    chanedit = rmfield(chanedit, 'plotrad');
    if ischar(plotrad) && ~isempty(str2num(plotrad)), plotrad = str2num(plotrad); end
    chaninfo.plotrad = plotrad;
end
if isfield(chanedit, 'shrink') && ~isempty(chanedit(1).shrink)
    shrinkorskirt = 1;
    if ~ischar(chanedit(1).shrink)
        plotrad = 0.5/(1-chanedit(1).shrink); % convert old values
    end
    chanedit = rmfield(chanedit, 'shrink');
    chaninfo.plotrad = plotrad;
end

% set non-existent fields to []
% -----------------------------
fields    = { 'labels' 'theta' 'radius' 'X'   'Y'   'Z'   'sph_theta' 'sph_phi' 'sph_radius' 'type' 'ref' 'urchan' };
fieldtype = { 'str'    'num'   'num'    'num' 'num' 'num' 'num'       'num'     'num'        'str'  'str' 'num'    };
check_newfields = true; %length(fieldnames(chanedit)) < length(fields);
if ~isempty(chanedit)
    for index = 1:length(fields)
        if check_newfields && ~isfield(chanedit, fields{index})
            % new field
            % ---------
            if strcmpi(fieldtype{index}, 'num')
                chanedit = setfield(chanedit, {1}, fields{index}, []);
            else
                for indchan = 1:length(chanedit)
                    chanedit = setfield(chanedit, {indchan}, fields{index}, '');
                end
            end
        else
            % existing fields
            % ---------------
            allvals = {chanedit.(fields{index})};
            if isnumeric(allvals{1}) && any(cellfun(@(x)~isempty(x) & all(isnan(x)), allvals))
                posNaNs = find(cellfun(@(x)~isempty(x) & all(isnan(x)), allvals));
                for iPos = 1:length(posNaNs)
                    chanedit = setfield(chanedit, {posNaNs(iPos)}, fields{index}, []);
                end
            end
            if strcmpi(fieldtype{index}, 'num')
                if ~all(cellfun('isclass',allvals,'double'))
                    numok = cellfun(@isnumeric, allvals);
                    if any(numok == 0)
                        for indConvert = find(numok == 0)
                            chanedit = setfield(chanedit, {indConvert}, fields{index}, []);
                        end
                    end
                end
            else
                strok = cellfun('isclass', allvals,'char');
                if strcmpi(fields{index}, 'labels'), prefix = 'E'; else prefix = ''; end
                if any(strok == 0)
                    for indConvert = find(strok == 0)
                        try
                              strval   = [ prefix num2str(getfield(chanedit, {indConvert}, fields{index})) ];
                              chanedit = setfield(chanedit, {indConvert}, fields{index}, strval);
                        catch
                            chanedit = setfield(chanedit, {indConvert}, fields{index}, '');
                        end
                    end
                end
            end
        end
    end
end

if ~isequal(fieldnames(chanedit)',fields)
    try
        chanedit = orderfields(chanedit, fields);
    catch, end
end


% check channel labels
if isfield(chanedit, 'labels')
    % prefix (EDF format specification)?
    if strfind([chanedit.labels], 'EEG') % `contains() is not back compatible
        chanprefixes = {'EEG-', 'EEG ', 'EEG'}; % order matters
        tmp = {chanedit.labels};
        if sum(~isnan(str2double( strrep(tmp, 'EEG', '')))) < 30 % more than 30 numerical channels, i.e., EEG001, do nothing
            disp('Detected/removing ''EEG'' prefix from channel labels')
            for idx = 1:length(chanprefixes)
                tmp = strrep(tmp, chanprefixes(idx), '');
            end
            [chanedit.labels] = deal(tmp{:});
        end
    end
        
    % duplicate labels?
    tmp = sort({chanedit.labels});
    if any(strcmp(tmp(1:end-1),tmp(2:end)))
        disp('Warning: some channels have the same label'); 
    end
    
    % empty labels?
     indEmpty = find(cellfun(@isempty, {chanedit.labels}));
    if ~isempty(indEmpty)
        tmpWarning = warning('backtrace'); 
        warning backtrace off;
        warning('channel labels should not be empty, creating unique labels');
        warning(tmpWarning); 
        for index = indEmpty
            chanedit(index).labels = sprintf('E%d', index);
        end
    end   

    % handle MEG
    if ~isfield(chaninfo, 'topoplot')
        if ~isempty(strfind([chanedit(1).labels], 'MLC11')) || (isfield(chanedit, 'type') && ~isempty(strfind(chanedit(1).type, 'meg')))
            disp('MEG data detected and topoplot options not set, so setting them in EEG.chaninfo')
            chaninfo.topoplot = { 'conv' 'on' 'headrad' 0.3 };
        end
    end
end

% remove fields
% -------------
if isfield(chanedit, 'sph_phi_besa'  ), chanedit = rmfield(chanedit, 'sph_phi_besa'); end
if isfield(chanedit, 'sph_theta_besa'), chanedit = rmfield(chanedit, 'sph_theta_besa'); end

% Check if some channels need conversion
% --------------------------------------
chanX        = cellfun('isempty',{ chanedit.X });
chanTheta    = cellfun('isempty',{ chanedit.theta });
chanSphTheta = cellfun('isempty',{ chanedit.sph_theta });
if any(~chanX & chanTheta) || any(~chanSphTheta & chanTheta) || any(~chanX & chanSphTheta)
    try
        % convert them all
        chanedit = convertlocs(chanedit,'auto');
    catch
        disp('eeg_checkchanlocs: Unable to convert electrode locations between coordinate systems');
    end
end

% reconstruct the chans structure
% -------------------------------
if complicated
    [chans, chaninfo.nodatchans] = getnodatchan( chanedit );
    if ~isfield(chaninfo, 'nodatchans'), chaninfo.nodatchans = []; end
    if isempty(chanedit)
        for iField = 1:length(fields)
            chanedit = setfield(chanedit, fields{iField}, []);
        end
    end
else
    chans = rmfield(chanedit,'datachan');
    chaninfo.nodatchans = [];
end

if processingEEGstruct
    tmpEEG.chanlocs = chans;
    tmpEEG.chaninfo = chaninfo;
    chans = tmpEEG;
end

% ---------------------------------------------
% separate data channels from non-data channels
% ---------------------------------------------
function [chans, fidsval] = getnodatchan(chans)
if isfield(chans,'datachan')
    [chans(cellfun('isempty',{chans.datachan})).datachan] = deal(0);
    fids = [chans.datachan] == 0;
    fidsval = chans(fids);    
    chans = rmfield(chans(~fids),'datachan');
else
    fids = [];
end

% ----------------------------------------
% fuse data channels and non-data channels
% ----------------------------------------
function [chans, chaninfo,complicated] = insertchans(chans, chaninfo, nchans)
if nargin < 3, nchans = length(chans); end
[chans.datachan] = deal(1);
complicated = false;        % whether we need complicated treatment of datachans & co further down the road.....

if isfield(chans,'type')
    mask = strcmpi({chans.type},'FID') | strcmpi({chans.type},'IGNORE');
    if any(mask)
        [chans(mask).datachan] = deal(0);
        complicated = true;
    end
end
if length(chans) > nchans && nchans ~= 0 % reference at the end of the structure
    chans(end).datachan = 0;
    complicated = true;
end
if isfield(chaninfo, 'nodatchans')
    if ~isempty(chaninfo.nodatchans) && isstruct(chaninfo.nodatchans)
        chanlen = length(chans);
        for index = 1:length(chaninfo.nodatchans)
            fields = fieldnames( chaninfo.nodatchans );
            ind = chanlen+index;
            for f = 1:length( fields )
                chans = setfield(chans, { ind }, fields{f}, getfield( chaninfo.nodatchans, { index },  fields{f}));
            end
            chans(ind).datachan = 0;
            complicated = true;
        end
        chaninfo = rmfield(chaninfo, 'nodatchans');
        
        % put these channels first
        % ------------------------
        % tmp = chans(chanlen+1:end);
        % chans(length(tmp)+1:end) = chans(1:end-length(tmp));
        % chans(1:length(tmp)) = tmp;
    end
end
