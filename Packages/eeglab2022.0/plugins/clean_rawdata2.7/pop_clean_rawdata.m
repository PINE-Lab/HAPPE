% pop_clean_rawdata(): Launches GUI to collect user inputs for clean_artifacts().
%                      ASR stands for artifact subspace reconstruction.
%                      To disable method(s), enter -1.
% Usage:
%   >>  EEG = pop_clean_rawdata(EEG);
%
% see also: clean_artifacts

% Author: Arnaud Delorme, Makoto Miyakoshi and Christian Kothe, SCCN,INC,UCSD
% History:
% 07/2019. Reprogrammed from Scratch (Arnaud Delorme)
% 07/31/2018 Makoto. Returns error if input data size is 3. 
% 04/26/2017 Makoto. Deletes existing EEG.etc.clean_channel/sample_mask. Try-catch to skip potential error in vis_artifact.
% 07/18/2014 ver 1.4 by Makoto and Christian. New channel removal method supported. str2num -> str2num due to str2num([a b]) == NaN.
% 11/08/2013 ver 1.3 by Makoto. Menu words changed. asr_process() line 168 bug fixed. 
% 10/07/2013 ver 1.2 by Makoto. Help implemented. History bug fixed.
% 07/16/2013 ver 1.1 by Makoto and Christian. Minor update for help and default values.
% 06/26/2013 ver 1.0 by Makoto. Created.

% Copyright (C) 2013, Arnaud Delorme, Makoto Miyakoshi and Christian Kothe, SCCN,INC,UCSD
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [EEG, com] = pop_clean_rawdata(EEG, varargin)

% Check input
com = '';
if size(EEG(1).data) == 3
    error('Input data must be continuous. This data seems epoched.')
end

fusechanrej = false;
if nargin < 2
    % Obtain user inputs.
    cb_filter = 'if get(gcbo, ''value''), set(findobj(gcbf, ''userdata'', ''filter''), ''enable'', ''on''); else set(findobj(gcbf, ''userdata'', ''filter''), ''enable'', ''off'', ''value'', 0); end';
    cb_chan   = 'if get(gcbo, ''value''), set(findobj(gcbf, ''userdata'', ''chan'')  , ''enable'', ''on''); else set(findobj(gcbf, ''userdata'', ''chan'')  , ''enable'', ''off'', ''value'', 0); end';
    cb_asr    = 'if get(gcbo, ''value''), set(findobj(gcbf, ''userdata'', ''asr'')   , ''enable'', ''on''); else set(findobj(gcbf, ''userdata'', ''asr'')   , ''enable'', ''off''); end';
    cb_rej    = 'if get(gcbo, ''value''), set(findobj(gcbf, ''userdata'', ''rej'')   , ''enable'', ''on''); else set(findobj(gcbf, ''userdata'', ''rej'')   , ''enable'', ''off''); end';
    cb_select1   = 'pop_chansel(get(gcbf, ''userdata''), ''field'', ''labels'', ''handle'', findobj(''parent'', gcbf, ''tag'', ''chanuse''));    set(findobj(''parent'', gcbf, ''tag'', ''chanuseflag''), ''value'', 1);    set(findobj(''parent'', gcbf, ''tag'', ''chanignoreflag''), ''value'', 0);';
    cb_select2   = 'pop_chansel(get(gcbf, ''userdata''), ''field'', ''labels'', ''handle'', findobj(''parent'', gcbf, ''tag'', ''chanignore'')); set(findobj(''parent'', gcbf, ''tag'', ''chanignoreflag''), ''value'', 1); set(findobj(''parent'', gcbf, ''tag'', ''chanuseflag''), ''value'', 0);';
    winsize   = max(0.5,1.5*EEG(1).nbchan/EEG(1).srate);
    uilist =    {...
        {'style' 'checkbox' 'string' 'Remove channel drift (data not already high-pass filtered)' 'fontweight' 'bold' 'tag' 'filter' 'callback' cb_filter} ...
        {} {'style' 'text' 'string' 'Linear filter (FIR) transition band [lo hi] in Hz            ' 'userdata' 'filter' 'enable' 'off' } ...
        {'style' 'edit' 'string' '0.25 0.75', 'enable' 'off' 'tag','filterfreqs', 'userdata' 'filter' 'tooltipstring', wordwrap('The first number is the frequency below which everything is removed, and the second number is the frequency above which everything is retained. There is a linear transition in between. For best performance of subsequent processing steps the upper frequency should be close to 1 or 2 Hz, but you can go lower if certain activities need to be retained.',80)} ...
        ...
        {} ...
        {'style' 'checkbox' 'string' 'Process/remove channels' 'fontweight' 'bold' 'tag' 'chanrm' 'callback' cb_chan 'value' 1 } ...
        ...
        {} { 'style' 'checkbox' 'string' 'Only consider these channels' 'tag' 'chanuseflag' 'userdata' 'chan' } ...
        { 'style' 'pushbutton' 'string' '...' 'userdata' 'chan' 'callback' cb_select1 } ...
        {'style' 'edit' 'string' '' 'userdata' 'chan'  'tag' 'chanuse'  } ...
        ...
        {} {'style' 'checkbox' 'string' 'Ignore these channels (ECG, EMG, ...)' 'userdata' 'chan' 'value' 0 'tag' 'chanignoreflag'  } ...
        {'style' 'pushbutton' 'string' '...', 'userdata' 'chan' 'callback' cb_select2 } ...
        {'style' 'edit' 'string' '', 'userdata' 'chan' 'tag' 'chanignore'  } ...
        ...
        {} {'style' 'checkbox' 'string' 'Remove channel if it is flat for more than (seconds)' 'tag' 'rmflat' 'userdata' 'chan' 'value' 1 } ...
        {'style' 'edit' 'string' '5', 'userdata' 'chan' 'tag' 'rmflatsec' 'tooltipstring', wordwrap('If a channel has a longer flatline than this, it will be removed. In seconds.',80)} ...
        ...
        {} {'style' 'checkbox' 'string' 'Max acceptable high-frequency noise std dev' 'value' 1 'tag' 'rmnoise' 'userdata' 'chan' } ...
        {'style' 'edit' 'string' '4',  'userdata' 'chan' 'tag' 'rmnoiseval' 'tooltipstring', wordwrap('If a channel has more line noise relative to its signal than this value, in standard deviations relative to the overall channel population, it will be removed.',80)} ...
        ...
        {} {'style' 'checkbox' 'string' 'Min acceptable correlation with nearby chans [0-1]' 'value' 1 'tag' 'rmcorr' 'userdata' 'chan'   } ...
        {'style' 'edit' 'string' '0.8', 'userdata' 'chan'  'tag' 'rmcorrval'  'tooltipstring', wordwrap('If a channel has lower correlation than this to an estimate of its activity based on other channels, and this applies to more than half of the recording, the channel will be removed. This method requires that channel locations are available and roughly correct; otherwise a fallback criterion will tried used using a default setting; you can customize the fallback method by directly calling clean_channels_nolocs in the command line.',80)} ...
        ...
        {} {'style' 'checkbox' 'string' 'Perform Artifact Subspace Reconstruction bad burst correction' 'fontweight' 'bold' 'value' 1 'tag' 'asr' 'callback' cb_asr } ...
        {} {'style' 'text' 'string' sprintf('Max acceptable %1.1f second window std dev', winsize) 'value' 1  'userdata' 'asr' } ...
        {'style' 'edit' 'string' '20', 'tag' 'asrstdval' 'userdata' 'asr' 'tooltipstring', wordwrap('Standard deviation cutoff for removal of bursts. Data portions whose variance is larger than this threshold relative to the calibration data are considered missing data and will be removed. The most aggressive value that can be used without losing much EEG is 3. A reasonably conservative value is 5, but some extreme EEG bursts (e.g., sleep spindles) can cross even 5. For new users it is recommended to at first visually inspect the difference between the original and cleaned data to get a sense of the removed content at various levels.',80)} ...
        {} {'style' 'checkbox' 'string' 'Use Riemanian distance metric (not Euclidean) - beta' 'userdata' 'asr' 'value' 0 'tag' 'distance' } {} ...
        {} {'style' 'checkbox' 'tag' 'asrrej' 'string' 'Remove bad data periods (instead of correcting them)' 'value' 1 'userdata' 'asr'} {} ...
        ...
        {} {'style' 'checkbox' 'string' 'Additional removal of bad data periods' 'fontweight' 'bold' 'value' 1 'tag' 'rejwin' 'callback' cb_rej } ...
        {} {'style' 'text' 'tag' 'asrwintext' 'string' 'Acceptable [min max] channel RMS range (+/- std dev)'  'userdata' 'rej'} ...
        {'style' 'edit' 'string' '-Inf 7','tag', 'rejwinval1', 'userdata' 'rej' 'tooltipstring', wordwrap('If a time window has a larger fraction of simultaneously corrupted channels than this (after the other cleaning attempts), it will be cut out of the data. This can happen if a time window was corrupted beyond the point where it could be recovered.',80)} ...
        {} {'style' 'text' 'tag' 'asrwintext' 'string' 'Maximum out-of-bound channels (%)'  'userdata' 'rej'} ...
        {'style' 'edit' 'string' '25','tag', 'rejwinval2', 'userdata' 'rej' 'tooltipstring', wordwrap('If a time window has a larger fraction of simultaneously corrupted channels than this (after the other cleaning attempts), it will be cut out of the data. This can happen if a time window was corrupted beyond the point where it could be recovered.',80)} ...
        ...
        {} {'style' 'checkbox' 'string' 'Pop up scrolling data window with rejected data highlighted' 'tag' 'vis' 'value' fastif(length(EEG) > 1, 0, 1) 'enable' fastif(length(EEG) > 1, 'off', 'on') }};

    % channel labels
    % --------------
    if ~isempty(EEG(1).chanlocs)
        tmpchanlocs = EEG(1).chanlocs;        
    else
        tmpchanlocs = [];
        for index = 1:EEG(1).nbchan
            tmpchanlocs(index).labels = int2str(index);
            tmpchanlocs(index).type = '';
        end
    end
    
    row   = [0.1 1 0.3];
    row4  = [0.1 0.8 0.2 0.3];
    row2  = [0.1 1.2 0.1];
    row3  = [0.9 0.2 0.3];
    geom =     { 1 row     1   1 row4 row4 row     row     row     1   1 row     row2 row2   1   1 row  row     1   1 };
    geomvert = [ 1 1       0.3 1 1    1    1       1       1       0.3 1 1       1    1      0.3 1 1    1       0.3 1 ];

    if length(EEG) > 1
        uilist = { uilist{1:6} {} {'style' 'checkbox' 'string' 'Fuse channel rejection for datasets with same subject and session' 'userdata' 'chan' 'fontweight' 'bold' 'tag' 'commonrej' 'value' 1} {} uilist{7:end} };
        geom = [ geom(1:4) { row2 } geom(5:end)];
        geomvert = [ geomvert(1:4) 1 geomvert(5:end) ];
    end

    [res,~,~,outs] = inputgui('title', 'pop_clean_rawdata()', 'geomvert', geomvert, 'geometry', geom, 'uilist',uilist, 'helpcom', 'pophelp(''clean_artifacts'');', 'userdata', tmpchanlocs);

    % Return error if no input.
    if isempty(res) return; end

    % process multiple datasets
    % -------------------------
    options = {};
    opt.FlatlineCriterion  = 'off';
    opt.ChannelCriterion   = 'off';
    opt.LineNoiseCriterion = 'off';
    opt.Highpass           = 'off';
    opt.BurstCriterion     = 'off';
    opt.WindowCriterion    = 'off';
    opt.BurstRejection     = 'off';
    opt.Distance           = 'Euclidian';
    
    if outs.filter, opt.Highpass = str2num(outs.filterfreqs); end
    
    if outs.chanrm
        if outs.chanignoreflag
            [ chaninds, chanlist ] = eeg_decodechan(EEG(1).chanlocs, outs.chanignore);
            if isempty(chanlist), chanlist = chaninds; end
            opt.channels_ignore = chanlist; 
        end
        if outs.chanuseflag
           [ chaninds, chanlist ] = eeg_decodechan(EEG(1).chanlocs, outs.chanuse);
           if isempty(chanlist), chanlist = chaninds; end
           opt.channels = chanlist;
        end
        if outs.rmflat, opt.FlatlineCriterion = str2num(outs.rmflatsec); end
        if outs.rmcorr, opt.ChannelCriterion  = str2num(outs.rmcorrval); end
        if outs.rmnoise, opt.LineNoiseCriterion = str2num(outs.rmnoiseval); end
    end
    
    if outs.asr
        opt.BurstCriterion = str2num(outs.asrstdval); 
        if outs.distance,  opt.Distance = 'Riemannian'; end
    end
    
    if outs.rejwin
        opt.WindowCriterionTolerances = str2num(outs.rejwinval1);
        opt.WindowCriterion = str2num(outs.rejwinval2)/100;
    end
    if outs.asrrej && ~strcmpi(opt.BurstCriterion, 'off')
        opt.BurstRejection = 'on';
    end
    if isfield(outs,'commonrej') && outs.commonrej
        opt.fusechanrej = outs.commonrej;
        fusechanrej = true;
    end

    % convert structure to cell
    options = fieldnames(opt);
    options(:,2) = struct2cell(opt);
    options = options';
    options = options(:)';
else
    options = varargin;
    for iOpt = 1:2:length(options)
        if strcmpi(options{iOpt}, 'fusechanrej')
            fusechanrej = true;
        end
    end
end

if length(EEG) > 1
    % process multiple datasets
    if nargin < 2
        [ EEG, com ] = eeg_eval( 'clean_artifacts', EEG, 'warning', 'on', 'params', options );
    else
        [ EEG, com ] = eeg_eval( 'clean_artifacts', EEG, 'params', options );
    end
    if fusechanrej
        if exist('pop_fusechanrej')
            EEG = pop_fusechanrej(EEG);
        else
            warning('Upgrade your version of EEGLAB to fuse channel rejection for datasets with same subject and session');
        end
    end
    return;
end

% Delete EEG.etc.clean_channel_mask and EEG.etc.clean_sample_mask if present.
if isfield(EEG.etc, 'clean_channel_mask')
    EEG.etc = rmfield(EEG.etc, 'clean_channel_mask');
    disp('EEG.etc.clean_channel_mask present: Deleted.')
end
if isfield(EEG.etc, 'clean_sample_mask')
    EEG.etc = rmfield(EEG.etc, 'clean_sample_mask');
    disp('EEG.etc.clean_sample_mask present: Deleted.')
end

cleanEEG = clean_artifacts(EEG, options{:});

% Apply Christian's function before and after comparison visualization.
if nargin < 2 && outs.vis == 1
    try
        vis_artifacts(cleanEEG,EEG);
    catch
        warning('vis_artifacts failed. Skipping visualization. Could be because of duplicate channel label.')
    end
end

% Update EEG.
EEG = cleanEEG;

% Output eegh.
com = sprintf('EEG = pop_clean_rawdata(EEG, %s);', vararg2str(options));

% Display the ending message.
disp('Done.')

function outtext = wordwrap(intext,nChars)
outtext = '';    
while ~isempty(intext)
    if length(intext) > nChars
        cutoff = nChars+find([intext(nChars:end) ' ']==' ',1)-1;
        outtext = [outtext intext(1:cutoff-1) '\n']; %#ok<*AGROW>
        intext = intext(cutoff+1:end);
    else 
        outtext = [outtext intext];
        intext = '';
    end
end
outtext = sprintf(outtext);
