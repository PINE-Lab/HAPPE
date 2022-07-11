% pop_loadbv() - load Brain Vision Data Exchange format dataset and
%                return EEGLAB EEG structure
%
% Usage:
%   >> [EEG, com] = pop_loadbv; % pop-up window mode
%   >> [EEG, com] = pop_loadbv(path, hdrfile);
%   >> [EEG, com] = pop_loadbv(path, hdrfile, srange);
%   >> [EEG, com] = pop_loadbv(path, hdrfile, [], chans);
%   >> [EEG, com] = pop_loadbv(path, hdrfile, srange, chans);
%
% Optional inputs:
%   path      - path to files
%   hdrfile   - name of Brain Vision vhdr-file (incl. extension)
%   srange    - scalar first sample to read (up to end of file) or
%               vector first and last sample to read (e.g., [7 42];
%               default: all)
%   chans     - vector channels channels to read (e.g., [1:2 4];
%               default: all)
%   metadata  - [true|false] when true, only read meta data. Default
%               false
%
% Outputs:
%   EEG       - EEGLAB EEG structure
%   com       - history string
%
% Note:
%   Import "Brain Vision Data Exchange" format files with this function.
%   Brain Vision Data Exchange files consist of a set of 3 files, a header
%   file (.vhdr), a marker file (.vmrk), and a data file. Export from
%   BrainVision Analyzer with "Generic Data" export. Select header and
%   marker file for export (text format; XML format is not yet supported).
%   Binary and text data formats, in both multiplexed and vectorized data
%   orientation are supported. Binary data formats offer higher precision
%   and faster file import.
%
% Author: Andreas Widmann & Arnaud Delorme, 2004-

% Copyright (C) 2004 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
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

% $Id: pop_loadbv.m 53 2010-05-22 21:57:38Z arnodelorme $
% Revision 1.5 2010/03/23 21:19:52 roy
% added some lines so that the function can deal with the space lines in the ASCII multiplexed data file

function [EEG, com] = pop_loadbv(path, hdrfile, srange, chans, metadata)

com = '';
EEG = [];
if nargin < 5
    metadata = false;
end
if nargin < 2
    [hdrfile, path] = uigetfile2({'*.vhdr' '*.ahdr'}, 'Select Brain Vision vhdr-file - pop_loadbv()');
    if hdrfile(1) == 0, return; end

    drawnow;
    uigeom = {[1 0.5] [1 0.5]};
    uilist = {{ 'style' 'text' 'string' 'Interval (samples; e.g., [7 42]; default: all):'} ...
              { 'style' 'edit' 'string' ''} ...
              { 'style' 'text' 'string' 'Channels (e.g., [1:2 4]; default: all):'} ...
              { 'style' 'edit' 'string' ''}};
    result = inputgui(uigeom, uilist, 'pophelp(''pop_loadbv'')', 'Load a Brain Vision Data Exchange format dataset');
    if isempty(result), return, end
    if ~isempty(result{1})
        srange = str2num(result{1});
    end
    if ~isempty(result{2})
        chans = str2num(result{2});
    end
end

% Header file
disp('pop_loadbv(): reading header file');
hdr = readbvconf(path, hdrfile);

% Common Infos
try
    EEG = eeg_emptyset;
catch
end
EEG.comments = ['Original file: ' hdr.commoninfos.datafile];
hdr.commoninfos.numberofchannels = str2double(hdr.commoninfos.numberofchannels);
EEG.srate = 1000000 / str2double(hdr.commoninfos.samplinginterval);

% Binary Infos
if strcmpi(hdr.commoninfos.dataformat, 'binary')
    switch lower(hdr.binaryinfos.binaryformat)
        case 'int_16',        binformat = 'int16'; bps = 2;
        case 'uint_16',       binformat = 'uint16'; bps = 2;
        case 'ieee_float_32', binformat = 'float32'; bps = 4;
        otherwise, error('Unsupported binary format');
    end
end

% Channel Infos
if ~exist('chans', 'var') || isempty(chans)
    chans = 1:hdr.commoninfos.numberofchannels;
    EEG.nbchan = hdr.commoninfos.numberofchannels;
else
    EEG.nbchan = length(chans);
end
if any(chans < 1) || any(chans > hdr.commoninfos.numberofchannels)
    error('chans out of available channel range');
end
if isfield(hdr, 'channelinfos')
    for chan = 1:length(chans)
        try
            [EEG.chanlocs(chan).labels, chanlocs(chan).ref, chanlocs(chan).scale, chanlocs(chan).unit] = strread(hdr.channelinfos{chans(chan)}, '%s%s%s%s', 1, 'delimiter', ',');
        catch % Octave compatible code below
            str  = hdr.channelinfos{chans(chan)};
            [EEG.chanlocs(chan).labels, str] = strtok(str, ',');
            [chanlocs(chan).ref,        str] = strtok(str, ',');
            [chanlocs(chan).scale,      str] = strtok(str, ',');
            [chanlocs(chan).unit,       str] = strtok(str, ',');
        end
        EEG.chanlocs(chan).labels = char(EEG.chanlocs(chan).labels);
        chanlocs(chan).scale = str2double(char(chanlocs(chan).scale));
%             chanlocs(chan).unit = native2unicode(double(char(chanlocs(chan).scale)), 'UTF-8');
%             EEG.chanlocs(chan).datachan = chans(chan);
    end
    if isempty([chanlocs.scale])
        chanlocs = rmfield(chanlocs, 'scale');
    end
end
%     [EEG.chanlocs.type] = deal([]);

% Coordinates
if isfield(hdr, 'coordinates')
    hdr.coordinates(end+1:length(chans)) = { [] };
    onenon0channel = 0;
    for chan = 1:length(chans)
        if ~isempty(hdr.coordinates{chans(chan)})
            if ismatlab
                [EEG.chanlocs(chan).sph_radius, theta, phi] = strread(hdr.coordinates{chans(chan)}, '%f%f%f', 'delimiter', ',');
            else
                str  = hdr.coordinates{chans(chan)};
                [EEG.chanlocs(chan).sph_radius, str] = strtok(str, ','); EEG.chanlocs(chan).sph_radius = str2num(EEG.chanlocs(chan).sph_radius);
                [theta,                         str] = strtok(str, ','); theta = str2num(theta);
                [phi,                           str] = strtok(str, ','); phi   = str2num(phi);
            end
            if EEG.chanlocs(chan).sph_radius == 0 && theta == 0 && phi == 0
                EEG.chanlocs(chan).sph_radius = [];
                EEG.chanlocs(chan).sph_theta = [];
                EEG.chanlocs(chan).sph_phi = [];
            else
                onenon0channel = 1;
                EEG.chanlocs(chan).sph_theta = phi - 90 * sign(theta);
                EEG.chanlocs(chan).sph_phi = -abs(theta) + 90;
            end
        end
    end
    try
        if onenon0channel
            [EEG.chanlocs, EEG.chaninfo] = pop_chanedit(EEG.chanlocs, 'convert', 'sph2topo');
            [EEG.chanlocs, EEG.chaninfo] = pop_chanedit(EEG.chanlocs, 'convert', 'sph2cart');
        end
    catch, end
end

% Open data file and find the number of data points
% -------------------------------------------------
disp('pop_loadbv(): reading EEG data');
[IN, message] = fopen(fullfile(path, hdr.commoninfos.datafile), 'r');
if IN == -1
    [IN, message] = fopen(fullfile(path, lower(hdr.commoninfos.datafile)));
    if IN == -1
        hdr.commoninfos.datafile = [ hdrfile(1:end-4) 'eeg' ];
        fprintf(2, 'The header file points to a binary file that does not exist\n');
        fprintf(2, 'trying to open a binary file with the same name as the header\nfile and a different extension...\n');
        [IN, message] = fopen(fullfile(path, hdr.commoninfos.datafile), 'r');
        if IN == -1
            [IN, message] = fopen(fullfile(path, lower(hdr.commoninfos.datafile)));
            if IN == -1
                error(message)
            end
        end
    end
end
if isfield( hdr.commoninfos, 'datapoints' ) && ~isempty( hdr.commoninfos.datapoints ) && isnumeric( str2double( hdr.commoninfos.datapoints ) ) && str2double( hdr.commoninfos.datapoints ) > 0
    hdr.commoninfos.datapoints = str2double(hdr.commoninfos.datapoints);
elseif strcmpi( hdr.commoninfos.dataformat, 'binary' )
    fseek( IN, 0, 'eof' );
    hdr.commoninfos.datapoints = ftell( IN ) / ( hdr.commoninfos.numberofchannels * bps );
    if round(hdr.commoninfos.datapoints) ~= hdr.commoninfos.datapoints
        fprintf(2, 'Binary file is truncated/incomplete\n');
        hdr.commoninfos.datapoints = floor(hdr.commoninfos.datapoints);
    end
    fseek( IN, 0, 'bof' );
else
    hdr.commoninfos.datapoints = NaN;
end

if ~strcmpi(hdr.commoninfos.dataformat, 'binary') % ASCII
    % tmppoint = hdr.commoninfos.datapoints;
    tmpchan = fscanf(IN, '%s', 1);
    
    % AW: Determination of number of datapoints will not work for files without chanlabels and/or multiplexed dataformat. Suggest trusting in header.
    % tmpdata = fscanf(IN, '%f', inf);
    % hdr.commoninfos.datapoints = length(tmpdata);
    % chanlabels = 1;
    if isnan(str2double(tmpchan)) 
        % hdr.commoninfos.datapoints = hdr.commoninfos.datapoints+1; 
        chanlabels = 1;
    else
        chanlabels = 0;
    end
end

% Sample range
if ~exist('srange', 'var') || isempty(srange)
    srange = [ 1 hdr.commoninfos.datapoints];
    EEG.pnts = hdr.commoninfos.datapoints;
elseif length(srange) == 1
    EEG.pnts = hdr.commoninfos.datapoints - srange(1) + 1;
else
    EEG.pnts = srange(2) - srange(1) + 1;
end
if any(srange < 1) || any(srange > hdr.commoninfos.datapoints)
    error('srange out of available data range');
end

% Read data
if metadata
    EEG.data = [];
elseif strcmpi(hdr.commoninfos.dataformat, 'binary')
    switch lower(hdr.commoninfos.dataorientation)
        case 'multiplexed'
            if EEG.nbchan == hdr.commoninfos.numberofchannels % Read all channels
                fseek(IN, (srange(1) - 1) * EEG.nbchan * bps, 'bof');
                EEG.data = fread(IN, [EEG.nbchan, EEG.pnts], [binformat '=>float32']);
            else % Read channel subset
                EEG.data = repmat(single(0), [EEG.nbchan, EEG.pnts]); % Preallocate memory
                for chan = 1:length(chans)
                    fseek(IN, (srange(1) - 1) * hdr.commoninfos.numberofchannels * bps + (chans(chan) - 1) * bps, 'bof');
                    EEG.data(chan, :) = fread(IN, [1, EEG.pnts], [binformat '=>float32'], (hdr.commoninfos.numberofchannels - 1) * bps);
                end
            end
        case 'vectorized'
            if isequal(EEG.pnts, hdr.commoninfos.datapoints) && EEG.nbchan == hdr.commoninfos.numberofchannels % Read entire file
                EEG.data = fread(IN, [EEG.pnts, EEG.nbchan], [binformat '=>float32']).';
            else % Read fraction of file
                EEG.data = repmat(single(0), [EEG.nbchan, EEG.pnts]); % Preallocate memory
                for chan = 1:length(chans)
                    fseek(IN, ((chans(chan) - 1) * hdr.commoninfos.datapoints + srange(1) - 1) * bps, 'bof');
                    EEG.data(chan, :) = fread(IN, [1, EEG.pnts], [binformat '=>float32']);
                end
            end
        otherwise
            error('Unsupported data orientation')
    end
else % ASCII data
%     disp('If this function does not work, export your data in binary format');
%     EEG.data = repmat(single(0), [EEG.nbchan, EEG.pnts]);
%     if strcmpi(lower(hdr.commoninfos.dataorientation), 'vectorized')
%         count = 1;
%         fseek(IN, 0, 'bof');
%         len = inf;
%         for chan = 1:hdr.commoninfos.numberofchannels
%             if chanlabels, tmpchan = fscanf(IN, '%s', 1); end;
%             tmpdata = fscanf(IN, '%f', len); len = length(tmpdata);
%             if ismember(chan, chans)
%                 EEG.data(count, :) = tmpdata(srange(1):srange(2))';
%                 count = count + 1;
%             end;
%         end;
%     elseif strcmpi(lower(hdr.commoninfos.dataorientation), 'multiplexed')
% %         fclose(IN);
% %         error('ASCII multiplexed reading not implemeted yet; export as a different format');
%         if EEG.nbchan == hdr.commoninfos.numberofchannels % Read all channels
%             tmpchan= fgetl(IN);
%             count = 1;
%             while ~feof(IN)
%                 tmpstr = fgetl(IN);
%                 if ~isempty(tmpstr)
%                     temp_ind = tmpstr==',';
%                     tmpstr(temp_ind) = '.';
%                     tmpdata = strread(tmpstr);
%                     EEG.data(:,count) = tmpdata';
%                     count = count + 1;
%                 end;
%             end;
%             EEG.pnts = count - 1;
%         else
%             
%         end;
%     end;

    % Rewritten by AW, 2013-10-02. Old version by Arno did no longer work. MATLAB changes?
    tmpdata = zeros([hdr.commoninfos.numberofchannels, hdr.commoninfos.datapoints], 'single');
    fseek(IN, 0, 'bof');
    switch lower(hdr.commoninfos.dataorientation)

        case 'vectorized'
            if chanlabels || (isfield(hdr, 'asciiinfos') && isfield(hdr.asciiinfos, 'decimalsymbol') && ~strcmp(hdr.asciiinfos.decimalsymbol, '.')) % Read line by line
 
                for iChan = 1:hdr.commoninfos.numberofchannels
                    tmpstr = fgetl(IN);
                    if chanlabels
                        [tmpchan, count, errmsg, nextindex] = sscanf(tmpstr, '%s', 1);
                        tmpstr = tmpstr(nextindex:end);
                    end
                    if isfield(hdr, 'asciiinfos') && isfield(hdr.asciiinfos, 'decimalsymbol') && ~strcmp(hdr.asciiinfos.decimalsymbol, '.')
                        tmpdata(iChan, :) = sscanf(regexprep(tmpstr, hdr.asciiinfos.decimalsymbol, '.'), '%f', inf);
                    else
                        tmpdata(iChan, :) = sscanf(tmpstr, '%f', inf);
                    end
                end
                
            else
                tmpdata = fscanf(IN, '%f', inf);
                if length(tmpdata) ~= hdr.commoninfos.numberofchannels*hdr.commoninfos.datapoints
                    error('File truncated: cannot be imported');
                end
                tmpdata = reshape(tmpdata, hdr.commoninfos.datapoints, hdr.commoninfos.numberofchannels)';
            end

        case 'multiplexed'
            if chanlabels
                tmpchan = fgetl(IN);
            end
            if isfield(hdr, 'asciiinfos') && isfield(hdr.asciiinfos, 'decimalsymbol') && ~strcmp(hdr.asciiinfos.decimalsymbol, '.')  % Read line by line
               for iPnt = 1:hdr.commoninfos.datapoints
                    tmpstr = fgetl(IN);
                    tmpdata(:, iPnt) = sscanf(regexprep(tmpstr, hdr.asciiinfos.decimalsymbol, '.'), '%f', inf);
                end
            else
                tmpdata = fscanf(IN, '%f', inf);
                if length(tmpdata) ~=  hdr.commoninfos.numberofchannels*hdr.commoninfos.datapoints
                    warning('File truncated');
                    hdr.commoninfos.datapoints = floor(length(tmpdata)/hdr.commoninfos.numberofchannels);
                    tmpdata = tmpdata(1:hdr.commoninfos.numberofchannels*hdr.commoninfos.datapoints);
                    srange(2) = min(srange(2), hdr.commoninfos.datapoints);
                end
                tmpdata = reshape(tmpdata, hdr.commoninfos.numberofchannels, hdr.commoninfos.datapoints);
            end

        otherwise
            error('Unknown data orientation')

    end

    EEG.data = tmpdata(chans, srange(1):srange(2));
end

fclose(IN);
EEG.trials = 1;
EEG.xmin   = 0;
EEG.xmax   = (EEG.pnts - 1) / EEG.srate;

% Convert to EEG.data to double for MATLAB < R14
if str2double(version('-release')) < 14
    EEG.data = double(EEG.data);
end

% Scale data
if exist('chanlocs', 'var') && isfield(chanlocs, 'scale') && ~isempty(EEG.data)
    disp('pop_loadbv(): scaling EEG data');
    for chan = 1:EEG.nbchan
        if ~isnan(chanlocs(chan).scale)
            EEG.data(chan, :) = EEG.data(chan, :) * chanlocs(chan).scale;
        end
    end
end

% Marker file
if isfield(hdr.commoninfos, 'markerfile')
    disp('pop_loadbv(): reading marker file');
    try 
        MRK = readbvconf(path, hdr.commoninfos.markerfile);
    catch
        MRK = readbvconf(path, [hdrfile(1:end-4) 'vmrk' ]);
    end
    if ~isequal(hdr.commoninfos.datafile, MRK.commoninfos.datafile)
        disp('pop_loadbv() warning: data files in header and marker files inconsistent.');
    end

    % Marker infos
    if isfield(MRK, 'markerinfos')
        EEG.event = parsebvmrk(MRK);

        % Correct event latencies by first sample offset
        tmpevent = EEG.event;
        for index = 1:length(EEG.event)
            tmpevent(index).latency = tmpevent(index).latency - srange(1) + 1;
        end
        EEG.event = tmpevent;

        % Remove unreferenced events
        EEG.event = EEG.event([tmpevent.latency] >= 1 & [tmpevent.latency] <= EEG.pnts);

        % Copy event structure to urevent structure
        EEG.urevent = rmfield(EEG.event, 'urevent');

        % find if boundaries at homogenous intervals
        % ------------------------------------------
        tmpevent = EEG.event;
        boundaries = strmatch('boundary', {tmpevent.type});
        boundlats = unique([tmpevent(boundaries).latency]);
        if (isfield(hdr.commoninfos, 'segmentationtype') && (strcmpi(hdr.commoninfos.segmentationtype, 'markerbased') || strcmpi(hdr.commoninfos.segmentationtype, 'fixtime'))) && length(boundaries) > 1 && length(unique(diff([boundlats EEG.pnts + 1]))) == 1
            EEG.trials = length(boundlats);
            EEG.pnts   = EEG.pnts / EEG.trials;
            EEG.event(boundaries) = [];

            % adding epoch field
            % ------------------
            tmpevent  = EEG.event;
            for index = 1:length(EEG.event)
                EEG.event(index).epoch = ceil(tmpevent(index).latency / EEG.pnts);
            end

            % finding minimum time
            % --------------------
            tles = strmatch('time 0', lower({tmpevent.code}))';
            if ~isempty(tles)
                for iTLE = tles(:)'
                    EEG.event(iTLE).type ='TLE';
                end
                EEG.xmin = -(tmpevent(tles(1)).latency - 1) / EEG.srate;
            end
        else
            for index = 1:length(boundaries)
                EEG.event(boundaries(index)).duration = NaN;
            end
        end
    end
end

EEG.ref = 'common';

try
    EEG = eeg_checkset(EEG);
catch
end

if nargout == 2
    com = sprintf('EEG = pop_loadbv(''%s'', ''%s'', %s, %s);', path, hdrfile, mat2str(srange), mat2str(chans));
end
