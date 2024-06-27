% readbvconf() - read Brain Vision Data Exchange format configuration 
%                file
%
% Usage:
%   >> CONF = readbvconf(pathname, filename);
%
% Inputs:
%   pathname  - path to file
%   filename  - filename
%
% Outputs:
%   CONF      - structure configuration
%
% Author: Andreas Widmann, University of Leipzig, 2007

% Copyright (C) 2007 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
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

% $Id: readbvconf.m 44 2009-11-12 02:00:56Z arnodelorme $

function CONF = readbvconf(pathname, filename)

if nargin < 2
    error('Not enough input arguments');
end
% 
% % Open and read file (old method, much slower)
% [IN, message] = fopen(fullfile(pathname, filename), 'r');
% if IN == -1
%     [IN, message] = fopen(fullfile(pathname, lower(filename)));
%     if IN == -1
%         error(message)
%     end;
% end
% raw={};
% while ~feof(IN)
%     raw = [raw; {fgetl(IN)}];
% end
% fclose(IN);
% 
% % Remove comments and empty lines
% raw(cellfun('isempty', raw) == true) = [];
% raw(strmatch(';', raw)) = [];

% Open and read file (automatically remove empty lines)
fid = fopen(fullfile(pathname, filename), 'r');
raw = textscan(fid, '%s', 'delimiter', '');
fclose(fid);
raw = raw{1};
raw(strmatch(';', raw)) = []; % remove comments

% Find sections
sectionArray = [strmatch('[', raw)' length(raw) + 1];
for iSection = 1:length(sectionArray) - 1

    % Convert section name
    tmpstr    = deblank(raw{sectionArray(iSection)});
    fieldName = lower(tmpstr(2:end-1));
    %fieldName = lower(char(strread(tmpstr(2:end), '[%s', 'delimiter', ']')));
    fieldName(isspace(fieldName) == true) = [];

    % Fill structure with parameter value pairs
    switch fieldName
        case {'commoninfos' 'binaryinfos' 'asciiinfos'}
            for line = sectionArray(iSection) + 1:sectionArray(iSection + 1) - 1
                splitArray = strfind(raw{line}, '=');
                fieldName2  = lower(raw{line}(1:splitArray(1) - 1));
                fieldName2(fieldName2 == ' ') = '_';
                fieldName2(fieldName2 == ':') = '_';
                CONF.(fieldName).(fieldName2) = raw{line}(splitArray(1) + 1:end);
            end
        case {'channelinfos' 'coordinates'}
            for line = sectionArray(iSection) + 1:sectionArray(iSection + 1) - 1
                splitArray = strfind(raw{line}, '=');
                CONF.(fieldName)(str2double(raw{line}(3:splitArray(1) - 1))) = {raw{line}(splitArray(1) + 1:end)};
            end
        case {'markerinfos'} % Allow discontinuity for markers (but not channelinfos and coordinates!)
            % try reading the whole section at once
            
            try
                fileData = raw(sectionArray(iSection) + 1:sectionArray(iSection + 1) - 1);
                fileData = sprintf('%s\n', fileData{:});
                newData = textscan(fileData', '%s', 'delimiter', '=');
                newData = newData{1};
                eventValues  = newData(2:2:end);
                eventMarkers = newData(1:2:end-1);
                eventMarkers = cellfun(@(x)str2double(x(3:end)), eventMarkers, 'uniformoutput', false);
                markerInfo = eventValues;
                markerInfo(:,2) = eventMarkers;
                CONF.( [ fieldName ] ) = markerInfo;
            catch
                % old method (much slower)
                for line = sectionArray(iSection) + 1:sectionArray(iSection + 1) - 1
                    splitArray = strfind(raw{line}, '=');
                    CONF.(fieldName)(line - sectionArray(iSection), :) = {raw{line}(splitArray(1) + 1:end) str2double(raw{line}(3:splitArray(1) - 1))};
                end
            end
            
            if ~isfield(CONF, fieldName)
                disp('No event found');
            else
                if any(cellfun(@isempty, CONF.(fieldName)(:,1)))
                    warning('Empty event(s).')
                end                        
                if ~all(1:size(CONF.(fieldName), 1) == [CONF.(fieldName){:, 2}])
                    warning('Marker number discontinuity.')
                end
            end
            
        case 'comment'
            CONF.(fieldName) = raw(sectionArray(iSection) + 1:sectionArray(iSection + 1) - 1);
        otherwise
            fprintf('Unrecognized entry: %s\n', fieldName);
    end
end

% Handle ahdr file type exceptions
[ ~, ~, ext ] = fileparts( fullfile(pathname, filename) );
if strcmp( ext, '.ahdr' )
    
    if ~isfield( CONF.commoninfos, 'numberofchannels' )
        error( 'Common infos field numberofchannels required.' )
    else
        ahdrChan = str2double( CONF.commoninfos.numberofchannels ) + 1;
        CONF.commoninfos.numberofchannels = num2str( ahdrChan );
    end
    
    CONF.channelinfos(ahdrChan) = { 'test,,1,mV' };
    
    CONF.coordinates(ahdrChan) = { '0,0,0' };
    
end

disp('Done.')
