function [EEG, rejectionWindows] = erplab_deleteTimeSegments(EEG, timeThresholdMS, beforeEventcodeBufferMS, afterEventcodeBufferMS, varargin)
% ERPLAB_DELETETIMESEGMENTSS Deletes data segments between if the length of time between any 2 consecutive event codes (string or number) 
% is greater than a user-specified threshold (msec)
%
% FORMAT:
%
%   EEG = erplab_deleteTimeSegments(EEG, timeThresholdMS, beforeEventcodeBufferMS, afterEventcodeBufferMS, ignoreUseEventcodes, ignoreUseType, ignoreBoundary, displayEEG);
%
%
% INPUT:
%
%   EEG                      - continuous EEG dataset (EEGLAB's EEG struct)
%   timeThresholdMS          - user-specified time threshold between event codes. 
%   beforeEventcodeBufferMS  - time buffer before event code in a block of trials, preserves this data before the event code
%   afterEventCodeBufferMS   - time buffer after event code in a block of trials, preserves this data after the event code
%
%
% OPTIONAL INPUT:
%
%   ignoreUseEventcodes      - (array) event code numbers to use or ignore (Default: [])
%   ignoreUseType            - (string) How to interpret the `ignoreUseEventcodes` array (Default: 'ignore')
%                              - 'ignore' - (default) look for time spec between all event codes EXCEPT for the listed eventcodes
%                              - 'use'    - look for time spec between these specific event codes
%   ignoreBoundary           - (0/1) If true (1), adds boundary events
%                               to list of event codes to ignore.
%                               (Default: 0)
%   displayEEG               - (true/false)  - (boolean) Display a plot of the EEG when finished
%
% OUTPUT:
%
%   EEG                      - (EEG-struct) continuous EEG dataset (EEGLAB's EEG struct)
%
%
% EXAMPLE: 
%
%   Delete data segments when there is greater than 3000 ms (3 secs) 
%   in between any consecutive event codes. Do not ignore any eventcodes. 
%   Display EEG plot at the end
%
%   EEG = erplab_deleteTimeSegments(EEG, 3000, 100, 200, [], 'ignore',0,true);   
%
%
%
%
% Requirements:
%   - none
%
% See also ...
%
%
% *** This function is part of ERPLAB Toolbox ***
% Author: Jason Arita
% Center for Mind and Brain
% University of California, Davis,
% Davis, CA
% 2009
% 
% Updated:
% Aaron Matthew Simmons, 2021

%b8d3721ed219e65100184c6b95db209bb8d3721ed219e65100184c6b95db209b
%
% ERPLAB Toolbox
% Copyright � 2007 The Regents of the University of California
% Created by Javier Lopez-Calderon and Steven Luck
% Center for Mind and Brain, University of California, Davis,
% javlopez@ucdavis.edu, sjluck@ucdavis.edu
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


%% Error check the input variables
if nargin<1
    help erplab_deleteTimeSegments
    return
elseif nargin<4
    error('ERPLAB:erplab_deleteTimeSegments: needs 4 inputs.')
elseif length(varargin) > 4                                     % only want 4 optional inputs at most
    error('ERPLAB:erplab_deleteTimeSegments:TooManyInputs', ...
        'requires at most 4 optional inputs');
else
    disp('Working...')
end

% Throw error if empty EEG-dataset passed into input
if(isempty(EEG.data));   error('Input EEG dataset is empty'); end


%% Handle optional variables
optargs = {[], 'ignore', 0, false}; % Default optional inputs

% Put defaults into the valuesToUse cell array,
% and overwrite the ones specified in varargin.
optargs(1:length(varargin)) = varargin;     % if vargin is empty, optargs keep their default values. If vargin is specified then it overwrites

% Place optional args into variable names
[ignoreUseEventcodes, ignoreUseType, ignoreBoundaryEvent, eegplotGUIFeedback] = optargs{:};


%% Convert all timing info to samples & Switch timeperiods 
%(startbuffer = before event code; endbuffer =after event code);
%AMS

timeThresholdSample     = round(timeThresholdMS     *(EEG.srate/1000));  % ms to samples
%startPeriodBufferSample = round(startEventcodeBufferMS *(EEG.srate/1000));  % ms to samples
startPeriodBufferSample = round(afterEventcodeBufferMS  *(EEG.srate/1000));  % ms to samples
%endPeriodBufferSample   = round(endEventcodeBufferMS   *(EEG.srate/1000));  % ms to samples
endPeriodBufferSample   = round(beforeEventcodeBufferMS   *(EEG.srate/1000));  % ms to samples


%% Set up WORKING event codes + IGNORED event codes

% Convert events to table
tb_events         = struct2table(EEG.event);

% Convert event codes to strings in order to use setdiff
if(isnumeric(tb_events.type))
    tb_events.type = arrayfun(@num2str, tb_events.type, 'UniformOutput', false);
elseif(iscell(tb_events.type))
    tb_events.type = cellfun(@num2str, tb_events.type, 'UniformOutput', false);
end
    
% Find set of all unique event codes
eventcodes_all    = unique(tb_events.type);



% Convert ignore/use event codes set to strings
if(~isempty(ignoreUseEventcodes))
    ignoreUseEventcodes = textscan(num2str(ignoreUseEventcodes), '%s');
    ignoreUseEventcodes = ignoreUseEventcodes{1};
end

switch lower(ignoreUseType)
    case 'use'
        analyzedEventcodes    = ignoreUseEventcodes;
    case 'ignore'     
        % remove ignored eventcodes from the set of all eventcodes to analyze
        analyzedEventcodes    = setdiff(eventcodes_all, ignoreUseEventcodes);
    otherwise
        error('Unrecognized value for `ignoreUseType`');
end

%ams: add ignoreboundary event
if ignoreBoundaryEvent == 1
    analyzedEventcodes = setdiff(analyzedEventcodes, 'boundary'); 
end

analyzedEventIndices  = ismember(tb_events.type, analyzedEventcodes);
analyzedSamples       = round(tb_events(analyzedEventIndices, :).latency)';                   % Convert event codes to samples


if analyzedSamples(1) ~= 1
    analyzedSamples = [1 analyzedSamples];          % add first time point index
end
if analyzedSamples(end) ~= EEG.pnts
    analyzedSamples = [analyzedSamples EEG.pnts];   % add first time point index
end

lastSample              = 1;
rejectionWindows        = zeros(length(analyzedSamples), 2); % [];


%% Find large segments between time samples
for ii=1:length(analyzedSamples)
    if abs(analyzedSamples(ii)-lastSample)>=timeThresholdSample
        t1          = lastSample;
        t2          = analyzedSamples(ii);
        
        % If at the beginning of the data array, don't add initial buffer
        if t1 == 1
            rejWin  = [t1 ...
                       t2 - endPeriodBufferSample];
        
        % If at the end of the data array, don't add end buffer
        elseif t2 == EEG.pnts
            rejWin  = [t1 + startPeriodBufferSample ...
                       t2 ];
        % else add time buffer to inital and end data points
        else
            rejWindowMin = t1 + startPeriodBufferSample;
            rejWindowMax = t2 - endPeriodBufferSample;
            
            % Test to ensure overlapping buffer windows do not delete
            % the overlapping time segment
            if rejWindowMin < rejWindowMax
                rejWin  = [rejWindowMin rejWindowMax];
            else
                rejWin = [];
            end
        end
        
        rejectionWindows = vertcat(rejectionWindows, rejWin); %#ok<AGROW>
        
    end
    lastSample = analyzedSamples(ii);
end

rejectionWindows(any(rejectionWindows==0,2),:) = []; % trim empty rows


%% Display input EEG plot to user with rejection windows marked
% Needs to occur before actually deleting the time segments in EEG
if eegplotGUIFeedback
    % Plot EEG data with to-be-rejected time windows
    windowChannelMatrix    = zeros(size(rejectionWindows,1),EEG.nbchan, 1);                                % do not mark any channel in EEGPLOT
    windowColorMatrix      = repmat([1 0 0], size(rejectionWindows,1),1);                                  % color matrix for EEGPLOT highlighting
    windowMatrix           = [rejectionWindows windowColorMatrix windowChannelMatrix];   % combined rejection window highlighting for EEGPLOT
    
%     assignin('base', 'rejectionWindows', rejectionWindows); % not sure why this is needed

    eegplotoptions = { ...
        'events',       EEG.event,          ...
        'srate',        EEG.srate,          ...
        'winlength',    20                  ...
        'winrej',       windowMatrix};


    % If channel labels exist then display labels instead of numbers
    if ~isempty(EEG.chanlocs)
        eegplotoptions = [ eegplotoptions {'eloc_file', EEG.chanlocs}];
    end
    
    % Run EEGPLOT
    %        'winrej'     - [start end R G B e1 e2 e3 ...] Matrix giving data periods to mark
    %                      for rejection, each row indicating a different period
    %                        [start end]    = period limits (in frames from beginning of data);
    %                        [R G B]        = specifies the marking color;
    %                        [e1 e2 e3 ...] = a (1,nchans) logical [0|1] vector giving
    %                           channels (1) to mark and (0) not mark for rejection.
    eegplot(EEG.data, eegplotoptions{:});
    
    
    fprintf('\n %g rejection segments marked.\n\n', size(rejectionWindows,1));
end



%% Via EEGLAB.EEG_EEGREJ, delete the rejected windows
rejectionWindowCount = size(rejectionWindows, 1);
if rejectionWindowCount < 1
    fprintf('\nNote: No large segments were found.\n')
else
    
    if(rejectionWindowCount > 1)
        rejectionWindows = joinclosesegments(rejectionWindows, [], 5);
    end
    
    EEG = eeg_eegrej( EEG, rejectionWindows);
      
end

%% Delete first boundary event code when it is the first sample.
if ischar(EEG.event(1).type)
    if strcmpi(EEG.event(1).type,'boundary') && EEG.event(1).latency<=1 % in sample
        EEG = pop_editeventvals(EEG,'delete',1);
    end
end


%% Warn if previously created EVENTLIST detected
if(isfield(EEG, 'EVENTLIST') && ~isempty(EEG.EVENTLIST))
    warning_txt = sprintf('Previously Created ERPLAB EVENTLIST Detected & Deleted \n _________________________________________________________________________\n\n This function changes your event codes, thus your prior eventlist is now obsolete and WILL BE DELETED. \n\n Remember to re-create a new ERPLAB EVENTLIST\n _________________________________________________________________________\n');
    warning(warning_txt); %#ok<SPWRN>
    
    % DELETE PRIOR EVENTLIST
    EEG.EVENTLIST = [];
end





