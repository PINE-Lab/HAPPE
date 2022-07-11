% eeg_regepochs() - Convert a continuous dataset into consecutive epochs of 
%                   a specified regular length by adding dummy events of type 
%                   and epoch the data around these events. Alternatively
%                   only insert events for extracting these epochs.
%                   May be used to split up continuous data for
%                   artifact rejection followed by ICA decomposition.
%                   The computed EEG.icaweights and EEG.icasphere matrices
%                   may then be exported to the continuous dataset and/or 
%                   to its epoched descendents.
% Usage:
%     >> EEGout = eeg_regepochs(EEG); % use defaults
%     >> EEGout = eeg_regepochs(EEG, 'key', val, ...); 
%
% Required input:
%     EEG                 - EEG continuous data structure (EEG.trials = 1)
%
% Optional inputs:
%     'recurrence' - [in s] the regular recurrence interval of the dummy
%                    events used as time-locking events for the 
%                    consecutive epochs {default: 1 s}
%     'limits'     - [minsec maxsec] latencies relative to the time-locking
%                    events to use as epoch boundaries. Stated epoch length 
%                    will be reduced by one data point to avoid overlaps 
%                    {default: [0 recurrence_interval]}
%     'rmbase'     - [NaN|latency] remove baseline (s). NaN does not remove
%                    baseline. 0 remove the pre-0 baseline. To
%                    remove the full epoch baseline, enter a value
%                    larger than the upper epoch limit. Default is 0.
%     'eventtype'  - [string] name for the event type. Default is 'X'
%     'extractepochs' - ['on'|'off'] extract data epochs ('on') or simply
%                    insert events ('off'). Default is 'on'.
%
% Outputs:
%     EEGout              - the input EEG structure separated into consecutive 
%                           epochs.
%
% See also: pop_editeventvals(), pop_epoch(), rmbase();
%
% Authors: Hilit Serby, Arnaud Delorme & Scott Makeig, SCCN/INC/UCSD, Sep 02, 2005

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, Sep 02, 2005, hilit@sccn.ucsd.edu
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

function EEG = eeg_regepochs(EEG, varargin)

if nargin < 1
    help eeg_regepochs;
    return;
end

if length(EEG) > 1
    EEG = pop_mergeset(EEG, [1:length(EEG)]);
end

% test input variables
% --------------------
if ~isstruct(EEG) || ~isfield(EEG,'event')
   error('first argument must be an EEG structure')
elseif EEG.trials > 1
   error('input dataset must be continuous data (1 epoch)');
end

if nargin > 1 && ~ischar(varargin{1})
    options = {};
    if nargin >= 2, options = { options{:} 'recurrence' varargin{1} }; end
    if nargin >= 3, options = { options{:} 'limits'     varargin{2} }; end
    if nargin >= 4, options = { options{:} 'rmbase'     varargin{3} }; end
else
    options = varargin;
end
g = finputcheck(options, { 'recurrence'    'real'  []  1;
                            'limits'        'real'  []  [0 1];
                            'rmbase'        'real'  []  0;
                            'eventtype'     'string' {} 'X';
                            'extractepochs' 'string' { 'on','off' } 'on' }, 'eeg_regepochs');
if ischar(g), error(g); end

if g.recurrence < 0 || g.recurrence > EEG.xmax
  error('recurrence_interval out of bounds');
end

if nargin < 3
  g.limits = [0 g.recurrence];
end

if length(g.limits) ~= 2 || g.limits(2) <= g.limits(1) 
   error('epoch limits must be a 2-vector [minsec maxsec]')
end

% calculate number of events to add
% ---------------------------------
bg = 0;        % beginning of data
en = EEG.xmax; % end of data in sec
nu = floor((EEG.xmax+1/EEG.srate)/g.recurrence); % number of type 'X' events to add and epoch on

% bg = EEG.event(1).latency/EEG.srate;   % time in sec of first event
% en = EEG.event(end).latency/EEG.srate; % time in sec of last event
% nu = length((bg+g.recurrence):g.recurrence:(en-g.recurrence)); % number of 'X' events, one every 'g.recurrence' sec

if nu < 1
  error('specified recurrence_interval too long')
end

% print info on commandline
% -------------------------
eplength = g.limits(2)-g.limits(1);
fprintf('The input dataset will be split into %d epochs of %g s\n',nu,eplength);
if (eplength-g.recurrence)/eplength*100 > 0
    fprintf('Epochs will overlap by %2.0f%%.\n',(eplength-g.recurrence)/eplength*100);
end

% insert events and urevents at the end of the current (ur)event tables
% ---------------------------------------------------------------------
fprintf('Inserting %d type ''%s'' events: ', nu, g.eventtype);
nevents = length(EEG.event);
% convert all event types to strings
for k = 1:nevents
    EEG.event(k).type = num2str(EEG.event(k).type);
end

nurevents = length(EEG.urevent);
for k = 1:nu
   if rem(k,40)
      fprintf('.')
   else
      fprintf('%d',k)
   end
   if k==40 || ( k>40 && ~rem(k-40,70))
     fprintf('\n');
   end

   EEG.event(nevents+k).type = g.eventtype;
   EEG.event(nevents+k).latency = g.recurrence*(k-1)*EEG.srate+1;

   EEG.urevent(nurevents+k).type = g.eventtype;
   EEG.urevent(nurevents+k).latency = g.recurrence*(k-1)*EEG.srate+1;
   EEG.event(nevents+k).urevent = nurevents+k;
end
fprintf('\n');

% sort the events based on their latency
% --------------------------------------
fprintf('Sorting the event table.\n');
EEG = pop_editeventvals( EEG, 'sort', {'latency' 0}); 

% split the dataset into epochs
% ------------------------------
if strcmpi(g.extractepochs, 'on')
    fprintf('Splitting the data into %d %2.2f-s epochs\n',nu,eplength); 
    setname = sprintf('%s - %g-s epochs', EEG.setname, g.recurrence);
    EEG = pop_epoch( EEG, { g.eventtype }, g.limits, 'newname', ...
                                      setname, 'epochinfo', 'yes');
    % baseline zero the epochs
    % ------------------------
    if ~isnan(g.rmbase) && g.limits(1) < g.rmbase
        fprintf('Removing the pre %3.2f second baseline mean of each epoch.\n', g.rmbase);
        EEG = pop_rmbase( EEG, [g.limits(1) g.rmbase]*1000);
    end
end;                              

