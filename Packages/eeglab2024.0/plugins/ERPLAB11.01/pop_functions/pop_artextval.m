% PURPOSE  : Marking epochs containing activity above an upper threshold and below a lower threshold
%
% FORMAT   :
%
% >> EEG = pop_artextval(EEG, parameters);
%
% INPUTS   :
%
% EEG           - input dataset
%
% The available parameters are as follows:
%
%        'Twindow' 	- time period (in ms) to apply this tool (start end). Example [-200 800]
%        'Threshold'    - range of amplitude (in uV). e.g  -100 100
%        'Channel' 	- channel(s) to search artifacts.
%        'LowPass' - Apply low pass filter at provided half-amplitude
%                           cutoff (FIR @ 26 filter order). 
%                           Default: -1/Do Not Apply
%        'Flag'         - flag value between 1 to 8 to be marked when an artifact is found.(1 value)
%        'Review'       - open a popup window for scrolling marked epochs.
%
% OUTPUTS  :
%
% EEG           - updated output dataset
%
% EXAMPLE  :
%
% EEG  = pop_artextval( EEG , 'Channel',  1:16, 'Flag',  1, 'Threshold', [ -100 100], 'Twindow', [ -200 798] );%
%
% See also pop_artblink pop_artderiv pop_artdiff pop_artflatline pop_artmwppth pop_artstep artifactmenuGUI.m markartifacts.m
%
% *** This function is part of ERPLAB Toolbox ***
% Author: Javier Lopez-Calderon
% Center for Mind and Brain
% University of California, Davis,
% Davis, CA
% 2009

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

function [EEG, com] = pop_artextval(EEG, varargin)
com = '';
if nargin<1
        help pop_artextval
        return
end
if isobject(EEG) % eegobj
        whenEEGisanObject % calls a script for showing an error window
        return
end
if nargin==1
        serror = erplab_eegscanner(EEG, 'pop_artextval', 2, 0, 1, 2, 1);
        if serror
                return
        end
        prompt = {'Test period (start end) [ms]','Voltage limits[uV] (e.g. -100 100):', 'Channel(s)'};
        dlg_title = 'Extreme Values';
        defx = {[EEG(1).xmin*1000 EEG(1).xmax*1000] [-100 100] [1:EEG(1).nbchan] 30 0 0};
        def = erpworkingmemory('pop_artextval');
        
        if isempty(def)
                def = defx;
        else
                if def{1}(1)<EEG(1).xmin*1000
                        def{1}(1) = single(EEG(1).xmin*1000);
                end
                if def{1}(2)>EEG(1).xmax*1000
                        def{1}(2) = single(EEG(1).xmax*1000);
                end
                def{3} = def{3}(ismember_bc2(def{3},1:EEG(1).nbchan));
        end
        try
                chanlabels = {EEG(1).chanlocs.labels};
        catch
                chanlabels = [];
        end
        
        %
        % Call GUI
        %
        answer = artifactmenuGUI(prompt, dlg_title, def, defx, chanlabels);
        
        if isempty(answer)
                disp('User selected Cancel')
                return
        end
        
        testwindow = answer{1};
        ampth      = answer{2};
        chanArray  = unique_bc2(answer{3}); % avoids repeated channels
        lpfilt      = answer{4};
        lpopt      = answer{5};
        flag       = answer{6};
        viewer    =  answer{end};
        
        if viewer
                viewstr = 'on';
        else
                viewstr = 'off';
        end
        if size(ampth,1)~=1 || size(ampth,2)~=2
                msgboxText{1} =  'ERROR, you have to specify 2 values for voltage limits. E.g. -100 100';
                title = 'ERPLAB: Flag input';
                errorfound(msgboxText, title);
                return
        end
        if ~isempty(find(flag<1 | flag>16, 1))
                msgboxText{1} =  'ERROR, flag cannot be greater than 16 nor lesser than 1';
                title = 'ERPLAB: Flag input';
                errorfound(msgboxText, title);
                return
        end
        erpworkingmemory('pop_artextval', {answer{1} answer{2} answer{3} answer{4} answer{5} answer{6}});
        if length(EEG)==1
                EEG.setname = [EEG.setname '_ar']; %suggest a new name
        end
        %
        % Somersault
        %
        [EEG, com] = pop_artextval(EEG, 'Twindow', testwindow, 'Threshold', ampth,...
                'Channel', chanArray, 'LowPass', lpfilt, 'Flag', flag, 'Review', viewstr, 'History', 'gui');
        return
end

%
% Parsing inputs
%
p = inputParser;
p.FunctionName  = mfilename;
p.CaseSensitive = false;
p.addRequired('EEG');

t1 = single(EEG(1).xmin*1000);
t2 = single(EEG(1).xmax*1000);
p.addParamValue('Twindow', [t1 t2], @isnumeric);
p.addParamValue('Threshold', 100, @isnumeric);
p.addParamValue('Channel', 1:EEG(1).nbchan, @isnumeric);
p.addParamValue('LowPass', -1, @isnumeric); 
p.addParamValue('Flag', 1, @isnumeric);
p.addParamValue('Review', 'off', @ischar); % to open a window with the marked epochs
p.addParamValue('History', 'script', @ischar); % history from scripting

p.parse(EEG, varargin{:});

testwindow =  p.Results.Twindow;
ampth      =  p.Results.Threshold;
chanArray  =  p.Results.Channel; % avoids repeated channels
lpval      =  p.Results.LowPass; 
flag       =  p.Results.Flag;

if strcmpi(p.Results.Review, 'on')% to open a window with the marked epochs
        eprev = 1;
else
        eprev = 0;
end

if ~isempty(find(chanArray<1 | chanArray>EEG(1).nbchan, 1))
        error('ERPLAB says: error at pop_artextval(). Channel indices cannot be greater than EEG.nbchan')
end
if ~isempty(find(flag<1 | flag>16, 1))
        error('ERPLAB says: error at pop_artextval(). Flag cannot be greater than 16 or lesser than 1')
end
if strcmpi(p.Results.History,'implicit')
        shist = 3; % implicit
elseif strcmpi(p.Results.History,'script')
        shist = 2; % script
elseif strcmpi(p.Results.History,'gui')
        shist = 1; % gui
else
        shist = 0; % off
end

%
% process multiple datasets. Updated August 23, 2013 JLC
%
if length(EEG) > 1
        options1 = {'Twindow', p.Results.Twindow, 'Threshold', p.Results.Threshold,...
                'Channel', p.Results.Channel, 'Flag', p.Results.Flag, 'Review', p.Results.Review, 'History', 'gui'};
        [ EEG, com ] = eeg_eval( 'pop_artextval', EEG, 'warning', 'on', 'params', options1);
        return;
end

fs       = EEG.srate;
nch      = length(chanArray);
ntrial   = EEG.trials;

[p1, p2, checkw] = window2sample(EEG, testwindow, fs,chanArray);

if checkw==1
        error('pop_artextval() error: time window cannot be larger than epoch.')
elseif checkw==2
        error('pop_artextval() error: too narrow time window')
end
if nch>EEG.nbchan
        error('Error: pop_artextval() number of tested channels cannot be greater than total.')
end
if isempty(EEG.reject.rejmanual)
        EEG.reject.rejmanual  = zeros(1,ntrial);
        EEG.reject.rejmanualE = zeros(EEG.nbchan, ntrial);
end

interARcounter = zeros(1,ntrial); % internal counter, for statistics
fprintf('channel #\n ');

%
% Tests RT info
%
isRT = 1; % there is RT info by default
if ~isfield(EEG.EVENTLIST.bdf, 'rt')
        isRT = 0;
else
        valid_rt = nnz(~cellfun(@isempty,{EEG.EVENTLIST.bdf.rt}));
        if valid_rt==0
                isRT = 0;
        end
end



%option to apply low-pass prior to artifact detection
if lpval > 1
    %if user elects to low pass data prior to AD, create EEG_lowfilt copy
    EEG_lowfilt = basicfilter(EEG, chanArray ,0, lpval, 26, 1, 0,[]);
    
end


for ch=1:nch
        
        fprintf('%g ',chanArray(ch));
        
        for i=1:ntrial
            
                if lpval >1
                     dataline = EEG_lowfilt.data(chanArray(ch), p1:p2 ,i);
                else
                     dataline = EEG.data(chanArray(ch), p1:p2 ,i);
                end
               
                criteria1 = max(dataline)>= ampth(2);
                criteria2 = min(dataline)<= ampth(1);
                if criteria1 || criteria2
                        interARcounter(i) = 1;      % internal counter, for statistics
                        
                        %
                        % subroutine
                        %
                        % flaf 1 is obligatory
                        [EEG errorm]= markartifacts(EEG, flag, chanArray, ch, i, isRT);
                        if errorm==1
                                error(['ERPLAB: There was not latency at the epoch ' num2str(i)])
                        elseif errorm==2
                                error('ERPLAB: invalid flag (0<=flag<=16)')
                        end
                end
        end
end
fprintf('\n');

% performance
perreject = nnz(interARcounter)/ntrial*100;
fprintf('pop_artextval() rejected a %.1f %% of total trials.\n', perreject);
fprintf('\n');
pop_summary_AR_eeg_detection(EEG, ''); % show table at the command window
EEG = eeg_checkset( EEG );

if eprev==1
        namefig = 'Extreme Values detection view';
        pop_plotepoch4erp(EEG, namefig)
end

skipfields = {'EEG', 'Review', 'History'};
fn  = fieldnames(p.Results);
com = sprintf( '%s  = pop_artextval( %s ', inputname(1), inputname(1));
for q=1:length(fn)
        fn2com = fn{q};
        if ~ismember_bc2(fn2com, skipfields)
                fn2res = p.Results.(fn2com);
                if ~isempty(fn2res)
                        if ischar(fn2res)
                                if ~strcmpi(fn2res,'off')
                                        com = sprintf( '%s, ''%s'', ''%s''', com, fn2com, fn2res);
                                end
                        else
                                if iscell(fn2res)
                                        fn2resstr = vect2colon(cell2mat(fn2res), 'Sort','on');
                                        fnformat = '{%s}';
                                else
                                        fn2resstr = vect2colon(fn2res, 'Sort','on');
                                        fnformat = '%s';
                                end
                                com = sprintf( ['%s, ''%s'', ' fnformat], com, fn2com, fn2resstr);
                        end
                end
        end
end
com = sprintf( '%s );', com);

% get history from script
switch shist
        case 1 % from GUI
                com = sprintf('%s %% GUI: %s', com, datestr(now));
                %fprintf('%%Equivalent command:\n%s\n\n', com);
                displayEquiComERP(com);
        case 2 % from script
                EEG = erphistory(EEG, [], com, 1);
        case 3
                % implicit
                %EEG = erphistory(EEG, [], com, 1);
                %fprintf('%%Equivalent command:\n%s\n\n', erpcom);
        otherwise %off or none
                com = '';
end

%
% Completion statement
%
msg2end
return