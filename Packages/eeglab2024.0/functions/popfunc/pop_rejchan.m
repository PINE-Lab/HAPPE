% POP_REJCHAN - reject artifacts channels in an EEG dataset using joint 
%                  probability of the recorded electrode.
%
% Usage:
%   >> pop_rejchan( INEEG ) % pop-up interactive window mode
%   >> [EEG, indelec, measure, com] = ...
%		= pop_rejchan( INEEG, 'key', 'val');
%
% Inputs:
%   INEEG      - input dataset
%
% Optional inputs:
%   'elec'     - [n1 n2 ...] electrode number(s) to take into 
%                consideration for rejection
%   'threshold' - [max] absolute threshold or activity probability 
%                 limit(s) (in std. dev.) if norm is 'on'.
%   'measure'  - ['prob'|'kurt'|'spec'] compute probability 'prob', kurtosis 'kurt'
%                or spectrum 'spec' for each channel. Default is 'kurt'.
%   'norm'     - ['on'|'off'] normalize measure above (using trimmed 
%                normalization as described in the function JOINTPROB
%                and REJKURT. Default is 'off'.
%   'precomp'  - [float array] use this array instead of computing the 'prob' 
%                or 'kurt' measures.
%   'freqrange' - [min max] frequency range for spectrum computation.
%                Default is 1 to sampling rate divided by 2. The average
%                of the log spectral power is computed over the frequency 
%                range of interest.
%
% Outputs:
%   OUTEEG    - output dataset with updated joint probability array
%   indelec   - indices of rejected electrodes
%   measure   - measure value for each electrode
%   com       - executed command
%
% Author: Arnaud Delorme, CERCO, UPS/CNRS, 2008-
%
% See also: JOINTPROB, REJKURT

% Copyright (C) 2008 Arnaud Delorme, CERCO, UPS/CNRS
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

function [EEG, indelec, measure, com] = pop_rejchan( EEG, varargin);

com = '';
indelec = [];
measure = [];
freq    = [];
if nargin < 1
   help pop_rejchan;
   return;
end;  

if nargin < 2

	% which set to save
	% -----------------
    cb_select = [ 'if get(gcbo, ''value'') == 3,' ...
                  '   set(findobj(gcbf, ''tag'', ''spec''), ''enable'', ''on'');' ...
                  'else,' ...
                  '   set(findobj(gcbf, ''tag'', ''spec''), ''enable'', ''off'');' ...
                  'end;' ];
    cb_norm = [ 'if get(gcbo, ''value''),' ...
                  '   set(findobj(gcbf, ''tag'', ''normlab''), ''string'', ''Z-score threshold [max] or [min max]'');' ...
                  'else,' ...
                  '   set(findobj(gcbf, ''tag'', ''normlab''), ''string'', ''Absolute threshold [max] or [min max]'');' ...
                  'end;' ];
    uilist = { { 'style' 'text' 'string' 'Electrode (number(s); Ex: 2 4 5)' } ...
               { 'style' 'edit' 'string' ['1:' int2str(EEG.nbchan)] } ...
               { 'style' 'text' 'string' 'Measure to use' } ...
               { 'style' 'popupmenu' 'string' 'Probability|Kurtosis|Spectrum' 'value' 2 'callback' cb_select } ...
               { 'style' 'text' 'string' 'Normalize measure (check=on)' } ...
               { 'style' 'checkbox' 'string' '' 'value' 1 'callback' cb_norm } { } ...
               { 'style' 'text' 'string' 'Z-score threshold [max] or [min max]' 'tag' 'normlab' } ...
               { 'style' 'edit' 'string' '5' } ...
               { 'style' 'text' 'string' 'Spectrum freq. range' 'enable' 'off' 'tag' 'spec' } ...
               { 'style' 'edit' 'string' num2str([1 EEG.srate/2])  'enable' 'off' 'tag' 'spec' } }; % 7/16/2014 Ramon
    geom = { [2 1.3] [2 1.3] [2 0.4 0.9] [2 1.3] [2 1.3] };
    result = inputgui( 'uilist', uilist, 'geometry', geom, 'title', 'Reject channel -- pop_rejchan()', ...
        'helpcom', 'pophelp(''pop_rejchan'')');
    if isempty(result), return; end
    
    options = { 'elec' eval( [ '[' result{1} ']' ] ) 'threshold' str2num(result{4}) };
    if result{3}, 
         options = { options{:} 'norm', 'on' }; 
    else options = { options{:} 'norm', 'off' }; 
    end
    
    if result{2} == 1
        options = { options{:} 'measure', 'prob' };
    elseif result{2} == 2 
        options = { options{:} 'measure', 'kurt' }; 
    else
        options = { options{:} 'measure', 'spec' };
        options = { options{:} 'freqrange', str2num(result{5})};
    end

else
    options = varargin;
end

opt = finputcheck( options, { 'norm'      'string'    { 'on';'off' }       'off';
                              'verbose'   'string'    { 'on';'off' }       'on';
                              'indexonly' 'string'    { 'on';'off' }       'off';
                              'measure'   'string'    { 'prob';'kurt';'spec';'std' }    'kurt';
                              'precomp'   'real'      []                   [];
                              'freqrange' 'real'      []                   [1 EEG.srate/2];
                              'elec'      'integer'   []                   [1:EEG.nbchan];
                              'threshold' 'real'   []                      400 }, 'pop_rejchan');
if ischar(opt), error(opt); end
if strcmpi(opt.verbose, 'on')
    myfprintf(1);
else
    myfprintf(0);
end

% compute the joint probability
% -----------------------------
if strcmpi(opt.norm, 'on')
    normval = 2;
else
    normval = 0;
end
if strcmpi(opt.measure, 'prob')
    myfprintf('Computing probability for channels...\n');
    [ measure, indelec ] = jointprob( reshape(EEG.data(opt.elec,:,:), length(opt.elec), size(EEG.data,2)*size(EEG.data,3)), opt.threshold, opt.precomp, normval);
elseif strcmpi(opt.measure, 'kurt')
    myfprintf('Computing kurtosis for channels...\n');
    [ measure, indelec ] = rejkurt( reshape(EEG.data(opt.elec,:,:), length(opt.elec), size(EEG.data,2)*size(EEG.data,3)), opt.threshold, opt.precomp, normval);
elseif strcmpi(opt.measure, 'std')
    % norm is ignored
    measure = std(EEG.data(:,:)');
    if length(opt.threshold) > 1
        indelec = measure < opt.threshold(1) | measure > opt.threshold(end);
    else
        indelec = measure > opt.threshold(1);
    end
elseif strcmpi(opt.measure, 'spec')
    myfprintf('Computing spectrum for channels...\n');
    if ~isempty(opt.precomp)
        measure = opt.precomp;
    else
        [measure, freq] = pop_spectopo(EEG, 1, [], 'EEG' , 'plot','off');
        measure = measure(opt.elec,:); % selecting channels

        % select frequency range
        if ~isempty(opt.freqrange)
            [~, fBeg] = min(abs(freq-opt.freqrange(1)));
            [~, fEnd] = min(abs(freq-opt.freqrange(2)));
            measure = measure(:, fBeg:fEnd);
        end
    end
    
    % consider that data below 20 db has been filtered and remove it
    indFiltered = find(mean(measure) < -20);
    if ~isempty(indFiltered) && indFiltered(1) > 11, 
        measure = measure(:,1:indFiltered(1)-10); 
        myfprintf('Removing spectrum data below -20dB (most likely filtered out)\n'); 
    end
    meanSpec = mean(measure);
    stdSpec  = std( measure);
    
%     for indChan = 1:size(measure,1)
%         if any(measure(indChan,:) > meanSpec+stdSpec*opt.threshold), indelec(indChan) = 1; end
%     end
    if strcmpi(opt.norm, 'on')
        measure1  = max(bsxfun(@rdivide, bsxfun(@minus, measure, meanSpec), stdSpec),[],2);
        if length(opt.threshold) > 1
            measure2 = min(bsxfun(@rdivide, bsxfun(@minus, measure, meanSpec), stdSpec),[],2);
            indelec = measure2 < opt.threshold(1) | measure1 > opt.threshold(end);
            myfprintf('Selecting minimum and maximum normalized power over the frequency range\n');
        else
            indelec = measure1 > opt.threshold(1);
            myfprintf('Selecting maximum normalized power over the frequency range\n');
        end
    else
        measure1 = max(measure,[],2);
        if length(opt.threshold) > 1
            measure2 = min(measure,[],2);
            indelec = measure2 < opt.threshold(1) | measure1 > opt.threshold(end);
            myfprintf('Selecting minimum and maximum power over the frequency range\n');
        else
            indelec = measure1 > opt.threshold(1);
            myfprintf('Selecting maximum power over the frequency range\n');
        end
    end
    measure = measure1;
end
colors = cell(1,length(opt.elec)); colors(:) = { 'k' };
colors(find(indelec)) = { 'r' }; colors = colors(end:-1:1);
myfprintf('%d electrodes labeled for rejection\n', length(find(indelec)));

% output variables
indelec = find(indelec)';
tmpchanlocs = EEG.chanlocs;
if ~isempty(EEG.chanlocs), tmplocs = EEG.chanlocs(opt.elec); tmpelec = { tmpchanlocs(opt.elec).labels }';
else                       tmplocs = []; tmpelec = mattocell([opt.elec]'); % tmpelec = mattocell([1:EEG.nbchan]');%Ramon on 8/7/2014
end
if exist('measure2', 'var')
    myfprintf('#\tElec.\t[min]\t[max]\n');
    tmpelec(:,3) = mattocell(measure2);
    tmpelec(:,4) = mattocell(measure);
else
    myfprintf('#\tElec.\tMeasure\n');
    tmpelec(:,3) = mattocell(measure);
end
tmpelec(:,2) = tmpelec(:,1);
tmpelec(:,1) = mattocell([1:length(measure)]');
for index = 1:size(tmpelec,1)
    if exist('measure2', 'var')
         myfprintf('%d\t%s\t%3.2f\t%3.2f', tmpelec{index,1}, tmpelec{index,2}, tmpelec{index,3}, tmpelec{index,4});
    elseif  ~isempty(EEG.chanlocs)
        myfprintf('%d\t%s\t%3.2f'       , tmpelec{index,1}, tmpelec{index,2}, tmpelec{index,3});
    else % Ramon on 8/7/2014
        myfprintf('%d\t%d\t%3.2f'       , tmpelec{index,1}, tmpelec{index,2}, tmpelec{index,3});
    end
    if any(indelec == index), myfprintf('\t*Bad*\n');
    else                      myfprintf('\n');
    end
end
if isempty(indelec), return; end

com = sprintf('EEG = pop_rejchan(EEG, %s);', vararg2str(options));
if nargin < 2
    tmpcom = [ 'EEGTMP = pop_select(EEG, ''nochannel'', [' num2str(opt.elec(indelec)) ']);' ];
    tmpcom = [ tmpcom ...
            'LASTCOM = ' vararg2str(com) ';' ...
            '[ALLEEG EEG CURRENTSET tmpcom] = pop_newset(ALLEEG, EEGTMP, CURRENTSET);' ...
            '   if ~isempty(tmpcom),' ... 
            '     EEG = eegh(LASTCOM, EEG);' ...
            '     eegh(tmpcom);' ...
            '     eeglab(''redraw'');' ...
            '  end; clear EEGTMP tmpcom;' ];
    eegplot(EEG.data(opt.elec,:,:), 'srate', EEG.srate, 'title', 'Scroll component activities -- eegplot()', ...
			 'limits', [EEG.xmin EEG.xmax]*1000, 'color', colors(end:-1:1), 'eloc_file', tmplocs, 'command', tmpcom);
else
    if ~strcmpi(opt.indexonly, 'on')
        EEG = pop_select(EEG, 'nochannel', opt.elec(indelec));
    end
end

function myfprintf(varargin)
persistent verbose;

if isnumeric(varargin{1})
    verbose = varargin{1};
    return;
end
if verbose
    fprintf(varargin{:});
end

return;
