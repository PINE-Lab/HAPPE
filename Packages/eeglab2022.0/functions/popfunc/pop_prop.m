% pop_prop() - plot the properties of a channel or of an independent
%              component. 
% Usage:
%   >> pop_prop( EEG);           % pops up a query window 
%   >> pop_prop( EEG, typecomp); % pops up a query window 
%   >> pop_prop( EEG, typecomp, chanorcomp, winhandle,spectopo_options);
%
% Inputs:
%   EEG        - EEGLAB dataset structure (see EEGGLOBAL)
%
% Optional inputs:
%   typecomp   - [0|1] 1 -> display channel properties 
%                0 -> component properties {default: 1 = channel}
%   chanorcomp - channel or component number[s] to display {default: 1}
%
%   winhandle  - if this parameter is present or non-NaN, buttons 
%                allowing the rejection of the component are drawn. 
%                If non-zero, this parameter is used to back-propagate
%                the color of the rejection button.
%   spectopo_options - [cell array] optional cell array of options for 
%                the spectopo() function. 
%                For example { 'freqrange' [2 50] }
% 
% Note: Type "set(gcf, ''renderer'', ''painter'')" before saving the figure 
% in postscript (epsc) or jpg format
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: pop_runica(), eeglab()

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% hidden parameter winhandle

% 01-25-02 reformated help & license -ad 
% 02-17-02 removed event index option -ad
% 03-17-02 debugging -ad & sm
% 03-18-02 text settings -ad & sm
% 03-18-02 added title -ad & sm

function com = pop_prop(EEG, typecomp, chanorcomp, winhandle, spec_opt)

com = '';
if nargin < 1
	help pop_prop;
	return;   
end
if nargin < 5
	spec_opt = {};
end
if nargin == 1
	typecomp = 1;    % defaults
        chanorcomp = 1;
end
if typecomp == 0 && isempty(EEG.icaweights)
   error('No ICA weights recorded for this dataset -- first run ICA on it');
end
chanoristr = '';
if nargin == 2
    commandchan = 'tmpchanlocs = EEG(1).chanlocs; [tmp tmpval] = pop_chansel({tmpchanlocs.labels}, ''withindex'', ''on'', ''selectionmode'', ''single''); set(findobj(gcbf, ''tag'', ''chan''), ''string'',tmpval); clear tmp tmpchanlocs tmpval'; 
    
    uitext = { { 'style' 'text' 'string' fastif(typecomp,'Channel index(ices) to plot:','Component index(ices) to plot:') } ...
               { 'style' 'edit' 'string' '1', 'tag', 'chan' } ...
               { 'style' 'pushbutton' 'string'  '...', 'enable' fastif(~isempty(EEG(1).chanlocs) && typecomp, 'on', 'off') 'callback' commandchan } ...               
               { 'style' 'text' 'string' 'Spectral options (see spectopo() help):' } ...
               { 'style' 'edit' 'string' '''freqrange'', [2 50]' } {} };
    uigeom = { [2 1 0.5 ] [2 1 0.5] };
    result = inputgui('geometry', uigeom, 'uilist', uitext, 'helpcom', 'pophelp(''pop_prop'');', ...
							 'title', fastif( typecomp, 'Channel properties - pop_prop()', 'Component properties - pop_prop()'));
	if size( result, 1 ) == 0 return; end
   
    chanoristr = result{1};
    chanorcomp   = eeg_decodechan(EEG.chanlocs, result{1} );
    spec_opt     = eval( [ '{' result{2} '}' ] );
end

if ischar(chanorcomp)
    chanoristr = chanorcomp;
    chanorcomp = eeg_decodechan(EEG.chanlocs, chanorcomp );
elseif isempty(chanoristr)
    chanoristr = int2str(chanorcomp);
end

% plotting several component properties
% -------------------------------------
if length(chanorcomp) > 1
    for index = chanorcomp
        pop_prop(EEG, typecomp, index, 0, spec_opt);  % call recursively for each chanorcomp
    end
	com = sprintf('pop_prop( EEG, %d, [%s], NaN, %s);',...
                  typecomp, int2str(chanorcomp), vararg2str( { spec_opt } ));
    return;
end

% should test for > number of components ??? -sm. 
% yes (max num components is not necessarily same as nbchan). -jri
if typecomp == 1
  maxChanorcomp = EEG.nbchan;
else
  maxChanorcomp = size(EEG.icaweights, 1);
end

if chanorcomp < 1 || chanorcomp > maxChanorcomp 
   error('Component index out of range');
end 

% assumed input is chanorcomp
% -------------------------
try 
    icadefs; 
catch
	BACKCOLOR = [0.8 0.8 0.8];
	GUIBUTTONCOLOR   = [0.8 0.8 0.8]; 
end
basename = [fastif(typecomp,'Channel ', 'Component ') chanoristr ];

fhandle = figure('name', ['pop_prop() - ' basename ' properties'], 'color', BACKCOLOR, 'numbertitle', 'off', 'visible', 'off');
pos     = get(fhandle,'Position');
set(fhandle,'Position', [pos(1) pos(2)-500+pos(4) 500 500], 'visible', 'on');
hh = axes('parent',fhandle);
pos      = get(hh,'position'); % plot relative to current axes
q = [pos(1) pos(2) 0 0];
s = [pos(3) pos(4) pos(3) pos(4)]./100;
axis(hh,'off');

% plotting topoplot
% -----------------
h = axes('parent',fhandle,'Units','Normalized', 'Position',[-10 60 40 42].*s+q);

%topoplot( EEG.icawinv(:,chanorcomp), EEG.chanlocs); axis square; 

if isfield(EEG.chanlocs, 'theta')
    if typecomp == 1 % plot single channel locations
        topoplot( chanorcomp, EEG.chanlocs, 'chaninfo', EEG.chaninfo, ...
                 'electrodes','off', 'style', 'blank', 'emarkersize1chan', 12); axis square;
    else             % plot component map
        topoplot( EEG.icawinv(:,chanorcomp), EEG.chanlocs, 'chaninfo', EEG.chaninfo, ...
                 'shading', 'interp', 'numcontour', 3); axis square;
    end
else
    axis(h,'off');
end
basename = [fastif(typecomp,'Channel ', 'IC') chanoristr ];
% title([ basename fastif(typecomp, ' location', ' map')], 'fontsize', 14); 
title(basename, 'fontsize', 14);

% plotting erpimage
% -----------------
hhh = axes('Parent', fhandle,'Units','Normalized', 'Position',[45 62 48 38].*s+q);
eeglab_options; 
if EEG.trials > 1
    % put title at top of erpimage
    axis(hhh,'off');
    hh = axes('Parent', fhandle,'Units','Normalized', 'Position',[45 62 48 38].*s+q);
    EEG.times = linspace(EEG.xmin, EEG.xmax, EEG.pnts);
    if EEG.trials < 6
      ei_smooth = 1;
    else
      ei_smooth = 3;
    end
    if typecomp == 1 % plot channel
         offset = nan_mean(EEG.data(chanorcomp,:));
         erp=nan_mean(squeeze(EEG.data(chanorcomp,:,:))')-offset;
         erp_limits=get_era_limits(erp);
         erpimage( EEG.data(chanorcomp,:)-offset, ones(1,EEG.trials)*10000, EEG.times*1000, ...
                       '', ei_smooth, 1, 'caxis', 2/3, 'cbar','erp','erp_vltg_ticks',erp_limits);   
    else % plot component
         icaacttmp  = eeg_getdatact(EEG, 'component', chanorcomp);
         offset     = nan_mean(icaacttmp(:));
         era        = nan_mean(squeeze(icaacttmp)')-offset;
         era_limits = get_era_limits(era);
         erpimage( icaacttmp-offset, ones(1,EEG.trials)*10000, EEG.times*1000, ...
                       '', ei_smooth, 1, 'caxis', 2/3, 'cbar','erp', 'yerplabel', '','erp_vltg_ticks',era_limits);   
    end
    axes(hhh);
    title(sprintf('%s activity \\fontsize{10}(global offset %3.3f)', basename, offset), 'fontsize', 14);
else
    % put title at top of erpimage
    EI_TITLE = 'Continous data';
    axis(hhh,'off');
    hh = axes('Parent', fhandle,'Units','Normalized', 'Position',[45 62 48 38].*s+q);
    ERPIMAGELINES = 200; % show 200-line erpimage
    while size(EEG.data,2) < ERPIMAGELINES*EEG.srate
       ERPIMAGELINES = 0.9 * ERPIMAGELINES;
    end
    ERPIMAGELINES = round(ERPIMAGELINES);
    if ERPIMAGELINES > 2   % give up if data too small
        if ERPIMAGELINES < 10
            ei_smooth = 1;
        else
            ei_smooth = 3;
        end
      erpimageframes = floor(size(EEG.data,2)/ERPIMAGELINES);
      erpimageframestot = erpimageframes*ERPIMAGELINES;
      eegtimes = linspace(0, erpimageframes-1, length(erpimageframes)); % 05/27/2014 Ramon: length(erpimageframes) by EEG.srate/1000  in eegtimes = linspace(0, erpimageframes-1, EEG.srate/1000);
      if typecomp == 1 % plot channel
           offset = nan_mean(EEG.data(chanorcomp,:));
           % Note: we don't need to worry about ERP limits, since ERPs
           % aren't visualized for continuous data
           erpimage( reshape(EEG.data(chanorcomp,1:erpimageframestot),erpimageframes,ERPIMAGELINES)-offset, ones(1,ERPIMAGELINES)*10000, eegtimes , ...
                         EI_TITLE, ei_smooth, 1, 'caxis', 2/3, 'cbar');  
      else % plot component
         icaacttmp = eeg_getdatact(EEG, 'component', chanorcomp);
         offset = nan_mean(icaacttmp(:));
         erpimage(reshape(icaacttmp(:,1:erpimageframestot),erpimageframes,ERPIMAGELINES)-offset,ones(1,ERPIMAGELINES)*10000, eegtimes , ...
                    EI_TITLE, ei_smooth, 1, 'caxis', 2/3, 'cbar','yerplabel', '');
      end
    else
            axis(hh,'off');
            text(0.1, 0.3, [ 'No erpimage plotted' 10 'for small continuous data']);
    end
    axes(hhh);
end

% plotting spectrum
% -----------------
if ~exist('winhandle')
    winhandle = NaN;
end
if ishandle(winhandle) 
	h = axes('Parent', fhandle,'units','normalized', 'position',[5 10 95 35].*s+q);
else
	h = axes('Parent', fhandle,'units','normalized', 'position',[5 0 95 40].*s+q);
end
%h = axes('units','normalized', 'position',[45 5 60 40].*s+q);
try
	eeglab_options; 
	if typecomp == 1
		[spectra, freqs] = spectopo( EEG.data(chanorcomp,:), EEG.pnts, EEG.srate, spec_opt{:} );
	else 
		if option_computeica  
			[spectra, freqs] = spectopo( EEG.icaact(chanorcomp,:), EEG.pnts, EEG.srate, 'mapnorm', EEG.icawinv(:,chanorcomp), spec_opt{:} );
        else
    		icaacttmp = (EEG.icaweights(chanorcomp,:)*EEG.icasphere)*reshape(EEG.data(EEG.icachansind,:,:), length(EEG.icachansind), EEG.trials*EEG.pnts); 
			[spectra, freqs] = spectopo( icaacttmp, EEG.pnts, EEG.srate, 'mapnorm', EEG.icawinv(:,chanorcomp), spec_opt{:} );
		end
	end
    % set up new limits
    % -----------------
    %freqslim = 50;
	%set(gca, 'xlim', [0 min(freqslim, EEG.srate/2)]);
    %spectra = spectra(find(freqs <= freqslim));
	%set(gca, 'ylim', [min(spectra) max(spectra)]);
    
	%tmpy = get(gca, 'ylim');
    %set(gca, 'ylim', [max(tmpy(1),-1) tmpy(2)]);
	set( get(h, 'ylabel'), 'string', 'Power 10*log_{10}(\muV^{2}/Hz)', 'fontsize', 14); 
	set( get(h, 'xlabel'), 'string', 'Frequency (Hz)', 'fontsize', 14); 
	title('Activity power spectrum', 'fontsize', 14); 
catch err
	axis off;
    text(0.1, 0.3, [ 'Error: no spectrum plotted' 10 err.message 10]);
%     lasterror
% 	text(0.1, 0.3, [ 'Error: no spectrum plotted' 10 ' make sure you have the ' 10 'signal processing toolbox']);
end
	
% display buttons
% ---------------
if ishandle(winhandle)
	COLREJ = '[1 0.6 0.6]';
	COLACC = '[0.75 1 0.75]';
	% CANCEL button
	% -------------
	h  = uicontrol(fhandle, 'Style', 'pushbutton', 'backgroundcolor', GUIBUTTONCOLOR, 'string', 'Cancel', 'Units','Normalized','Position',[-10 -10 15 6].*s+q, 'callback', 'close(gcf);');

	% VALUE button
	% -------------
	hval  = uicontrol(fhandle, 'Style', 'pushbutton', 'backgroundcolor', GUIBUTTONCOLOR, 'string', 'Values', 'Units','Normalized', 'Position', [15 -10 15 6].*s+q);

	% REJECT button
	% -------------
    if ~isempty(EEG.reject.gcompreject)
    	status = EEG.reject.gcompreject(chanorcomp);
    else
        status = 0;
    end
	hr = uicontrol(fhandle, 'Style', 'pushbutton', 'backgroundcolor', eval(fastif(status,COLREJ,COLACC)), ...
				'string', fastif(status, 'REJECT', 'ACCEPT'), 'Units','Normalized', 'Position', [40 -10 15 6].*s+q, 'userdata', status, 'tag', 'rejstatus');
	command = [ 'set(gcbo, ''userdata'', ~get(gcbo, ''userdata''));' ...
				'if get(gcbo, ''userdata''),' ...
				'     set( gcbo, ''backgroundcolor'',' COLREJ ', ''string'', ''REJECT'');' ...
				'else ' ...
				'     set( gcbo, ''backgroundcolor'',' COLACC ', ''string'', ''ACCEPT'');' ...
				'end;' ];					
	set(hr, 'callback', command); 

	% HELP button
	% -------------
	h  = uicontrol(fhandle, 'Style', 'pushbutton', 'backgroundcolor', GUIBUTTONCOLOR, 'string', 'HELP', 'Units','Normalized', 'Position', [65 -10 15 6].*s+q, 'callback', 'pophelp(''pop_prop'');');

	% OK button
	% ---------
 	command = [ 'global EEG;' ...
 				'tmpstatus = get( findobj(''parent'', gcbf, ''tag'', ''rejstatus''), ''userdata'');' ...
 				'EEG.reject.gcompreject(' num2str(chanorcomp) ') = tmpstatus;' ]; 
	if winhandle ~= 0
        if VERS < 8.04
            command = [ command ...
                sprintf('if tmpstatus set(%3.15f, ''backgroundcolor'', %s); else set(%3.15f, ''backgroundcolor'', %s); end;', ...
                winhandle, COLREJ, winhandle, COLACC)];
        elseif VERS >= 8.04
            command = [ command ...
                sprintf('if tmpstatus set(findobj(''Tag'', ''%s''), ''backgroundcolor'', %s); else set(findobj(''Tag'',''%s''), ''backgroundcolor'', %s); end;', ...
                winhandle.Tag, COLREJ, winhandle.Tag, COLACC)];
        end
    end			
	command = [ command 'close(gcf); clear tmpstatus' ];
	h  = uicontrol(fhandle, 'Style', 'pushbutton', 'string', 'OK', 'backgroundcolor', GUIBUTTONCOLOR, 'Units','Normalized', 'Position',[90 -10 15 6].*s+q, 'callback', command);

	% draw the figure for statistical values
	% --------------------------------------
	index = num2str( chanorcomp );
	command = [ ...
		'figure(''MenuBar'', ''none'', ''name'', ''Statistics of the component'', ''numbertitle'', ''off'');' ...
		'' ...
		'pos = get(gcf,''Position'');' ...
		'set(gcf,''Position'', [pos(1) pos(2) 340 340]);' ...
		'pos = get(gca,''position'');' ...
		'q = [pos(1) pos(2) 0 0];' ...
		's = [pos(3) pos(4) pos(3) pos(4)]./100;' ...
		'axis off;' ...
		''  ...
		'txt1 = sprintf(''(\n' ...
						'Entropy of component activity\t\t%2.2f\n' ...
					    '> Rejection threshold \t\t%2.2f\n\n' ...
					    ' AND                 \t\t\t----\n\n' ...
					    'Kurtosis of component activity\t\t%2.2f\n' ...
					    '> Rejection threshold \t\t%2.2f\n\n' ...
					    ') OR                  \t\t\t----\n\n' ...
					    'Kurtosis distribution \t\t\t%2.2f\n' ...
					    '> Rejection threshold\t\t\t%2.2f\n\n' ...
					    '\n' ...
					    'Current thresholds sujest to %s the component\n\n' ...
					    '(after manually accepting/rejecting the component, you may recalibrate thresholds for future automatic rejection on other datasets)'',' ...
						'EEG.stats.compenta(' index '), EEG.reject.threshentropy, EEG.stats.compkurta(' index '), ' ...
						'EEG.reject.threshkurtact, EEG.stats.compkurtdist(' index '), EEG.reject.threshkurtdist, fastif(EEG.reject.gcompreject(' index '), ''REJECT'', ''ACCEPT''));' ...
		'' ...				
		'uicontrol(gcf, ''Units'',''Normalized'', ''Position'',[-11 4 117 100].*s+q, ''Style'', ''frame'' );' ...
		'uicontrol(gcf, ''Units'',''Normalized'', ''Position'',[-5 5 100 95].*s+q, ''String'', txt1, ''Style'',''text'', ''HorizontalAlignment'', ''left'' );' ...
		'h = uicontrol(gcf, ''Style'', ''pushbutton'', ''string'', ''Close'', ''Units'',''Normalized'', ''Position'', [35 -10 25 10].*s+q, ''callback'', ''close(gcf);'');' ...
		'clear txt1 q s h pos;' ];
	set( hval, 'callback', command); 
	if isempty( EEG.stats.compenta )
		set(hval, 'enable', 'off');
	end
	
	com = sprintf('pop_prop( EEG, %d, %d, 0, %s);', typecomp, chanorcomp, vararg2str( { spec_opt } ) );
else
	com = sprintf('pop_prop( EEG, %d, %d, NaN, %s);', typecomp, chanorcomp, vararg2str( { spec_opt } ) );
end
disp('Note: Type "set(gcf, ''renderer'', ''painter'')" before saving the figure in postscript (epsc) or jpg format');


return;

function out = nan_mean(in)

    nans = find(isnan(in));
    in(nans) = 0;
    sums = sum(in);
    nonnans = ones(size(in));
    nonnans(nans) = 0;
    nonnans = sum(nonnans);
    nononnans = find(nonnans==0);
    nonnans(nononnans) = 1;
    out = sum(in)./nonnans;
    out(nononnans) = NaN;

    
function era_limits=get_era_limits(era)
%function era_limits=get_era_limits(era)
%
% Returns the minimum and maximum value of an event-related
% activation/potential waveform (after rounding according to the order of
% magnitude of the ERA/ERP)
%
% Inputs:
% era - [vector] Event related activation or potential
%
% Output:
% era_limits - [min max] minimum and maximum value of an event-related
% activation/potential waveform (after rounding according to the order of
% magnitude of the ERA/ERP)

mn=min(era);
mx=max(era);
mn=orderofmag(mn)*round(mn/orderofmag(mn));
mx=orderofmag(mx)*round(mx/orderofmag(mx));
era_limits=[mn mx];


function ord=orderofmag(val)
%function ord=orderofmag(val)
%
% Returns the order of magnitude of the value of 'val' in multiples of 10
% (e.g., 10^-1, 10^0, 10^1, 10^2, etc ...)
% used for computing erpimage trial axis tick labels as an alternative for
% plotting sorting variable

val=abs(val);
if val>=1
    ord=1;
    val=floor(val/10);
    while val>=1
        ord=ord*10;
        val=floor(val/10);
    end
    return;
else
    ord=1/10;
    val=val*10;
    while val<1
        ord=ord/10;
        val=val*10;
    end
    return;
end

