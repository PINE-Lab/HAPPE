% pop_comments() - edit comments
%
% Usage:
%   >> newcomments = pop_comments( oldcomments);
%   >> newcomments = pop_comments( oldcomments, title, newcomments, concat);
%
% Inputs:
%   oldcomments - old comments (string or cell array of strings)
%   title       - optional window title (string)
%   newcomments - new comments (string or cell array of strings)
%                 to assign (during commandline calls only)
%   concat      - [0|1] 1 concatenate the newcomments to the old one.
%                 Default is 0.
%
% Outputs:
%   newcomments - new comments, string
%
% Note: if new comments are given as input, there are simply
%       converted and returned by the function; otherwise a
%       window pops up.
%
% Example
%  EEG.comments = pop_comments( { 'This is the first line.' ' ' ...
%               'This is the third line.' }, 'Editing');
% EEG.comments = pop_comments(EEG.comments,'','This is the fourth line.",1);
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeglab()

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

% 01-25-02 reformated help & license -ad 
% 03-16-02 text interface editing -sm & ad 

function [newcomments, com] = pop_comments( comments, plottitle, newcomments, concat );

com = '';
if exist('comments') ~=1, comments = '';
elseif iscell(comments), comments = strvcat(comments{:}); 
end

% remove trailing blanks and make multiline
comments = strmultiline( comments, 53);

if nargin < 3
    newcomments = comments;
	try, icadefs;
	catch,
		BACKCOLOR  =  [.8 .8 .8];     
		GUIBUTTONCOLOR   = [.8 .8 .8];    
	end
	fig = figure('menubar', 'none', 'tag', 'comment', 'color', BACKCOLOR, 'userdata', 0, ...
		   'numbertitle', 'off', 'name', 'Read/Enter comments -- pop_comments()');
	pos = get(gca,'position'); % plot relative to current axes
	q = [pos(1) pos(2) 0 0];
	s = [pos(3) pos(4) pos(3) pos(4)]./100;
	if exist('plottitle') ~=1, plottitle = ''; end
	
	h = title(plottitle);
	set(h, 'fontname','Helvetica','fontweight', 'bold', 'interpreter', 'none');
	
	axis off;

	% create the buttons
	% ------------------
  	uicontrol('Parent',fig, ...
  	'Units','Normalized', ...
	'Position', [0 -5 20 10].*s+q, ...
	'backgroundcolor', GUIBUTTONCOLOR, ...
	'string','CANCEL', 'callback', 'close(findobj(''tag'', ''comment''));' );
		
  	uicontrol('Parent',fig, ...
  	'Units','Normalized', ...
	'Position', [80 -5 20 10].*s+q, ...
	'backgroundcolor', GUIBUTTONCOLOR, ...
	'string','SAVE', 'callback', ...
		[ 'set(gcbf, ''userdata'', ' ...
		'get(findobj(''parent'', gcbf, ''tag'', ''edit''), ''string''));' ]);

	%hh = text( q(1), 100*s(2)+q(2), comments, 'tag', 'edit');
	%set( hh, 'editing', 'on', 'verticalalignment', 'top');

    %hh = uicontrol('Parent',fig, ...
  	%'Units','Normalized', ...
  	%'style', 'text', ...
	%'Position', [0 100 105 5].*s+q, ...
	%'string', 'Warning: each blank line must contain at least a ''space'' character', ...
	%'horizontalalignment', 'left', ...
    %'backgroundcolor', BACKCOLOR );

    hh = uicontrol('Parent',fig, ...
  	'Units','Normalized', ...
  	'style', 'edit', ...
  	'tag', 'edit', ... 
	'Position', [0 10 105 85].*s+q, ...
	'string', comments, ...
	'backgroundcolor', [ 1 1 1], ...
	'horizontalalignment', 'left', ...
	'max', 3, ...
	'fontsize', 12);

    % Try to use 'courier' since it has constant character size
    try
      lf = listfonts; % not compatible with Octave
      tmppos = strmatch('Courier', lf);
      if ~isempty(tmppos)
          set(hh, 'fontname', lf{tmppos(1)}, 'fontsize', 10);
      end
    catch, end; 
    
    waitfor(fig, 'userdata');

    % find return mode
    if ~ishghandle(fig), return; end
    tmp = get(fig, 'userdata');
    if ~isempty(tmp) && ischar(tmp)    
        newcomments = tmp; % ok button
    else return;
    end

	close(fig);
else
    if iscell(newcomments)
        newcomments = strvcat(newcomments{:});
    end
    if nargin > 3 && concat == 1
        newcomments = strvcat(comments, newcomments);
    end
    return;
end

I = find( comments(:) == '''');
comments(I) = ' ';  
if nargout > 1
        if ~strcmp( comments, newcomments)
          allsame = 1;
            for index = 1:size(comments, 1)
                if ~strcmp(comments(index,:), newcomments(index,:)), allsame = 0; end
            end
        else
            allsame = 0;
        end
        if allsame && ~isempty(comments)
             com =sprintf('EEG.comments = pop_comments(EEG.comments, '''', %s, 1);', vararg2str(newcomments(index+1:end,:)));
        else 
            com =sprintf('EEG.comments = pop_comments('''', '''', %s);', vararg2str(newcomments));     
        end
end
return;
