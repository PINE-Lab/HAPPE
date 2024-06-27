% INPUTGUI - A comprehensive gui automatic builder. This function helps
%              to create GUI very quickly without bothering about the 
%              positions of the elements. After creating a geometry, 
%              elements just place themselves in the predefined 
%              locations. It is especially useful for figures in which 
%              you intend to put text buttons and descriptions.
%
% Usage:
%   >> [ outparam ] = inputgui( 'key1', 'val1', 'key2', 'val2', ... );
%   >> [ outparam userdat strhalt outstruct tags] = ...
%             inputgui( 'key1', 'val1', 'key2', 'val2', ... );
% 
% Inputs:
%   'geom'       - cell array of cell array of integer vector. Each cell
%                  array defines the coordinate of a given input in the 
%                  following manner: { nb_row nb_col [x_topcorner y_topcorner]
%                  [x_bottomcorner y_bottomcorner] };
%   'geometry'   - cell array describing horizontal geometry. This corresponds 
%                  to the supergui function input 'geomhoriz'
%   'geomvert'   - vertical geometry argument, this argument is passed on to
%                  the supergui function
%   'uilist'     - list of uicontrol lists describing elements properties
%                  { { ui1 }, { ui2 }... }, { 'uiX' } being GUI matlab 
%                  uicontrol arguments such as { 'style', 'radiobutton', 
%                  'String', 'hello' }. See Matlab function UICONTROL for details.
%   'helpcom'    - optional help command 
%   'helpbut'    - text for help button
%   'title'      - optional figure title
%   'userdata'   - optional userdata input for the figure
%   'mode'       - ['normal'|'noclose'|'plot' fignumber]. Either wait for
%                  user to press OK or CANCEL ('normal'), return without
%                  closing window input ('noclose'), only draw the gui ('plot')
%                  or process an existing window which number is given as 
%                  input (fignumber). Default is 'normal'.
%   'eval'       - [string] command to evaluate at the end of the creation 
%                  of the GUI but before waiting for user input. 
%   'screenpos'  - see supergui.m help message.
%   'skipline'   - ['on'|'off'] skip a row before the "OK" and "Cancel"
%                  button. Default is 'on'.
%
% Output:
%   outparam   - list of outputs. The function scans all lines and
%                add up an output for each interactive uicontrol, i.e
%                edit box, radio button, checkbox and listbox.
%   userdat    - 'userdata' value of the figure.
%   strhalt    - the function returns when the 'userdata' field of the
%                button with the tag 'ok' is modified. This returns the
%                new value of this field.
%   outstruct  - returns outputs as a structure (only tagged ui controls
%                are considered). The field name of the structure is
%                the tag of the ui and contain the ui value or string.
%   instruct   - returns inputs provided in the same format as 'outstruct'
%                This allow to compare in/outputs more easy.
%   tags       - uicontrols by tags
%
% Note: the function also adds three buttons at the bottom of each 
%       interactive windows: 'CANCEL', 'HELP' (if callback command
%       is provided) and 'OK'.
%
% Example:
%   res = inputgui('geometry', { 1 1 }, 'uilist', ...
%                         { { 'style' 'text' 'string' 'Enter a value' } ...
%                           { 'style' 'edit' 'string' '' } });
%
%   res = inputgui('geom', { {2 1 [0 0] [1 1]} {2 1 [1 0] [1 1]} }, 'uilist', ...
%                         { { 'style' 'text' 'string' 'Enter a value' } ...
%                           { 'style' 'edit' 'string' '' } });
%
% Author: Arnaud Delorme, CNL / Salk Institute, La Jolla, 1 Feb 2002
%
% See also: SUPERGUI, EEGLAB

% Copyright (C) Arnaud Delorme, CNL/Salk Institute, 27 Jan 2002, arno@salk.edu
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

function [result, userdat, strhalt, resstruct, instruct, alltags] = inputgui( varargin);

result = [];
userdat = [];
strhalt = [];
resstruct = [];
if nargin < 2
   help inputgui;
   return;
end	

% decoding input and backward compatibility
% -----------------------------------------
if ischar(varargin{1})
    options = varargin;
else
    options = { 'geometry' 'uilist' 'helpcom' 'title' 'userdata' 'mode' 'geomvert' };
    options = { options{1:length(varargin)}; varargin{:} };
    options = options(:)';
end

% checking inputs
% ---------------
g = finputcheck(options, { 'geom'     'cell'                []      {}; ...
                           'geometry' {'cell','integer'}    []      []; ...
                           'uilist'   'cell'                []      {}; ...
                           'cancel'   'string'              []      'Cancel'; ...
                           'helpcom'  { 'string','cell' }   { [] [] }      ''; ...
                           'title'    'string'              []      ''; ...
                           'eval'     'string'              []      ''; ...
                           'helpbut'  'string'              []      'Help'; ...
                           'skipline' 'string'              { 'on' 'off' } 'on'; ...
                           'addbuttons' 'string'            { 'on' 'off' } 'on'; ...
                           'userdata' ''                    []      []; ...
                           'getresult' 'real'               []      []; ...
                           'minwidth'  'real'               []      200; ...
                           'screenpos' ''                   []      []; ...
                           'mode'     ''                    []      'normal'; ...
                           'geomvert' 'real'                []       [] ...
                          }, 'inputgui');
if ischar(g), error(g); end

if isempty(g.getresult)
    if ischar(g.mode)
        fig = figure('visible', 'off');
        set(fig, 'name', g.title);
        set(fig, 'userdata', g.userdata);
        if ~iscell( g.geometry )
            oldgeom = g.geometry;
            g.geometry = {};
            for row = 1:length(oldgeom)
                g.geometry = { g.geometry{:} ones(1, oldgeom(row)) };
            end
        end
        
        % skip a line
        if strcmpi(g.skipline, 'on')  
            g.geometry = { g.geometry{:} [1] };
            if ~isempty(g.geom)
                for ind = 1:length(g.geom)
                    g.geom{ind}{2} = g.geom{ind}{2}+1; % add one row
                end
                g.geom = { g.geom{:} {1 g.geom{1}{2} [0 g.geom{1}{2}-2] [1 1] } };
            end
            g.uilist   = { g.uilist{:}, {} };
        end
        
        % add buttons
        if strcmpi(g.addbuttons, 'on')
            g.geometry = { g.geometry{:} [1 1 1 1] };
            if ~isempty(g.geom)
                for ind = 1:length(g.geom)
                    g.geom{ind}{2} = g.geom{ind}{2}+1; % add one row
                end
                g.geom = { g.geom{:} ...
                      {4 g.geom{1}{2} [0 g.geom{1}{2}-1] [1 1] }, ... 
                      {4 g.geom{1}{2} [1 g.geom{1}{2}-1] [1 1] }, ... 
                      {4 g.geom{1}{2} [2 g.geom{1}{2}-1] [1 1] }, ...
                      {4 g.geom{1}{2} [3 g.geom{1}{2}-1] [1 1] } };
            end
            if ~isempty(g.helpcom)
                if ~iscell(g.helpcom)
                    g.uilist = { g.uilist{:}, { 'width' 80 'align' 'left' 'Style', 'pushbutton', 'string', g.helpbut, 'tag', 'help', 'callback', g.helpcom } {} };
                else
                    g.uilist = { g.uilist{:}, { 'width' 80 'align' 'left' 'Style', 'pushbutton', 'string', 'Help gui', 'callback', g.helpcom{1} } };
                    g.uilist = { g.uilist{:}, { 'width' 80 'align' 'left' 'Style', 'pushbutton', 'string', 'More help', 'callback', g.helpcom{2} } };
                end
            else
                g.uilist = { g.uilist{:}, {} {} };
            end
            g.uilist = { g.uilist{:}, { 'width' 80 'align' 'right' 'Style', 'pushbutton', 'string', g.cancel, 'tag' 'cancel' 'callback', 'close(gcbf)' } };
            g.uilist = { g.uilist{:}, { 'width' 80 'align' 'right' 'stickto' 'on' 'Style', 'pushbutton', 'tag', 'ok', 'string', 'OK', 'callback', 'set(gcbo, ''userdata'', ''retuninginputui'');' } };
        end
        
        % add the three buttons (CANCEL HELP OK) at the bottom of the GUI
        % ---------------------------------------------------------------
        if ~isempty(g.geom)
            [~, ~, allobj, alltags] = supergui( 'fig', fig, 'minwidth', g.minwidth, 'geom', g.geom, 'uilist', g.uilist, 'screenpos', g.screenpos);
        elseif isempty(g.geomvert)
            [~, ~, allobj, alltags] = supergui( 'fig', fig, 'minwidth', g.minwidth, 'geomhoriz', g.geometry, 'uilist', g.uilist, 'screenpos', g.screenpos);
        else
            if strcmpi(g.skipline, 'on'),  g.geomvert = [g.geomvert(:)' 1]; end
            if strcmpi(g.addbuttons, 'on'),g.geomvert = [g.geomvert(:)' 1]; end
            [~, ~, allobj, alltags] = supergui( 'fig', fig, 'minwidth', g.minwidth, 'geomhoriz', g.geometry, 'uilist', g.uilist, 'screenpos', g.screenpos, 'geomvert', g.geomvert(:)' );
        end
    else 
        fig = g.mode;
        set(findobj('parent', fig, 'tag', 'ok'), 'userdata', []);
        allobj = findobj('parent',fig);
        allobj = allobj(end:-1:1);
    end

    % evaluate command before waiting?
    % --------------------------------
    if ~isempty(g.eval), eval(g.eval); end
    instruct = outstruct(allobj); % Getting default values in the GUI. 

    % create figure and wait for return
    % ---------------------------------
    if ischar(g.mode) && (strcmpi(g.mode, 'plot') || strcmpi(g.mode, 'return') )
        if strcmpi(g.mode, 'plot')
           return; % only plot and returns
        end
    else 
        waitfor( findobj('parent', fig, 'tag', 'ok'), 'userdata');
    end
else
    fig = g.getresult;
    allobj = findobj('parent',fig);
    allobj = allobj(end:-1:1);
end

result    = {};
userdat   = [];
strhalt   = '';
resstruct = [];

if ~(ishandle(fig)), return; end % Check if figure still exist
 
% output parameters
% -----------------
strhalt = get(findobj('parent', fig, 'tag', 'ok'), 'userdata');
[resstruct,result] = outstruct(allobj); % Output parameters  
userdat = get(fig, 'userdata');

if isempty(g.getresult) && ischar(g.mode) && ( strcmp(g.mode, 'normal') || strcmp(g.mode, 'return') )
	close(fig);
end
drawnow; % for windows

% function for gui res (deprecated)
% --------------------
% function g = myguihandles(fig, g)
% 	h = findobj('parent', fig);
%         if ~isempty(get(h(index), 'tag'))
% 			try, 
% 				switch get(h(index), 'style')
% 				 case 'edit', g = setfield(g, get(h(index), 'tag'), get(h(index), 'string'));
% 				 case { 'value' 'radio' 'checkbox' 'listbox' 'popupmenu' 'radiobutton'  }, ...
% 					  g = setfield(g, get(h(index), 'tag'), get(h(index), 'value'));
%                 end
% 			catch, end
% 		end

function [resstructout, resultout] = outstruct(allobj)
counter   = 1;
resultout    = {};
resstructout = [];

for index=1:length(allobj)
    if iscell(allobj), currentobj = allobj{index};
    else               currentobj = allobj(index);
    end
    if isnumeric(currentobj) || ~isprop(currentobj,'GetPropertySpecification') % To allow new object handles
        try
            objstyle = get(currentobj, 'style');
            switch lower( objstyle )
                case { 'listbox', 'checkbox', 'radiobutton' 'popupmenu' 'radio' }
                    resultout{counter} = get( currentobj, 'value');
                    if ~isempty(get(currentobj, 'tag')), 
                        try
                            resstructout = setfield(resstructout, get(currentobj, 'tag'), resultout{counter}); 
                        catch
                            fprintf('Warning: tag "%" may not be use as field in output structure', get(currentobj, 'tag'));
                        end
                    end
                    counter = counter+1;
                case 'edit'
                    resultout{counter} = get( currentobj, 'string');
                    if ~isempty(get(currentobj, 'tag')), 
                        try
                            resstructout = setfield(resstructout, get(currentobj, 'tag'), resultout{counter});                         
                        catch
                            fprintf('Warning: tag "%" may not be use as field in output structure', get(currentobj, 'tag'));
                        end
                    end
                    counter = counter+1;
            end
        catch, end
    else
        ps              = currentobj.GetPropertySpecification;
        resultout{counter} = arg_tovals(ps,false);
        count = 1;
        while isfield(resstructout, ['propgrid' int2str(count)])
            count = count + 1;
        end
        resstructout = setfield(resstructout, ['propgrid' int2str(count)], arg_tovals(ps,false));
    end
end
