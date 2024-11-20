% eegplugin_legacy() - Plugin to import legacy EGI data formats
%
% Usage:
%   >> eegplugin_legacy(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer] eeglab figure.
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Author: Arnaud Delorme, UCSD, 2018
%
% See also: eeglab()

% Copyright (C) 2018 Arnaud Delorme, UCSD
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
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1.07  USA

function vers = eegplugin_egilegacy(fig, trystrs, catchstrs)

vers = 'egilegacy1.0';
if nargin < 3
    error('eegplugin_egilegacy requires 3 arguments');
end;

% find tools menu
% ---------------
neuro_m = findobj(fig, 'tag', 'import data');

% command to check that the '.source' is present in the EEG structure
% -------------------------------------------------------------------
cb_readegi     = [ 'try, [EEG LASTCOM] = pop_readegi;'      catchstrs.new_and_hist ];
cb_readsegegi  = [ 'try, [EEG LASTCOM] = pop_readsegegi;'   catchstrs.new_and_hist ];
cb_readegiepo  = [ 'try, [EEG LASTCOM] = pop_importegimat;' catchstrs.new_and_hist ];

uimenu( neuro_m, 'Label', 'From Netstation binary simple file'    , 'CallBack', cb_readegi,    'Separator', 'on');
uimenu( neuro_m, 'Label', 'From Multiple seg. Netstation files'   , 'CallBack', cb_readsegegi);
uimenu( neuro_m, 'Label', 'From Netstation Matlab files'          , 'CallBack', cb_readegiepo);
