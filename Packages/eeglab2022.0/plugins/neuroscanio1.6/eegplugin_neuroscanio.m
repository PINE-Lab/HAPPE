% eegplugin_neuroscanio() - Plugin to import Neuroscan EEG data format
%
% Usage:
%   >> eegplugin_neuroscanio(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer] eeglab figure.
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Author: Arnaud Delorme, UCSD, 2018
%
% See also: eeglab()

% Copyright (C) 2019 Arnaud Delorme, UCSD
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

function vers = eegplugin_neuroscanio(fig, trystrs, catchstrs)

vers = 'neuroscanio1.6';
if nargin < 3
    error('eegplugin_neuroscanio requires 3 arguments');
end

% find tools menu
% ---------------
neuro_m = findobj(fig, 'tag', 'import data');
event_m = findobj(fig, 'tag', 'import event');
epoch_m = findobj(fig, 'tag', 'import epoch');

% command to check that the '.source' is present in the EEG structure
% -------------------------------------------------------------------
cb_loadcnt     = [ 'try, [EEG LASTCOM] = pop_loadcnt;'      catchstrs.new_and_hist ]; 
cb_loadeeg     = [ 'try, [EEG LASTCOM] = pop_loadeeg;'      catchstrs.new_and_hist ]; 
cb_loaddat     = [ trystrs.check_epoch '[EEG LASTCOM]= pop_loaddat(EEG);'    catchstrs.store_and_hist ]; 
cb_importev2   = [ trystrs.check_data  '[EEG LASTCOM]= pop_importev2(EEG);'  catchstrs.store_and_hist ]; 

uimenu( neuro_m, 'Label', 'From Neuroscan .CNT file', 'CallBack', cb_loadcnt,    'Separator', 'on');
uimenu( neuro_m, 'Label', 'From Neuroscan .EEG file', 'CallBack', cb_loadeeg);

uimenu( epoch_m, 'Label', 'From Neuroscan .DAT file', 'CallBack', cb_loaddat);
uimenu( event_m, 'Label', 'From Neuroscan .ev2 file', 'CallBack', cb_importev2);


