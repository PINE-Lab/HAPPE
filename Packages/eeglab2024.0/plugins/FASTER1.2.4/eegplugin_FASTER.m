% eegplugin_FASTER() - EEGLAB plugin for using FASTER processing on EEG datasets
%
% Usage:
%   >> eegplugin_FASTER(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer]  EEGLAB figure
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Author: Hugh Nolan, Robert Whelan, Richard Reilly, Trinity College Dublin, 2010

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2010 Hugh Nolan, Robert Whelan, Richard Reilly, Trinity College Dublin, nolanhu@tcd.ie
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

function vers = eegplugin_FASTER(fig, trystrs, catchstrs)

    vers = 'FASTERv1.2.4';
    if nargin < 3
        error('eegplugin_FASTER requires 3 arguments');
    end

    % add folder to path
    % -----------------------
    if ~exist('FASTER_GUI')
        p = which('FASTER_GUI');
        p = p(1:findstr(p,'FASTER_GUI.m')-1);
        addpath(p);
    end

    % find import data menu
    % ---------------------
    menu = findobj(fig, 'Label', 'Tools');

    % menu callbacks
    % --------------
    com_open_GUI = [trystrs.no_check 'EEG=FASTER_GUI(1);' catchstrs.add_to_hist];

    % create menus if necessary
    % -------------------------
    submenu = uimenu( menu, 'Label', 'Process with FASTER', 'CallBack', com_open_GUI, 'Separator', 'on', 'tag', 'FASTER', 'ForegroundColor',[1,0.1,0.2]);
%     set(submenu,'enable','on');
%     c=get(menu,'Children');
%     i_faster=find(c==findobj(fig,'tag','FASTER'));
%     i_study=find(c==findobj(fig,'Label','Load existing study'));
%     c_faster=c(i_faster);
%     if (i_faster > i_study)
%        c(i_study+1:i_faster)=c(i_study:i_faster-1);
%        c(i_study)=c_faster;
%     else
%        c(i_faster:i_study-1)=c(i_faster+1:i_study);
%        c(i_study)=c_faster;
%     end
%     set(menu,'Children',c);
