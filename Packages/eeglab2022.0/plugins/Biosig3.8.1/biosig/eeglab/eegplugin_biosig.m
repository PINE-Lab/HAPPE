% eegplugin_biosig() - EEGLAB plugin for importing data using BIOSIG Matlab toolbox
%
% Usage:
%   >> eegplugin_biosig(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer]  EEGLAB figure
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks. 
%
% Create a plugin:
%   For more information on how to create an EEGLAB plugin see the
%   help message of eegplugin_besa() or visit http://www.sccn.ucsd.edu/eeglab/contrib.html
%
% Author: Arnaud Delorme, SCCN/INC, 25 Nov. 2003

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2003 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function vers = eegplugin_biosig(fig, trystrs, catchstrs)
    
  vers = 'biosig (rev.>=0.07)';
  if nargin < 3
      error('eegplugin_eepimport requires 3 arguments');
  end;
    
  % add besa folder to path
  % -----------------------
  if exist('sload') ~=2
      p = which('eegplugin_biosig');
      p = p(1:findstr(p,'eegplugin_biosig.m')-1);
      delim = p(end);
      if exist('sload') ~= 2  % detect BIOSIG
          addpath([ p delim '..' delim 't200' ] );
      end;
  end;
  
  % test if plugin already present
  % ------------------------------
  if ~isempty(findobj(fig, 'tag', 'biosig'))
      vers = 'biosig (ERROR: plugin already present)';
  end;
  
  % find import data menu
  % ---------------------
  menu = findobj(fig, 'tag', 'import data');
  
  % menu callbacks
  % --------------
  combdf = [ trystrs.no_check '[EEGTMP LASTCOM] = pop_readbdf;' catchstrs.new_non_empty ]; 
  comedf = [ trystrs.no_check '[EEGTMP LASTCOM] = pop_readedf;' catchstrs.new_non_empty ]; 
  combio = [ trystrs.no_check '[EEGTMP LASTCOM] = pop_biosig;'  catchstrs.new_non_empty ]; 
  
  % create menus if necessary
  % -------------------------
  if exist('sload') == 2  % detect BIOSIG
      pathsload   = which('sload');   pathsload   = pathsload(1:end-7);
      pathloadeeg = which('loadeeg'); pathloadeeg = pathloadeeg(1:end-9);
      if strcmpi(pathsload, pathloadeeg)
          pathloadeeg = which('loadeeg');
          disp('Warning: Deleting BIOSIG function ''loadeeg'' which is in conflict');
          disp('         with the EEGLAB ''loadeeg'' function (Neuroscan)');
          disp('         Note: the ''loadeeg'' BIOSIG function does not perform any processing');
          disp('               and is just an alias to the BIOSIG ''sload'' function');
          try, delete(pathloadeeg);
          catch, 
              disp([ 'Warning: cannot delete file ' pathloadeeg ] );
          end;
      end;
      uimenu( menu, 'Label', 'From Biosemi .BDF file using BIOSIG', 'CallBack', combdf, 'Separator', 'on'); 
      uimenu( menu, 'Label', 'From .EDF file using BIOSIG'        , 'CallBack', comedf, 'tag', 'biosig'); 
      uimenu( menu, 'Label', 'From other formats using BIOSIG'    , 'CallBack', combio); 
  end;
