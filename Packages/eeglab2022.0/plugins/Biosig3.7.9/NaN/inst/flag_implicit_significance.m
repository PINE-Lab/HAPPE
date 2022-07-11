function alpha=flag_implicit_significance(i)
% The use of FLAG_IMPLICIT_SIGNIFICANCE is in experimental state. 
% flag_implicit_significance might even become obsolete.
%
% FLAG_IMPLICIT_SIGNIFICANCE sets and gets default alpha (level) of any significance test
% The default alpha-level is stored in the global variable FLAG_implicit_significance
% The idea is that the significance must not be assigned explicitely. 
% This might yield more readable code. 
%
% Choose alpha low enough, because in alpha*100% of the cases, you will 
% reject the Null hypothesis just by change. For this reason, the default
% alpha is 0.01. 
% 
%   flag_implicit_significance(0.01) 
%	sets the alpha-level for the significance test
% 
% alpha = flag_implicit_significance()
% 	gets default alpha
%
% flag_implicit_significance(alpha)
% 	sets default alpha-level
%
% alpha = flag_implicit_significance(alpha)
%	gets and sets alpha 
%
% features:
% - compatible to Matlab and Octave
%
% see also: CORRCOEF, PARTCORRCOEF

%	$Id$
%	Copyright (C) 2000-2002,2009,2010 by Alois Schloegl <alois.schloegl@gmail.com>	
%       This function is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/
%
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; If not, see <http://www.gnu.org/licenses/>.


persistent FLAG_implicit_significance;
DEFAULT_ALPHA = 0.01;

%%% check whether FLAG was already defined 
if ~exist('FLAG_implicit_significance','var'),
	FLAG_implicit_significance = DEFAULT_ALPHA;  % default value 
end;
if isempty(FLAG_implicit_significance),
	FLAG_implicit_significance = DEFAULT_ALPHA;  % default value 
end;

if nargin>0,
        fprintf(2,'Warning: flag_implicit_significance is in an experimental state\n');
        fprintf(2,'It might become obsolete.\n');
        FLAG_implicit_significance = i; 
end;

alpha = FLAG_implicit_significance;
