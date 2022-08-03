function [D, M] = medAbsDev(X, DIM)
% medAbsDev calculates the median absolute deviation 
%
% Usage:  D = medAbsDev(X, DIM)  
%    or:  [D, M] = medAbsDev(X, DIM)
% Input:  X  : data
%         DIM: dimension along which mad should be calculated (1=columns, 2=rows) 
%               (optional, default=first dimension with more than 1 element
% Output: D  : median absolute deviations
%         M  : medians (optional)


%	Copyright (C) 2003 Patrick Houweling%	Copyright (C) 2009 Alois Schloegl	
%	$Id$
%       This function is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/

%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.


% input checks
if any(size(X)==0),
  	return; 
end;

if nargin<2,
	M = median(X);
else
        M = median(X, DIM);
end;

% median absolute deviation: median of absolute deviations to median
D = median(abs(X - repmat(M, size(X)./size(M))), DIM);