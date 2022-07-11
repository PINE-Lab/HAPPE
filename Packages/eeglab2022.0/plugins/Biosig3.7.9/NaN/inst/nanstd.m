function [y] = nanstd(x,FLAG,DIM)
% NANSTD same as STD but ignores NaN's. 
% NANSTD is OBSOLETE; use NaN/STD instead. NANSTD is included 
%    to fix a bug in alternative implementations and to 
%    provide some compatibility. 
%
% Y = nanstd(x, FLAG, [,DIM])
% 
% x     data
% FLAG  0: [default] normalizes with (N-1), N = sample size
% FLAG  1: normalizes with N, N = sample size
% DIM	dimension
%	1 sum of columns
%	2 sum of rows
%	default or []: first DIMENSION with more than 1 element
% Y	resulting standard deviation
% 
% see also: SUM, SUMSKIPNAN, NANSUM, STD

%    $Id$
%    Copyright (C) 2000-2003,2006,2008,2009,2010 by Alois Schloegl <alois.schloegl@gmail.com>	
%    This is part of the NaN-toolbox. For more details see
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/

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

if nargin<2,
	FLAG = 0; 
end;
	
if nargin<3,
	DIM = []; 
end;
if isempty(FLAG), 
        FLAG = 0; 
end;
if isempty(DIM), 
        DIM = find(size(x)>1,1);
        if isempty(DIM), DIM=1; end;
end;

[y,n,ssq] = sumskipnan(x,DIM);
if all(ssq(:).*n(:) > 2*(y(:).^2)),
	%% rounding error is neglectable 
	y = ssq - y.*y./n;
else
	%% rounding error is not neglectable 
	[y,n] = sumskipnan(center(x,DIM).^2,DIM);
end; 

if (FLAG==1)
        y = sqrt(y./n);	% normalize with N
else
	% default method
        y = sqrt(y./max(n-1,0));	% normalize with N-1
end;


%!assert(nanstd(0),NaN)

