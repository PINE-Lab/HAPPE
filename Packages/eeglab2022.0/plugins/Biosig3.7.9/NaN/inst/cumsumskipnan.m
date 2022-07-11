function [x] = cumsumskipnan(x, DIM)
% CUMSUMSKIPNAN  Cumulative sum while skiping NaN's. 
% If DIM is omitted, it defaults to the first non-singleton dimension.
% 
% Y = cumsumskipnan(x [,DIM])
% 
% x	input data 	
% DIM	dimension (default: [])
% y	resulting sum
%
% see also: CUMSUM, SUMSKIPNAN


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

%	$Id$
%    	Copyright (C) 2011 by Alois Schloegl <alois.schloegl@gmail.com>	
%       This function is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/


%%% TODO: implement as mex-function

i = isnan(x);
x(i) = 0;

if nargin==2,
	x = cumsum(x,DIM);
	x(i) = NaN;
elseif nargin==1,
	x = cumsum(x);
	x(i) = NaN;
else
	help cumsumskipnan
end; 	



