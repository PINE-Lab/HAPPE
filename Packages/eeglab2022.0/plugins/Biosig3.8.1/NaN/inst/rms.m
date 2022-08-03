function o=rms(x,DIM,W)
% RMS calculates the root mean square
%   can deal with complex data. 
%
% y = rms(x,DIM,W)
%
% DIM	dimension
%	1 STD of columns
%	2 STD of rows
% 	N STD of  N-th dimension 
%	default or []: first DIMENSION, with more than 1 element
% W	weights to compute weighted s.d. (default: [])
%	if W=[], all weights are 1. 
%	number of elements in W must match size(x,DIM) 
%
% y	estimated standard deviation
%
% features:
% - can deal with NaN's (missing values)
% - weighting of data 
% - dimension argument also in Octave
% - compatible to Matlab and Octave
%
% see also: SUMSKIPNAN, MEAN


%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; If not, see <http://www.gnu.org/licenses/>.


%	Copyright (C) 2000-2019 by Alois Schloegl <alois.schloegl@gmail.com>
%       This function is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/


if nargin<2,
	[o,N,ssq] = sumskipnan(x);
elseif nargin<3
	[o,N,ssq] = sumskipnan(x,DIM);
else
	[o,N,ssq] = sumskipnan(x,DIM,W);
end;

o = sqrt(ssq./N);
   
%!assert(rms([1,5,1,NaN]) == 3)
%!assert(rms([1,5,1,-Inf,NaN]) == Inf)

