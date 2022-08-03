function [y] = harmmean(x,DIM,W)
% HARMMEAN calculates the harmonic mean of data elements. 
% The harmonic mean is the inverse of the mean of the inverse elements.
% 
% 	y = harmmean(x [,DIM [,W]]) is the same as 
% 	y = mean(x,'H' [,DIM [,W]]) 
%
% DIM	dimension
%	1 STD of columns
%	2 STD of rows
%	default or []: first DIMENSION, with more than 1 element
% W	weights to compute weighted mean (default: [])
%	if W=[], all weights are 1. 
%	number of elements in W must match size(x,DIM) 
%
% features:
% - can deal with NaN's (missing values)
% - weighting of data 
% - dimension argument also in Octave
% - compatible to Matlab and Octave
%
% see also: SUMSKIPNAN, MEAN, GEOMEAN
%

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


%    Copyright (C) 2000-2019 by Alois Schloegl <alois.schloegl@gmail.com>
%    This is part of the NaN-toolbox. For more details see
%    http://pub.ist.ac.at/~schloegl/matlab/NaN/


if nargin<2
        DIM=min(find(size(x)>1));
        if isempty(DIM), DIM=1; end;
end;
if nargin<3
	W = [];
end;

[y, n] = sumskipnan(1./x,DIM,W);
y = n./y;

%!assert(harmmean([1,1/4,1,NaN]) == 0.5)
%!assert(harmmean([1,2,0,NaN]) == 0)
%!assert(harmmean([1,2,Inf,NaN]) == 2)
%!assert(harmmean([1,2,-Inf,NaN]) == 2)

