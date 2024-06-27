function [SE,M]=sem(x,DIM, W)
% SEM calculates the standard error of the mean
% 
% [SE,M] = SEM(x [, DIM [,W]])
%   calculates the standard error (SE) in dimension DIM
%   the default DIM is the first non-single dimension
%   M returns the mean. 
%   Can deal with complex data, too. 
%
% DIM	dimension
%	1: SEM of columns
%	2: SEM of rows
% 	N: SEM of  N-th dimension 
%	default or []: first DIMENSION, with more than 1 element
% W	weights to compute weighted mean and s.d. (default: [])
%	if W=[], all weights are 1. 
%	number of elements in W must match size(x,DIM) 
%
% features:
% - can deal with NaN's (missing values)
% - weighting of data 
% - dimension argument 
% - compatible to Matlab and Octave
%
% see also: SUMSKIPNAN, MEAN, VAR, STD

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

%	Copyright (C) 2000-2019 by Alois Schloegl <alois.schloegl@gmail.com>	
%       This function is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/


if nargin>2,
	[S,N,SSQ] = sumskipnan(x,DIM,W);
elseif nargin>1,
	[S,N,SSQ] = sumskipnan(x,DIM);
else    
	[S,N,SSQ] = sumskipnan(x);
end

M  = S./N;
SE = (SSQ.*N - real(S).^2 - imag(S).^2)./(N.*N.*(N-1)); 
SE(SE<=0) = 0; 	% prevent negative value caused by round-off error  
SE = sqrt(real(SE));

%!assert(sem([1,4,1,NaN])    == 1)
%!assert(sem([1,4,1,NaN]',1) == 1)
%!assert(sem([1,4,1,NaN] ,2) == 1)

