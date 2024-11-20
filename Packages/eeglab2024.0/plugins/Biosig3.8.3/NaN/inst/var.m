function y=var(x,opt,DIM,W)
% VAR calculates the variance.
% 
% y = var(x [, opt[, DIM]])
%   calculates the variance in dimension DIM
%   the default DIM is the first non-single dimension
%
% opt   0: normalizes with N-1 [default]
%	1: normalizes with N 
% DIM	dimension
%	1: VAR of columns
%	2: VAR of rows
% 	N: VAR of  N-th dimension 
%	default or []: first DIMENSION, with more than 1 element
% W	weights to compute weighted variance (default: [])
%	if W=[], all weights are 1. 
%	number of elements in W must match size(x,DIM) 
% 
% usage: 
%	var(x)	
%	var(x, opt, DIM)	
%	var(x, [], DIM)	
%	var(x, W, DIM)
%	var(x, opt, DIM, W)	
%
% features:
% - can deal with NaN's (missing values)
% - weighting of data 
% - dimension argument 
% - compatible to Matlab and Octave
%
% see also: MEANSQ, SUMSQ, SUMSKIPNAN, MEAN, RMS, STD,

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
%	Copyright (C) 2000-2003,2006,2009,2010 by Alois Schloegl <alois.schloegl@gmail.com>
%       This is part of the NaN-toolbox for Octave and Matlab 
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/

if nargin<3,
	DIM = []; 
end; 

if nargin==1,
	W = [];
	opt = [];

elseif any(nargin==[2,3]) 
	if (numel(opt)<2),
		W = [];
	else
		W = opt;
		opt = []; 
	end;
elseif (nargin==4) && (numel(opt)<2) && (numel(DIM)<2),
	;
else 
	fprintf(1,'Error VAR: incorrect usage\n');	
	help var; 
	return;
end; 	

if isempty(opt),
	opt = 0;  
end; 	

if isempty(DIM), 
        DIM = find(size(x)>1,1);
        if isempty(DIM), DIM=1; end;
end;

[y,n,ssq] = sumskipnan(x,DIM,W);
if all(ssq(:).*n(:) > 2*(y(:).^2)),
	%% rounding error is neglectable 
	y = ssq - y.*y./n;
else
	%% rounding error is not neglectable
	szx = size(x);
	szy = size(y);
	if length(szy)<length(szx);
        	szy(length(szy)+1:length(szx)) = 1;
	end;
	[y,n] = sumskipnan((x-repmat(y./n,szx./szy)).^2,DIM,W);
end; 

if (opt~=1)
    	n = max(n-1,0);			% in case of n=0 and n=1, the (biased) variance, STD and STE are INF
end;
y = y./n;	% normalize

