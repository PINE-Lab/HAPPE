function [i,m,s] = zscore(i,OPT, DIM, W)
% ZSCORE removes the mean and normalizes data 
% to a variance of 1. Can be used for pre-whitening of data, too. 
%
% [z,mu, sigma] = zscore(x [,OPT [, DIM])
%   z   z-score of x along dimension DIM
%   sigma is the inverse of the standard deviation
%   mu   is the mean of x
%
% The data x can be reconstucted with 
%     x = z*diag(sigma) + repmat(m, size(z)./size(m))  
%     z = x*diag(1./sigma) - repmat(m.*v, size(z)./size(m))  
%
% OPT   option 
%	0:  normalizes with N-1 [default] when computing sigma
%		provides the square root of best unbiased estimator of the variance [1]
%	1:  normalizes with N, when computing sigma
%		this provides the square root of the second moment around the mean
% 	otherwise: 
%               best unbiased estimator of the standard deviation (see [1])      
%
% DIM	dimension
%	1: STATS of columns
%	2: STATS of rows
%	default or []: first DIMENSION, with more than 1 element
%
% see also: SUMSKIPNAN, MEAN, STD, DETREND
%
% REFERENCE(S):
% [1] http://mathworld.wolfram.com/z-Score.html

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

%	Copyright (C) 2000-2019 by Alois Schloegl <alois.schloegl@ist.ac.at>	
%       This function is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/


if any(size(i)==0); return; end;

if nargin<2
        OPT=[]; 
end
if nargin<3
        DIM=[]; 
end
if nargin<4
        W = []; 
end
if ~isempty(OPT) && ~any(OPT==[0,1])
	error('OPT must be 0, 1 or empty.')
end
if isempty(DIM), 
        DIM=min(find(size(i)>1));
        if isempty(DIM), DIM=1; end;
end;


% pre-whitening
[S,N,SSQ] = sumskipnan(i, DIM, W);
m = S./N; 
i = i-repmat(m, size(i)./size(m));  % remove mean
s = std (i, OPT, DIM, W);
s(s==0)=1;
i = i ./ repmat(s,size(i)./size(s)); % scale to var=1

%!assert(isequaln(zscore([1,1,-1,-1,0,NaN]),[1,1,-1,-1,0,NaN]))

