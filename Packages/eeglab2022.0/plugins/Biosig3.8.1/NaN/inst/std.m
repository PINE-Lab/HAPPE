function [o,v]=std(x,opt,DIM,W)
% STD calculates the standard deviation.
% 
% [y,v] = std(x [, opt[, DIM [, W]]])
% 
% opt   option 
%	0:  normalizes with N-1 [default]
%		provides the square root of best unbiased estimator of the variance
%	1:  normalizes with N, 
%		this provides the square root of the second moment around the mean
% 	otherwise: 
%               best unbiased estimator of the standard deviation (see [1])      
%
% DIM	dimension
% 	N STD of  N-th dimension 
%	default or []: first DIMENSION, with more than 1 element
% W	weights to compute weighted s.d. (default: [])
%	if W=[], all weights are 1. 
%	number of elements in W must match size(x,DIM) 
%
% y	estimated standard deviation
%
% features:
% - provides an unbiased estimation of the S.D. 
% - can deal with NaN's (missing values)
% - weighting of data 
% - dimension argument also in Octave
% - compatible to Matlab and Octave
%
% see also: RMS, SUMSKIPNAN, MEAN, VAR, MEANSQ,
%
%
% References(s):
% [1] http://mathworld.wolfram.com/StandardDeviationDistribution.html


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

if nargin<4,
	W = []; 
end;
if nargin<3,
	DIM = []; 
end;
if isempty(DIM), 
        DIM = find(size(x)>1,1);
        if isempty(DIM), DIM=1; end;
end;


[y,n,ssq] = sumskipnan(x,DIM,W);
if all(ssq(:).*n(:) > 2*(y(:).^2))
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


if nargin<2,
        opt = 0;
end;
if isempty(opt),
        opt = 0;
end;


if opt==0, 
        % square root if the best unbiased estimator of the variance 
        ib = inf;
        o  = sqrt(y./max(n-1,0));	% normalize
        
elseif opt==1, 
	ib = NaN;        
        o  = sqrt(y./n);

else
        % best unbiased estimator of the mean
        if exist('unique','file'), 
		% usually only a few n's differ
                [N,tmp,tix] = unique(n(:));	% compress n and calculate ib(n)
        	ib = sqrt(N/2).*gamma((N-1)./2)./gamma(N./2);	%inverse b(n) [1]
	        ib = ib(reshape(tix,size(y)));	% expand ib to correct size
                
        elseif exist('histo3','file'), 
		% usually only a few n's differ
                [N,tix] = histo3(n(:)); N = N.X;
                ib = sqrt(N/2).*gamma((N-1)./2)./gamma(N./2);	%inverse b(n) [1]
	        ib = ib(reshape(tix,size(y)));	% expand ib to correct size
                
        else	% gamma is called prod(size(n)) times 
                ib = sqrt(n/2).*gamma((n-1)./2)./gamma(n./2);	%inverse b(n) [1]
        end;	
        ib = reshape(ib,size(y));
        o  = sqrt(y./n).*ib;
end;

if nargout>1,
	v = y.*((max(n-1,0)./(n.*n))-1./(n.*ib.*ib)); % variance of the estimated S.D. ??? needs further checks
end;


