function R = skewness(i,DIM)
% SKEWNESS estimates the skewness 
%
% y = skewness(x,DIM)
%   calculates skewness of x in dimension DIM
%
% DIM	dimension
%	1: STATS of columns
%	2: STATS of rows
%	default or []: first DIMENSION, with more than 1 element
%
% features:
% - can deal with NaN's (missing values)
% - dimension argument 
% - compatible to Matlab and Octave
%
% see also: SUMSKIPNAN, STATISTIC
%
% REFERENCE(S):
% http://mathworld.wolfram.com/

%    $Id$
%    Copyright (C) 2000-2003,2010 by Alois Schloegl <alois.schloegl@gmail.com>
%    This function is part of the NaN-toolbox
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



% check input arguments 

if nargin==1,
        DIM = find(size(i)>1,1);
        if isempty(DIM), DIM=1; end;
end;

[R.SUM,R.N,R.SSQ] = sumskipnan(i,DIM);	% sum

R.MEAN 	= R.SUM./R.N;			% mean 
R.SSQ0	= R.SSQ - real(R.SUM).*real(R.MEAN) - imag(R.SUM).*imag(R.MEAN); % sum square with mean removed

%if flag_implicit_unbiased_estim;    %% ------- unbiased estimates ----------- 
    n1 	= max(R.N-1,0);			% in case of n=0 and n=1, the (biased) variance, STD and SEM are INF
%else
%    n1	= R.N;
%end;

R.VAR	= R.SSQ0./n1;	     		% variance (unbiased) 
R.STD	= sqrt(R.VAR);		     	% standard deviation

i       = i - repmat(R.MEAN,size(i)./size(R.MEAN));
R.CM3 	= sumskipnan(i.^3,DIM)./n1;
%R.CM4 	= sumskipnan(i.^4,DIM)./n1;

R = R.CM3./(R.STD.^3);
%R = R.CM4./(R.VAR.^2)-3;
