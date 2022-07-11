function M=moment(i,p,opt,DIM)
% MOMENT estimates the p-th moment 
% 
% M = moment(x, p [,opt] [,DIM])
% M = moment(H, p [,opt])
%   calculates p-th central moment from data x in dimension DIM
%	of from Histogram H
%
% p	moment of order p
% opt   'ac': absolute 'a' and/or central ('c') moment
%	DEFAULT: '' raw moments are estimated
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
% see also: STD, VAR, SKEWNESS, KURTOSIS, STATISTIC, 
%
% REFERENCE(S):
% http://mathworld.wolfram.com/Moment.html

%    $Id$
%    Copyright (C) 2000-2002,2010 by Alois Schloegl <alois.schloegl@gmail.com>	
%    This functions is part of the NaN-toolbox
%    http://pub.ist.ac.at/~schloegl/matlab/NaN/

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

if nargin==2,
        DIM=[];
	opt=[];        
elseif nargin==3,
        DIM=[];
elseif nargin==4,
        
else
        fprintf('Error MOMENT: invalid number of arguments\n');
        return;
end;

if p<=0;
        fprintf('Error MOMENT: invalid model order p=%f\n',p);
        return;
end;

if isnumeric(opt) || ~isnumeric(DIM),
        tmp = DIM;
        DIM = opt;
	opt = tmp;        
end;
if isempty(opt), 
	opt='r';
end;
if isempty(DIM), 
        DIM = find(size(i)>1,1);
        if isempty(DIM), DIM=1; end;
end;

N = nan;        
if isstruct(i),
    if isfield(i,'HISTOGRAM'),
	sz = size(i.H)./size(i.X);
	X  = repmat(i.X,sz);
	if any(opt=='c'),
		N = sumskipnan(i.H,1);	% N
	        N = max(N-1,0);		% for unbiased estimation 
		S = sumskipnan(i.H.*X,1);	% sum
		X = X - repmat(S./N, size(X)./size(S)); % remove mean
	end;
	if any(opt=='a'),
		X = abs(X);
	end;
    	[M,n] = sumskipnan(X.^p.*i.H,1);
    else
	warning('invalid datatype')		
    end;
else
	if any(opt=='c'),
	    	[S,N] = sumskipnan(i,DIM);	% gemerate N and SUM
	        N = max(N-1,0);			% for unbiased estimation
		i = i - repmat(S./N, size(i)./size(S)); % remove mean
	end;
	if any(opt=='a'),
		i = abs(i);	
	end;
	[M,n] = sumskipnan(i.^p,DIM);
end;

if isnan(N), N=n; end; 
M = M./N;
