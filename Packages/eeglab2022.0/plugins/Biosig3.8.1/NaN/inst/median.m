function [y]=median(x,DIM)
% MEDIAN data elements, 
% [y]=median(x [,DIM])
%
% DIM	dimension
%	1: median of columns
%	2: median of rows
% 	N: median of  N-th dimension 
%	default or []: first DIMENSION, with more than 1 element
%
% features:
% - can deal with NaN's (missing values)
% - accepts dimension argument like in Matlab in Octave, too. 
% - compatible to Matlab and Octave 
%
% see also: SUMSKIPNAN
 
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

%	$Id$
%	Copyright (C) 2000-2003,2009 by Alois Schloegl <alois.schloegl@gmail.com>	
%       This function is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/

global FLAG_NANS_OCCURED;

% check dimension of x
sz=size(x);

% find the dimension for median
if nargin<2,
        DIM=min(find(sz>1));
        if isempty(DIM), DIM=1; end;
end;

if DIM>length(sz),
        sz = [sz,ones(1,DIM-length(sz))];
end;

D1 = prod(sz(1:DIM-1));
D2 = sz(DIM);
D3 = prod(sz(DIM+1:length(sz)));
D0 = [sz(1:DIM-1),1,sz(DIM+1:length(sz))];
y  = repmat(nan,D0);
flag_MexKthElement = exist('kth_element','file')==3;

for k = 0:D1-1,
for l = 0:D3-1,
        xi = k + l * D1*sz(DIM) + 1 ;
        xo = k + l * D1 + 1;
        t = x(xi+(0:sz(DIM)-1)*D1);
        t = t(~isnan(t));
        n = length(t);

        if n==0,
                y(xo) = nan;
        elseif flag_MexKthElement, 
	        if (D1==1) t = t+0.0; end;	% make sure a real copy (not just a reference to x) is used
        	flag_KthE = 0; % fast kth_element can be used, because t does not contain any NaN and there is need to care about in-place sorting
                if ~rem(n,2),
                        y(xo) = sum( kth_element( double(t), n/2 + [0,1], flag_KthE) ) / 2;
                elseif rem(n,2),
                        y(xo) = kth_element(double(t), (n+1)/2, flag_KthE);
                end;
        else
                t = sort(t);
                if ~rem(n,2),
                        y(xo) = (t(n/2) + t(n/2+1)) / 2;
                elseif rem(n,2),
                        y(xo) = t((n+1)/2);
                end;
        end

        if (n<D2) 
	        FLAG_NANS_OCCURED = 1; 
        end; 
end;
end;


%!assert(median([1,NaN,3,inf,-inf]),2)
%!assert(median([1,NaN,3,inf,4,-inf]),3)

