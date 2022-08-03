function y=trimean(x,DIM)
% TRIMEAN yields the weighted mean of the median and the quartiles
%    m = TRIMEAN(y).
%
% The trimean is  m = (Q1+2*MED+Q3)/4
%    with quartile Q1 and Q3 and median MED   
%
% N-dimensional data is supported
% 
% REFERENCES:
% [1] http://mathworld.wolfram.com/Trimean.html


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

%	Copyright (C) 1996-2019 by Alois Schloegl <alois.schloegl@gmail.com>	
%       This function is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/

global FLAG_NANS_OCCURED;

% check dimension
sz=size(x);

% find the dimension
if nargin==1,
        DIM = find(size(x)>1,1);
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
q  = repmat(nan,3,1);
for k = 0:D1-1,
for l = 0:D3-1,
        xi = k + l * D1*sz(DIM) + 1 ;
        xo = k + l * D1 + 1;
        t = x(xi+(0:sz(DIM)-1)*D1);
        t = sort(t(~isnan(t)));
        t = t(:);
	n = length(t); 
        if (n<D2) 
	        FLAG_NANS_OCCURED = 1; 
        end; 
        
        % q = flix(t,x); 	% The following find the quartiles and median.
        			% INTERP1 is not an alternative since it fails for n<2;
        x  = n*[0.25;0.50;0.75] + [0.75;0.50;0.25]; 
	d  = x - floor(x);	% distance to next sample	 

        t  = t(:);
	ix = ~logical(d);     	% find integer indices
	q(ix) = t(x(ix)); 	% put integer indices
	ix = ~ix;	     	% find non-integer indices
	q(ix) = t(floor(x(ix))).*(1-d(ix)) + t(ceil(x(ix))).*d(ix);  
        
        y(xo) = (q(1) + 2*q(2) + q(3))/4;
end;
end;

%!assert(trimean([1:5,-Inf,Inf,NaN]) == 3)
%!assert(trimean([1:5,-Inf,Inf,NaN]',1) == 3)
%!assert(trimean([1:5,-Inf,Inf,NaN],2) == 3)

