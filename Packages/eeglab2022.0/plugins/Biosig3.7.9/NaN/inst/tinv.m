function y = tinv(x,n)
% TINV returns inverse cumulative function of the student distribution
%
% x = tinv(p,v);
%
% Computes the quantile (inverse of the CDF) of a the student
%    cumulative distribution with mean m and standard deviation s
% p,v must be matrices of same size, or any one can be a scalar. 
%
% see also: TPDF, TCDF, NORMPDF, NORMCDF, NORMINV 

% Reference(s):

%	$Id$
%	Copyright (C) 2000-2003,2009 by Alois Schloegl <alois.schloegl@gmail.com>	
%       This function is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/

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


% allocate output memory and check size of arguments
if all(size(x)==1)
	x = repmat(x,size(n));
elseif all(size(n)==1)
	n = repmat(n,size(x));
elseif all(size(x)==size(n))
	;	%% OK, do nothing
else
	error('size of input arguments must be equal or scalar')
end; 	

y = norminv(x); % do special cases, like x<=0, x>=1, isnan(x), n > 10000;
y(~(n>0)) = NaN; 

ix = find(~isnan(x) & (n>0) & (n<10000));
if ~isempty(ix)
        y(ix) = (sign(x(ix) - 1/2).*sqrt(n(ix)./betainv(2*min(x(ix), 1-x(ix)), n(ix)/2, 1/2) - n(ix)));
end;

y = reshape(y,size(x));

%!assert(tinv(NaN,4),NaN)
