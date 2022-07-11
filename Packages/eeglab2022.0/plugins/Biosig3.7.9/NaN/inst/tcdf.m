function p = tcdf(x,n)
% TCDF returns student cumulative distribtion function
%
% cdf = tcdf(x,DF);
%
% Computes the CDF of the students distribution 
%    with DF degrees of freedom 
% x,DF must be matrices of same size, or any one can be a scalar. 
%
% see also: NORMCDF, TPDF, TINV 

% Reference(s):

%	Copyright (C) 2000-2003,2009 by Alois Schloegl <alois.schloegl@gmail.com>	
%    	This is part of the NaN-toolbox. For more details see
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


% check size of arguments
if nargin~=2
	error('invalid number of input arguments')
elseif ~isreal(x) || ~isreal(n)
	error('input arguments must be real')
elseif all(size(x)==1)
	x = repmat(x,size(n));
elseif all(size(n)==1)
	n = repmat(n,size(x));
elseif all(size(x)==size(n))
	;	%% OK, do nothing
else
	error('size of input arguments must be equal or scalar')
end; 	
		
% allocate memory
p = zeros(size(x));
p((x==Inf) & (n>0)) = 1;

% workaround for invalid arguments in BETAINC
ix   = isnan(x) | ~(n>0);
p(ix)= NaN;

ix    = (x > -Inf) & (x < Inf) & (n > 0);
p(ix) = betainc (n(ix) ./ (n(ix) + x(ix).^2), n(ix)/2, 1/2) / 2;

ix    = ix & (x>0);
p(ix) = 1 - p(ix);

% shape output
p = reshape(p,size(x));

%!assert(tcdf(NaN,4),NaN)
%!assert(tcdf(inf,3),1)
