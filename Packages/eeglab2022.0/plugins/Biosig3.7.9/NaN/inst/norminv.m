function x = norminv(p,m,s)
% NORMINV returns inverse cumulative function of the normal distribution
%
% x = norminv(p,m,s);
%
% Computes the quantile (inverse of the CDF) of a the normal 
%    cumulative distribution with mean m and standard deviation s
%    default: m=0; s=1;
% p,m,s must be matrices of same size, or any one can be a scalar. 
%
% see also: NORMPDF, NORMCDF 

% Reference(s):

%    $Id$
%    Copyright (C) 2000-2003,2010,2011 by Alois Schloegl <alois.schloegl@gmail.com>	
%    This function is part of the NaN-toolbox
%    http://pub.ist.ac.at/~schloegl/matlab/NaN/

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

if nargin==1,
        m=0; s=1;
elseif nargin==2,
        s=1;
end;        

% allocate output memory and check size of arguments
x = sqrt(2)*erfinv(2*p - 1).*s + m;  % if this line causes an error, input arguments do not fit. 

x((p>1) | (p<0) | isnan(p) | isnan(m) | isnan(s) | (s<0)) = nan;

k = (s==0) & ~isnan(m);		% temporary variable, reduces number of tests.

x((p==0) & k) = -inf;

x((p==1) & k) = +inf;

k = (p>0) & (p<1) & k;
if numel(m)==1,
        x(k) = m;
else
        x(k) = m(k);
end;


%!assert(sum(~isnan(norminv([-inf,-.2,0,.2,.5,1,2,inf,nan],2,0))),4)


