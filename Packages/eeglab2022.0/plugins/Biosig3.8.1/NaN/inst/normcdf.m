function p = normcdf(x,m,s)
% NORMCDF returns normal cumulative distribtion function
%
% cdf = normcdf(x,m,s);
%
% Computes the CDF of a the normal distribution 
%    with mean m and standard deviation s
%    default: m=0; s=1;
% x,m,s must be matrices of same size, or any one can be a scalar. 
%
% see also: NORMPDF, NORMINV 

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
z = (x-m)./s;	  % if this line causes an error, input arguments do not fit. 

p = erfc(z/-sqrt(2))/2;

z = (s==0);
p((x<m) & z) = 0;

p((x==m)& z) = 0.5;

p((x>m) & z) = 1;

p(isnan(x) | isnan(m) | isnan(s) | (s<0)) = nan;

%!assert(sum(isnan(normcdf([-inf,-2,-1,-.5,0,.5,1,2,3,inf,nan]',2,0))),1)






