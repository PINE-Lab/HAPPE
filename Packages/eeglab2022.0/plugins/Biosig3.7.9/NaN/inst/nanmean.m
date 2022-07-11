function [o] = nanmean(i,DIM)
% NANMEAN same as SUM but ignores NaN's. 
% NANMEAN is OBSOLETE; use MEAN instead. NANMEAN is included 
%    to provide backward compatibility 
%
% Y = nanmean(x [,DIM])
% 
% DIM	dimension
%	1 sum of columns
%	2 sum of rows
%	default or []: first DIMENSION with more than 1 element
% Y	resulting mean
%
% 
% see also: MEAN, SUMSKIPNAN, NANSUM 

%	$Id$
%    	Copyright (C) 2009 by Alois Schloegl <alois.schloegl@gmail.com>	
%    	This is part of the NaN-toolbox. For more details see
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/
%
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


if nargin>1
        [o,n] = sumskipnan(i,DIM);
else
        [o,n] = sumskipnan(i);
end;
o=o./n;
