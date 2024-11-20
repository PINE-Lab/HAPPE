function Q=iqr(Y,DIM)
% IQR calculates the interquartile range  
%  Missing values (encoded as NaN) are ignored. 
%
%  Q = iqr(Y)
%  Q = iqr(Y,DIM)
%     returns the IQR along dimension DIM of sample array Y.
%
%  Q = iqr(HIS)
%     returns the IQR from the histogram HIS. 
%     HIS must be a HISTOGRAM struct as defined in HISTO2 or HISTO3.
%
% see also: MAD, RANGE, HISTO2, HISTO3, PERCENTILE, QUANTILE


%	Copyright (C) 2009-2019 by Alois Schloegl <alois.schloegl@gmail.com>	
%       This function is part of the NaN-toolbox
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

if nargin<2,
        DIM = [];
end;
if isempty(DIM),
        DIM = min(find(size(Y)>1));
        if isempty(DIM), DIM = 1; end;
end;


if nargin<1,
	help iqr
        
else
	Q = quantile(Y,[1,3]/4,DIM); 
	Q = diff(Q,[],DIM); 
end;

%!assert(iqr([1:5,NaN]) == 2)
%!assert(iqr([1:5,NaN]',1) == 2)
%!assert(iqr([1:5,NaN],2) == 2)

