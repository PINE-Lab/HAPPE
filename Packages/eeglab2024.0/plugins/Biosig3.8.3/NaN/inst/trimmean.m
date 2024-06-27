function Q=trimmean(Y,p,DIM)
% TRIMMEAN calculates the trimmed mean by removing the fraction of p/2 upper and 
% p/2 lower samples. Missing values (encoded as NaN) are ignored and not taken into account. 
% The same number from the upper and lower values are removed, and is compatible to various
% spreadsheet programs including GNumeric [1], LibreOffice, OpenOffice and MS Excel.
%
%  Q = trimmean(Y,p)
%  Q = trimmean(Y,p,DIM)
%     returns the TRIMMEAN along dimension DIM of sample array Y.
%  If p is a vector, the TRIMMEAN for each p is computed. 
%
% see also: MAD, RANGE, HISTO2, HISTO3, PERCENTILE, QUANTILE
%
% References:
% [1] http://www.fifi.org/doc/gnumeric-doc/html/C/gnumeric-trimmean.html


%	Copyright (C) 2009-2019 by Alois Schloegl <alois.schloegl@ist.ac.at>
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

if nargin<3,
        DIM = [];
end;
if isempty(DIM),
        DIM = find(size(Y)>1,1);
        if isempty(DIM), DIM = 1; end;
end;

if nargin<2,
	help trimmean
        
else
	sz = size(Y);
	if DIM > length(sz),
	        sz = [sz,ones(1,DIM-length(sz))];
	end;

	D1 = prod(sz(1:DIM-1));
	D2 = length(p);
	D3 = prod(sz(DIM+1:length(sz)));
	Q  = repmat(nan,[sz(1:DIM-1),D2,sz(DIM+1:length(sz))]);
	for k = 0:D1-1,
		for l = 0:D3-1,
		        xi = k + l * D1*sz(DIM) + 1 ;
			xo = k + l * D1*D2;
		        t  = Y(xi:D1:xi+D1*sz(DIM)-1);
		        t  = sort(t(~isnan(t)));
		        N  = length(t); 
			for m=1:D2,
				n  = floor(N*p(m)/2);
			        f  = sum(t(1+n:N-n))/(N-2*n);
				Q(xo + 1 + (m-1)*D1) = f;
			end; 
		end;
	end;
end;

%!assert(trimmean([11.4, 17.3, 21.3, 25.9, 40.1],.2),23.2)

