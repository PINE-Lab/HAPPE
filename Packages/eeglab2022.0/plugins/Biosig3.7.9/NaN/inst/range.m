function Q=range(Y,DIM)
% RANGE calculates the range of Y 
%  Missing values (encoded as NaN) are ignored. 
%
%  Q = range(Y)
%  Q = range(Y,DIM)
%     returns the range along dimension DIM of sample array Y.
%
%  Q = range(HIS)
%     returns the RANGE from the histogram HIS.
%     HIS must be a HISTOGRAM struct as defined in HISTO2 or HISTO3.
%
% see also: IQR, MAD, HISTO2, HISTO3, PERCENTILE, QUANTILE

%    Copyright (C) 2009-2019 by Alois Schloegl <alois.schloegl@gmail.com>
%    This function is part of the NaN-toolbox
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

if nargin<2,
        DIM = [];
end;
if isempty(DIM),
        DIM = find(size(Y)>1,1);
        if isempty(DIM), DIM = 1; end;
end;


if nargin<1,
	help range

else
        SW = isstruct(Y);
        if SW, SW = isfield(Y,'datatype'); end;
        if SW, SW = strcmp(Y.datatype,'HISTOGRAM'); end;
        if SW,
		Q = repmat(NaN,1,size(Y.H,2)); 
		for k=1:size(Y.H,2);
			t = Y.X(find(Y.H(:,k)>0),min(size(Y.X,2),k));
			Q(1,k) = max(t)-min(t);
		end;
        elseif isnumeric(Y) && nargin==1,
		Q = max(Y) - min(Y); 
        elseif isnumeric(Y) && nargin==2,
		Q = max(Y,[],DIM) - min(Y,[],DIM); 
	else 
		help range
	end; 
end;


%!assert(range([1:5,NaN]) == 4)

