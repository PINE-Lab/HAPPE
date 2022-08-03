function [R, tix] = histo4(Y, W)
% HISTO4 calculates histogram of multidimensional data samples 
%   and supports data compression
%
% R = HISTO4(Y)
% R = HISTO4(Y, W)
%	Y    data: on sample per row, each sample has with size(Y,2) elements 
%	W    weights of each sample (default: [])
%	     W = [] indicates that each sample has equal weight	
% 	R is a struct with these fields: 
%       R.X  are the bin-values 
%       R.H  is the frequency of occurence of value X (weighted with W)
%  	R.N  are the total number of samples (or sum of W)
%
% HISTO4 might be useful for data compression, because
% [R,tix] = histo4(Y) 
%     	is the compression step
% R.X(tix,:) 
%  	is the decompression step
%
% The effort (in memory and speed) for compression is O(n*log(n))
% The effort (in memory and speed) for decompression is only O(n)
% 
% see also: HISTO, HISTO2, HISTO3, HISTO4
%

%	Copyright (C) 1996-2019 by Alois Schloegl <alois.schloegl@ist.ac.at>	
%    	This is part of the NaN-toolbox 
%	https://octave.sourceforge.io/nan/index.html
%	https://pub.ist.ac.at/~schloegl/matlab/NaN/
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.


%%%%% check input arguments %%%%%
[yr, yc] = size(Y);
if nargin<2, 
	W = []; 
end; 
if ~isempty(W) && (yr ~= numel(W)),
	error('number of rows of Y does not match number of elements in W');
end; 
R.datatype = 'HISTOGRAM';
if isempty(Y)
	R.N = 0; 
	R.X = zeros(size(Y));
	R.H = [];
	return
end

%%%%% identify all possible X's and generate overall Histogram %%%%%
[Y, idx] = sortrows(Y);

d  = diff(Y,[],1);
ix = any( (~isnan(d) & (d~=0) ) | diff(isnan(Y),[],1), 2);

tmp = [find(ix); yr];
R.X = Y(tmp,:);
if isempty(W)
	R.H = [tmp(1); diff(tmp)];
	R.N = yr;
else
	W   = cumsum(W(idx));
	R.H = [W(tmp(1)); diff(W(tmp))];
	R.N = W(end);
end; 

%%%%% generate inverse index %%%%%
if nargout>1,
        tix = cumsum([1;ix]);	% rank 
        cc  = 1;
        tmp = sum(ix);
	if tmp < 2^8;
                tix = uint8(tix);
                cc = 8/1;
        elseif tmp < 2^16;
                tix = uint16(tix);
                cc = 8/2;
        elseif tmp < 2^32;
                tix = uint32(tix);
                cc = 8/4;
        end;
	[tmp, idx] = sort(idx);        % inverse index
        tix = tix(idx);		% inverse sort rank
        
        R.compressionratio = (prod(size(R.X)) + yr/cc) / (yr*yc);
        R.tix = tix;
end;

%!assert(getfield(histo4([]),'N'), 0)
%!assert(getfield(histo4(1),'N'), 1)
%!assert(getfield(histo4([1;1]),'H'), 2)
%!assert(getfield(histo4([repmat(NaN,4,2),[1;1;1;3]]),'N')==4)
