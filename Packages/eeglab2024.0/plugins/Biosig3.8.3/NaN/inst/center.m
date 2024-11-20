function [i,S] = center(i,DIM,W)
% CENTER removes the mean 
%
% [z,mu] = center(x,DIM,W)
%   removes mean x along dimension DIM
%
% x	input data 
% DIM	dimension
%	1: column
%	2: row
%	default or []: first DIMENSION, with more than 1 element
% W	weights to computed weighted mean (default: [], all weights = 1)
%	numel(W) must be equal to size(x,DIM)
%
% features:
% - can deal with NaN's (missing values)
% - weighting of data 
% - dimension argument 
% - compatible to Matlab and Octave
%
% see also: SUMSKIPNAN, MEAN, STD, DETREND, ZSCORE
%
% REFERENCE(S):

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


%    Copyright (C) 2000-2019 by Alois Schloegl <alois.schloegl@gmail.com>
%    This is part of the NaN-toolbox. For more details see
%    http://pub.ist.ac.at/~schloegl/matlab/NaN/


if any(size(i)==0); return; end;

if nargin<3,
	W = [];
end; 
if nargin>1,
        [S,N] = sumskipnan(i,DIM,W);
else
        [S,N] = sumskipnan(i,[],W);
end;

S     = S./N;
szi = size(i);
szs = size(S);
if length(szs)<length(szi);
        szs(length(szs)+1:length(szi)) = 1;
end;
i = i - repmat(S,szi./szs);		% remove mean

%!assert(abs(mean(center([5;NaN;3;NaN;NaN;-1]))) < 2*eps)
