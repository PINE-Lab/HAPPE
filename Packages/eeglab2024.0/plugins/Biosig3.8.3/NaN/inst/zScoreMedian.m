function Z = zScoreMedian(X, DIM)
% zScoreMedian removes the median and standardizes by the 1.483*median absolute deviation
%
% Usage:  Z = zScoreMedian(X, DIM)
% Input:  X  : data
%         DIM: dimension along which z-score should be calculated (1=columns, 2=rows) 
%              (optional, default=first dimension with more than 1 element
% Output: Z  : z-scores

%	Copyright (C) 2003 Patrick Houweling
%	Copyright (C) 2009 Alois Schloegl
%	$Id$
%       This function is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/

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


% input checks
if any(size(X)==0), 
	return; 
end;

% robust moment estimators: 
% - mean: median
% - standard deviation: 1.483 * median absolute deviation (medAbsDev);
%   the factor 1.483 is the ratio of the standard deviation of a normal random variable to its MAD.
if nargin<2,
        [D, M] = medAbsDev(X);
else
        [D, M] = medAbsDev(X, DIM);
end;

% z-score: subtract M and divide by 1.483*D
Z = (X - repmat(M, size(X)./size(M))) ./ repmat(1.483*D, size(X)./size(D));
