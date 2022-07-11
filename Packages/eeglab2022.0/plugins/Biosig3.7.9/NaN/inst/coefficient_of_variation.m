function cv=coefficient_of_variation(i,DIM)
% COEFFICIENT_OF_VARIATION returns STD(X)/MEAN(X)
% 
% cv=coefficient_of_variation(x [,DIM])
%  cv=std(x)/mean(x) 
%
% see also: SUMSKIPNAN, MEAN, STD
%
%   REFERENCE(S):
%   http://mathworld.wolfram.com/VariationCoefficient.html

%	$Id$
%	Copyright (C) 1997-2003 by Alois Schloegl <alois.schloegl@gmail.com>	
%       This function is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/

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


if nargin<2,
        DIM = find(size(i)>1,1);
        if isempty(DIM), DIM=1; end;
end;

[S,N,SSQ] = sumskipnan(i,DIM);

% sqrt((SSQ-S.*S./N)./max(N-1,0))/(S./N);    % = std(i)/mean(i)

cv = sqrt(SSQ.*N./(S.*S)-1);

%if flag_implicit_unbiased_estim,
	cv = cv.*sqrt(N./max(N-1,0));
%end;
