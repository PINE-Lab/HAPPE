function r = rankcorr(X,Y)
% RANKCORR calculated the rank correlation coefficient.
% This function is replaced by CORRCOEF. 
% Significance test and confidence intervals can be obtained from CORRCOEF, too. 
%
% R = CORRCOEF(X, [Y, ] 'Rank');
%
% The rank correlation   r = corrcoef(ranks(x)). 
% is often confused with Spearman's rank correlation.  
% Spearman's correlation is defined as 
%   r(x,y) = 1-6*sum((ranks(x)-ranks(y)).^2)/(N*(N*N-1))
% The results are different. Here, the former version is implemented. 
%
% see also: CORRCOEF, SPEARMAN, RANKS
%
% REFERENCES:
% [1] http://mathworld.wolfram.com/SpearmanRankCorrelationCoefficient.html
% [2] http://mathworld.wolfram.com/CorrelationCoefficient.html

%    $Id$
%    Copyright (C) 2000-2003 by Alois Schloegl <alois.schloegl@gmail.com>	
%    This function is part of the NaN-toolbox
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


% warning('RANKCORR might become obsolete; use CORRCOEF(ranks(x)) or CORRCOEF(...,''Rank'') instead');

if nargin < 2
        r = corrcoef(ranks(X));
else
        r = corrcoef(ranks(X),ranks(Y));
end