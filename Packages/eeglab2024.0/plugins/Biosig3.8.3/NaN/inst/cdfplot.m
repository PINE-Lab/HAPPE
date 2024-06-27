function [h,stats] = cdfplot(X, varargin)
% CDFPLOT plots empirical commulative distribution function
%
%   cdfplot(X)
%   cdfplot(X, FMT)
%   cdfplot(X, PROPERTY, VALUE,...)
%   h = cdfplot(...)
%   [h,stats] = cdfplot(X)
%
%  X contains the data vector
% 	(matrix data is currently changed to a vector, this might change in future) 
%  FMT,PROPERTY,VALUE 
%	are used for formating; see HELP PLOT for more details 
%  h 	graphics handle to the cdf curve
%  stats 
%	a struct containing various summary statistics including
%	mean, std, median, min, max.
%
% see also: ecdf, median, statistics, hist2res, plot
%
% References: 

%       $Id$
%       Copyright (C) 2009,2010 by Alois Schloegl <alois.schloegl@gmail.com>
%       This function is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
% of the  License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street - Fifth Floor, Boston, MA 02110-1301, USA.


his = histo(X(:));
cdf = cumsum(his.H,1) ./ sum(his.H,1);
ix1 = ceil ([1:2*size(his.X,1)]'/2);  
ix2 = floor([2:2*size(his.X,1)]'/2);  
hh  = plot (his.X(ix1), [0; cdf(ix2)], varargin{:});

if nargout>0,
        h = hh; 
end;
if nargout>1,
        stats = hist2res(his);
        stats.median = quantile(his,.5); 
end;


