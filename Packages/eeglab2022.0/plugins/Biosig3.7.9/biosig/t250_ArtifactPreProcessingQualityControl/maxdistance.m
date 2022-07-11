function [d,ix] = maxdistance(s,Q1,Q2)
% MAXDISTANCE finds the point with the largest distances from a line 
%   connecting the two endpoints in 2-D space. 
%
%  The function is used to find the onset of the a burst in function
%  BURSTANALYSIS, but can be used for other examples too.
%
% usage: 
%	[d,ix] = maxdistance(...)
%	...    = maxdistance(s)
%	...    = maxdistance(s, Q1)
%	...    = maxdistance(s, Q1, Q2)
%
%	d  normal distance of each point in s
%	ix index of point with largest distance
%
%
% see also: BURSTANALYSIS
%
% REFERENCES:
% [1] http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
% [2] http://www.mathworks.com/matlabcentral/newsreader/view_thread/164048
% [3] http://www.mathworks.de/support/solutions/en/data/1-1BYSR/

%    Copyright (C) 2013 by Alois Schloegl <alois.schloegl@gmail.com>
%    This is part of the BIOSIG-toolbox http://biosig.sf.net/
%
%    BioSig is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    BioSig is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with BioSig.  If not, see <http://www.gnu.org/licenses/>.


if any(size(s))==1,
	s = [[1:length(s)]',s(:)];
end; 
if (nargin<2)
	Q1=s(1,:);
end; 
if (nargin<3)
	Q2=s(end,:);
end; 


P  = [s(:,1)-s(1,1), s(:,2)-s(1,2)];	
dQ = [s(end,1)-s(1,1), s(end,2)-s(1,2)];
D  = norm(dQ);

d  = zeros(size(s,1),1);
for k = 1:length(s)
 	d(k) = abs(det([dQ; P(k,:)]))/D;
end;
[tmp,ix] = max(d); 


