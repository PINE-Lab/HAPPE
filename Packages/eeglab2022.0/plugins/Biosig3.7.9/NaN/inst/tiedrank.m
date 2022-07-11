function R=tiedrank(X,flag1,flag2)
% TIEDRANK compute rank of samples, the mean value is used in case of ties
%  this function is just a wrapper for RANKS, and provided for compatibility 
%  with the statistics toolbox of matlab(tm)
% 
%    R = tiedrank(X)
%	computes the rank R of vector X
%    
% see also: RANKS


%	Copyright (C) 2009,2010,2017 by Alois Schloegl <alois.schloegl@gmail.com>	
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


if nargin>3,
	error('more than 3 input argument is currently not supported ')
end; 	
if nargin<2,
	flag1=0;
end;
if nargin<3,
	flag2=0;
end;

if nargout>2,
	warning('more than 1 output argument is currently not supported ')
end; 	

if nargin<2,
        DIM = [];
end;
if isempty(DIM),
        DIM = find(size(X)>1,1);
        if isempty(DIM), DIM = 1; end;
end
if (DIM<1), DIM = 1; end; %% Hack, because min([])=0 for FreeMat v3.5

R = ranks(X,DIM); 

