function [y,z] = teager(x,UC,A,cini)
% TEAGER computes the teager-kaiser operaterion [1]
%   
%  y = teager(x)
%  
%
% see also: 
%
% REFERENCE(S):
% [1] J. K. Kaiser, “On a simple algorithm to calculate the energy of
% a signal,” Proc. IEEE ICASSP 90, vol 1, pp. 381-384, 1990.
%

%	$Id$
%	Copyright (C) 2009 by Alois Schloegl <alois.schloegl@gmail.com>
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

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


y = x.*x; 

y(2:end-1)=y(2:end-1)-x(1:end-2).*x(3:end);



