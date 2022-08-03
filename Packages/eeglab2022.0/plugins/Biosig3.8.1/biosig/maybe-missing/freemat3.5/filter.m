%% Copyright (C) 2008 Alois Schloegl
%% $Id: filter.m,v 1.1 2008-01-19 21:55:45 schloegl Exp $
%% This function is part of BioSig http://biosig.sf.net 
%% This function is needed, if your system does not provide it (e.g. FreeMat v3.5)
%%
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2 of the License, or
%% (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

%% [Y,Z]=filter(B,A,X,Z)
%%    The z-transform of Y,B,A and X correspond in the following way. 
%%    Y(z) = B(z)/A(z)*X(z)   

function [Y,Z]=filter(B,A,X,Z)
	
	if nargin<3
		usage('filter(B,A,X)');
		return; 
	elseif nargin<4
		Z = []; 
	end; 

	sx = size(X); 	
	if all(sx>1)
		fprintf('X must be a vector');
		return; 
	end; 	
	[Y,Z]=mvfilter(B(:).',A(:).',X(:).',Z);
	Y = reshape(Y,sx); 
