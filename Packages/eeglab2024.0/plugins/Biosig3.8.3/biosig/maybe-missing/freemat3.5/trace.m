%% Copyright (C) 2008 Alois Schloegl 
%%
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 3 of the License, or
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

%% trace(f)
%%    trace of matrix f 

function y = trace(x);
  	if ndims(x)>2,
  		error('input argument must be a matrix');
	elseif isempty(x)
		y = 0;
  	elseif any(size(x)==1)
    		y = x(1);
  	else
    		y = sum (diag (x));
  	end;
  	
  	

