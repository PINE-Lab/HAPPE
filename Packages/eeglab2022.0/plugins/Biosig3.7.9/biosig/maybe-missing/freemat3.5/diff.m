%% Copyright (C) 2008 Alois Schloegl 
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

%% diff(x [,p[,DIM]])
%%    returns length argument 
%% 
%% diff, needed for FreeMat 3.5

function y = diff(x,p,DIM);

	if nargin<2, p = 1; end;
	if isempty(p) p=1; end; 
	if (p~=1) 
		fprintf(1,'Error: Only 1st order difference (p=1) is supported');
		return;
	end; 	
	
	if nargin<3,
		DIM = []; 
	end;
	if isempty(DIM), 
	        DIM = min(find(size(x)>1));
	        if isempty(DIM), DIM=1; end;
	        if (DIM<1) DIM = 1; 
	end;

	if (DIM==1)
		y= x(2:end,:)-x(1:end-1,:);
	elseif (DIM==2)
		y= x(:,2:end)-x(:,1:end-1);
	end;  			
end


