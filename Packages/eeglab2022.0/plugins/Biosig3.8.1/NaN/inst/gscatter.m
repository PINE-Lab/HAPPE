function [h] = gscatter(x,y,group,clr,sym,siz,doleg,xname,yname)
% GSCATTER scatter plot of groups 
%
%  gscatter(x,y,group)
%  gscatter(x,y,group,clr,sym,siz)
%  gscatter(x,y,group,clr,sym,siz,doleg)
%  gscatter(x,y,group,clr,sym,siz,doleg,xname,yname)
%  h = gscatter(...) 
%
%  x,y, group: 	vectors with equal length 
%  clf: 	color vector, default 'bgrcmyk'
%  sym:		symbol, default '.'
%  siz: 	size of Marker
%  doleg:  'on' (default) shows legend, 'off' turns of legend 
%  xname, yname: name of axis
%
%
% see also: ecdf, cdfplot
%
% References: 

%	Copyright (C) 2009,2019 by Alois Schloegl <alois.schloegl@gmail.com>
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

[b,i,j] = unique(group); 

if nargin<3
	help gscatter; 
	error('invalid number of arguments;')
end;
if nargin<4
	clr = [];
end
if nargin<5
	sym = [];
end
if nargin<6
	siz = [];
end
if nargin<7
	doleg = [];
end
if nargin<8
	xname = [];
end
if nargin<9
	yname = [];
end; 
if isempty(clr), clr='bgrcmyk'; end; 
if isempty(sym), sym='.'; end; 
if isempty(doleg), doleg='on'; end; 

for k=1:length(b);
	%ix = find(k==j);
	c = clr(mod(k-1,length(clr))+1); 
	s = sym(mod(k-1,length(sym))+1); 
	hh(k) = plot(x(k==j),y(k==j),[c,s]);
	if ~isempty(siz)	
		z = siz(mod(k-1,length(siz))+1); 
		set(hh(k),'MarkerSize',z);
	end	
	hold on; 
end; 
hold off; 

if strcmpi(doleg,'off')
	if isnumeric(b)
		b={num2str(b(:))};
	end;	
	legend(b,'box','off');
end;
if ~isempty(xname)
	xlabel(xname); 
end; 
if ~isempty(yname)
	ylabel(yname); 
end; 
	
if nargout>0,
	h = hh; 
end;

