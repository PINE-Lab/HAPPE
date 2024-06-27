% TEST_XPTOPEN tests XPTOPEN  

%	$Id$
%	Copyright (C) 2010 by Alois Schloegl <alois.schloegl@gmail.com>	
%       This function is part of the NaN-toolbox
%       http://biosig-consulting.com/matlab/NaN/

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

%x.c = [-1000,-2,-1,0,1,2,NaN,10,100,1000,10000,1e6,1e7,1e8]';
%y.Y = [1,2,NaN,1]'+10;

if 1, 
X.a = [-2,-0,NaN,10,444,-pi]';%,100,1000,10000,1e6,1e7,1e8]';
X.d = [1,2,NaN,1,Inf,-Inf]'+10;
X.b = {'a','B',' ','*','Z','zzz'}';

fn  = 'test.xpt';
Y   = xptopen(fn,'w',X)
Z   = xptopen(fn,'r')


end;

fn = {'buy','humid','prdsale'};
for k1 = 1:length(fn);
	X = xptopen(fn{k1},'r');
	xptopen([fn{k1},'.xpt'],'w',X);
	f = fieldnames(X);

	fid = fopen([fn{k1},'.csv'],'w');
	for k1=1:length(f)
		if k1>1, fprintf(fid,';'); end; 	
		fprintf(fid,'%s',f{k1});
	end; 	
	fprintf(fid,'\n');

	for k2=1:length(X.(f{1}));
	for k1=1:length(f)
		if k1>1, fprintf(fid,';'); end; 	
		v = X.(f{k1})(k2);
		if isnumeric(v)
		        if strcmp(f{k1},'DATE'),
		        	fprintf(fid,'%s',datestr(v + datenum([1960,1,1]),1));
		        elseif strcmp(f{k1},'MONTH'),
		        	fprintf(fid,'%s',datestr(v + datenum([1960,1,1]),3));
		        elseif v==ceil(v),
				fprintf(fid,'%i',v);
			else
				fprintf(fid,'%f',v);
			end	
		elseif iscell(v) && ischar(v{1})	
			fprintf(fid,'%s',v{1});
		else
			fprintf(fid,'--');
		end;	
	end; 	
	fprintf(fid,'\n');
	end; 	

	fclose(fid);	
end; 

