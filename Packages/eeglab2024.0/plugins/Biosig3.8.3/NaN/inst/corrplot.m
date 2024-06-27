function RES = corrplot(data, varargin)
% CORRPLOT displays the correlation plot 
%
%   corrplot(data)
%   corrplot(data,'type',TYPE)
%   [R,PValue,H] = corrplot(data,Name,Value)
% 
% Input: 
%   data    
%   TYPE:  'Pearson' (default), 'Kendall', 'Spearman'
% 
% 
%	Copyright (C) 2021 by Alois Schloegl <alois.schloegl@gmail.com>	
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

Mode=[];
alpha=0.05;

k=1; 
while k<=nargin,
	if strcmpi(varargin{k},'type')
		Mode=varargin{k+1};
		k=k+1;
	elseif strcmpi(varargin{k},'alpha')
		alpha=varargin{k+1};
		
	end

	k=k+1;
end

[nr,nc]=size(data); 

for k1=1:nc
for k2=1:nc
	subplot(nc,nc,k1*nc+k2-nc)
	plot(data(:,k1),data(:,k2),'d')
end
end

R = corrcoef(data);


