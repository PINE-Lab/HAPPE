function [x,z] = nanfilter1uc(uc,x,z);
% NANFILTER1UC is an adaptive filter for data with missing values encoded as NaN. 
%       
%      [Y,Z] = nanfilter1uc(uc,X [, Z]);  
%
% if X contains no missing data, NANFILTER behaves like FILTER(uc,[1,uc-1],X[,Z]).
%
% see also: FILTER, NANFILTER, SUMSKIPNAN

%	$Id$
%	Copyright (C) 2010,2011 by Alois Schloegl <alois.schloegl@gmail.com>	
%       This function is part of the NaN-toolbox available at 
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/ and 
%	http://octave.svn.sourceforge.net/viewvc/octave/trunk/octave-forge/extra/NaN/inst/

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
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA


na = 2; %length(A);
nb = 2; %length(B);
if any(size(x)==1)
	nc = 1; 
else	
	nc = size(x,2);
end; 

acN = zeros(1,nc);
if nargin<3,
	z = zeros(1,nc);
end;
acc = NaN(1,nc);
for k = 1:size(x,1),
	ix = isnan(x(k,:));
	acN = acN.*ix+1;
	UC1 = ((1-uc).^acN);
        acc(~ix) = (1-UC1(~ix)) .* x(k,~ix) + z(~ix);  % / A{1};
	ix = isnan(acc);
	acc(ix) = x(k,ix); 
	z   = (1-uc) * acc;
        x(k,:) = acc;
end;

