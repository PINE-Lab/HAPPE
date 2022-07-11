function [Y,Z] = nanfilter(B,A,X,z);
% NANFILTER is able to filter data with missing values encoded as NaN. 
%       
%      [Y,Z] = nanfilter(B,A,X [, Z]);  
%
% If X contains no missing data, NANFILTER should behave like FILTER. 
% NaN-values are handled gracefully. 
%
% WARNING: missing values can introduce aliasing - causing unintended results.
%    Moreover, the behavior of bandpass and highpass filters in case of missing values 
%    is not fully understood, and might contain some pitfalls.  
%
% see also: FILTER, SUMSKIPNAN, NANFFT, NANCONV, NANFILTER1UC

%	$Id$
%	Copyright (C) 2005,2011 by Alois Schloegl <alois.schloegl@gmail.com>		
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


warning('NANFILTER is experimental. For more details see HELP NANFILTER');

na = length(A);
nb = length(B);
if any(size(X)==1)
	nc = 1; 
else	
	nc = size(X,2);
end; 

if nargin<4,
        [t,Z.S] = filter(B,A,zeros(na+nb,nc));
        [t,Z.N] = filter(B,A,zeros(na+nb,nc));
elseif isnumeric(z),
        Z.S = z; 
        [t, Z.N] = filter(B, A, zeros(na+nb,nc));
elseif isstruct(z),
        Z = z; 
end;

NX        = isnan(X);
X(NX)     = 0;

[Y , Z.S] = filter(B, A,   X, Z.S);
[NY, Z.N] = filter(B, A, ~NX, Z.N);
Y = (sum(B)/sum(A)) * Y./NY;

