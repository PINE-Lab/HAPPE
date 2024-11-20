function [C,N,c] = nanconv(X,Y,arg3)
% NANCONV computes the convolution for data with missing values. 
%  X and Y can contain missing values encoded with NaN.
%  NaN's are skipped, NaN do not result in a NaN output. 
%  The output gives NaN only if there are insufficient input data
%
% [...] = NANCONV(X,Y);
%      calculates 2-dim convolution between X and Y
% [C]   = NANCONV(X,Y);
%
% WARNING: missing values can introduce aliasing - causing unintended results.
%    Moreover, the behavior of bandpass and highpass filters in case of missing values 
%    is not fully understood, and might contain some pitfalls.  
%
% see also: CONV, NANCONV2, NANFFT, NANFILTER

%	$Id: conv2nan.m 6973 2010-02-28 20:19:12Z schloegl $
%	Copyright (C) 2000-2005,2010,2011 by Alois Schloegl <alois.schloegl@gmail.com>	
%       This function is part of the NaN-toolbox
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
%    along with this program; If not, see <http://www.gnu.org/licenses/>.


warning('NANCONV is experimental. For more details see HELP NANCONV');


if nargin~=2,
        fprintf(2,'Error NANCONV2: incorrect number of input arguments\n');
end;

m = isnan(X);
n = isnan(Y);

X(m) = 0;
Y(n) = 0;

C = conv(X,Y);         % 2-dim convolution
N = conv(real(~m),real(~n));     % normalization term
c = conv(ones(size(X)),ones(size(Y))); % correction of normalization

if nargout==1,
        C = C.*c./N;
elseif nargout==2,
        N = N./c;
end;

