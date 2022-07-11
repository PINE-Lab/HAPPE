function [C,N,LAGS] = xcovf(X,Y,MAXLAG,SCALEOPT)
% XCOVF generates cross-covariance function. 
% XCOVF is the same as XCORR except 
%   X and Y can contain missing values encoded with NaN.
%   NaN's are skipped, NaN do not result in a NaN output. 
%   The output gives NaN only if there are insufficient input data
%
% [C,N,LAGS] = xcovf(X,MAXLAG,SCALEOPT);
%      calculates the (auto-)correlation function of X
% [C,N,LAGS] = xcovf(X,Y,MAXLAG,SCALEOPT);
%      calculates the crosscorrelation function between X and Y
%
%  SCALEOPT   [character string] specifies the type of scaling applied
%          to the correlation vector (or matrix). is one of:
%    'none'      return the unscaled correlation, R,
%    'biased'    return the biased average, R/N, 
%    'unbiased'  return the unbiassed average, R(k)/(N-|k|), 
%    'coeff'     return the correlation coefficient, R/(rms(x).rms(y)),
%          where "k" is the lag, and "N" is the length of X.
%          If omitted, the default value is "none".
%          If Y is supplied but does not have the ame length as X,
%          scale must be "none".
%
%
% see also: COVM, XCORR

%	Copyright (C) 2005,2010,2011 by Alois Schloegl <alois.schloegl@gmail.com>	
%       This function is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/

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

if nargin<2,
        Y = [];
        MAXLAG = [];
        SCALEOPT = 'none';
elseif ischar(Y),
        MAXLAG = Y;
        SCALEOPT=MAXLAG;
        Y=[];
elseif all(size(Y)==1),
        if nargin<3
                SCALEOPT = 'none';
        else
                SCALEOPT = MAXLAG;
        end;
        MAXLAG = Y; 
        Y = [];
end;

if 0,

elseif isempty(Y) && isempty(MAXLAG)
        NX = isnan(X);
        X(NX) = 0;
        [C,LAGS] = xcorr(X,'none');
        [N,LAGS] = xcorr(1-NX,'none');
elseif ~isempty(Y) && isempty(MAXLAG)
        NX = isnan(X);
        NY = isnan(Y);
        X(NX) = 0;
        Y(NY) = 0;
        [C,LAGS] = xcorr(X,Y,'none');
        [N,LAGS] = xcorr(1-NX,1-NY,'none');
elseif isempty(Y) && ~isempty(MAXLAG)
        NX = isnan(X);
        X(NX) = 0;
        [C,LAGS] = xcorr(X,MAXLAG,'none');
        [N,LAGS] = xcorr(1-NX,MAXLAG,'none');
elseif ~isempty(Y) && ~isempty(MAXLAG)
        NX = isnan(X);
        NY = isnan(Y);
        X(NX) = 0;
        Y(NY) = 0;
        [C,LAGS] = xcorr(X,Y,MAXLAG,'none');
        [N,LAGS] = xcorr(1-NX,1-NY,MAXLAG,'none');
end;

if 0,

elseif strcmp(SCALEOPT,'none')
	% done

elseif strcmp(SCALEOPT,'coeff')
	ix = find(LAGS==0);
	if ~any(size(X)==1), %% ~isvector(X)	
		c  = C(ix,1:size(X,2)+1:end);	%% diagonal elements
		v  = c.^-0.5; % sqrt(1./c(:));
		v  = v'*v;
		C  = C.*repmat(v(:).',size(C,1),1);
	elseif isempty(Y)
		C = C/C(ix);
	else 
		C = C/sqrt(sumsq(X)*sumsq(Y));
	end;		

elseif strcmp(SCALEOPT,'biased')
	C = C./repmat(max(N),size(C,1),1);

elseif strcmp(SCALEOPT,'unbiased')
	C = C./(repmat(max(N),size(C,1),1)-repmat(LAGS,1,size(C,2)));

else
        warning('invalid SCALEOPT - not supported');
end;

