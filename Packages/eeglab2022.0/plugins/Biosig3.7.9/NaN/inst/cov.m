function CC = cov(X,Y,Mode)
% COV covariance matrix
% X and Y can contain missing values encoded with NaN.
% NaN's are skipped, NaN do not result in a NaN output. 
% The output gives NaN only if there are insufficient input data
% The mean is removed from the data. 
% 
% Remark: for data contains missing values, the resulting 
% matrix might not be positiv definite, and its elements have magnitudes
% larger than one. This ill-behavior is more likely for small sample 
% sizes, but there is no garantee that the result "behaves well" for larger
% sample sizes. If you want the a "well behaved" result (i.e. positive 
% definiteness and magnitude of elements not larger than 1), use CORRCOEF. 
% However, COV is faster than CORRCOEF and might be good enough in some cases.
%
% C = COV(X [,Mode]);
%      calculates the (auto-)correlation matrix of X
% C = COV(X,Y [,Mode]);
%      calculates the crosscorrelation between X and Y. 
%      C(i,j) is the correlation between the i-th and jth 
%      column of X and Y, respectively. 
%   NOTE: Octave and Matlab have (in some special cases) incompatible implemenations. 
%       This implementation follows Octave. If the result could be ambigous or  
%       incompatible, a warning will be presented in Matlab. To avoid this warning use: 
%       a) use COV([X(:),Y(:)]) if you want the traditional Matlab result. 
%       b) use C = COV([X,Y]), C = C(1:size(X,2),size(X,2)+1:size(C,2)); if you want to be compatible with this software.  
%
% Mode = 0 [default] scales C by (N-1)
% Mode = 1 scales C by N. 
%
% see also: COVM, COR, CORRCOEF, SUMSKIPNAN
%
% REFERENCES:
% http://mathworld.wolfram.com/Covariance.html

%	$Id$
%	Copyright (C) 2000-2003,2005,2009,2011,2012 by Alois Schloegl <alois.schloegl@ist.ac.at>	
%       This function is part of the NaN-toolbox 
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/

%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; If not, see <http://www.gnu.org/licenses/>.


if nargin==1
        Mode = 0;
        Y = [];
elseif nargin==2,
	% if all(size(Y)==1) & any(Y==[0,1]); 	% This is not compatible with octave 
	% short-circuit evaluation is required  
	% but for compatibility to matlab, && is avoided  
	SW = all(size(Y)==1);
	if SW, SW = any(Y==[0,1]); end;
	if SW,
		Mode = Y;
                Y = [];
        else
                Mode = 0;        
	end;
elseif nargin==3, 
		        
else
	fprintf(2,'Error COV: invalid number of arguments\n');
end;

if ~exist('OCTAVE_VERSION','builtin') && ~isempty(Y) && (size(X,2)+size(Y,2)~=2), 	
        % COV in Matlab is differently defined than COV in Octave. 
        % For compatibility reasons, this branch reflects the difference. 
        fprintf(2,'Warning NaN/COV: This kind of use of COV is discouraged because it produces different results for Matlab and Octave. \n');
        fprintf(2,'  (a) the traditional Matlab result can be obtained with:  C = COV([X(:),Y(:)]).\n');
        fprintf(2,'  (b) the traditional Octave result can be obtained with:  C = COV([X,Y]); C = C(1:size(X,2),size(X,2)+1:size(C,2)).\n');

        if numel(Y)~=numel(X),
                error('The lengths of X and Y must match.');
        end;
        X = [X(:),Y(:)];
        Y = [];
end;

if isempty(Y)
	CC = covm(X,['D',int2str(Mode>0)]);	
else        
        CC = covm(X,Y,['D',int2str(Mode>0)]);	
end;

