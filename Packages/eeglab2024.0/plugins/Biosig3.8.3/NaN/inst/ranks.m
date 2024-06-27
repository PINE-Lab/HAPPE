function r = ranks(X,DIM,Mode)
% RANKS gives the rank of each element in a vector.
% This program uses an advanced algorithm with averge effort O(m.n.log(n)) 
% NaN in the input yields NaN in the output.
% 
% r = ranks(X[,DIM])
%   if X is a vector, return the vector of ranks of X adjusted for ties.
%   if X is matrix, the rank is calculated along dimension DIM. 
%   if DIM is zero or empty, the lowest dimension with more then 1 element is used. 
% r = ranks(X,DIM,'traditional')
%   implements the traditional algorithm with O(n^2) computational 
%   and O(n^2) memory effort
% r = ranks(X,DIM,'mtraditional')
%   implements the traditional algorithm with O(n^2) computational 
%   and O(n) memory effort
% r = ranks(X,DIM,'advanced   ')
%   implements an advanced algorithm with O(n*log(n)) computational 
%   and O(n.log(n)) memory effort
% r = ranks(X,DIM,'advanced-ties')
%   implements an advanced algorithm with O(n*log(n)) computational 
%   and O(n.log(n)) memory effort
%   but without correction for ties 
%   This is the fastest algorithm 
%
% see also: CORRCOEF, SPEARMAN, RANKCORR
%
% REFERENCES:
% --


%    $Id$
%    Copyright (C) 2000-2002,2005,2010,2013 by Alois Schloegl <alois.schloegl@gmail.com>	
%    This script is part of the NaN-toolbox
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

% Features:
% + is fast, uses an efficient algorithm for the rank correlation
% + computational effort is O(n.log(n)) instead of O(n^2)
% + memory effort is O(n.log(n)), instead of O(n^2). 
%     Now, the ranks of 8000 elements can be easily calculated
% + NaNs in the input yield NaN in the output 
% + compatible with Octave and Matlab
% + traditional method is also implemented for comparison. 


if nargin<2, DIM = 0; end;
if ischar(DIM),
	Mode= DIM; 
	DIM = 0; 
elseif (nargin<3), 
	Mode = '';
end; 
if isempty(Mode),
	Mode='advanced   '; 
end;

sz_orig = size (X);
X = squeeze (X); %remove singleton dimensions for convenience
nd = ndims (X);
if (~DIM)
	 DIM = 1;
end
if DIM > 1 %shift the array so that the dimension to sort over is first
  perm = [DIM 1:(DIM-1) (DIM+1):nd];
  X = permute (X, perm);
end
if nd > 2  %convert X to 2-D if it has >2 dimensions
  sz = size(X);
  N = sz(1);
  M = prod(sz(2:end));
  X = reshape(X, N, M);
else
  [N,M] = size(X);
end

if strcmp(Mode(1:min(11,length(Mode))),'traditional'), % traditional, needs O(m.n^2)
% this method was originally implemented by: KH <Kurt.Hornik@ci.tuwien.ac.at>
% Comment of KH: This code is rather ugly, but is there an easy way to get the ranks adjusted for ties from sort?

r = zeros(size(X));
        for i = 1:M;
                p = X(:, i(ones(1,N)));
                r(:,i) = (sum (p < p') + (sum (p == p') + 1) / 2)';
        end;
        % r(r<1)=NaN;
        
elseif strcmp(Mode(1:min(12,length(Mode))),'mtraditional'), 
        % + memory effort is lower
        
	r = zeros(size(X));
        for k = 1:N;
        for i = 1:M;
                r(k,i) = (sum (X(:,i) < X(k,i)) + (sum (X(:,i)  == X(k,i)) + 1) / 2);
        end;
        end;
        % r(r<1)=NaN;
        
elseif strcmp(Mode(1:min(13,length(Mode))),'advanced-ties'), % advanced
        % + uses sorting, hence needs only O(m.n.log(n)) computations         
        % - does not fix ties

        r = zeros(size(X));
        [sX, ix] = sort(X,1); 
	for k=1:M,
	        [tmp,r(:,k)] = sort(ix(:,k),1);	    % r yields the rank of each element 	
	end;        
	r(isnan(X)) = nan;

        
elseif strcmp(Mode(1:min(8,length(Mode))),'advanced'), % advanced
        % + uses sorting, hence needs only O(m.n.log(n)) computations         
        
        % [tmp,ix] = sort([X,Y]);     
        % [tmp,r] = sort(ix);     % r yields rank. 
        % but because sort does not work accordingly for cell arrays, 
        % and DIM argument not supported by Octave 
        % and DIM argument does not work for cell-arrays in Matlab
        % we sort each column separately:
        
        r = zeros(size(X));
        n = N;
        for k = 1:M,
                [sX,ix] = sort(X(:,k)); 
                [tmp,r(:,k)] = sort(ix);	    % r yields the rank of each element 	
                
                % identify multiple occurences (not sure if this important, but implemented to be compatible with traditional version)
                if isnumeric(X)
                        n=sum(~isnan(X(:,k)));
                end;
                x = [0;find(sX~=[sX(2:N);n])];    % for this reason, cells are not implemented yet.   
                d = find(diff(x)>1);
                
                % correct rank of multiple occurring elements
                for l = 1:length(d),
                        t = (x(d(l))+1:x(d(l)+1))';
                        r(ix(t),k) = mean(t);
                end;
        end;
        r(isnan(X)) = nan;
        
elseif strcmp(Mode,'=='), 
% the results of both algorithms are compared for testing.    
%
% if the Mode-argument is omitted, both methods are applied and 
% the results are compared. Once the advanced algorithm is confirmed, 
% it will become the default Mode. 

        r  = ranks(X,'advanced   ');
        r(isnan(r)) = 1/2;
        
        if N>100,
	        r1 = ranks(X,'mtraditional');  % Memory effort is lower 
        else
                r1 = ranks(X,'traditional');
        end;
        if ~all(all(r==r1)),
                fprintf(2,'WARNING RANKS: advanced algorithm does not agree with traditional one\n Please report to <alois.schloegl@gmail.com>\n');
                r = r1;
        end;
        r(isnan(X)) = nan;
end;

%reshape r to match the input X
if nd > 2
  r = reshape (r, sz);
end
if (DIM > 1)
	r = ipermute (r, perm);
end
r = reshape (r, sz_orig); %restore any singleton dimensions


%!shared z, r1, r2
%! z = magic (4);
%! r1 = [4   1   1   4; 2   3   3   2; 3   2   2   3; 1   4   4   1];
%! r2 = [4   1   2   3; 1   4   3   2; 3   2   1   4; 2   3   4   1];
%!assert (ranks(z), r1);
%!assert (ranks(z, 2), r2);
%! z = nan(2, 2, 2);  
%! z(:, :, 1) = [1 2; 3 4];
%! z(:, :, 2) = [4 3; 2 1];
%! r1 = cat(3, [1 1; 2 2], [2 2; 1 1]);
%! r2 = cat(3, [1 2; 1 2], [2 1; 2 1]);
%!assert (ranks(z), r1);
%!assert (ranks(z, 2), r2);
%!assert (ranks(z, 3), r1);
