function [H,p,ks2stat,D] = kstest2(x,y,varargin);
% KSTEST2 computes the two-sampleKolmogorov-Smirnov.
%
% Usage: 
%    [H,p,ks2stat,D] = kstest2(x,y);
%    [...] = kstest2(x, y, [, 'alpha', alpha] [, 'tail', '>'] );
%
% Input: 
%    x, y input vectors for comparison
%    X  matrix whos colums are pairwise compared, such
%
% Output:
%    H   1: statistical significance (p < alpha)
%    D   maximum absolute difference between sample data
%        D(k,l) is the m.a.d. from X(:,k) and X(:,l)
%    df  is the degree-of freedom 
%        df(k,l) = n(k)*n(l)/(n(k)+n(l)) with n samples of corresponding 
%        column X. 
%    p   p-value, it's also a matrix where
%        pval(k,l) is the p-value from column k and l
% 
% see also:
%    kolmogorov_smirnov

%    Copyright (C) 2019,2020 by Alois Schloegl <alois.schloegl@gmail.com>
%    	This is part of the NaN-toolbox 
%	https://octave.sourceforge.io/nan/index.html
%	https://pub.ist.ac.at/~schloegl/matlab/NaN/
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.


% default
alpha=0.05;
tail='unequal';

k=2;
while k < length(varargin)
      if strcmpi(varargin{k}, 'tail')  
           tail = varargin{k+1};
           k    = k+1;
      elseif strcmpi(varargin{k}, 'alpha')
           alpha = varargin{k+1};
           k    = k+1;
      else 
           error(sprintf('argument %d not supported - ignored', k))
      end
      k=k+1;
end

[D, ks2stat, p, df] = kolmogorov_smirnov(x(:),y(:), 'tail', tail);
H = p < alpha; 

%! assert(kstest2([1:5]',[1:5]'+5))




