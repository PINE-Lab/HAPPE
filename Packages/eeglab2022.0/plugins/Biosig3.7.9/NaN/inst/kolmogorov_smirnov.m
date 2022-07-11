function [D,ks,p,df] = kolmogorov_smirnov(x,varargin);
% KOLMOGOROV_SMIRNOV computes the two-sample Kolmogorov-Smirnov test for 
%   each pair columns. If data size does not match, the data can be 
%   filled up with not-a-number (NaN).
%
% Usage: 
%    [D,ks,pval,df] = kolmogorov_smirnov(x,y);
%    [...] = kolmogorov_smirnov(X);
%    [...] = kolmogorov_smirnov(..,'Tail',tail);
%
% Input: 
%    x,y input vectors for comparison
%    X   data matrix, each column represents a sample distribution
%        in case, the number of samples do not match, the matrix can
%        can be filled up with not-a-number (NaN) values. 
%    tail: 'unequal' (default), 'larger', 'smaller'
%
% Output:
%    D   maximum absolute difference between sample data
%        D(k,l) is the m.a.d. from X(:,k) and X(:,l) 
%
%    df  is the degree-of freedom 
%        df(k,l) = n(k)*n(l)/(n(k)+n(l)) with n samples of corresponding 
%        column X. 
%    pval  p-value, it's also a matrix where 
%        pval(k,l) is the p-value from column k and l. 
%
%  NOTE: For M data sets (M-columns of X), there are M*(M-1)/2 tests, 
%    you might want to consider some correction for multiple 
%    comparision. 
% 
% see also:
%    kstest2, statistic/kolmogorov_smirnov_test_2


%	Copyright (C) 2019 by Alois Schloegl <alois.schloegl@ist.ac.at>	
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
tail='unequal';

if nargin<2
        X=x;
        IX=1:size(X,2);
        IY=IX;
else
    k=1; 
    if isnumeric(varargin{k})
        y = varargin{k};
        n1=size(x,1);
        n2=size(y,1);
        if n1==n2;
                X=[x,y];
        elseif n1<n2
                X=[[x;repmat(NaN,n2-n1,size(x,2))],y];
        else
                X=[x,[y;repmat(NaN,n1-n2,size(y,2))]];
        end
        k=2;
        IX=1:size(x,2);
        IY=[1:size(y,2)]+size(x,2);
    else
        X=x;
        IX=1:size(X,2);
        IY=IX;
    end
    if strcmpi(varargin{k},'tail')
        tail = varargin{k+1};
    end
  
end
sz=size(X);

HIS=histo3(X);

% compute compulative distribution function
cdf=cumsum(HIS.H,1)./repmat(HIS.N,size(HIS.H,1),1);

% maxium absolute differences between all pairs of columns
[ix,iy] = meshgrid(IX,IY);
if strcmp(tail, 'unequal')
        D  = reshape(max(abs(cdf(:,ix)-cdf(:,iy))),length(IX),length(IY));
elseif strcmp(tail, 'larger')
        D  = reshape(max(cdf(:,ix)-cdf(:,iy)),length(IX),length(IY));
elseif strcmp(tail, 'smaller')
        D  = reshape(-min(cdf(:,ix)-cdf(:,iy)),length(IX),length(IY));
end

% degree of freedom 
df = HIS.N(ix).*HIS.N(iy)./(HIS.N(ix)+HIS.N(iy));
ks = sqrt(df).*D;
p  = exp(-2.*ks.*ks);



