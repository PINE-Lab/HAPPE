function LMOM = lmom(data,P,opt)
% LMOM estimates the L-Moments [1,2] from a sample distribution  
%   and might be a useful density estimation [1,3].
%   LMOM is equivalent to samlmo.F from dataplot [4].
%
% Usage: 
%   XMOM = lmom(X,P)
%   XMOM = lmom(X,P,'ratios')
%
%   X	input data, NaN's are ignored
%   P   maximum order, L moments 1:P are estimated
%   option: default 'false', 
%           'ratios': compute L-moment ratios
%   XMOM  vector of L-Moments from 1:P
% 	in case option='ratios', XMOM(3:P) will 
% 	return the L-moment rations (i.e. scaled L-moments).
%
% The current implementation is tested only on data sets up to 1000 samples
% and P=10. The algorithm has not been analyzed with respect to accuracy and
% computational efficiency. Eventually, this implementation should be
% compared also to samlmu.F from dataplot [4], which is also used in [5]. 
%
% References: 
% [1] Hosking (1990), L-MOMENTS: ANALYSIS AND ESTIMATION OF DISTRIBUTIONS, 
%        J. R. Statist. Soc. B (1990), 52,No. 1,pp. 105-124
% [2] https://en.wikipedia.org/wiki/L-moment
% [3] https://en.wikipedia.org/wiki/Density_estimation
% [4] Hosking, function samlmo.F from https://github.com/usnistgov/dataplot
% [5] 'lmom'-package for R, available from 
%     https://www.rdocumentation.org/packages/lmom/versions/2.8/topics/lmom-package

% Copyright (C) 2019,2020 by Alois Schlögl <alois.schloegl@gmail.com>
%       This function is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/

%% This program is free software: you can redistribute it and/or
%% modify it under the terms of the GNU General Public License as
%% published by the Free Software Foundation, either version 3 of the
%% License, or (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful, but
%% WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%% General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; see the file COPYING.  If not, see
%% <http://www.gnu.org/licenses/>.


if nargin<3
	opt=0;
end
opt = strcmp(opt,'ratios');

data(isnan(data))=[];
F = [0:length(data)]'/length(data);
u = sort(data);
u = u([1,1:end]);

% TODO: one might do this more efficiently 
p=repmat(NaN,P,P);
for r = 1:P
for k = 1:r
	p(r,k)=bincoeff(r,k)*bincoeff(r+k,k)*(-1)^(r-k);
end
end

for k = 1:P,
	xi(k) = trapz(F, u .* F.^(k-1));
	if k==1,
		LMOM(k) = xi(k);
	else 
		LMOM(k) = xi(1:k) * [(-1).^(k-1); p(k-1, 1:k-1)'];
	end
end
if (opt && (P>2))
	LMOM(3:P) = LMOM(3:P)/LMOM(2);
end

return

% lambda(1) = trapz(F, [u] ); 
% lambda(2) = trapz(F, [u] .* (2*F-1)) ;
% lambda(3) = trapz(F, [u] .* (6*F.^2 - 6*F + 1)) ;
% lambda(4) = trapz(F, [u] .* (20*F.^3 - 30*F.^2 + 12*F - 1)) ;
%


