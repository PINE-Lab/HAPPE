function p = fishers_exact_test(a,b,c,d)
% FISHERS_EXACT_TEST implements Fisher's exact test for the analysis of 
% contincency tables e.g. "Lady tasting tea" experiment [1-6].
%
% Usage: 
%      p = fishers_exact_test(H) 
%      p = fishers_exact_test(a,b,c,d) 
% 
% with H being a 2x2 matrix representing a contincency table H = [[a,b];[c,d]]
% and p is the resulting p-value. The implementation provides exact results,
% when (1) the symbolic toolbox (with vpa) is loaded, or (2) for small sample
% sizes. In the latter case, the result might be subject to the limited accuracy of
% floating point numbers for large sample sizes (a warning might be shown);
% in the case, the symbolic toolbox should be loaded.
%
% References:
% [1] https://en.wikipedia.org/wiki/Fisher%27s_exact_test
% [2] https://en.wikipedia.org/wiki/Lady_tasting_tea
% [3] Fisher, R. A. (1922). "On the interpretation of χ2 from contingency 
%     tables, and the calculation of P". 
%     Journal of the Royal Statistical Society. 85 (1): 87–94. 
%     doi:10.2307/2340521. JSTOR 2340521.
% [4] Fisher, R.A. (1954). Statistical Methods for Research Workers. 
%     Oliver and Boyd. ISBN 0-05-002170-2.
% [5] Agresti, Alan (1992). "A Survey of Exact Inference for Contingency Tables". 
%     Statistical Science. 7 (1): 131–153. 
%     CiteSeerX 10.1.1.296.874. doi:10.1214/ss/1177011454. JSTOR 2246001.
% [6] Fisher, Sir Ronald A. (1956) [The Design of Experiments (1935)]. 
%     "Mathematics of a Lady Tasting Tea". In James Roy Newman (ed.). 
%     The World of Mathematics, volume 3. Courier Dover Publications. 
%     ISBN 978-0-486-41151-4.

% Copyright (C) 2019 Alois Schloegl <alois.schloegl@gmail.com>
% This is part of the BIOSIG-toolbox http://biosig.sf.net/
%
% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% version 3 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Library General Public License for more details.
%
% You should have received a copy of the GNU Library General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc., 59 Temple Place - Suite 330,
% Boston, MA  02111-1307, USA.

if (nargin==1) && isequal(size(a),[2,2]),
	H = a;
	a=H(1,1);
	b=H(1,2);
	c=H(2,1);
	d=H(2,2);
elseif (nargin==4) && isscalar(a) && isscalar(b) && isscalar(c) && isscalar(d)
	H=[[a,b];[c,d]];
else
	error('invalid input argument')
end

try
	% use symbolic package if available
	a = vpa(a);
	b = vpa(b);
	c = vpa(c);
	d = vpa(d);
end

u = nchoosek(a+b,a);
v = nchoosek(c+d,c);
w = nchoosek(a+b+c+d,a+c);

if strcmp(lastwarn(),'nchoosek: possible loss of precision')
	printf('It is recommended to load the symbolic package, and re-run fishers_exact_test.\n')
end

p = u * v / w;

%!assert((double(fishers_exact_test(1,1,1,1))-2/3)<eps)
%!assert(abs(double(fishers_exact_test(10,1,1,10))-1.715261003186700e-04) < eps)


