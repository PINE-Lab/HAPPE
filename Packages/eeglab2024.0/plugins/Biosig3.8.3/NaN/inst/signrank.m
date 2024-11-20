function [pval, h, stats] = signrank (x, m, alpha, tail, DIM)
% SIGNRANK - Wilcoxon signed-rank test
%     The Wilcoxon signed-rank test is a non-parametric statistical hypothesis
%     test used to compare two related samples whether their population median
%     ranks differ [1-4]. SIGNRANK treads NaNs as "Missing values" and ignores these.
%     The current implementation has been validated against wilcox.test of R (v4.0.4).
%
% pval = signrank(X,m [, ...])
%	tests Null-hypothesis that median of x is m.
%       size(m,DIM) must be one, all other dimensions must match those of X
% pval = signrank(x,y [, ...])
% 	size of x and size of y must match, it is tested whether the
%	difference x-y is significantly different to m=0;
% pval = signrank(x,y,alpha)
% pval = signrank(x,y,alpha,tail)
% pval = signrank(x,y,alpha,tail,DIM)
% [pval,H,stats] = signrank(...)
%
%     H=1 indicates a rejection of the Null-hypothesis at a significance
%     level of alpha (default alpha = 0.05).
%
%     With the optional argument string TAIL, the alternative of interest
%     can be selected.  If TAIL is '!=' or '<>' or 'both', the null is tested
%     against the two-sided Alternative `mean (X) ~= mean (Y)'.  If TAIL
%     is '>' or 'right', the one-sided Alternative `mean (X) > mean (Y)' is used.
%     Similarly for '<' or 'left', the one-sided Alternative `mean (X) < mean
%     (Y)' is used.  The default is the two-sided case.
%
%     H returns whether the Null-Hypotheses must be rejected.
%     The p-value of the test is returned in PVAL.
%
%     signrank works on the first non-singleton dimension of X and Y, or on DIM.
%
%     If no output argument is given, the p-value of the test is
%     displayed.
%
% Reference(s):
% [1] Glenn A Walker, (2002)
%     Common Statistical Methods for Clinical Research (with SAS examples), 2nd edition
%     Chapter 12 The Wilcoxon Signed-Rank Test.
% [2] https://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test
% [3] https://math.stackexchange.com/questions/1414794/wilcoxon-signed-rank-test
% [4] https://de.wikipedia.org/wiki/Wilcoxon-Vorzeichen-Rang-Test

%    Copyright (C) 2010,2019,2021 by Alois Schloegl <alois.schloegl@gmail.com>
%    This function is part of the NaN-toolbox
%    http://pub.ist.ac.at/~schloegl/matlab/NaN/

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

  if ((nargin < 2) || (nargin > 5) || nargout > 4)
        print_usage ;
  end

  if (nargin == 2)
    alt  = '~=';
  end
  if (nargin < 3) || isempty(alpha)
    alpha = .05;
  end

  if (nargin < 4) || isempty(tail)
    tail = '~=';
  end
  if (~isnumeric(alpha))
    error ('signrank: alpha must be numeric');
  end
  if (~ischar(tail))
    error ('signrank: tail must be a string');
  end
  if nargin<5,
       DIM = find(size(x)>1,1);
  end;
  if isempty(DIM), DIM=1; end;

  szx = size(x);
  szm = size(m);
  szx(DIM) = 1;
  szm(DIM) = 1;
  if size(m,DIM)==1
  	x = x-repmat(m,size(x)./szx);
  elseif size(x,DIM) == size(m,DIM)
  	x = x-m;
  	m = zeros(szm);
  else
        error ('signrank: dimension of X and Y do not fit');
  end

  % Algorithm using "signed rank procedure" [2]
  x(x==0)=NaN;
  Nr=sum(~isnan(x) & (x~=0), DIM);

  % Order and Rank the data using "modified ranks" [2]
  Rix = tiedrank(abs(x), DIM);

  % Glenn A Walker [3,4]: compute correction for ties; 
  sz  = size(x);
  P   = [DIM,1:DIM-1,DIM+1:length(sz)];
  HIS = histo2(reshape(permute(Rix,P),sz(DIM),prod(sz([1:DIM-1,DIM+1:end]))));
  m   = HIS.H;
  m(HIS.H<=1) = NaN;
  C   = ipermute(reshape(sumskipnan(m.*(m-1).*(m+1),1),[1,sz([1:DIM-1,DIM+1:end])]),P);

  W   = sumskipnan(sign(x).*Rix,DIM);
  z   = W./sqrt(Nr.*(Nr+1).*(2*Nr+1)/6 - C/12);
  p2  = 2*abs(W)./(Nr.*Nr);

  % https://math.stackexchange.com/questions/1414794/wilcoxon-signed-rank-test
  Tplus  = sumskipnan((x>0) .* Rix, DIM);
  Tminus = sumskipnan((x<0) .* Rix, DIM);

  stats.signedrank  = min(Tplus, Tminus);
  #
  # simple formula is:
  #     z = (Tplus - Nr.*(Nr+1)/4)./sqrt(Nr.*(Nr+1).*(2*Nr+1)/24);
  # but for approximation for N<60, and ties (see [4])
  #     z = (abs(Tplus - Nr.*(Nr+1)/4)-0.5)./sqrt(Nr.*(Nr+1).*(2*Nr+1)/24 - C/48);
  V = sqrt(Nr.*(Nr+1).*(2*Nr+1)./24 - C/48);

  stats.z = z;
  % effective size [2] 
  stats.r = 2*(Tplus-Tminus)./(Nr.*(Nr+1));

  % see also NaN/ttest
  if (strcmp (tail, '~=') || strcmp (tail, '!=') || strcmp (tail, '<>')) || strcmp(tail,'both'),
    z    = -(abs(Tplus  - Nr.*(Nr+1)/4) - .5) ./ V;
    cdf  = normcdf(z);
    pval = 2 * min (cdf, 1 - cdf);
  elseif strcmp (tail, '>') || strcmp(tail,'right'),
    z    = -(Tplus  - Nr.*(Nr+1)/4 - .5) ./ V;
    pval = normcdf(z);
  elseif strcmp (tail, '<') || strcmp(tail,'left'),
    z    = -(Tminus - Nr.*(Nr+1)/4 - .5) ./ V;
    pval = normcdf(z);
  else
    error ('signrank: option %s not recognized', tail);
  end

  h = pval < alpha;
  if (nargout == 0)
    fprintf(1,'  pval: %g\n', pval);
  end

%!test
%! x = [15,8;10,3;6,7;5,13;10,2;15,12;7,14;5,8;8,13;12,3;4,9;13,3;8,10;10,2;11,4;13,7;6,1;6,11;,9,3; 5,5;10,2;9,8;11,5;8,8];
%! [p,h,stats] = signrank( x(:,1), x(:,2) );
%! assert ( abs(p - 0.045350) < 1e-6)
%!
%!test
%! [p,h,stats]=signrank([1:15]',8);
%! assert ( h==0 )
%!
%!test
%! [p,h,stats]=signrank([17.6, 20.6, 22.2, 15.3, 20.9, 21.0, 18.9, 18.9, 18.9, 18.2]',25);
%! assert ( abs(p - 0.005793) < 1e-5)
%!
%!test
%! [p,h,stats]=signrank([17.6, 20.6, 22.2, 15.3, 20.9, 21.0, 18.9, 18.9, 18.9, 18.2]',25,[],'<');
%! assert ( abs(p - 2.8965e-03) < 1e-6)
%!
%!test
%! [p,h,stats]=signrank([17.6, 20.6, 22.2, 15.3, 20.9, 21.0, 18.9, 18.9, 18.9, 18.2]', 25, [], '>');
%! assert ( abs(p - 0.9979) < 1e-5)
%!
%!test
%! [p,h,stats]=signrank([17.6, 20.6, 22.2, 15.3, 20.9, 21.0, 18.9, 18.9, 18.9, 18.2]', 19);
%! assert ( abs(p - 0.838) < 1e-4)
%!
%!test
%! [p,h,stats]=signrank([17.6, 20.6, 22.2, 15.3, 20.9, 21.0, 18.9, 18.9, 18.9, 18.2]', 19, [], '<');
%! assert ( abs(p - 0.6204) < 1e-4)
%!
%!test
%! [p,h,stats]=signrank([17.6, 20.6, 22.2, 15.3, 20.9, 21.0, 18.9, 18.9, 18.9, 18.2]', 19, [], '>');
%! assert ( abs(p - 0.419) < 1e-4)
%!
%!test
%! [p,h,stats]=signrank([0, 2, 3, 4, 6, 7, 8, 9, 11, 14, 15, 17, -18]+1,1);
%! assert ( abs(p - 0.037633) < 1e-6)
%!
%!test
%! [p,h,stats]=signrank([0, 2, 3, 4, 6, 7, 8, 9, 11, 14, 15, 17, -18],0);
%! assert ( abs(p - 0.037633) < 1e-6)
%!
%!test
%! [p,h,stats]=signrank([   2, 3, 4, 6, 7, 8, 9, 11, 14, 15, 17, -18],0);
%! assert ( abs(p - 0.037633) < 1e-6)
%!
%!test
%! [p,h,stats]=signrank([-1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,-13],0);
%! assert ( abs(p - 0.030276) < 1e-6)
%!
%!test
%! [p,h,stats]=signrank([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, -12],0);
%! assert ( abs(p - 0.037633) < 1e-6)
%!
%!test
%! x = [125 	110;115	122;130	125;140	120;140	140;115	124;140	123;125	137;140	135;135	145];
%! [p,h,stats] = signrank( x(:,1), x(:,2) );
%%! assert ( abs(p - 0.6113) < 1e-4)	% this does not match [2]
%! assert ( abs(p - 0.6353) < 1e-4)	% matches the result in R
%!
%!test
%! x = [4	1;3	2; 2	3; 5	0; 5	0; 3	2];
%! [p,h,stats] = signrank( x(:,1), x(:,2) );
%!
%!test
%! SZ = [100,5,3]; X=round(rand(SZ)*10); Y=round(rand(SZ)*10);
%! p1 = signrank(X,Y,1,[],1);
%! assert(size(p1,2)==SZ(2) && size(p1,3)==SZ(3))
%! p2 = signrank(X,Y,1,[],2);
%! assert(size(p2,1)==SZ(1) && size(p2,3)==SZ(3))
%! p3 = signrank(X,Y,1,[],3);
%! assert(size(p3,2)==SZ(2) && size(p3,1)==SZ(1))

