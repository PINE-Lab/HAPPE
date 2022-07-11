function [pval, h, stats] = signrank (x, m, alpha, tail, DIM)
% SIGNRANK - Wilcoxon signed-rank test
%     The Wilcoxon signed-rank test is a non-parametric statistical hypothesis 
%     test used to compare two related samples whether their population median 
%     ranks differ [1-3]. SIGNRANK treads NaNs as "Missing values" and ignores these. 
%     Octave's statistical package has also wilcoxon_test, however, this works only 
%     for data with N>25 samples, signrank is based on the works [1-3] and can 
%     be used also for smaller sample sizes.
%
% pval = signrank(x,m)
%	tests Null-hypothesis that median of x is m. 		
% pval = signrank(x,y)
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
%     signrank works on the first non-singleton dimension or on DIM. 
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

%       Copyright (C) 2010,2019 by Alois Schloegl <alois.schloegl@ist.ac.at>	
%       This function is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/

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
  if (~ ischar (tail))
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
  	;
  elseif size(x,DIM) == size(m,DIM)
  	x = x-m;
  	m = zeros(szm);
  else
    error ('signrank: dimension of X and Y do not fit');
  end 	  
  
  % Algorithm according to 
  % https://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test

  % Step 2: exclude 0, and get N_r
  x(abs(x)==0)=NaN; 
  Nr=sum(~isnan(x),DIM);

  % Step 3 and 4: Order and Rank the data
  Rix = tiedrank(abs(x));

  % Glenn A Walker [3]
  % compute correction for ties; 
  sz  = size(x);
  P   = [DIM,1:DIM-1,DIM+1:length(sz)];
  HIS = histo2(reshape(permute(Rix,P),sz(DIM),sz([1:DIM-1,DIM+1:end])));	
  m   = HIS.H; 
  m(HIS.H<=1) = NaN;
  C   = ipermute(reshape(sumskipnan(m.*(m-1).*(m+1),1),[1,sz([1:DIM-1,DIM+1:end])]),P);

  %Step 5:
  W = sum(sign(x).*Rix,DIM);
  % z = W./sqrt(Nr.*(Nr+1).*(2*Nr+1)./(6*(Nr-1)));

  % https://math.stackexchange.com/questions/1414794/wilcoxon-signed-rank-test
  Tplus  = sumskipnan((x>0).*Rix, DIM);
  Tminus = sumskipnan((x<0).*Rix, DIM);

  stats.z = (max(Tplus,Tminus)-Nr.*(Nr+1)/4)./sqrt(Nr.*(Nr+1).*(2*Nr+1)./24);
  stats.signedrank = max(Tplus,Tminus);
  
  S   = (Tplus - Tminus) / 2;
  V   = (Nr.* (Nr+1).*(2*Nr+1) - C/2) / 24;
  t   = S .* sqrt(max(Nr-1,0)) ./ sqrt(Nr.*V - S.*S);
  cdf = tcdf(t, Nr);

  % see also NaN/ttest  
  if (strcmp (tail, '~=') || strcmp (tail, '!=') || strcmp (tail, '<>')) || strcmp(tail,'both'),
    pval = 2 * min (cdf, 1 - cdf);
  elseif strcmp (tail, '>') || strcmp(tail,'right'),
    pval = 1 - cdf;
  elseif strcmp (tail, '<') || strcmp(tail,'left'),
    pval = cdf;
  else
    error ('signrank: option %s not recognized', tail);
  end

  h = pval < alpha;	
  if (nargout == 0)
    fprintf(1,'  pval: %g\n', pval);
  end
  stats.t=t;	

%!test
%! % example from [3] 
%! x = [15,8;10,3;6,7;5,13;10,2;15,12;7,14;5,8;8,13;12,3;4,9;13,3;8,10;10,2;11,4;13,7;6,1;6,11;,9,3; 5,5;10,2;9,8;11,5;8,8]; 
%! [p,h,stats] = signrank( x(:,1), x(:,2) );
%! assert ( abs(stats.t - 2.184) < .01)

