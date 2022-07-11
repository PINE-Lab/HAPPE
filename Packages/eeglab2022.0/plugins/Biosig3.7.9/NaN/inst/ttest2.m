function [h, pval, ci, stats, df] = ttest2 (x, y, alpha, tail, vartype, DIM)
% TTEST2 (unpaired) t-test
%     For two samples x and y from normal distributions with unknown
%     means and unknown equal variances, perform a two-sample t-test of
%     the null hypothesis of equal means.  Under the null, the test
%     statistic T follows a Student distribution with DF degrees of
%     freedom.
%
%     TTEST2 treads NaNs as "Missing values" and ignores these. 
%
% H = ttest2(x,y)
% H = ttest2(x,y,alpha)
% H = ttest2(x,y,alpha,tail)
% H = ttest2(x,y,alpha,tail,vartype)
% H = ttest2(x,y,alpha,tail,vartype,DIM)
% [H,PVAL] = ttest2(...)
% [h,p,ci,stats] = ttest2(...)
%
%     H=1 indicates a rejection of the Null-hypothesis at a significance 
%     level of alpha (default alpha = 0.05).	 
% 
%     With the optional argument string TAIL, the Alternative of interest
%     can be selected.  If TAIL is '!=' or '<>' or 'both', the null is tested
%     against the two-sided Alternative `mean (X) ~= mean (Y)'.  If TAIL
%     is '>' or 'right', the one-sided Alternative `mean (X) > mean (Y)' is used.
%     Similarly for '<' or 'left', the one-sided Alternative `mean (X) < mean
%     (Y)' is used.  The default is the two-sided case.
% 
%     vartype support only 'equal' (default value); the value 'unequal' is not supported. 
%
%     H returns whether the Null-Hypotheses must be rejected. 
%     The p-value of the test is returned in PVAL. 
% 
%     TTEST2 works on the first non-singleton dimension or on DIM. 
%
%     If no output argument is given, the p-value of the test is
%     displayed.
%

%%% not supported yet 
% [h,p,ci] = ttest2(...)
% [h,p,ci,stats] = ttest2(...)

%       $Id$
%       Copyright (C) 1995, 1996, 1997, 1998, 2000, 2002, 2005, 2006, 2007
%               Kurt Hornik
%       Copyright (C) 2010 by Alois Schloegl <alois.schloegl@gmail.com>	
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

  if ((nargin < 2) || (nargin > 6) || nargout > 4)
        print_usage ;
  end

  if (nargin < 3) || isempty(alpha)
    alpha = .05;
  end

  if (nargin < 4) || isempty(tail)
    tail = '~=';
  end
  if (~ ischar (tail))
    error ('ttest2: tail must be a string');
  end
  if (nargin < 5) || isempty(vartype)
    vartype = 'equal';
  end
  if ~strcmp(vartype,'equal')
    error ('test: vartype not supported')
  end	
  if nargin<6,
       DIM = find(size(x)>1,1);
  end;
  if isempty(DIM), DIM=1; end;

  szx = size(x); 
  szy = size(y);	
  szy(DIM) = 1;	  
  szx(DIM) = 1;
  
  if (any(szx-szy))
    error ('ttest2: dimension of X and Y do not fit');
  end

  [SX, NX] = sumskipnan(x, DIM);
  [SY, NY] = sumskipnan(y, DIM);
  stats.df = NX + NY - 2;
  MX = SX ./ NX;
  MY = SY ./ NY;

  if any(size(x)==0) || any(size(y)==0)
  	v = NaN;
  else	
  	v = sumsq(x-repmat(MX,size(x)./size(MX))) + sumsq(y-repmat(MY,size(y)./size(MY)));
  end; 	
  stats.sd    = sqrt(v/stats.df);
  stats.tstat = (MX - MY) .* sqrt ((NX .* NY .* stats.df) ./ (v .* (NX + NY)));
  cdf         = tcdf (stats.tstat, stats.df);

  if (strcmp (tail, '~=') || strcmp (tail, '!=') || strcmp (tail, '<>')) || strcmp(tail,'both'),
    pval = 2 * min (cdf, 1 - cdf);
  elseif strcmp (tail, '>') || strcmp(tail,'right'),
    pval = 1 - cdf;
  elseif strcmp (tail, '<') || strcmp(tail,'left'),
    pval = cdf;
  else
    error ('ttest2: option %s not recognized', tail);
  end

  h = pval < alpha;	

  if (nargout == 0)
    fprintf(1,'  pval: %g\n', pval);
  end

