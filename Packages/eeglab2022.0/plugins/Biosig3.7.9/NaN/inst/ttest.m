function [h, pval, ci, stats] = ttest (x, m, varargin)
% TTEST (paired) t-test
%     For a sample X from a normal distribution with unknown mean and
%     variance, perform a t-test of the null hypothesis `mean (X) == M'.
%     Under the null, the test statistic T follows a Student
%     distribution with `DF = length (X) - 1' degrees of freedom.
%
%     TTEST treads NaNs as "Missing values" and ignores these. 
%
% H = ttest(x,m)
%	tests Null-hypothesis that mean of x is m. 		
% H = ttest(x,y)
% 	size of x and size of y must match, it is tested whether the 
%	difference x-y is significantly different to m=0; 
% H = ttest(x,y,alpha)
% H = ttest(x,y,alpha,tail)
% H = ttest(x,y,alpha,tail,DIM)
% [H,PVAL] = ttest(...)
% [H,PVAL,CI] = ttest(...)
% [H,PVAL,CI,stats] = ttest(...)
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
%     TTEST works on the first non-singleton dimension or on DIM. 
%
%     If no output argument is given, the p-value of the test is
%     displayed.
%

%       Copyright (C) 2014 Tony Richardson
%       Copyright (C) 2010,2020 by Alois Schloegl <alois.schloegl@gmail.com>
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

  % Set default arguments
  my_default = 0;
  alpha = 0.05;
  tail  = 'both';

  % Find the first non-singleton dimension of x
  DIM = min(find(size(x)~=1));
  if isempty(DIM), DIM = 1; end

  if (nargin == 1)
    m = my_default;
  end

  i = 1;
  while ( i <= length(varargin) )
    switch lower(varargin{i})
      case 'alpha'
        i = i + 1;
        alpha = varargin{i};
      case 'tail'
        i = i + 1;
        tail = varargin{i};
      case 'dim'
        i = i + 1;
        DIM = varargin{i};
      otherwise
        error('Invalid Name argument.',[]);
    end
    i = i + 1;
  end

  if ~isa(tail, 'char')
    error('tail argument to ttest must be a string\n',[]);
  end

  if any(and(~isscalar(m),size(x)~=size(m)))
    error('Arrays in paired test must be the same size.');
  end

  % Set default values if arguments are present but empty
  if isempty(m)
    m = my_default;
  end

  % This adjustment allows everything else to remain the
  % same for both the one-sample t test and paired tests.
  x = x - m;

  szx = size(x); 
  szm = size(m);	
  szx(DIM) = 1;	  
  szm(DIM) = 1;
  
  [S, N] = sumskipnan(x, DIM);
  x_bar = S./N;
  stats.df = N - 1;
  stats.sd = std (x, 0, DIM);
  x_bar_std = stats.sd./sqrt(N);
  tval = (x_bar)./x_bar_std;
  stats.tstat = tval;

  if (strcmp (tail, '~=') || strcmp (tail, '!=') || strcmp (tail, '<>')) || strcmp(tail,'both'),
    pval = 2*(1 - tcdf(abs(tval), N-1));
    tcrit = -tinv(alpha/2,N-1);
    ci = [x_bar-tcrit.*x_bar_std; x_bar+tcrit.*x_bar_std] + m;
  elseif strcmp (tail, '>') || strcmp(tail,'right'),
    pval = tcdf(tval, N-1);
    tcrit = -tinv(alpha, N-1);
    ci = [m+x_bar-tcrit.*x_bar_std; inf*ones(size(x_bar))];
  elseif strcmp (tail, '<') || strcmp(tail,'left'),
    pval = tcdf(tval, N-1);
    tcrit = -tinv(alpha,N-1);
    ci = [-inf*ones(size(x_bar)); m+x_bar+tcrit.*x_bar_std];
  else
    error ('ttest: option %s not recognized', tail);
  end

  % Reshape the ci array to match MATLAB shaping
  if and(isscalar(x_bar), DIM==2)
    ci = ci(:)';
  elseif size(x_bar,2)<size(x_bar,1)
    ci = reshape(ci(:),length(x_bar),2);
  end

  h = double(pval < alpha);
  if (nargout == 0)
    fprintf(1,'  pval: %g\n', pval);
  end

%!test
%! x = [8:0.1:12,repmat(NaN,1,40)];
%! [h, pval, ci] = ttest (x, 10);
%! assert (h, 0)
%! assert (pval, 1, 10*eps)
%! assert (ci, [9.6219 10.3781], 1E-5)
%! [h, pval, ci0] = ttest (x, 0);
%! assert (h, 1)
%! assert (pval, 0)
%! assert (ci0, ci, 2e-15)
%! [h, pval, ci] = ttest (x, 10, "tail", "right", "dim", 2, "alpha", 0.05);
%! assert (h, 0)
%! assert (pval, 0.5, 10*eps)
%! assert (ci, [9.68498 Inf], 1E-5)
