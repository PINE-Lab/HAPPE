%% Copyright (C) 1995, 1996, 1997  Kurt Hornik
%%               2006, 2011 Alois Schloegl 
%%
%% This file is part of Biosig.
%%
%% Biosig is free software; you can redistribute it and/or modify it
%% under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 3, or (at your option)
%% any later version.
%%
%% Biosig is distributed in the hope that it will be useful, but
%% WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%% General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with Octave; see the file COPYING.  If not, write to the Free
%% Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%% 02110-1301, USA.

%% -*- texinfo -*-
%% @deftypefn {Function File} {[@var{pval}, @var{z}] =} wilcoxon_test (@var{x}, @var{y}, @var{alt})
%% For two matched-pair sample vectors @var{x} and @var{y}, perform a
%% Wilcoxon signed-rank test of the null hypothesis PROB (@var{x} >
%% @var{y}) == 1/2.  Under the null, the test statistic @var{z}
%% approximately follows a standard normal distribution.
%%
%% With the optional argument string @var{alt}, the alternative of
%% interest can be selected.  If @var{alt} is @code{'~='} or
%% @code{'<>'}, the null is tested against the two-sided alternative
%% PROB (@var{x} > @var{y}) ~= 1/2.  If alt is @code{'>'}, the one-sided
%% alternative PROB (@var{x} > @var{y}) > 1/2 is considered.  Similarly
%% for @code{'<'}, the one-sided alternative PROB (@var{x} > @var{y}) <
%% 1/2 is considered.  The default is the two-sided case.
%% If @var{x} and @var{y} are matrices (must have same size), the test
%% is applied to each column.
%%
%% The p-value and z-score of the test are returned in @var{pval} and
%% @var{z}, resp.
%%
%% If no output argument is given, the p-value of the test is displayed.
%% @end deftypefn

%% Author: KH <Kurt.Hornik@wu-wien.ac.at>
%% Description: Wilcoxon signed-rank test

function [pval, z] = wilcoxon_test (x, y, alt)

  if ((nargin < 2) || (nargin > 3))
    help wilcoxon_test;
  end

  d = x - y;
  d(d==0) = NaN;	
  n = sum(~isnan(d),1);
  r = ranks (abs (d), 1);
  r(d <= 0) = NaN; 
  z = sumskipnan(r,1);	

  z = ((z - n .* (n + 1) / 4) ./ sqrt (n .* (n + 1) .* (2 * n + 1) / 24));
  z(n<=0) = 0;	

  cdf = normcdf(z,0,1);

  if (nargin == 2)
    alt = '~=';
  end

  if (strcmp (alt, '~=') || strcmp (alt, '<>') || (alt==0))
    pval = 2 * min (cdf, 1 - cdf);
  elseif (strcmp (alt, '>') || (alt==1))
    pval = 1 - cdf;
  elseif (strcmp (alt, '<') || (alt==-1))
    pval = cdf;
  else
    error ('wilcoxon_test: option %s not recognized', alt);
  end

  if (nargout == 0)
    printf ('  pval: %g\n', pval);
  end

end
