%% Copyright (C) 1995, 1996, 1997  Kurt Hornik
%%
%% This file is part of Octave.
%%
%% Octave is free software; you can redistribute it and/or modify it
%% under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2, or (at your option)
%% any later version.
%%
%% Octave is distributed in the hope that it will be useful, but
%% WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%% General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with Octave; see the file COPYING.  If not, write to the Free
%% Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%% 02110-1301, USA.

%% -*- texinfo -*-
%% @deftypefn {Function File} {} binocdf (@var{x}, @var{n}, @var{p})
%% For each element of @var{x}, compute the CDF at @var{x} of the
%% binomial distribution with parameters @var{n} and @var{p}.
%% @end deftypefn

%% Author: KH <Kurt.Hornik@wu-wien.ac.at>
%% Description: CDF of the binomial distribution

function cdf = binocdf (x, n, p)

  if (nargin ~= 3)
    error;
  end

  if (~isscalar (n) || ~isscalar (p))
    [retval, x, n, p] = common_size (x, n, p);
    if (retval > 0)
      error ('binocdf: x, n and p must be of common size or scalar');
    end
  end

  sz = size (x);
  cdf = zeros (sz);

  k = find (isnan (x) | ~(n >= 0) | (n ~= round (n)) | ~(p >= 0) | ~(p <= 1));
  if (any (k))
    cdf(k) = NaN;
  end

  k = find ((x >= n) & (n >= 0) & (n == round (n)) & (p >= 0) & (p <= 1));
  if (any (k))
    cdf(k) = 1;
  end

  k = find ((x >= 0) & (x < n) & (n == round (n)) & (p >= 0) & (p <= 1));
  if (any (k))
    tmp = floor (x(k));
    if (isscalar (n) && isscalar (p))
      cdf(k) = 1 - betainc (p, tmp + 1, n - tmp);
    else
      cdf(k) = 1 - betainc (p(k), tmp + 1, n(k) - tmp);
    end
  end

end
