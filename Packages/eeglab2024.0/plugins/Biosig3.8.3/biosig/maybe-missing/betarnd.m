## Copyright (C) 1995, 1996, 1997  Kurt Hornik
##
## This file is part of Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, write to the Free
## Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.

## -*- texinfo -*-
## @deftypefn {Function File} {} betarnd (@var{a}, @var{b}, @var{r}, @var{c})
## @deftypefnx {Function File} {} betarnd (@var{a}, @var{b}, @var{sz})
## Return an @var{r} by @var{c} or @code{size (@var{sz})} matrix of 
## random samples from the Beta distribution with parameters @var{a} and
## @var{b}.  Both @var{a} and @var{b} must be scalar or of size @var{r}
##  by @var{c}.
##
## If @var{r} and @var{c} are omitted, the size of the result matrix is
## the common size of @var{a} and @var{b}.
## @end deftypefn

## Author: KH <Kurt.Hornik@wu-wien.ac.at>
## Description: Random deviates from the Beta distribution

function rnd = betarnd (a, b, r, c)

  if (nargin > 1)
    if (!isscalar(a) || !isscalar(b)) 
      [retval, a, b] = common_size (a, b);
      if (retval > 0)
	error ("betarnd: a and b must be of common size or scalar");
     end
   end
 end

  if (nargin == 4)
    if (! (isscalar (r) && (r > 0) && (r == round (r))))
      error ("betarnd: r must be a positive integer");
   end
    if (! (isscalar (c) && (c > 0) && (c == round (c))))
      error ("betarnd: c must be a positive integer");
   end
    sz = [r, c];

    if (any (size (a) != 1)
	&& (length (size (a)) != length (sz) || any (size (a) != sz)))
      error ("betarnd: a and b must be scalar or of size [r,c]");
   end
  elseif (nargin == 3)
    if (isscalar (r) && (r > 0))
      sz = [r, r];
    elseif (isvector(r) && all (r > 0))
      sz = r(:)';
    else
      error ("betarnd: r must be a postive integer or vector");
   end

    if (any (size (a) != 1)
	&& (length (size (a)) != length (sz) || any (size (a) != sz)))
      error ("betarnd: a and b must be scalar or of size sz");
   end
  elseif (nargin == 2)
    sz = size(a);
  else
    print_usage ();
 end

  if (isscalar(a) && isscalar(b))
    if (find (!(a > 0) | !(a < Inf) | !(b > 0) | !(b < Inf)))
      rnd = NaN * ones (sz);
    else
      rnd = betainv (rand(sz), a, b);
   end
  else
    rnd = zeros (sz);

    k = find (!(a > 0) | !(a < Inf) | !(b > 0) | !(b < Inf));
    if (any (k))
      rnd(k) = NaN * ones (size (k));
   end

    k = find ((a > 0) & (a < Inf) & (b > 0) & (b < Inf));
    if (any (k))
      rnd(k) = betainv (rand (size (k)), a(k), b(k));
   end
 end

endfunction
