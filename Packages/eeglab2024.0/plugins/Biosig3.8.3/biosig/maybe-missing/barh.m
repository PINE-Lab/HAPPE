## Copyright (C) 1996, 1997 John W. Eaton
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
## @deftypefn {Function File} {} barh (@var{x}, @var{y})
## Given two vectors of x-y data, @code{bar} produces a horizontal bar graph.
##
## If only one argument is given, it is taken as a vector of y-values
## and the x coordinates are taken to be the indices of the elements.
##
## If two output arguments are specified, the data are generated but
## not plotted.  For example,
##
## @example
## barh (x, y);
## @end example
##
## @noindent
## and
##
## @example
## [xb, yb] = barh (x, y);
## plot (xb, yb);
## @end example
##
## @noindent
## are equivalent.
## @seealso{bar, plot, semilogx, semilogy, loglog, polar, mesh, contour,
## stairs, replot, xlabel, ylabel, title}
## @end deftypefn

## Author: jwe

function [xb, yb] = bar (x, y)

  if (nargin == 1)
    if (isvector (x))
      len = 3 * length (x) + 1;
      tmp_xb = tmp_yb = zeros (len, 1);
      tmp_xb(1) = 0.5;
      tmp_yb(1) = 0;
      k = 1;
      for i = 2:3:len
        tmp_xb(i) = k-0.5;
        tmp_xb(i+1) = k+0.5;
        tmp_xb(i+2) = k+0.5;
        tmp_yb(i) = x(k);
        tmp_yb(i+1) = x(k);
        tmp_yb(i+2) = 0.0;
        k++;
      end
    else
      error ("bar: argument must be a vector");
   end
  elseif (nargin == 2)
    if (isvector (x) && isvector (y))
      xlen = length (x);
      ylen = length (y);
      if (xlen == ylen)
        len = 3 * xlen + 1;
        tmp_xb = tmp_yb = zeros (len, 1);
        cutoff = zeros (1, xlen);
        for i = 1:xlen-1
          cutoff(i) = (x(i) + x(i+1)) / 2.0;
        end
        delta_p = cutoff(1) - x(1);
        delta_m = delta_p;
        tmp_xb(1) = x(1) - delta_m;
        tmp_yb(1) = 0.0;
        k = 1;
        for i = 2:3:len
          tmp_xb(i) = tmp_xb(i-1);
          tmp_xb(i+1) = x(k) + delta_p;
          tmp_xb(i+2) = tmp_xb(i+1);
          tmp_yb(i) = y(k);
          tmp_yb(i+1) = y(k);
          tmp_yb(i+2) = 0.0;
          if (k < xlen)
            if (x(k+1) < x(k))
              error ("barh: x vector values must be in ascending order");
           end
            delta_m = x(k+1) - cutoff(k);
            k++;
            if (k < xlen)
              delta_p = cutoff(k) - x(k);
            else
              delta_p = delta_m;
           end
         end
        end
      else
        error ("barh: arguments must be the same length");
     end
    else
      error ("barh: arguments must be vectors");
   end
  else
    usage ("[xb, yb] = barh (x, y)");
 end

  if (nargout == 0)
    % modified by Alois Schloegl
    plot (tmp_yb, tmp_xb);
  else
    xb = tmp_xb;
    yb = tmp_yb;
 end

endfunction
