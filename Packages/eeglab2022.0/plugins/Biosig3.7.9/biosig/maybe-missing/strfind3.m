## Copyright (C) 2004 by Alois Schloegl
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

## -*- texinfo -*-
## @deftypefn {Function File} {@var{idx} =} strfind (@var{str}, @var{pattern})
## @deftypefnx {Function File} {@var{idx} =} strfind (@var{cellstr}, @var{pattern})
## Search for @var{pattern} in the string @var{str} and return the
## starting index of every such occurrence in the vector @var{idx}.
## If there is no such occurrence, or if @var{pattern} is longer
## than @var{str}, then @var{idx} is the empty array @code{[]}.
##
## If the cell array of strings @var{cellstr} is specified instead of the
## string @var{str}, then @var{idx} is a cell array of vectors, as specified
## above.
## @end deftypefn
##
## @seealso{findstr, strmatch, strcmp, strncmp, strcmpi, strncmpi}

## Author: alois schloegl <a dot schloegl at ieee dot org>
## Created: 1 November 2004
## Adapted-By: wpyh

function idx = strfind (text, pattern)

  lp = length (pattern);

  if (nargin != 2)
    usage ("idx = strfind (text, pattern)");
  elseif (! ischar (pattern))
    error ("strfind: pattern must be a string value");
 end

  if (ischar (text))
    text_was_str = true;
    text = {text};
    sz = [1 1];
    n = 1;
  elseif (iscellstr (text))
    text_was_str = false;
    sz = size (text);
    n = numel (text);
  else
    error ("strfind: text must be a string or cell array of strings");
 end

  idx = cell (sz);
  for i = 1:n
    tmp_text = text{i};
    tmp_idx = 1:(length (tmp_text) - lp + 1);
    k = 0;
    while ((k < lp) && ! isempty (tmp_idx))
      tmp_idx = tmp_idx(tmp_text(tmp_idx + k) == pattern(++k));
    endwhile
    idx{i} = tmp_idx;
  end

  if (text_was_str)
    idx = idx{:};
 end

endfunction

# tests taken from examples in Matlab's documentation
%!shared S, cstr
%! S = "Find the starting indices of the pattern string";
%! cstr = {"How much wood would a woodchuck chuck";
%!         "if a woodchuck could chuck wood?"};
%!assert(strfind(S,"in"),[2,15,19,45]);
%!assert(isempty(strfind(S,"In")));
%!assert(strfind(S," "),[5,9,18,26,29,33,41]);
%!assert(strfind(cstr,"wood"),{[10,23];[6,28]});
# demo, also from examples in Matlab's documentation
%!demo
%! cstr = {"How much wood would a woodchuck chuck";
%!         "if a woodchuck could chuck wood?"};
%! idx = strfind (cstr, "wood");
%! idx{:,:}
