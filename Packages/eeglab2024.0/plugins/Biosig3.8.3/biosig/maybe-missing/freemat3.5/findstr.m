%% Copyright (C) 1996 Kurt Hornik
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
%% @deftypefn {Function File} {} findstr (@var{s}, @var{t}, @var{overlap})
%% Return the vector of all positions in the longer of the two strings
%% @var{s} and @var{t} where an occurrence of the shorter of the two starts.
%% If the optional argument @var{overlap} is nonzero, the returned vector
%% can include overlapping positions (this is the default).  For example,
%%
%% @example
%% findstr ("ababab", "a")
%%      @result{} [ 1, 3, 5 ]
%% findstr ("abababa", "aba", 0)
%%      @result{} [ 1, 5 ]
%% @end example
%% @end deftypefn

%% Note that this implementation swaps the strings if second one is longer
%% than the first, so try to put the longer one first.
%%
%% Author: Kurt Hornik <Kurt.Hornik@wu-wien.ac.at>
%% Adapted-By: jwe

function v = findstr (s, t, overlap)

  if (nargin < 2 || nargin > 3)
    usage ('findstr (s, t, overlap)');
  end

  if (all (size (s) > 1) || all (size (t) > 1))
    error ('findstr: arguments must have only one non-singleton dimension');
  end

  if (nargin == 2)
    overlap = 1;
  end

  %% Make S be the longer string.
  if (length (s) < length (t))
    tmp = s;
    s = t;
    t = tmp;
  end
  
  l_s = length (s);
  l_t = length (t);
  
  if (l_t == 0)
    %% zero length target: return empty set
    v = [];
    
  elseif (l_t == 1)
    %% length one target: simple find
    v = find (s == t);
    
  elseif (l_t == 2)
    %% length two target: find first at i and second at i+1
    v = find (s(1:l_s-1) == t(1) & s(2:l_s) == t(2));
    
  else
    %% length three or more: match the first three by find then go through
    %% the much smaller list to determine which of them are real matches
    limit = l_s - l_t + 1;
    v = find (s(1:limit) == t(1) & s(2:limit+1) == t(2) & s (3:limit+2) == t(3));
  end

  %% Need to search the index vector if our find was too short
  %% (target length > 3), or if we don't allow overlaps.  Note though
  %% that there cannot be any overlaps if the first character in the
  %% target is different from the remaining characters in the target,
  %% so a single character, two different characters, or first character
  %% different from the second two don't need to be searched.
  if ((l_t >= 3) || (~overlap && l_t > 1 && any (t(1) == t(2:l_t))))
  
    %% force strings to be both row vectors or both column vectors
    if (all (size (s) ~= size (t)))
      t = t.';
    end
    
    %% determine which ones to keep
    keep = zeros (size (v));
    ind = 0:l_t-1;
    if (overlap)
      for idx = 1:length (v)
	keep(idx) = all (s(v(idx) + ind) == t);
      end
    else
      next = 1; %% first possible position for next non-overlapping match
      for idx = 1:length (v)
	if (v(idx) >= next && s(v(idx) + ind) == t)
	  keep(idx) = 1;
	  next = v(idx) + l_t; %% skip to the next possible match position
	else
	  keep(idx) = 0;
	end
      end
    end
    if (~ isempty (v))
      v = v(find (keep));
    end
  end

  if (isempty (v))
    v = [];
  end

  %% Always return a column vector, because that's what the old one did
  if (size(v,1) > 1) 
    v = v.';
  end

end
