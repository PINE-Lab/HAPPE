%% Copyright (C) 2008 Alois Schloegl
%% $Id: tf2zp.m,v 1.1 2008-01-19 22:08:37 schloegl Exp $
%% This function is part of BioSig http://biosig.sf.net 
%% Originally, it was part of Octave. It was modified for the use with FreeMat
%%
%% Octave is free software; you can redistribute it and/or modify it
%% under the terms of the GNU General Public License as published by the
%% Free Software Foundation; either version 2, or (at your option) any
%% later version.
%%
%% Octave is distributed in the hope that it will be useful, but WITHOUT
%% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
%% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
%% for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with Octave; see the file COPYING.  If not, write to the Free
%% Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%% 02110-1301 USA.

%% -*- texinfo -*-
%% @deftypefn {Function File} {[@var{zer}, @var{pol}, @var{k}] =} tf2zp (@var{num}, @var{den})
%% Converts transfer functions to poles-and-zero representations.
%%
%% Returns the zeros and poles of the @acronym{SISO} system defined 
%% by @var{num}/@var{den}.
%% @var{k} is a gain associated with the system zeros.
%% @end deftypefn

%% Author: A. S. Hodel <a.s.hodel@eng.auburn.edu>
%% Adapted for FreeMat by: AS <alois.schloegl@gmail.com> Jan 2008

function [zer, pol, k] = tf2zp (num, den)

  if (nargin == 2)
    if (length (den) > 1)
      pol = roots (den);
    else
      pol = [];
    end

    if (length (num) > 1)
      zer = roots (num);
    else
      zer = [];
    end
  else
    print_usage ();
  end

  k = num(1) / den(1);

