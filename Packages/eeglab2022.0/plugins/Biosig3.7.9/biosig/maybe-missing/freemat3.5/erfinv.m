%% Copyright (C) 1994, 1995, 1996  Kurt Hornik
%% Copyright (C) 2008 Alois Schloegl
%% $Id: erfinv.m,v 1.3 2008-01-23 09:34:40 schloegl Exp $
%% This function is part of BioSig http://biosig.sf.net 
%% Originally, it was part of Octave. It was modified for the use with FreeMat
%%
%% BioSig is free software; you can redistribute it and/or modify it
%% under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 3, or (at your option)
%% any later version.

%% ERFINV returns the inverse of the error function ERF. 
%%   y = erfinv(x)
%%

function [y, iterations] = erfinv (x)

  if (nargin ~= 1)
    usage ('erfinv (x)');
  end

  maxit = 100;
  tol = 2*eps;

  iterations = 0;

  sz = size (x);

  y = zeros (size(x));

  y(~(abs(x)<= 1)) = NaN;  %% x<1, x>1, x==NaN
  y(x == -1) = -Inf;
  y(x == +1) = +Inf;

  i = find ((x > -1) & (x < 1));
  if any(i)
    s = sqrt (pi) / 2;
    z = sqrt (-log (1 - abs (x(i)))) .* sign (x(i));
    while (any (abs (erf (z) - x(i)) > tol * abs (x(i))))
      z = z - (erf (z) - x(i)) .* exp (z.^2) * s;
      iterations = iterations+1;      
      if (iterations > maxit)
        warning ('erfinv: iteration limit exceeded');
        break;
     end
   end
    y(i) = z;
  end
  y = reshape (y, sz);

end
