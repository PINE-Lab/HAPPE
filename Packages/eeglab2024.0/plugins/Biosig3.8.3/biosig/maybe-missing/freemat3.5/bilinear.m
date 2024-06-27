%% Copyright (C) 2008 Alois Schloegl
%% $Id: bilinear.m,v 1.1 2008-01-19 21:55:45 schloegl Exp $
%% This function is part of BioSig http://biosig.sf.net 
%% Originally, it was part of Octave. It was modified for the use with FreeMat
%%
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2 of the License, or
%% (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

%% usage: [Zz, Zp, Zg] = bilinear(Sz, Sp, Sg, T)
%%        [Zb, Za] = bilinear(Sb, Sa, T)
%%
%% Transform a s-plane filter specification into a z-plane
%% specification. Filters can be specified in either zero-pole-gain or
%% transfer function form. The input form does not have to match the
%% output form. T is the sampling frequency represented in the z plane.
%%
%% Theory: Given a piecewise flat filter design, you can transform it
%% from the s-plane to the z-plane while maintaining the band edges by
%% means of the bilinear transform.  This maps the left hand side of the
%% s-plane into the interior of the unit circle.  The mapping is highly
%% non-linear, so you must design your filter with band edges in the
%% s-plane positioned at 2/T tan(w*T/2) so that they will be positioned
%% at w after the bilinear transform is complete.
%%
%% The following table summarizes the transformation:
%%
%% Transform         Zero at x                  Pole at x
%% ----------------  -------------------------  ------------------------
%% Bilinear          zero: (2+xT)/(2-xT)        pole: (2+xT)/(2-xT)
%%      2 z-1        pole: -1                   zero: -1
%% S -> - ---        gain: (2-xT)/T             gain: (2-xT)/T
%%      T z+1
%% ----------------  -------------------------  ------------------------
%%
%% With tedious algebra, you can derive the above formulae yourself by
%% substituting the transform for S into H(S)=S-x for a zero at x or
%% H(S)=1/(S-x) for a pole at x, and converting the result into the
%% form:
%%
%%    H(Z)=g prod(Z-Xi)/prod(Z-Xj)
%%
%% Please note that a pole and a zero at the same place exactly cancel.
%% This is significant since the bilinear transform creates numerous
%% extra poles and zeros, most of which cancel. Those which do not
%% cancel have a "fill-in" effect, extending the shorter of the sets to
%% have the same number of as the longer of the sets of poles and zeros
%% (or at least split the difference in the case of the band pass
%% filter). There may be other opportunistic cancellations but I will
%% not check for them.
%%
%% Also note that any pole on the unit circle or beyond will result in
%% an unstable filter.  Because of cancellation, this will only happen
%% if the number of poles is smaller than the number of zeros.  The
%% analytic design methods all yield more poles than zeros, so this will
%% not be a problem.
%%
%% References:
%%
%% Proakis & Manolakis (1992). Digital Signal Processing. New York:
%% Macmillan Publishing Company.

%% Author: Paul Kienzle <pkienzle@users.sf.net>
%% Adapted for FreeMat by: AS <alois.schloegl@gmail.com> Jan 2008

function [Zz, Zp, Zg] = bilinear(Sz, Sp, Sg, T)

  if nargin==3
    T = Sg;
    [Sz, Sp, Sg] = tf2zp(Sz, Sp);
  elseif nargin~=4
    usage('[Zz, Zp, Zg]=bilinear(Sz,Sp,Sg,T) or [Zb, Za]=blinear(Sb,Sa,T)');
  end;

  p = length(Sp);
  z = length(Sz);
  if z > p || p==0
    error('bilinear: must have at least as many poles as zeros in s-plane');
  end

%% ----------------  -------------------------  ------------------------
%% Bilinear          zero: (2+xT)/(2-xT)        pole: (2+xT)/(2-xT)
%%      2 z-1        pole: -1                   zero: -1
%% S -> - ---        gain: (2-xT)/T             gain: (2-xT)/T
%%      T z+1
%% ----------------  -------------------------  ------------------------
  Zg = real(Sg * prod((2-Sz*T)/T) / prod((2-Sp*T)/T));
  Zp = (2+Sp*T)./(2-Sp*T);
  if isempty(Sz)
    Zz = -ones(size(Zp));
  else
    Zz = [(2+Sz*T)./(2-Sz*T)];
    Zz = postpad(Zz, p, -1);
  end

  if nargout==2, [Zz, Zp] = zp2tf(Zz, Zp, Zg); end;

end
