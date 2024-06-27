%% Copyright (C) 2008 Alois Schloegl
%% $Id: isfinite.m,v 1.1 2008-02-01 21:14:50 schloegl Exp $
%% This function is part of BioSig http://biosig.sf.net 
%% Originally, it was part of Octave. It was modified for the use with FreeMat
%%
%% BioSig is free software; you can redistribute it and/or modify it
%% under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 3, or (at your option)
%% any later version.

%% ISFINITE returns the inverse of the error function ERF. 
%%   y = isfinite(x)
%%

function [y] = isfinite(x)
	y = (x < inf) & (x > -inf);
end
