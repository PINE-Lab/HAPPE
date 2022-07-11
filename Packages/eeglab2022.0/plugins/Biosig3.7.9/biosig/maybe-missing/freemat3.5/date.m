%% DATE returns current date
%%   y = date
%%

%% Copyright (C) 2009 Alois Schloegl
%% $Id$
%% This function is part of BioSig http://biosig.sf.net 
%% Originally, it was part of Octave. It was modified for the use with FreeMat
%%
%% BioSig is free software; you can redistribute it and/or modify it
%% under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 3, or (at your option)
%% any later version.


function [A] = date()
	t = clock;
	month = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};	
	A = sprintf('%02i-%3s-%04i',t(3),month{t(2)},t(1));
end
