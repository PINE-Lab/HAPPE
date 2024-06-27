function [N] = numel(X)
%% NUMEL returns number of elements
%%

%% This program is free software; you can redistribute it and/or
%% modify it under the terms of the GNU General Public License
%% as published by the Free Software Foundation; either version 2
%% of the License, or (at your option) any later version.
%% 
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%% 
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

%%	$Id: numel.m,v 1.1 2007-07-19 15:50:32 schloegl Exp $
%%	Copyright (C) 2005 by Alois Schloegl <alois.schloegl@gmail.com>
%%      This function is part of BIOSIG http://biosig.sf.net/

N = prod(size(X));







