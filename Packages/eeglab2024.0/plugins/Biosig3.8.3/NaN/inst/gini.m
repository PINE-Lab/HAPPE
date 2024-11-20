function G = gini(data)
% GINI computes the gini-coefficient [1] using by
%   computing the L-moments [2]. 
%
% USAGE: 
%   G = gini(data)
% 
% 
% 
% References: 
% [1] https://en.wikipedia.org/wiki/Gini_coefficient
% [2] https://en.wikipedia.org/wiki/L-moment 

% Copyright (C) 2019,2020 by Alois Schlögl <alois.schloegl@gmail.com>
%       This function is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/

%% This program is free software: you can redistribute it and/or
%% modify it under the terms of the GNU General Public License as
%% published by the Free Software Foundation, either version 3 of the
%% License, or (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful, but
%% WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%% General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; see the file COPYING.  If not, see
%% <http://www.gnu.org/licenses/>.


L = lmom(data,2);
G = L(2)/L(1);

