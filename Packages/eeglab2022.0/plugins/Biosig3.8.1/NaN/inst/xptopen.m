% XPTOPEN read of several file formats and writing of the SAS Transport Format (*.xpt)
%
%        X = xptopen(filename)
%        X = xptopen(filename,'r')
%                read filename and return variables in struct X
%        Supported are ARFF, SAS-XPT and STATA files.
%
%        X = xptopen(filename,'w',X)
%                save fields of struct X in filename.
%
%        The fields of X must be column vectors of equal length.
%        Each vector is either a numeric vector or a cell array of strings.
%
% The SAS-XPT format stores Date/Time as numeric value counting the number of days since 1960-01-01.

%    Copyright (C) 2015 by Alois Schloegl <alois.schloegl@gmail.com>     
%    This is part of the NaN-toolbox. For more details see
%    https://pub.ist.ac.at/~schloegl/matlab/NaN/
%
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; If not, see <http://www.gnu.org/licenses/>.
