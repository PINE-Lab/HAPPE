function [RC,efinal] = ac2rc(AC);
% converts the autocorrelation function into reflection coefficients 
% [RC,r0] = ac2rc(r)
%
% see also ACOVF ACORF AR2RC RC2AR DURLEV AC2POLY, POLY2RC, RC2POLY, RC2AC, AC2RC, POLY2AC
% 

%       $Id$
%       Copyright (C) 1998-2002,2008 by Alois Schloegl <alois.schloegl@gmail.com>
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.



if all(size(AC))>1,
        fprintf(2,'Error AC2RC: "r" must be a vector\n');
        return;
end;

mfilename='AC2RC';
if ~exist('durlev','file')
        fprintf(2,'Error %s: DURLEV.M not found. \n Download TSA toolbox from http://pub.ist.ac.at/~schloegl/matlab/tsa/\n',mfilename);
        return;
end;

[AR,RC,PE] = durlev(AC(:).');
RC=-RC;
efinal=AC(1);
