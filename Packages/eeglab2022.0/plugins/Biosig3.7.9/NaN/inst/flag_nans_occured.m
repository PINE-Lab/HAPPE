function [flag]=flag_nans_occured()
% FLAG_NANS_OCCURED checks whether the last call(s) to sumskipnan or covm 
% contained any not-a-numbers in the input argument. Because many other 
% functions like mean, std, etc. are also using sumskipnan, 
% also these functions can be checked for NaN's in the input data. 
% 
% A call to FLAG_NANS_OCCURED() resets also the flag whether NaN's occured. 
% Only sumskipnan or covm can set the flag again. 
%
% see also: SUMSKIPNAN, COVM

%	$Id$
%	Copyright (C) 2009 by Alois Schloegl <alois.schloegl@gmail.com>	
%       This function is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/

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

global FLAG_NANS_OCCURED;

%%% check whether FLAG was already defined 
if isempty(FLAG_NANS_OCCURED),
	FLAG_NANS_OCCURED = logical(0);  % default value 
end;

flag = FLAG_NANS_OCCURED;		% return value 

FLAG_NANS_OCCURED = logical(0);		% reset flag

return;
