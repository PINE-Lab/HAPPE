function FLAG = flag_accuracy_level(i)
% FLAG_ACCURACY_LEVEL sets and gets accuracy level 
%   used in SUMSKIPNAN_MEX and COVM_MEX
%   The error margin of the naive summation is N*eps (N is the number of samples),
%   the error margin is only 2*eps if Kahan's summation is used [1].    
%
%	0: maximum speed [default]
%	   accuracy of double (64bit) with naive summation (error = N*2^-52) 
%	1: accuracy of extended (80bit) with naive summation (error = N*2^-64) 
%	2: accuracy of double (64bit) with Kahan summation (error = 2^-52) 
%	3: accuracy of extended (80bit) with Kahan summation  (error = 2^-64)  
%
%   Please note, level 3 might be equally accurate but slower than 1 or 2 on
%   some platforms. In order to determine what is good for you, you might want
%   to run ACCTEST. 
%
% FLAG = flag_accuracy_level()
% 	gets current level
% flag_accuracy_level(FLAG) 
% 	sets accuracy level  
% 
% see also: ACCTEST
% 
% Reference:
% [1] David Goldberg, 
%       What Every Computer Scientist Should Know About Floating-Point Arithmetic
%       ACM Computing Surveys, Vol 23, No 1, March 1991. 


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
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA

%	$Id$
% 	Copyright (C) 2009 by Alois Schloegl <alois.schloegl@gmail.com>
%       This function is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/


persistent FLAG_ACCURACY_LEVEL;

%% if strcmp(version,'3.6'), FLAG_ACCURACY_LEVEL=1; end;	%% hack for the use with Freemat3.6

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set the default accuracy level for your platform, ACCTEST might help to determine the optimum for your platform.  
%% If you use Matlab, use level 0 or 2; 1 and 3 are much slower but do not show a better accuracy  
%% Octave seems to be able to use all 4 levels, were the differences of accuracy between succeeding levels become smaller   
DEFAULT_ACCURACY_LEVEL = 0;	%% maximum speed, accuracy sufficient for most needs.
%% DEFAULT_ACCURACY_LEVEL = 2;	%% slower, but better accuracy for: AMDx64 Opteron, Phenom, Intel Pentium
%% DEFAULT_ACCURACY_LEVEL = 1;	%% slower, but better accuracy for: Octave on Intel Atom (no improvement with Matlab, just slower) 
%% DEFAULT_ACCURACY_LEVEL = 3;	%% similar accuracy than 1 or 2 (depending on platform) but even slower. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% set DEFAULT value of FLAG
if isempty(FLAG_ACCURACY_LEVEL),
	FLAG_ACCURACY_LEVEL = DEFAULT_ACCURACY_LEVEL;
end;

if nargin>0,
	if (i>3), i=3; end;
	if (i<0), i=0; end;
	FLAG_ACCURACY_LEVEL = double(i); 
end;
FLAG = FLAG_ACCURACY_LEVEL;

