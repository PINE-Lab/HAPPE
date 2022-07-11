function [y]=mean(x,DIM,opt,W)
% MEAN calculates the mean of data elements. 
% 
%  y = mean(x [,DIM] [,opt] [, W])
%
% DIM	dimension
%	1 MEAN of columns
%	2 MEAN of rows
% 	N MEAN of  N-th dimension 
%	default or []: first DIMENSION, with more than 1 element
%
% opt	options 
%	'A' arithmetic mean
%	'G' geometric mean
%	'H' harmonic mean
%
% W	weights to compute weighted mean (default: [])
%	if W=[], all weights are 1. 
%	number of elements in W must match size(x,DIM) 
%
% usage: 
%	mean(x)
%	mean(x,DIM)
%	mean(x,opt)
%	mean(x,opt,DIM)
%	mean(x,DIM,opt)
%	mean(x,DIM,W)
%	mean(x,DIM,opt,W); '
%
% features:
% - can deal with NaN's (missing values)
% - weighting of data 
% - dimension argument also in Octave
% - compatible to Matlab and Octave
%
% see also: SUMSKIPNAN, MEAN, GEOMEAN, HARMMEAN
%

%	$Id$
%	Copyright (C) 2000-2004,2008,2009,2011 by Alois Schloegl <alois.schloegl@gmail.com>	
%    	This is part of the NaN-toolbox. For more details see
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/
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

if nargin==1,
	%------ case:  mean(x)
	W = []; 
        DIM=[]; 
        opt='a';
elseif (nargin==2)
	W = []; 
        %if ~isnumeric(DIM), %>=65;%abs('A'), 
        if (DIM>64) %abs('A'), 
		%------ case:  mean(x,opt)
                opt=DIM;
                DIM=[]; 
        elseif (DIM>length(size(x)))
		y=x; 
		return; 	
        else
		%------ case:  mean(x,DIM)
                opt='a';
        end;	
elseif (nargin == 3), 
        if isnumeric(DIM) && isnumeric(opt)
		%------ case:  mean(x,DIM,W)
		W = opt; 
		opt='a';
	elseif (DIM>64) %abs('A'), 
		%------ case:  mean(x,opt,DIM)
        	%if ~isnumeric(DIM), %>=65;%abs('A'), 
                tmp=opt;
                opt=DIM;
                DIM=tmp;
                W = []; 
        elseif (DIM>length(size(x)))
		y=x; 
		return; 	
        else 
        	%------ case:  mean(x,DIM,opt)
        	W = [];
        end;
elseif nargin==4,
		%------ case: mean(x,DIM,opt,W)
	; 
else
	help mean 
%	fprintf(1,'usage: mean(x) or mean(x,DIM) or mean(x,opt,DIM) or mean(x,DIM,opt) or mean(x,DIM,W) or mean(x,DIM,opt,W); '
end;

if isempty(opt)
	opt = 'A';
elseif any(opt=='aAgGhH')
	opt = upper(opt); % eliminate old version 
else 
	error('Error MEAN: invalid opt argument');
end; 

if  (opt == 'A')
	[y, n] = sumskipnan(x,DIM,W);
        y = y./n;
elseif (opt == 'G')
	[y, n] = sumskipnan(log(x),DIM,W);
    	y = exp (y./n);
elseif (opt == 'H')
	[y, n] = sumskipnan(1./x,DIM,W);
    	y = n./y;
else
    	fprintf (2,'mean:  option `%s` not recognized', opt);
end 

%!assert(mean([1,NaN],1),[1,NaN])
%!assert(mean([1,NaN],2),1)
%!assert(mean([+inf,-inf]),NaN)
%!assert(mean([+0,-0],'h'),NaN)
%!assert(mean([1,4,NaN],'g'),2)
      
      

