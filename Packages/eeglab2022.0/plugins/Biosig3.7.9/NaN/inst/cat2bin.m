function [B,BLab]=cat2bin(D, Label, MODE)
% CAT2BIN converts categorical into binary data 
%   each category of each column in D is converted into a logical column
% 
%   B = cat2bin(C); 
%   [B,BinLabel] = cat2bin(C,Label); 
%   [B,BinLabel] = cat2bin(C,Label,MODE)
%
%  C        categorical data 
%  B        binary data 
%  Label    description of each column in C
%  BinLabel description of each column in B
%  MODE     default [], ignores NaN
%           'notIgnoreNAN' includes binary column for NaN 
%           'IgnoreZeros'  zeros do not get a separate category 
%           'IgnoreZeros+NaN' zeros and NaN are ignored
%
%  example: 
%     cat2bin([1;2;5;1;5]) results in 
%             1     0     0
%             0     1     0
%             0     0     1
%             1     0     0
%             0     0     1

%	$Id$
%	Copyright (C) 2009 by Alois Schloegl <alois.schloegl@gmail.com>
%       This function is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
% of the  License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street - Fifth Floor, Boston, MA 02110-1301, USA.

if nargin<3,
         MODE = []; 
end; 

% convert data 
B = []; 

c     = 0; 
k1    = 0; 
BLab  = [];
for m = 1:size(D,2) 
        h = histo_mex(D(:,m));
        x = h.X(h.H>0);
        if strcmpi(MODE,'notIgnoreNaN')
                ;
        elseif strcmpi(MODE,'IgnoreZeros')
                x = x(x~=0);
        elseif strcmpi(MODE,'IgnoreZeros+NaN')
                x = x((x~=0) & (x==x));
        else 
                x = x(x==x);
        end; 
        for k = 1:size(D,1),
                if ~isnan(D(k,m))
                        B(k, c + find(D(k,m)==x)) = 1;
                elseif isnan(x(end)),
                        B(k, c + length(x)) = 1;
                end;
        end;

        c = c + length(x);
        if nargout>1,
                for k = 1:length(x),
                        k1 = k1+1;
                        if isempty(Label)
                                BLab{k1} = ['#',int2str(m),':',int2str(x(k))];
                        else
                                BLab{k1} = [Label{m},':',int2str(x(k))];
                        end;
                end;
        end;
end;


%!assert(cat2bin([1;2;5;1;5]),[1,0,0;0,1,0;0,0,1;1,0,0;0,0,1])

