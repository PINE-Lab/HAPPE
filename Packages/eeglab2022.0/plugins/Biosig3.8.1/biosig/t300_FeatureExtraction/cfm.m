function c = cfm(s,arg2,arg3,arg4,mode)
% CFM - cerebral function monitor 
%       c = CFM(filename)
%       c = CFM(s,Fs)
%       c = CFM(s,HDR)

% INPUT:
%    s          raw data, one channel per column
%    Fs         sampling rate
%    HDR        header information  
%    filename   filename 
% 
% OUTPUT:
%    c  s is filtered (2-12Hz), rectified, and downsampled to (1/min)


%	$Id$
%	Copyright (C) 2009 by Alois Schloegl <alois.schloegl@gmail.com>
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

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
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


if ischar(s) 
        [s,HDR]=sload(s);
else
        if size(s,1)<size(s,2),
                warning('signal data must be ordered in columns')
        end;
        if isscalar(arg2)
               HDR.SampleRate = arg2;
        end;        
end;

W = 60; % window length is 60 s
B = fir1(HDR.SampleRate,[2,12]/HDR.SampleRate*2);
c = rs(abs(filter(B,1,s)),W*HDR.SampleRate,1); 
% bp = log(filter( ones(W*HDR.SampleRate,1), W*HDR.SampleRate, abs(filter(B,1,tmp)) ))];


