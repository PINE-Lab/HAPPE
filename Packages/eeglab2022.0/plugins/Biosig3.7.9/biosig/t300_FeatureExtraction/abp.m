function pos = abp(x,fs); 
% ABP is Arterial blood pressure (ABP) pulse detector 
%   The algorithm is based on [1]  
%
%  pos = abp(x,fs);
% 
%
%  x	blood pressure signal 
%  fs   sampling rate
%  pos  (tentative) position of blood pressure onset. 
%
%
% Reference(s):
% [1] Zong, W.; Heldt, T.; Moody, G.B.; Mark, R.G.
%     An open-source algorithm to detect onset of arterial blood pressure pulses
%     Computers in Cardiology, 2003
%     Volume , Issue , 21-24 Sept. 2003 Page(s): 259 - 262

%  Copyright (C) 2009 by Alois Schloegl <alois.schloegl@gmail.com>
%  This is part of the BIOSIG-toolbox http://biosig.sf.net/

% BioSig is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% Version 3 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Library General Public License for more details.
%
% You should have received a copy of the GNU Library General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc., 59 Temple Place - Suite 330,
% Boston, MA  02111-1307, USA.


% Lowpass filtering : 
% (this is equivalent to a haar wavelet)
n = 0.020*fs; 
h1 = ones(1,n); 
%h2 = cumsum([h1,-h1]); 
%y = filter(h2,1,x); 
dy = filter(h1,1,x); % this computes already the y(k)-y(k-1)

% slope sum function:
du = dy;
du(du<0)=0; 
h3 = ones(1,0.128*fs);
z  = filter(h3,1,du); 

% decision rule 
warning('algorithm is not completed, decision rule not implemented yet'); 

th = quantile(z,.9);	% fixed threshold (unlike in [1]) 
pos = gettrigger(z,th);


