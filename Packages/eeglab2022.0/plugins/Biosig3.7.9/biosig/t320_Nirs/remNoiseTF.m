function corrSignal=remNoiseTF(signal,noise,fs,windowlength)
%
% remNoiseTF removes respiration an blood pressure related noise from
% [(de)oxy-Hb] signals by using a transfer function model [1,2]. For a detailed
% description of the model see [3].
%  
%   [corrSignal]=remNoiseTF(signal,noise,fs,shift)
%
% Input:
%   signal     ... signal with noise (either [oxy-Hb] or [deoxy-Hb])
%   noise      ... noise signal from a different source (respiration, BPdia, HR, ...    
%   fs         ... sampling frequency
% 
% Optional input parameter: 
%   windowlength ... length of one segment for the correction, default=240 seconds
%
% Output:
%   corrSignal:  ... corrected signal without influence of noise
%
%[1] Priestley, M.B., 1981. Spectral Analysis and Time Series. Vol. 1 and 2.
%    Academic Press, London, pp. 671.
%[2] Wei, W.W.S., 1990. Time Series Analysis; Univariate and Multivariate Methods. Addison Wesley, New York, pp. 289.
%[3] Florian G, Stancak A, Pfurtscheller G. Cardiac response induced by
%    voluntary selfpaced finger movement. International Journal of Psychophysiology, 28: 273-283, 1998.


% Copyright (C) 2012 by the Institute for Knowledge Discovery, Graz University of Technology 
% Rupert Ortner and Guenther Bauernfeind <g.bauernfeind@tugraz.at>
% WWW: http://bci.tugraz.at/
% $Id: remNoiseTF.m v0.1 2012-02-13 10:00:00 ISD$
%
% This software was implemented in the framework of the Styrian government project "Einfluss von
% Herz-Kreislauf-Parametern auf das Nah-Infrarot-Spektroskopie (NIRS) Signal" (A3-22.N-13/2009-8) 
% and is part of the BIOSIG-toolbox http://biosig.sf.net/
%
% LICENSE:
% 
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.

if nargin==3
    windowlength=240; %240
end


global mmax;
mmax=15;     %maximum filter order in function Korr_p (for details see [3])

seq=(fs*windowlength);
partn=ceil(length(signal)/seq);
Corr=zeros(length(signal),1); %pre-allocation of the correction term
endact=0;                      % end of actual part of signal
for i=1:partn        
    onset=(i-1)*seq+1;
    ending=i*seq;
    
    if ending > length(signal)
        ending=length(signal);
    end 
    
    if onset==1
        noise_e=[ones(mmax,1)*noise(1);noise(onset:ending)];
    else
        noise_e=noise(onset-mmax:ending);     
    end
    
    %Calculation of the correction term using partTF
    Corr_p=partTF(noise_e,signal(onset:ending));
    Corr(endact+1:endact+length(Corr_p))=Corr_p;
    endact=endact+length(Corr_p); 
end

endsignal=length(signal);
Corr(end+1:endsignal)=Corr(end);


corrSignal=signal-Corr;      % Corrected Signal  
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function Corr_p=partTF(noise_e,signal_p)
% partTF is the implementation of the transfer function equations from [3].
%
% Input
%     noise_e       ...  noise
%     signal_p      ...  signal
%    
%
% Output:
%      Corr_p       ... Correction term

% Copyright (C) 2012 by the Institute for Knowledge Discovery, Graz University of Technology 
% Rupert Ortner and Guenther Bauernfeind <g.bauernfeind@tugraz.at>
% WWW: http://bci.tugraz.at/
% $Id: partTF.m v0.1 2012-02-13 10:00:00 ISD$
%
% This software was implemented in the framework of the Styrian government project "Einfluss von
% Herz-Kreislauf-Parametern auf das Nah-Infrarot-Spektroskopie (NIRS) Signal" (A3-22.N-13/2009-8) 
% and is part of the BIOSIG-toolbox http://biosig.sf.net/
%
% LICENSE:
% 
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.


global mmax;
noise_p=noise_e(mmax+1:end);
mmin=5;
g_yy=xcov(noise_p,mmax,'biased');
g_xy = xcov(signal_p,noise_p,mmax,'biased');%cross-cov. of signal and noise 
g_xx0 = xcov(signal_p,0,'biased');          %Autocov. of signal_p

index=mmax+1;   %tau=0
lambda=zeros(1,mmax-mmin+1);
gu=cell(mmax,1);
for m=mmin:mmax
    g_xy_pj=g_xy(index:index+m);         %Gamma_xy from 0 to m
    g_yy_uj=g_yy(index:index+m);         %Gamma_yy from 0 to m
    G=toeplitz(g_yy_uj);
    gu{m}=inv(G)*g_xy_pj;          
    Snn=g_xx0-gu{m}'*g_xy_pj; 
    lambda(m-mmin+1)=length(noise_p)*log(Snn)+2*(m+1);
end
[minm,minind]=min(lambda);       %minimizing lambda
mind=mmin+minind-1;        
b=gu{mind};

S=filter(b,1,noise_e);
Corr_p=S(mmax+1:end);