function [HDR] = sysdetect(signal,fs)
% SYSDETECT - detection of SYS_BP_points
%
%   [HDR] = sysdetect(signal,fs)
%
% Input:
%   signal          ... BP signal data 
%   fs              ... Sampling frequency 
%
% Output:
%   HDR.EVENT.pos   ... fiducially points of systolic BP
%

% Copyright (C) 2012 by the Institute for Knowledge Discovery, Graz University of Technology 
% Guenther Bauernfeind <g.bauernfeind@tugraz.at>
% WWW: http://bci.tugraz.at/
% $Id: sysdetect.m v0.1 2012-02-13 10:00:00 ISD$
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



%Remove of the trend
signal=detrend(signal,1);

%Reduction of noise, LP filtering 10 Hz
[N,Wn]=buttord(10/(fs/2),30/(fs/2),3,60); 
[b,a]=butter(N,Wn,'low');
signal=filtfilt(b,a,signal);

%smooth
signal=smooth(signal);
%1st derivation
d_signal=diff(signal);
%tresholding
TH=std(d_signal);
k=find(d_signal>TH); 
%seperation
sep=fs/3; %Puls lies below 180 bpm
% Detection
kpp(1)=k(1);
kpn=[];
for i=1:length(k)-1
    if k(i+1)-k(i)>sep    %To wave must be separated at least by sep
        kpp=[kpp,k(i+1)]; %The point before sep is the last one of the wave 
        kpn=[kpn,k(i)];   %the next one is the starting point of the next 
                          %wave
    end
end
kpn=[kpn,k(length(k))];   %the last point of k is the last one of the last
                          %wave

for i=1:length(kpp)
    [m,n(i)]=max(d_signal(kpp(i):kpn(i)));
    pos_detect(i)=kpp(i)+n(i)-1;
end
if pos_detect(1)<=40
    pos_detect=pos_detect(2:end);
end
if pos_detect(end)>=length(signal)-40
    pos_detect=pos_detect(1:end-1);
end
for i=1:length(pos_detect)
    [m1,n1(i)]=max(signal(pos_detect(i):pos_detect(i)+40));
    pos(i)=pos_detect(i)+n1(i)-1;
end

HDR.EVENT.POS=pos;

 
