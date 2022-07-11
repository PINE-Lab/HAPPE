function linBP = calcBPlin(signal,fs,typ)
% calcBPlin calculates the linear interpolated BPdia or BPsys signal and is partially based
% on the calcHR function written by Clemens Brunner.
%
%
% [linBP] = calcBPlin(signal,fs,typ)
%
% Input:
%   signal   ... Raw continuous BP signal
%   fs       ... Sampling frequency
%   typ      ... [1] systolic BP
%            ... [2] diastolic BP
%
%
% Output:
%   linBP    ... linear interpolated BPdia or BPsys signal
%

% Copyright (C) 2012 by the Institute for Knowledge Discovery, Graz University of Technology 
% Guenther Bauernfeind <g.bauernfeind@tugraz.at>
% WWW: http://bci.tugraz.at/
% $Id: calcBPlin.m v0.1 2012-02-13 10:00:00 ISD$
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


% Calculate BPsys or BPdia

if typ==1
    h_BP = sysdetect(signal, fs);
end
if typ==2
    h_BP = diadetect(signal, fs);
end
time = h_BP.EVENT.POS;

timesec=time./fs;


bp=signal(time); 

%%%%%%linear interpolation%%%%%%
X=timesec;                              %old time axis
XI = (timesec(1):1/fs:timesec(end))';   %new time axis for interpolation
linBP = interp1(X,bp,XI,'linear');      %HR linear interpolated
linBP=[linBP(1)*ones(time(1)-1,1);linBP;linBP(end)*ones(length(signal)-time(end),1)];

