function [cleanSignal]=adaptPulseremove(dirtySignal,Noise,fs)
% adaptPulseremove removes the influence of the Pulse-Signal from the
% [(de)oxy-Hb] signal
%
% [cleanSignal]=adaptPulseremove(dirtySignal,Noise,fs)
%
% Input:
%   dirtySignal     ... dirty signal (either [oxy-Hb] or [deoxy-Hb]) with pulse influence
%   Noise           ... artificial Noise signal (either BP signal or signal from a
%                       fingerPulse sensor
%   fs              ... Sampling frequency    
%
% Output:
%   cleanSignal     ... clean signal without influence


% Copyright (C) 2012 by the Institute for Knowledge Discovery, Graz University of Technology 
% Guenther Bauernfeind <g.bauernfeind@tugraz.at>
% WWW: http://bci.tugraz.at/
% $Id: adaptPulseremove.m v0.1 2012-02-13 10:00:00 ISD$
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

% Pre-processing (band-pass filtering) of the Noise-Signal to get only the Pulse influence
    WnPulse=[0.8 1.8];                    %Borders of the filter 
    N=200;
    b = fir1(N,WnPulse/(fs/2),'bandpass');
    Noise=filtfilt(b,1,Noise);

    Noise=Noise/max(Noise); %Normalise 

% Definition/calculation of the unknown influence (using LMS algorithm)
    calcPulseNoise = calcInfluence(dirtySignal, Noise); 

    
% Post processing (high pass filtering) to avoid removements below 0.8 Hz
    Wn=0.8;                    
    N=200;
    b = fir1(N,Wn/(fs/2),'high');
    calcPulseNoise=filtfilt(b,1,calcPulseNoise);
    
% Calculation of the cleand signal
    cleanSignal=dirtySignal-calcPulseNoise;



end

