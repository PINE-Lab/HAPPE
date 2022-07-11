function r=calcNIRSspectra(signal,fs)
% calcNIRSspectra calculates the Spectrum of the NIRS signal
%
% r=calcNIRSspectra(signal,fs)
%
% Input:
%   signal  ... NIRS signal (either [oxy-Hb] or [deoxy-Hb] 
%   fs      ... Sampling frequency    
%
% Output:
%   r       ... structure containing the spectrum.


% Copyright (C) 2012 by the Institute for Knowledge Discovery, Graz University of Technology 
% Guenther Bauernfeind <g.bauernfeind@tugraz.at>
% WWW: http://bci.tugraz.at/
% $Id: calcNIRSspectra.m v0.1 2012-02-13 10:00:00 ISD$
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

% Used setting
windowtyp= 'Hanning';   % Window typ
window_length = 100;    % Windowlength in [s]
window_overlap = 50;    % Define Overlap in [s]
fft_length     = 200;



    if strcmp(windowtyp, 'Hamming') == 1
        window = hamming(window_length*fs);
    elseif strcmp(windowtyp, 'Hanning') == 1
        window = hann(window_length*fs);
    elseif strcmp(windowtyp, 'Triang') == 1
        window = triang(window_length*fs);
    elseif strcmp(windowtyp, 'Bartlett') == 1
        window = bartlett(window_length*fs);
    end



 
[p, f] = psd(signal-mean(signal), fft_length*fs, fs, window, window_overlap*fs,'linear');
    

r{1}.p = p;
r{1}.f = f;
r{1}.fs = fs;




end
