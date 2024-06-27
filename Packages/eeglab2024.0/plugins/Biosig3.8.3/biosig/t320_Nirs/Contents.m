% 
% T320: Removal of Physiolgical Artifacts from NIRS
% ---------------------------------------------
%	adaptPulsremove.m     removes the Pulse influence from NIRS signals using calcInfluence.m
% 	calcInfluence.m       calculates the unknown influence of the Pulse-Noise using the LMS algorithm
% 	remNoiseTF.m          removes respiration and MayerWave artefacts by the transfer function approach
%   remNoiseICA.m         removes respiration and MayerWave artefacts by ICA approach
% 	remNoiseCAR.m         removes respiration and MayerWave artefacts by CAR approach
%
% Additional files
% -------------------
%   calcNIRSspectra.m     calculates the Spectrum of the NIRS signal
%   calcBPlin.m           calculates the linear interpolated BPdia or BPsys signal using
%   sysdetect.m           calculates fiducially points of systolic BP
%   diadetect.m           calculates fiducially points of diatolic BP
%   Illustration_multichannel_spectra.m   uses the calcNIRSspectra function to calculate and illustrates the [(de)oxy-Hb] spectra
%                                       of each used NIRS channel from a 3*11 measurement grid



 
% Copyright (C) 2012 by Guenther Bauernfeind <g.bauernfeind@tugraz.at>
% WWW: http://bci.tugraz.at/bauernfeind/
% $Id: Contents.m  2012-02-13 10:00:00 bauernfeind $
%
% This software was written in the framework of the Styrian government project "Einfluss von
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

