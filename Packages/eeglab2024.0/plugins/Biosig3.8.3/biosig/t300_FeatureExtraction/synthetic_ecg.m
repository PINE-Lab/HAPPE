function ecg = synthetic_ecg(heart_rate,heart_peak,durration,srate)
% function ecg = synthetic_ecg(heart_rate,heart_peak,durration,srate)
% heart_rate is the desired frequency in beats per minute, and heart_peak
% describes the R- amplitude of the generated ECG in microvolts.
% "durration" is the durration of the recording in seconds, and finally 
% srate is the sampling rate.

%	$Revision: 1.1 $
%	$Id$
%	Copyright (c) 2004 by George Townsend
%	townsend@auc.ca	

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% version 2 of the License, or (at your option) any later version.
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

heart_rate=heart_rate/60;
len=durration*srate;

% The following information was produced from the analysis of a real two minute
% long ECG recording. The recording was detrended, and notch filtered, and
% all resulting heartbeats were superimposed and averaged. The result was
% then duplicated to produce a string of a hundred identical heartbeats
% which were then fft'd. Only the fundamental frequency and integral
% harmonics were present in the final result, and the 30 element vectors
% below capture the results in the following format:
% heart_ampl=[ fundamental 2ndHarmonic 3rdHarmonic ... 30thHarmonic ]
% heart_phase=[ fund_phase 2ndHarm_phase ... 30thHarmonic_phase ]

	heart_ampl=[0.84 0.73 0.49 0.93 1.00 0.71 0.75 0.74 0.73 0.74 0.70 0.69 ...
			    0.64 0.56 0.51 0.47 0.43 0.38 0.33 0.27 0.23 0.19 0.17 0.14 ...
			    0.11 0.09 0.07 0.06 0.04 0.04 0.03 0.02 0.01 0.01 0.01];
	heart_phase=[-2.37 -0.35  1.39 -1.40  1.38 -1.84  0.89 -2.76  0.17  3.00 ...
			     -0.47  2.40 -1.00  1.83 -1.64  1.15 -2.29  0.55 -2.88 -0.06 ...
			      2.71 -0.80  1.98 -1.48  1.29 -2.24  0.52 -3.02 -0.27  2.47 ...
			     -1.03 1.69 -1.87   0.88 -2.67];
% Produce a time vector for the ECG     
t=1:len;
t=t-1;
t=t/srate;
	
% Now produce the appropriate amplitude and phase for each harmonic then scale the result    
ecg = heart_ampl * cos(2*pi*heart_rate*(1:length(heart_ampl)).'*t + heart_phase(ones(1,length(t)),:).');
ecg = heart_peak * ecg / (max(ecg)-min(ecg));    
