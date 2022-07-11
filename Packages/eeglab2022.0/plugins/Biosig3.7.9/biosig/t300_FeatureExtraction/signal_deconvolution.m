function d = signal_deconvolution(r,t,fs,highpass,lowpass)
% SIGNAL_DECONVOLUTION deconvolves some raw data with some given template in order
%   to improve the detection of miniature epsc's, and ipsc's. 
%
% d = SIGNAL_DECONVOLUTION(raw,template,samplerate,highpass,lowpass)
% ... SIGNAL_DECONVOLUTION(raw,template,samplerate,[highpass,lowpass])
% ... SIGNAL_DECONVOLUTION(raw,template,samplerate,[lowpass, highpass])
%
% INPUT:
%    raw: raw data (a Nx1 data vector)
%    template: template (a Mx1 data vector)
%	it is assumed that the template starts immidiately with the first sample 		
%    samplerate: sampling rate in Hz
%    highpass: edge frequency of highpass filter in Hz, default 0.1 Hz.
%    lowpass: edge frequency of lowpass filter in Hz, default 100 Hz. 
% 	
% Output: 
%    d: detection trace 
%
% see also: get_local_maxima_above_threshold
%
% Reference(s): 
%  [1] A. Pernía-Andrade, S.P. Goswami, Y. Stickler, U. Fröbe, A. Schlögl, and P. Jonas (2012)
%     A deconvolution-based method with high sensitivity and temporal resolution for 
%     detection of spontaneous synaptic currents in vitro and in vivo.
%     Biophysical Journal Volume 103 October 2012 1–11.

%  Copyright (C) 2012,2013 by Alois Schloegl, IST Austria <alois.schloegl@ist.ac.at>
%  This is part of the BIOSIG-toolbox http://biosig.sf.net/


%% check filter settings - input arguments
if nargin>4
    if numel(highpass)==2,
	B = highpass;
    else
	B = [lowpass, highpass]; 
    end;
    B = [min(B),max(B)];
else
    B = [];
end;

%% transform into frequency domain 
H = fft(t,size(r,1));
R = fft(r); 

%% compute deconvolution in frequency domain
D = R./H;

%% filter in frequency domain
f = [0:size(r,1)-1] * fs / size(r,1);
if isempty(B),
	;
elseif 0,
    %% rectangular window in interval [B(1),B(2)]
    D(f < B(1) | B(2) < f) = 0;
    D = D*2;
elseif 1,
    %% rectangular window in intervals [B(1),B(2)] and [fs-B(2),fs-B(1)]
    D( f<B(1) | ( B(2) < f & f < (fs-B(2)) ) | (fs-B(1)) < f ) = 0;
else
    %% Gaussian window
    w = 1/sqrt(2*pi*B(2)/fs) * exp (-0.5*min([f;fs-f]/B(2),[],1).^2);
    w( f<B(1) | fs-B(1) < f ) = 0; 
    D = fs*w(:).*D;
end; 


%% convert from frequency domain into time domain. 
d = real(ifft(D));
	

