function bp = bandpower(s,arg2,arg3,arg4,mode)
% BANDPOWER calculation
%       bp = bandpower(s, Fs, bands, smoothing, mode)
%
% INPUT:
%    s          raw data, one channel per column
%    Fs         sampling rate
%    bands      each row has two elements indicating the lower and upper frequency
%               default value is [10,12;16,24] indicating two bands of
%               10-12 and 16-24 Hz.
%    smoothing  length of smoothing window in seconds. The default value is 1 [s]
% 		for mode==6 (Hilbert transform) this parameter is ignored.
%    mode       mode == 1 uses FIR filter and log10
%               mode == 2 uses Butterworth IIR filter and log10
%               mode == 3 udes FIR filter and ln
%               mode == 4 uses Butterworth IIR filter and ln
%               mode == 5 uses FFT filter and ln
%               mode == 6 uses Hilbert transform, and returns log(envelope^2)
%               the default value is mode == 4 (Butterworth filter and ln)
%
% OUTPUT:
%    bp is log(bandpower) of s
%       the order of the features is
%       [f1(#1), f1(#2),...,f1(#n), f2(#1),...,f2(#n),...,fm(#1),...,fm(#n)]
%       First, the first frequency band of all channels, is followed by
%       the the second band of all channels, until the last
%       last f-band of all channels
%

%	Copyright (C) 2007,2014,2016 by Alois Schloegl <alois.schloegl@ist.ac.at>
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
        error('signal data must be ordered in columns')
    end;
    HDR.SampleRate = arg2;
end;

if nargin<3
    F = [10,12;16,24];
else
    F = arg3;
end;
if nargin<4
    W = 1;
else
    W = arg4;
end;
if nargin<5
    mode = 4;
end;

bp = [];
tmp = double(s); tmp(isnan(tmp))=0;
if mode == 1  % FIR and log10
    for k=1:size(F,1),
        B  = fir1(HDR.SampleRate,F(k,:)/HDR.SampleRate*2);
        bp = [bp,log10(filter(ones(W*HDR.SampleRate,1),W*HDR.SampleRate,filter(B,1,tmp).^2 ))];
    end;
elseif mode == 2  % Butterworth and log10
    for k=1:size(F,1),
        [B,A] = butter(4,F(k,:)/HDR.SampleRate*2);
        bp = [bp,log10(filter(ones(W*HDR.SampleRate,1),W*HDR.SampleRate,filter(B,A,tmp).^2 ))];
    end;
elseif mode == 3  % FIR and ln
    for k=1:size(F,1),
        B  = fir1(HDR.SampleRate,F(k,:)/HDR.SampleRate*2);
        bp = [bp,log(filter(ones(W*HDR.SampleRate,1),W*HDR.SampleRate,filter(B,1,tmp).^2 ))];
    end;
elseif mode == 4  % Butterworth and ln
    for k=1:size(F,1),
        [B,A] = butter(4,F(k,:)/HDR.SampleRate*2);
        bp = [bp,log(filter(ones(W*HDR.SampleRate,1),W*HDR.SampleRate,filter(B,A,tmp).^2 ))];
    end;
elseif mode == 5  % FFT
    S = fft(tmp);
    f = [0:size(S,1)-1]*HDR.SampleRate/size(S,1);
    for k=1:size(F,1),
	SS = zeros(size(S));
	ix = (F(k,1) <= f) & (f <= F(k,2));
	SS(ix) = S(ix);
	bp = [bp, log(abs(ifft(SS)).^2) ];
    end;
elseif mode == 6  % Hilbert
    fts = fft(tmp,[],1);        % Fourier Transform
    f = HDR.SampleRate * [0:size(tmp,1)-1]'/size(tmp,1);
    for k=1:size(F,1),
	FTS = fts;
	FTS(f < F(k,1) | F(k,2) < f, :) = 0;
	HTS = 2*ifft(FTS,[],1);		% factor of 2 because only half of the spectrum is considered
	bp  = [bp, 2*log(abs(HTS))];	% factor of 2 because power (not amplitude) is returned
    end
end;
