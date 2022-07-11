function resamp = eeg_autocorr_fftw(EEG, pct_data)

if exist('resample') ~= 2
     error('Signal Processing Toolbox absent: You need the signal processing toolbox to run ICLabel');  
end

if ~exist('pct_data', 'var') || isempty(pct_data)
    pct_data = 100;
end

nfft = 2^nextpow2(2*EEG.pnts-1);

% calc autocorrelation
fftw('planner', 'hybrid');
ac = zeros(size(EEG.icaact, 1), nfft);
for it = 1:size(EEG.icaact, 1)
    X = fft(EEG.icaact(it, :, :), nfft, 2);
    ac(it, :) = mean(abs(X).^2, 3);
end
ac = ifft(ac, [], 2);
if EEG.pnts < EEG.srate
    ac = [ac(:, 1:EEG.pnts, :) zeros(size(ac, 1), EEG.srate - EEG.pnts + 1)];
else
    ac = ac(:, 1:EEG.srate + 1, :);
end

% normalize by 0-tap autocorrelation
ac = bsxfun(@rdivide, ac(:, 1:EEG.srate + 1, :), ac(:, 1)); 

% resample to 1 second at 100 samples/sec
resamp = resample(ac', 100, EEG.srate)';
resamp(:, 1) = [];