% BIOSIG/T300 contains signal processing functions for Matlab/Octave
%     It requires also the TSA-toolbox [1].
%
% Currently the following methods are supported:  
%   Multivariate Autoregressive Analysis   
%       baccala2001
%       mvar
%       mvfilter
%	mvfreqz
%   Time-varying Autoregressive spectral estimation  
%       aar, amarma
%       tfar
%	tvaar	wrapper for AAR estimators
%   Time-varying Multivariate Autoregressive Analysis
%       mvar
%       mvaar
%       tfmvar
%   EEG 
%       lumped model
%       kemp's feedback loop model
%	evoked potential 
%	arspectrum
%	cfm - cerebral function monitor
%   EMG analysis
%       Paynter
%   ECG analysis
%	qrsdetect
%       ecg_wave_analysis (not complete yet)
%	berger
%       ecgbcorr
%       qrscorr
%       tvaar
%	heartratevariability
%   SaO2 - Oxygen saturation
%	desatur
%   Bloodpressure
%	abp - arterial blood pressure pulse detection 
%   Others:
%	Brainrate (including SEF90, SEF95)
%       Hjorth
%       Bandpower
%	Wackermann
%	Hurst (Hurst coefficients)
% 
%
% REFERENCES: 
%  [1]  A. Schloegl, Time Series Analysis toolbox for Matlab and Octave. 1996-2004.
%     available online: http://www.dpmi.tu-graz.ac.at/~schloegl/matlab/tsa/download.html
%
%

%	$Id$
%	Copyright (C) 1997-2004,2007 by Alois Schloegl <alois.schloegl@gmail.com>
%	This is part of the BIOSIG project http://biosig.sf.net/

