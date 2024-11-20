% APPLY_GENERAL_CLASSIFIER shows how to apply a precomputed classifier
% to new data sets described in [1] and applied in [2]. It is assumed that
% that classifier has been obtained with DEMO_MODY, such that the same
% preprocessing method is used.
%
% In order to run this demo, the following software packages
% are needed:
% - Octave (or Matlab) and the Signal Processing toolbox or package
% - Biosig [3]
% - NaN-toolbox [4]
% - SigViewer [5] - for scoring the data and visualization of the
%    resulting detections
% - Data (recordings and scorings from [1,2] are available [6]
%
% Copyright (C) 2016-2023 Alois Schlögl, Institue of Science and Technology Austria (ISTA)
%
% References:
% [1] Zhang X, Schlögl A, Vandael D, Jonas P,
%     MOD: A novel machine-learning optimal-filtering method for accurate and efficient
%        detection of subthreshold synaptic events in vivo.
%     Journal of Neuroscience Methods, 2021.
%     doi:10.1016/j.jneumeth.2021.109125
% [2] Xiaomin Zhang, Alois Schlögl, Peter Jonas
%     Selective Routing of Spatial Information Flow from Input to Output in Hippocampal Granule Cells, Neuron, 2020.
%     https://doi.org/10.1016/j.neuron.2020.07.006
%     http://www.sciencedirect.com/science/article/pii/S0896627320305237
% [3] Biosig - an Free and Ppen source software library for biomedical signal processing
%     http://biosig.sourceforge.net/index.html
% [4] The NaN-toolbox: A statistics and machine learning toolbox for Octave and Matlab®
%     for data with and w/o MISSING VALUES encoded as NaN's.
%     https://pub.ist.ac.at/~schloegl/matlab/NaN/
% [5] SigViewer (v0.51, or v0.65)
%     https://git.ist.ac.at/alois.schloegl/sigviewer.git
%     https://pub.ist.ac.at/~schloegl/software/sigviewer/sigviewer-0.6.4-win64.exe
% [6] https://pub.ist.ac.at/~schloegl/software/mod/

% This script tries to do all the analysis at once using
%   - raw data
%   - scorings from both exports
%   - of all 12 data sets
%
% This allows for an efficient computation as several processing steps
%  (e.g. loading, or preprocessing, or feature extraction) are done only once
% while answering different questions (implemented by different XV procedures)
%
% Things considered:
% - time lag between two export scores
% - time lag of OUT vs E1 and E2
%
% Copyright (C) 2018,2012,2022,2023 Alois Schlögl, IST Austria


if exist('OCTAVE_VERSION','builtin')
        pkg load statistics
        pkg load nan
        pkg load tsa
        pkg load signal
        pkg load biosig

        close all
        graphics_toolkit gnuplot

elseif exist('/fs3/home/schloegl/src/octave-NaN/','dir')
        addpath /fs3/home/schloegl/src/octave-tsa/inst
        addpath /fs3/home/schloegl/src/octave-tsa/src
        addpath /fs3/home/schloegl/src/octave-NaN/inst
        addpath /fs3/home/schloegl/src/octave-NaN/src

elseif exist('/fs3/group/jonasgrp/','dir')
        addpath /fs3/group/jonasgrp/Software/Matlab/tsa/inst
        addpath /fs3/group/jonasgrp/Software/Matlab/tsa/src
        addpath /fs3/group/jonasgrp/Software/Matlab/NaN/inst
        addpath /fs3/group/jonasgrp/Software/Matlab/NaN/src
end

BIOSIG_MATLAB_PATH=getenv('BIOSIG_MATLAB_PATH')
if ispc()
    BIOSIG_MATLAB_PATH='C:\Users\sjamrich.ISTA\Desktop\biosig-code\biosig4matlab';
    mypp=pwd;
    cd(BIOSIG_MATLAB_PATH)
    install
    cd(mypp)
else
    pp=strsplit(BIOSIG_MATLAB_PATH,filesep);
    if strcmp(pp{end-2},'share'),
	    pp{end-2}='libexec';
	    pp{1}=filesep;
	    addpath(fullfile(pp{:}));
    end
end
if ~exist('hdr2ascii','file')
        addpath(fullfile(BIOSIG_MATLAB_PATH, 't200_FileAccess'));
end;
addpath(fullfile(BIOSIG_MATLAB_PATH, 'doc'));
if ~exist('mexSLOAD','file')
        addpath(fullfile(BIOSIG_MATLAB_PATH, '../biosig4c++/mex'));
end;
if ~exist('signal_deconvolution','file')
        addpath(fullfile(BIOSIG_MATLAB_PATH, 't300_FeatureExtraction'))
end;
if ~exist('rs','file')
        addpath(fullfile(BIOSIG_MATLAB_PATH, 't250_ArtifactPreProcessingQualityControl'))
end;
if ~exist('auc','file')
        addpath(fullfile(BIOSIG_MATLAB_PATH, 't490_EvaluationCriteria'));
end;
if ~exist('auc','file')
        addpath(fullfile(BIOSIG_MATLAB_PATH, 't490_EvaluationCriteria'));
end;
if ~exist('mod_optimal_detection_filter','file')
        addpath(fullfile(BIOSIG_MATLAB_PATH, 't340_MachineLearningOptimalDetection'));
end;


%% check prerequisites
if ~exist('sumskipnan','file') || ~exist('roc','file')
	fprintf(2,'Warning: sumskipnan() is missing - installing NaN-toolbox\n');
end
if isnan(mean([1,3,NaN]))
	fprintf(2,'Warning: mean() from NaN-toolbox is missing - install NaN-toolbox\n');
end
tmp=roc(randn(100,1),[1:100]'<50);
if ~isfield(tmp,'H_kappa')
	error('roc() from NaN-toolbox is missing - install NaN-toolbox');
end
if sum(isnan(detrend([1,NaN,3])))>1
	fprintf(2,'Warning: detrend() from NaN-toolbox is missing - install NaN-toolbox\n');
end
if ~exist('mexSLOAD','file')
        error('mexSLOAD is missing - mexBiosig or Biosig need to be installed');
end
if ~exist('detect_spikes_bursts','file')
        fprintf(2,'Error: detect_spikes_bursts() is missing - Biosig need to be installed\n');
end


DATADIR  = 'data';
DATAFILES2 = {
        '21-Feb-2018_001.dat'
        '23-Nov-2017_002.dat'
        'XM160126-02.dat'
        'XM160218-02.dat'
        'XM_2016-12-14_001.dat'
        'XM_2017-03-28.dat'
};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loading the data and do some basic quality checks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('outdir','var')
	outdir='output011'
end

load(fullfile(outdir,'GeneralClassifier.mat'), 'RES_C', 'HIGHPASS', 'LOWPASS', 'Fs', 'WINLEN', 'Mode', '-mat');
% HIGHPASS = 0;
% LOWPASS  = 1000;
% Fs       = 25000;	% target sampling rate [Hz]
MAXLAG = length(RES_C{1}.CLASSIFIER.A)-1;
% WINLEN   = 4; 	% window size [ms] for scoring trace
% Mode ; 		% defines preprocessing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Apply MOD to whole (including unscored) data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for E = 1 %:2,
	RES = RES_C{E};
	for k1 = 1:length(DATAFILES2);
		datafile = fullfile(DATADIR, DATAFILES2{k1});

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%         load data
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		chan     = 1; %find(strcmp('Adc0',HDR.Label));
		[data,HDR,rawdata] = sload4mod(datafile,chan,Fs,Mode,MAXLAG,[-0.150]);
		S = data(:,chan);
		[tmp_p,tmp_f,tmp_e]=fileparts(datafile);

		selpos = sort([1;HDR.EVENT.POS(HDR.EVENT.TYP==hex2dec('7ffe')); size(S,1)+1]);
		Stmp = S;
		%%% remove RS-pulse at 19.850-19.900, blank until end of sweep
		[ix,iy]=meshgrid(selpos-1,-round(-0.150*Fs):0);
		ix=ix(:)+iy(:);
		ix = ix( (0 < ix) & (ix <= size(S,1)) );
		Stmp(ix)=NaN;

		% the output parameters from the training steps are
		% RES.CLASSIFIER.A ; 			% filter coefficients of the Wiener filter
		TH    = RES.CLASSIFIER.THRESHOLD;	% detection threshold - based on maxKappa
		delay = RES.CLASSIFIER.delay;		% time shift of the filter

		% the classifier is applied to the whole data set.
		D     = filter(RES.CLASSIFIER.A,1,Stmp);

		% without smoothing
		WinHann = 1;
		D    = [D(delay+1:end);repmat(NaN,delay,1)];
		pos  = get_local_maxima_above_threshold(D,TH,5,0.001*min(WINLEN/2,1)*Fs);
		% [tmp, pos] = findpeaks(D, 'MinPeakHeight', TH, 'MinPeakDistance', 0.001*max(WINLEN,1)*Fs);

		% save results in a GDF file for visualization with SigViewer
		H0 = [];
		[tmp_p,tmp_f,tmp_e] = fileparts(HDR.FileName);
		H0.FileName = sprintf('%s.%s.w%d.e%d.h%d.%s.gdf',fullfile(outdir,tmp_f), 'mod',round(WINLEN*1000),E,WinHann,Mode);
		H2 = minidet2gdf([S, D, D>TH], [], pos, WINLEN, Fs, H0, HDR);
	end
end


