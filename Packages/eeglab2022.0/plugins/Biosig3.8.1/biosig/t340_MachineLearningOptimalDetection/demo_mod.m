% demo_mod shows the essential steps of using MOD for  
% detecting EPSP/EPSCs as described in [1] and applied in [2]. 
% 
% In order to run this demo, the following software packages 
% are needed: 
% - Octave (or Matlab) and the Signal Processing toolbox or package
% - Biosig [3]
% - NaN-toolbox [4] 
% - SigViewer [5] - for scoring the data and visualization of the 
%    resulting detections
%
% Copyright (C) 2016-2021 Alois Schlögl, IST Austria
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
%     https://github.com/cbrnr/sigviewer
%     https://git.ist.ac.at/alois.schloegl/sigviewer.git 
%     https://pub.ist.ac.at/~schloegl/software/sigviewer/sigviewer-0.6.4-win64.exe


if exist('OCTAVE_VERSION','builtin')
	try
		pkg load biosig
	end
	try
		pkg load mexbiosig
	end
 	pkg load nan
	pkg load tsa
	pkg load signal 
end
%% check prerequisites
if ~exist('sumskipnan','file') || ~exist('roc','file')
	addpath(fullfile(pwd,'NaN','src'));
        addpath(fullfile(pwd,'NaN','inst'));
	fprintf(2,'Warning: sumskipnan() is missing - installing NaN-toolbox\n');
end
if isnan(mean([1,3,NaN]))
	addpath(fullfile(pwd,'NaN','src'));
        addpath(fullfile(pwd,'NaN','inst'));
	fprintf(2,'Warning: mean() from NaN-toolbox is missing - install NaN-toolbox\n');
end
tmp=roc(randn(100,1),[1:100]'<50);
if ~isfield(tmp,'H_kappa')
	error('roc() from NaN-toolbox is missing - install NaN-toolbox');
end
if sum(isnan(detrend([1,NaN,3])))>1
	addpath(fullfile(pwd,'NaN','src'));
        addpath(fullfile(pwd,'NaN','inst'));
	fprintf(2,'Warning: detrend() from NaN-toolbox is missing - install NaN-toolbox\n');
end
if ~exist('mexSLOAD','file')
        error('mexSLOAD is missing - mexBiosig or Biosig need to be installed');
end
if ~exist('detect_spikes_bursts','file')
        fprintf(2,'Error: detect_spikes_bursts() is missing - Biosig need to be installed\n');
end


% read channel of raw data, the example data is available from here:
%    https://pub.ist.ac.at/~schloegl/software/mod/
filename='data/XM160126-02.dat'; chan=1;
if ~exist(filename,'file'),
	error(sprintf('file %s not available - you can download the example file from https://pub.ist.ac.at/~schloegl/software/mod/data/'));
end
[s,HDR]=mexSLOAD(filename,chan);
HDR.SampleRate=round(HDR.SampleRate);
HDR.EVENT.SampleRate=round(HDR.EVENT.SampleRate);
% clean up data, everything before this marker can be ignored. 
if isfield(HDR.EVENT,'CodeDesc')
	typ   = find(strcmp(HDR.EVENT.CodeDesc,'IV-1000'));
	if length(typ)>0,
		START = max(HDR.EVENT.POS(HDR.EVENT.TYP==typ));
		s(1:START) = NaN;
	end
end

% load the two segments of the scoring from Export E1 
[e11,EVT]=mexSLOAD('data/XM160126-02_E1S1.evt');
event_pos1 = EVT.EVENT.POS(EVT.EVENT.TYP==10); 
[e12,EVT]=mexSLOAD('data/XM160126-02_E1S2.evt');
event_pos2 = EVT.EVENT.POS(EVT.EVENT.TYP==10); 


% Here are some standard parameters we have used in [1,2].
if HDR.SampleRate < 50000,
	Fs=HDR.SampleRate;
else
	Fs=25000;
end
TemplateLength=0.040*Fs; % filterlength
Twin=0.004;	% window length in seconds [s]
WINLEN=round(Twin*Fs);


% load the two segment of the scoring from Export E2 
%   we ignore this here for the sake of simplicity 
% [e2,HDR] =mexSLOAD('data/XM160126-02_E2S12.evt')


% identify action potentials 
AP = [];
if exist('detect_spikes_bursts','file')
	AP = detect_spikes_bursts(filename, chan);
end

% downsample to 25 kHz 
if (HDR.SampleRate ~= Fs)
	DIV = round(HDR.SampleRate/Fs);
	s = rs(s,DIV,1);

	% adapt event positions to new sampling rate 
	HDR.EVENT.POS = round((HDR.EVENT.POS-1)*Fs/HDR.EVENT.SampleRate)+1;
	if isfield(HDR.EVENT,'DUR');
        	HDR.EVENT.DUR = round(HDR.EVENT.DUR*Fs/HDR.EVENT.SampleRate);
	end
	HDR.SampleRate = Fs; 
	HDR.EVENT.SampleRate = Fs;

	event_pos1 = round(event_pos1/DIV);
	event_pos2 = round(event_pos2/DIV);

	if isfield(AP,'EVENT')
		AP.EVENT.POS = round((AP.EVENT.POS-1)*Fs/AP.EVENT.SampleRate)+1;
		if isfield(HDR.EVENT,'DUR');
			AP.EVENT.DUR = round(AP.EVENT.DUR*Fs/AP.EVENT.SampleRate);
		end
		AP.EVENT.SampleRate = Fs;
	end
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       remove AP's and spikelets and replace with blanks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(AP)
	apix = find(AP.EVENT.TYP==hex2dec('0201'));
	[ix,iy]=meshgrid(apix,[0:TemplateLength]-round(TemplateLength/4));
	ix=ix(:)+iy(:); ix((ix<1) | (ix>length(s)))=[];
	s(ix)=NaN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       generate scoring trace 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TiS1=[min(event_pos1);max(event_pos1)];
TiS2=[min(event_pos2);max(event_pos2)];
C = repmat(NaN,size(s));  % set to NaN for unscored periods
C(TiS1(1):TiS1(2))=0;	% scored periods are now initialized
C(TiS2(1):TiS2(2))=0;
[ix,iy]=meshgrid([event_pos1(:);event_pos2(:)],round(-WINLEN/2):round(WINLEN/2));
ix=ix(:)+iy(:); ix(ix<1)=[]; ix(ix>length(s))=[];
C(ix) = 1;		% set the events with window size WINLEN in the scoring trace

% remove linear trend 
S = detrend(s);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       cross-validation (XV)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1, 	% disable or enable XV
	% Grouping for S1-S2 cross-validation (XV)
	G1 = zeros(size(s)); 
	G1(TiS1(1):TiS1(2))=1;
	G1(TiS2(1):TiS2(2))=2;

	% Grouping for A1B2-A2B1 cross-validation(XV)
	G2 = zeros(size(s)); 
	ti3=[sum(TiS1),sum(TiS2)]/2;
	G2(TiS1(1):ti3(1))=1;
	G2(ti3(1):TiS1(2))=2;
	G2(TiS2(1):ti3(2))=2;
	G2(ti3(2):TiS2(2))=1;

	% XV: S1->S2 cross-validiation
	ixtrain=find(G1==1);
	ixtest=find(G1==2);
	RES12S=mod_optimal_detection_filter(S, C, TemplateLength, ixtrain, ixtest);
	% XV: S2->S1 cross-validiation
	ixtrain=find(G1==2);
	ixtest=find(G1==1);
	RES21S=mod_optimal_detection_filter(S, C, TemplateLength, ixtrain, ixtest);

	% In [1,2], we combined the output from the two testing sets, and 
	% and applied the ROC analysis to this trace. 
	% for the sake of simplicity, this is omitted here. 
	fprintf(1,'AUC: (XV:S1->S2):\t %g\n', RES12S.test.AUC);
	fprintf(1,'AUC: (XV:S2->S1):\t %g\n', RES21S.test.AUC);
	fprintf(1,'Kappa (XV:S1-S2)    : %.3g/%.3g\n',kappa(RES12S.test.CM).kappa,kappa(RES21S.test.CM).kappa);

	% XV: A1B2->A2B1 cross-validiation
	ixtrain=find(G2==1);
	ixtest=find(G2==2);
	RES1A2B=mod_optimal_detection_filter(S, C, TemplateLength, ixtrain, ixtest);
	% XV: A2B1->A1B2 cross-validiation
	ixtrain=find(G2==2);
	ixtest=find(G2==1);
	RES2A1B=mod_optimal_detection_filter(S, C, TemplateLength, ixtrain, ixtest);

	% In [1,2], we combined the output from the two testing sets, and 
	% and applied the ROC analysis to this trace. 
	% for the sake of simplicity, this is omitted here. 
	fprintf(1,'AUC: (XV:A1B2->A2B1):\t %g\n', RES1A2B.test.AUC);
	fprintf(1,'AUC: (XV:A2B1->A1B2):\t %g\n', RES2A1B.test.AUC);
	fprintf(1,'Kappa (XV:A1B2-A2B1): %.3g/%.3g\n',kappa(RES1A2B.test.CM).kappa,kappa(RES2A1B.test.CM).kappa);

	% XV: LOOM
	% this is not implemented here, essentially the idea is to 
	% concatenate the data from multiple cells, and use the data
	% of m-1 cells for training, and the data from the remaining cell for
	% testing.  
end

% w/o XV: scored data is used for training
ixtrain=find(~isnan(C)); % use the scored data for training 	
% avoid edge issues when computing correlation functions
ixtrain=ixtrain((2*TemplateLength<ixtrain) & (ixtrain < (length(S)-2*TemplateLength)));
ixtest=[];
RES=mod_optimal_detection_filter(S, C, TemplateLength, ixtrain, ixtest);
% the output parameters from the training steps are
% RES.CLASSIFIER.A ; 			% filter coefficients of the Wiener filter 
TH    = RES.CLASSIFIER.THRESHOLD;	% detection threshold - based on maxKappa
delay = RES.CLASSIFIER.delay;		% time shift of the filter

% the classifier is applied to the whole data set. 
D     = filter(RES.CLASSIFIER.A,1,S);

% without smoothing 
WinHann = 1; 
D    = [D(delay+1:end);repmat(NaN,delay,1)];
pos  = get_local_maxima_above_threshold(D,TH,4,0.001*Fs);    


% save results in a GDF file for visualization with SigViewer 
mkdir('output');
H0 = [];
[tmp_p,tmp_f,tmp_e] = fileparts(HDR.FileName);
H0.FileName = sprintf('output/%s.%s.w%d.e1.h%d.gdf',tmp_f, 'mod',round(Twin*1000),WinHann);
H2 = minidet2gdf([s, D, D>TH], [], pos, WINLEN, Fs, H0, AP);


% apply a smoothing filter as discribed in [1] - this might not be needed 
D     = filter(RES.CLASSIFIER.A,1,S);
WinHann = 13; 
A    = hann(WinHann); A=A/sum(A);
nix  = isnan(D); D(nix) = mean(D);
D    = filtfilt(A,1,D);
D(nix)=NaN;
D    = [D(delay+1:end);repmat(NaN,delay,1)];
pos  = get_local_maxima_above_threshold(D,TH,4,0.001*Fs);

% save results in a GDF file for visualization with SigViewer 
H0 = [];
[tmp_p,tmp_f,tmp_e] = fileparts(HDR.FileName);
H0.FileName = sprintf('output/%s.%s.w%d.e1.h%d.gdf', tmp_f, 'mod', round(Twin*1000), WinHann);
H2 = minidet2gdf([s, D, D>TH], [], pos, WINLEN, Fs, H0, AP);

