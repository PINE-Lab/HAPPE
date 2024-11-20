% DEMO_MODY shows the essential steps of using MOD for
% detecting EPSP/EPSCs as described in [1] and applied in [2].
%
% DEMO_MODY replaces DEMO_MODX, and uses SLOAD4MOD to load the data, and
% apply preprocessing methods. It is an extended version of the simpler
% version of DEMO_MOD and adds the LOOM-crossvalidation procedure and
% handling the scoring from two experts.
% DEMO_MODY is also useful to investigate the influence of different
% preprocessing methods.
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
% - delay time correction for pj025 (corrected on Apr 30th)
%
% Copyright (C) 2018,2012,2022 Alois Schlögl, IST Austria


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

BIOSIG_MATLAB_PATH=getenv('BIOSIG_MATLAB_PATH');
if ispc()
    BIOSIG_MATLAB_PATH='C:\Users\sjamrich.ISTA\Desktop\biosig-code\biosig4matlab';
    mypp=pwd;
    cd(BIOSIG_MATLAB_PATH);
    install
    cd(mypp)
elseif exist('BIOSIG_MATLAB_PATH','var')
    pp=strsplit(BIOSIG_MATLAB_PATH,filesep);
    if strcmp(pp{end-2},'share'),
	    pp{end-2}='libexec';
	    pp{1}=filesep;
	    addpath(fullfile(pp{:}));
    end
end

if ~exist('OCTAVE_VERSION','builtin')
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
		addpath(fullfile('C:\Users\sjamrich.ISTA\Desktop\','NaN','src'));
		addpath(fullfile('C:\Users\sjamrich.ISTA\Desktop\','NaN','inst'));
		fprintf(2,'Warning: sumskipnan() is missing - installing NaN-toolbox\n');
	end
	if isnan(mean([1,3,NaN]))
		addpath(fullfile('C:\Users\sjamrich.ISTA\Desktop\','NaN','src'));
		addpath(fullfile('C:\Users\sjamrich.ISTA\Desktop\','NaN','inst'));
		fprintf(2,'Warning: mean() from NaN-toolbox is missing - install NaN-toolbox\n');
	end
	tmp=roc(randn(100,1),[1:100]'<50);
	if ~isfield(tmp,'H_kappa')
		error('roc() from NaN-toolbox is missing - install NaN-toolbox');
	end
	if sum(isnan(detrend([1,NaN,3])))>1
		addpath(fullfile('C:\Users\sjamrich.ISTA\Desktop\','NaN','src')); hex2dec('0211')

		addpath(fullfile('C:\Users\sjamrich.ISTA\Desktop\','NaN','inst'));
		fprintf(2,'Warning: detrend() from NaN-toolbox is missing - install NaN-toolbox\n');
	end
	if ~exist('mexSLOAD','file')
		error('mexSLOAD is missing - mexBiosig or Biosig need to be installed');
	end
	if ~exist('detect_spikes_bursts','file')
		fprintf(2,'Error: detect_spikes_bursts() is missing - Biosig need to be installed\n');
	end
end

if ~exist('WINLEN','var')
	WINLEN=4; 	% window size [ms] for scoring trace
end;
MAXLAG=1000;	% filter length [samples]
Fs = 25000;	% target sampling rate [Hz]




%%% setup loops for manual patternsearch
ftype='gauss2';
for HP = 1; [1,3,.3]; [10, .1, .001, 1, .01]; %[30, 3, .3, 0.003, 0.01, 0.03]; [0, 0.0001,0.00003,0.00001,.001,.003,.01,.03,.1,.3,1]; % [1,.3,.1,.03,.01]; [1,2,3,5,10,30]
for LP = 5000; [7000,10000]; [-1,1000,2000,500,5000]; [1:2:9,1,2:2:10]*1000;
for WINLEN = 3; [2,3,4]; [0.60:0.08:1.0]; [0.80, .6,.4,.2,1,2,4,.64,.68,.72,.76,.80,.84,.88,.92,.96,1.0]; %.6,.4,.2,1,2,4]; [4,3,2,1];
for MAXLAG = [400]; 800:200:1200; 400:200:1000; [100:100:1000];

datapath  = 'data';
DATAFILES = {
        '21-Feb-2018_001.dat'
        '23-Nov-2017_002.dat'
        'XM160126-02.dat'
        'XM160218-02.dat'
        'XM_2016-12-14_001.dat'
        'XM_2017-03-28.dat'
};

% the results from different preprocessing methiods
Mode=sprintf('%s_%.6f_%d',ftype,HP,LP);
outdir=sprintf('output051_%sHP%.6f_LP%d_WINLEN%.2f_MAXLAG%d',ftype,HP,LP,WINLEN,MAXLAG)
mkdir(outdir);
if exist(fullfile(outdir,'GeneralClassifier.mat'),'file'), continue; end

fid0=fopen(fullfile(outdir,'summary.txt'),'w');
fid1=fopen(fullfile(outdir,'summary001.txt'),'w');
fid2=fopen(fullfile(outdir,'summary003.txt'),'a');
fprintf(fid2,'dataset#\tDuration(S1,E1)\t#PSP(S1,E1)\tEventRate(S1,E1)\tmean(log(IEI11))\tstd(log(IEI11))\t Duration(S2,E1)\t#PSP(S2,E1)\tEventRate(S2,E1)\tmean(log(IEI21))\tstd(log(IEI21))\t Duration(S1,E2)\t#PSP(S1,E2)\tEventRate(S1,E2)\tmean(log(IEI12))\tstd(log(IEI12))\t Duration(S2,E2)\t#PSP(S2,E2)\tEventRate(S2,E2)\tmean(log(IEI22))\tstd(log(IEI22))\n');

NAN_BLOCK=repmat(NaN,3*Fs,1);
X.C1 = NAN_BLOCK;
X.C2 = NAN_BLOCK;
X.G0 = NAN_BLOCK;        % cell number
X.G1 = NAN_BLOCK;        % 1: david's data 2: xiaomin's data, 3: CA1 - not used in this demo
X.G2 = NAN_BLOCK;        % S1 - S2
X.G3 = NAN_BLOCK;        % A1B2 - A2B1
X.S  = NAN_BLOCK;        % data
X.isnan = NAN_BLOCK;     % data

RES.lag=repmat(NaN,length(DATAFILES),1);
RES.CM=repmat(NaN,[length(DATAFILES),2,2]);
RES.kappa=repmat(NaN,length(DATAFILES),1);
RES.kappa1=repmat(NaN,length(DATAFILES),1);
RES.kappa_sd=repmat(NaN,length(DATAFILES),1);
RES.kappa1_sd=repmat(NaN,length(DATAFILES),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loading the data and do some basic quality checks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k1=1:length(DATAFILES);
        CHAN = 1;
        datafile = fullfile(datapath, DATAFILES{k1});
        G1(k1) = ceil(k1/length(DATAFILES));
        G2(k1) = mod(k1-1,length(DATAFILES))+1;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         load scoring
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        TupLimit   = Inf;
        [p,f,e]    = fileparts(datafile);
        evtfile2   = fullfile(p,[f,'_s2.evt']);
        evtfile3   = fullfile(p,[f,'_s3.evt']);
        evtfile01  = fullfile(p,[f,'_david.evt']);
        evtfile01b = '';

        if isempty(dir(evtfile01));
                evtfile01 = fullfile(p,[f,'.evt']);
        end
        if isempty(dir(evtfile01));
                evtfile01 = fullfile(p,[f(1:end-4),'.evt']);
        end
        if isempty(dir(evtfile01));
                evtfile01 = fullfile(p,[f(1:end-4),'.evt']);
        end
        if isempty(dir(evtfile01));
                evtfile01 = fullfile(p,[f,'_VD.evt']);
        end
        if isempty(dir(evtfile01));
                evtfile01 = fullfile(p,[f,'_VD_S2.evt']);
                evtfile01b = fullfile(p,[f,'_VD_S3.evt']);
        end
        if isempty(dir(evtfile01b));
	        evtfile01b = fullfile(p,[f,'_missingpart_david.evt']);
        end

        if isempty(dir(evtfile2));
                evtfile2 = fullfile(p,[f(1:end-4),'_s2.evt']);
        end
        if isempty(dir(evtfile3));
                evtfile3 = fullfile(p,[f(1:end-4),'_s3.evt']);
        end

        if isempty(dir(evtfile01));
                GROUP(k1)=0;
                fprintf(2,'scoring file1 %s missing\n',evtfile01);
                continue
        end;
        if isempty(dir(evtfile2));
                GROUP(k1)=0;
                fprintf(2,'scoring file2 %s missing\n',evtfile2);
                continue
        end;
        if isempty(dir(evtfile3));
                GROUP(k1)=0;
                fprintf(2,'scoring file3 %s missing\n',evtfile3);
                continue
        end;

        [t,EVT01] = mexSLOAD(evtfile01);
        [t,EVT2]  = mexSLOAD(evtfile2);
        [t,EVT3]  = mexSLOAD(evtfile3);

        tix01 = sort(EVT01.EVENT.POS((EVT01.EVENT.TYP==10) | (EVT01.EVENT.TYP==11) ) / EVT01.EVENT.SampleRate);
        if exist(evtfile01b,'file');
                [t,EVT01b] = mexSLOAD(evtfile01b);
                tix01b = sort(EVT01b.EVENT.POS((EVT01b.EVENT.TYP==10) | (EVT01b.EVENT.TYP==11) )/EVT01b.EVENT.SampleRate);
		if strfind(evtfile01b,'missingpart_david')
	                tix01 = [tix01(tix01<200); tix01b];
		else
	                tix01 = [tix01; tix01b];
		end;
        end
        tix2  = sort(EVT2.EVENT.POS(EVT2.EVENT.TYP==10)/EVT2.EVENT.SampleRate);
        tix3  = sort(EVT3.EVENT.POS(EVT3.EVENT.TYP==10)/EVT3.EVENT.SampleRate);
        if (strcmp(f,'XM_2016-12-14_001') || strcmp(f,'XM_2017-03-28'))
                fprintf(2,'Warning: remove %d events from file %s\n',sum(tix2<563),evtfile2);
                tix2=tix2(tix2>563);
        end
        tix02 = [tix2; tix3];
        if isempty(tix01)
                fprintf(2,'scorings from 1st expert are missing\n');
        end
        if isempty(tix02)
                fprintf(2,'scorings from 2nd expert are missing\n');
        end

        tt1   = (min(tix01)+max(tix01))/2;
        tt2   = (min(tix02)+max(tix02))/2;
        tix11 = sort(tix01(tix01 < tt1));
        tix21 = sort(tix01(tix01 > tt1));
        tix12 = sort(tix02(tix02 < tt2));

        tix22 = sort(tix02(tix02 > tt2));


        % fprintf(fid0,'== %s\n\t%f\t%f\t%f\t%f\t%f\n\t%f\t%f\t%f\t%f\t%f\n',datafile,[min(tix11),max(tix11),tt1,min(tix21),max(tix21),min(tix12),max(tix12),tt2,min(tix22),max(tix22)]);
        % fprintf(fid1,'== %s\n\t%f\t%f\t%f\t%f\t%f\n\t%f\t%f\t%f\t%f\t%f\n',[f,e],[min(tix11),max(tix11),tt1,min(tix21),max(tix21),min(tix12),max(tix12),tt2,min(tix22),max(tix22)]);

        % fprintf(fid0,'E1:\t%s\n\t%s\nE2:\t%s\n\t%s\n',evtfile01,evtfile01b,evtfile2,evtfile3);
        % fprintf(fid0,'Duration [s]\t%f\t%f\n\t\t%f\t%f\n',tix11(end)-tix11(1),tix21(end)-tix21(1),tix12(end)-tix12(1),tix22(end)-tix22(1));
        % fprintf(fid0,'NoOfPSPs\t%d\t%d\n\t\t%d\t%d\n',length(tix11),length(tix21),length(tix12),length(tix22));
        % fprintf(fid0,'EventRate [Hz]\t%f\t%f\n\t\t%f\t%f\n',length(tix11)/(tix11(end)-tix11(1)),length(tix21)/(tix21(end)-tix21(1)),length(tix12)/(tix12(end)-tix12(1)),length(tix22)/(tix22(end)-tix22(1)));

        [tmp,f2,e2] = fileparts(evtfile2);
        % fprintf(fid1,'E1:\t%s\n\t%s\nE2:\t%s\n\t%s\n',[f1,e1],[f1b,e1b],[f2,e2],[f3,e3]);
        %fprintf(fid1,'Duration [s]\t%f\t%f\n\t\t%f\t%f\n',tix11(end)-tix11(1),tix21(end)-tix21(1),tix12(end)-tix12(1),tix22(end)-tix22(1));
        % fprintf(fid1,'NoOfPSPs\t%d\t%d\n\t\t%d\t%d\n',length(tix11),length(tix21),length(tix12),length(tix22));
        % fprintf(fid1,'EventRate [Hz]\t%f\t%f\n\t\t%f\t%f\n',length(tix11)/(tix11(end)-tix11(1)),length(tix21)/(tix21(end)-tix21(1)),length(tix12)/(tix12(end)-tix12(1)),length(tix22)/(tix22(end)-tix22(1)));

	iei11=diff(tix11);iei11=iei11(isfinite(iei11));
	iei21=diff(tix21);iei21=iei21(isfinite(iei21));
	iei12=diff(tix12);iei12=iei12(isfinite(iei12));
	iei22=diff(tix22);iei22=iei22(isfinite(iei22));

        % fprintf(fid2,'%i\t%.1f\t%d\t%.2f\t%.3f\t%.3f\t %.1f\t%d\t%.2f\t%.3f\t%.3f\t %.1f\t%d\t%.2f\t%.3f\t%.3f\t %.1f\t%d\t%.2f\t%.3f\t%.3f\n',k1, [ [tix11(end)-tix11(1), tix21(end)-tix21(1), tix12(end)-tix12(1), tix22(end)-tix22(1)]; [length(tix11), length(tix21), length(tix12), length(tix22)]; [length(tix11)/(tix11(end)-tix11(1)), length(tix21)/(tix21(end)-tix21(1)), length(tix12)/(tix12(end)-tix12(1)), length(tix22)/(tix22(end)-tix22(1))]; exp([mean(log(iei11)), mean(log(iei21)), mean(log(iei12)), mean(log(iei22)) ]); [std(log(iei11)), std(log(iei21)), std(log(iei12)), std(log(iei22)) ] ] );

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         load data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [HDR]    = mexSOPEN(datafile);
        chan     = find(strcmp('Vmon-1',HDR.Label));
        if isempty(chan) chan=1; end
        TemplateLength  = MAXLAG;
        BlankingInterval = [];
	[S,HDR,data] = sload4mod(datafile,chan,Fs,Mode,TemplateLength,BlankingInterval);

        tix11 = round(tix11*Fs);
        tix21 = round(tix21*Fs);
        tix12 = round(tix12*Fs);
        tix22 = round(tix22*Fs);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         run EPSP detection, blank AP's
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %assert(isequal(HDR.PhysDim(chan),{'V'}))

        % use scorings from both expoerts and use only common intervals for the analysis
        % t2 indicates the first scoring segment, t3 indicates the 2nd segment
        dt = round(median(diff([tix12;tix22])));
        tt = max([tix12(1)-dt,1]) : min([tix22(end)+dt,length(S)]);
        t2 = max([tix12(1)-dt,1]) : min([tix12(end)+dt,length(S)]);
        t3 = max([tix22(1)-dt,1]) : min([tix22(end)+dt,length(S)]);


if isempty(t3) || isempty(t2)
        fprintf(fid0,'# ERROR: t2 or t3 is empty %d/%d\n',length(t2),length(t3));
        continue
end

        s = S;
        % S = data(:,chan);
        s(1:tt(1)-1)=NaN; s(tt(end)+1:end)=NaN;
        % S(1:tt(1)-1)=NaN; S(tt(end)+1:end)=NaN;

        tix  = [t2(1):t2(end),t3(1):t3(end)]';
        tix1 = [t2(1):t2(end)]';
        tix2 = [t3(1):t3(end)]';
        % s(tix(isnan(s(tix)))) = mean(s(tix));

         %%% Classes - scoring trace
        CC.WINLEN = WINLEN/1000;
        win = round(CC.WINLEN*Fs);
        c1  = repmat(NaN,size(s));
        c2  = repmat(NaN,size(s));
        c1(tix) = 0;
        c2(tix) = 0;
        for ix = [tix11;tix21]',
                c1( max(1,ix-floor(win/2)):min(ix+floor(win/2),length(c1)) ) = 1;
        end;
        for ix = [tix12;tix22]',
                c2( max(1,ix-floor(win/2)):min(ix+floor(win/2),length(c2)) ) = 1;
        end;

        T1=-.050*Fs;
        T2= .050*Fs;
        if 0
                ; %noop

        elseif 1, %  store clean data to GDF so that we can avoid doing the
		  %  cleanup again when applying the classifier to the whole data set
                H1.FileName=fullfile('data',[HDR.FILE.Name,'.gdf']);
                H1.TYPE='GDF';
                H1.VERSION=3;
                H1.NS=1;
                H1.SampleRate=Fs;
                H1.NRec=HDR.NRec*HDR.SPR;
                H1.SPR=1;
                H1.PhysMax = HDR.PhysMax(chan);
                H1.PhysMin = HDR.PhysMin(chan);
                H1.DigMax  = HDR.DigMax(chan);
                H1.DigMin  = HDR.DigMin(chan);
                H1.PhysDimCode = HDR.PhysDimCode(chan);
                H1.Label = HDR.Label(chan);
                H1.Cal  = HDR.Cal(chan);
                H1.Off  = HDR.Off(chan);
                H1.GDFTYP  = 16; HDR.GDFTYP(chan);
                H1.FLAG.UCAL=0;
                if isfield(H1,'AS'); H1=rmfield(H1,'AS'); end
                H1=sopen(H1,'w'); H1=swrite(H1,S); H1=sclose(H1);
        end

        % identify time offset "toff" between experts,
        %   and compute interscorer agreement using Cohen's kappa
        %    before and after correction for time shift
        [xcf,lag]     = xcorr(c1(tix),c2(tix),0.02*Fs,'coeff');
        [tmp,tmpix]   = max(xcf);
        toff1         = tmpix - find(lag==0);                % that's the time offset between scorers
        toff          = lag(tmpix);                % that's the time offset between scorers

        RES.lag(k1,1) = toff;

        if 0
            % inter-scorer agreement
            HIS=histo3(c1(tix)+2*c2(tix));
            H=reshape(HIS.H,2,2);
            [RES.kappa(k1),RES.kappa_sd(k1),h,z(k1,1)] = kappa(H);
            RES.CM(k1,:,:) = h;

            % inter-scorer agreement with correction for time offset
            HIS=histo3(c1(tix+toff)+2*c2(tix));
            H=reshape(HIS.H,2,2);
            [RES.kappa1(k1),RES.kappa1_sd(k1),h,z1(k1,1)] = kappa(H);
            RES.CM1(k1,:,:) = h;

            subplot(5,4,k1);
            plot(lag/Fs,xcf,'-',lag(tmpix)/Fs,xcf(tmpix),'ro');grid on
            title(sprintf('lag=%.2f ms, kappa=%.3f (%.3f)', 1000*toff/Fs, RES.kappa(k1), RES.kappa1(k1) ))
            set(gca,'box','off')

            fprintf(fid0,'# Add delay of %.2f ms to scoring of E1 (%.2f, %d)\n', toff1*1000/Fs, toff*1000/Fs, toff1==toff);
        end
        % generate data with labels for further analysis
        X.C1    = [X.C1;    c1(tix1+toff);                 NAN_BLOCK; c1(tix2+toff);                 NAN_BLOCK]; % scoring Expert 1, compensate for offset
        X.C2    = [X.C2;    c2(tix1);                      NAN_BLOCK; c2(tix2);                      NAN_BLOCK]; % scoring Expert 2
        X.G0    = [X.G0;    repmat(k1,length(tix1),1);     NAN_BLOCK; repmat(k1,length(tix2),1);     NAN_BLOCK]; % cell number
        X.G1    = [X.G1;    repmat(G1(k1),length(tix1),1); NAN_BLOCK; repmat(G1(k1),length(tix2),1); NAN_BLOCK]; % 1: david's data 2: dentate 3: ca1
        X.G2    = [X.G2;    repmat(1,length(tix1),1);      NAN_BLOCK; repmat(2,length(tix2),1);      NAN_BLOCK]; % S1 - S2

        A1B2    = [1:length(tix1)]'*2 < length(tix1);
        A2B1    = [1:length(tix2)]'*2 > length(tix2);
        X.G3    = [X.G3;    A1B2;                          NAN_BLOCK; A2B1;                          NAN_BLOCK]; % A1B2 - A2B1

        X.S     = [X.S;     s(tix1);                       NAN_BLOCK; s(tix2);                       NAN_BLOCK];
        X.isnan = [X.isnan; isnan(S(tix1));             NAN_BLOCK; isnan(S(tix2));             	NAN_BLOCK];
end

ResTable1=[];
ResTable2=[];
for k1=[]; 1:length(DATAFILES);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%       within-cell cross-validation (XV)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% XV: S1->S2 cross-validiation
	ixtrain=find(X.G2==1 & X.G0==k1);
	ixtest=find(X.G2==2 & X.G0==k1);
	RES12S_C1=mod_optimal_detection_filter(X.S, X.C1, TemplateLength, ixtrain, ixtest);
	RES12S_C2=mod_optimal_detection_filter(X.S, X.C2, TemplateLength, ixtrain, ixtest);
	% XV: S2->S1 cross-validiation
	ixtrain=find(X.G2==2 & X.G0==k1);
	ixtest=find(X.G2==1 & X.G0==k1);
	RES21S_C1=mod_optimal_detection_filter(X.S, X.C1, TemplateLength, ixtrain, ixtest);
	RES21S_C2=mod_optimal_detection_filter(X.S, X.C2, TemplateLength, ixtrain, ixtest);

	% In [1,2], we combined the output from the two testing sets, and
	% and applied the ROC analysis to this trace.
	% for the sake of simplicity, this is omitted here.

	% XV: A1B2->A2B1 cross-validiation
	ixtrain=find(X.G3==1 & X.G0==k1);
	ixtest=find(X.G3==0 & X.G0==k1);
	RES1A2B_C1=mod_optimal_detection_filter(X.S, X.C1, TemplateLength, ixtrain, ixtest);
	RES1A2B_C2=mod_optimal_detection_filter(X.S, X.C2, TemplateLength, ixtrain, ixtest);
	% XV: A2B1->A1B2 cross-validiation
	ixtrain=find(X.G3==0 & X.G0==k1);
	ixtest=find(X.G3==1 & X.G0==k1);
	RES2A1B_C1=mod_optimal_detection_filter(X.S, X.C1, TemplateLength, ixtrain, ixtest);
	RES2A1B_C2=mod_optimal_detection_filter(X.S, X.C2, TemplateLength, ixtrain, ixtest);

	% In [1,2], we combined the output from the two testing sets, and
	% and applied the ROC analysis to this trace.
	% for the sake of simplicity, this is omitted here.

	fprintf(fid0,'\n#####\nResults from export scoring 1 of %s \n',DATAFILES{k1});
	fprintf(fid0,'AUC: (XV:S1->S2):\t %g\n', RES12S_C1.test.AUC);
	fprintf(fid0,'AUC: (XV:S2->S1):\t %g\n', RES21S_C1.test.AUC);
	fprintf(fid0,'Kappa (XV:S1-S2)    : %.3g/%.3g\n',kappa(RES12S_C1.test.CM).kappa,kappa(RES21S_C1.test.CM).kappa);
	fprintf(fid0,'AUC: (XV:A1B2->A2B1):\t %g\n', RES1A2B_C1.test.AUC);
	fprintf(fid0,'AUC: (XV:A2B1->A1B2):\t %g\n', RES2A1B_C1.test.AUC);
	fprintf(fid0,'Kappa (XV:A1B2-A2B1): %.3g/%.3g\n',kappa(RES1A2B_C1.test.CM).kappa,kappa(RES2A1B_C1.test.CM).kappa);

	fprintf(fid0,'\n#####\nResults from export scoring 2 of %s\n',DATAFILES{k1});
	fprintf(fid0,'AUC: (XV:S1->S2):\t %g\n', RES12S_C2.test.AUC);
	fprintf(fid0,'AUC: (XV:S2->S1):\t %g\n', RES21S_C2.test.AUC);
	fprintf(fid0,'Kappa (XV:S1-S2)    : %.3g/%.3g\n',kappa(RES12S_C2.test.CM).kappa,kappa(RES21S_C2.test.CM).kappa);
	fprintf(fid0,'AUC: (XV:A1B2->A2B1):\t %g\n', RES1A2B_C2.test.AUC);
	fprintf(fid0,'AUC: (XV:A2B1->A1B2):\t %g\n', RES2A1B_C2.test.AUC);
	fprintf(fid0,'Kappa (XV:A1B2-A2B1): %.3g/%.3g\n',kappa(RES1A2B_C2.test.CM).kappa,kappa(RES2A1B_C2.test.CM).kappa);

    ResTable1(k1,1:4)=[RES12S_C1.test.AUC, RES21S_C1.test.AUC, RES1A2B_C1.test.AUC, RES2A1B_C1.test.AUC];
    ResTable2(k1,1:4)=[RES12S_C2.test.AUC, RES21S_C2.test.AUC, RES1A2B_C2.test.AUC, RES2A1B_C2.test.AUC];
end

for k=1:max(X.G0), % length(DATAFILES),
	% XV: LOOM
	% this is not implemented here, essentially the idea is to
	% concatenate the data from multiple cells, and use the data
	% of m-1 cells for training, and the data from the remaining cell for
	% testing.

	ixtrain=find(X.G0~=k & ~isnan(X.G0));
	ixtest=find(X.G0==k);
	RES_LOOM_C1=mod_optimal_detection_filter(X.S, X.C1, TemplateLength, ixtrain, ixtest);
	RES_LOOM_C2=mod_optimal_detection_filter(X.S, X.C2, TemplateLength, ixtrain, ixtest);

	fprintf(fid0,'\n\n####################\nResults (XV:LOOM) from export scorings E1 and E2 for each cell [AUC]\n');
	fprintf(fid0,'cell#%d:\t AUC: %g\t%g\n', k, RES_LOOM_C1.test.AUC, RES_LOOM_C2.test.AUC);
	ResTable3(k,1:2) = [RES_LOOM_C1.test.AUC, RES_LOOM_C2.test.AUC]
end

% general classifier obtained from all scored data, one classifier from each of the export scorings.
ixtrain  = find(~isnan(X.G0));
ixtest   = [];
RES_C{1} = mod_optimal_detection_filter(X.S, X.C1, TemplateLength, ixtrain, ixtest);
RES_C{2} = mod_optimal_detection_filter(X.S, X.C2, TemplateLength, ixtrain, ixtest);
RES_C{1}.output = [];
RES_C{2}.output = [];
save(fullfile(outdir,'GeneralClassifier.mat'),'RES_C','Fs','WINLEN','Mode','ResTable3','ResTable2','ResTable1', '-mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Apply MOD to whole (including unscored) data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for E = 1:2,
RES = RES_C{E};
for k1 = 1:length(DATAFILES);
        datafile = fullfile(datapath, DATAFILES{k1});
	[HDR.FILE.Path,HDR.FILE.Name,HDR.FILE.Ext]=fileparts(datafile);
	% we use here the previously cleaned and stored GDF data
	% alternatively, the original data can be loaded, and the clean up procedure
	% (i.e. artifact remove, AP blanking, downsampling, etc) need to be applied here, too
	[data,HDR,rawdata] = sload4mod(datafile,chan,Fs,Mode,TemplateLength);

	% the output parameters from the training steps are
	% RES.CLASSIFIER.A ; 			% filter coefficients of the Wiener filter
	TH    = RES.CLASSIFIER.THRESHOLD;	% detection threshold - based on maxKappa
	delay = RES.CLASSIFIER.delay;		% time shift of the filter

	% the classifier is applied to the whole data set.
	D     = filter(RES.CLASSIFIER.A,1,data);

	% without smoothing
	WinHann = 1;
	D    = [D(delay+1:end);repmat(NaN,delay,1)];
	pos  = get_local_maxima_above_threshold(D,TH,5,0.001*min(WINLEN/2,1)*Fs);

	% save results in a GDF file for visualization with SigViewer
	H0 = [];
	[tmp_p,tmp_f,tmp_e] = fileparts(HDR.FileName);
	H0.FileName = sprintf('%s.%s.w%d.e%d.h%d.%s.train.gdf',fullfile(outdir,tmp_f), 'mod',round(WINLEN*1000),E,WinHann,Mode);
	H2 = minidet2gdf([data, D, D>TH], [], pos, WINLEN, Fs, H0, HDR);

	% apply a smoothing filter as discribed in [1] - this might not be needed
	D    = filter(RES.CLASSIFIER.A,1,data);
	WinHann = 13;
	A    = hann(WinHann); A=A/sum(A);
	nix  = isnan(D); D(nix) = mean(D);
	D    = filtfilt(A,1,D);
	D(nix) = NaN;
	D    = [D(delay+1:end);repmat(NaN,delay,1)];
	pos  = get_local_maxima_above_threshold(D,TH,5,0.001*min(WINLEN/2,1)*Fs);

	% save results in a GDF file for visualization with SigViewer
	H0 = [];
	[tmp_p,tmp_f,tmp_e] = fileparts(HDR.FileName);
	H0.FileName = sprintf('%s.%s.w%d.e%d.h%d.%s.train.gdf', fullfile(outdir, tmp_f), 'mod', round(WINLEN*1000), E, WinHann, Mode);
	H2 = minidet2gdf([data, D, D>TH], [], pos, WINLEN, Fs, H0, HDR);

end
end

apply_general_classifier;

end;
end;
end;
end;
