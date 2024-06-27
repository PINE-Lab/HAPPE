% DEMO_MODX shows the essential steps of using MOD for  
% detecting EPSP/EPSCs as described in [1] and applied in [2]. 
% 
% DEMO_MODX is an extended version of the simpler version of DEMO_MOD for  
% and add the LOOM-crossvalidation procedure and handling the scoring
% from two experts. 
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
% Copyright (C) 2016-2022 Alois Schlögl, Institue of Science and Technology Austria (ISTA)
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
BIOSIG_PATH=getenv('BIOSIG_MATLAB_PATH');
pp=strsplit(BIOSIG_PATH,filesep);
if strcmp(pp{end-2},'share'), 
	pp{end-2}='libexec'; 
	pp{1}=filesep;	
	addpath(fullfile(pp{:}));
end
if ~exist('hdr2ascii','file')
        addpath(fullfile(BIOSIG_PATH, 't200_FileAccess'));
end;
addpath(fullfile(BIOSIG_PATH, 'doc'));
if ~exist('mexSLOAD','file')
        addpath(fullfile(BIOSIG_PATH, '../biosig4c++/mex'));
end; 
if ~exist('signal_deconvolution','file')
        addpath(fullfile(BIOSIG_PATH, 't300_FeatureExtraction'))
end; 
if ~exist('rs','file')
        addpath(fullfile(BIOSIG_PATH, 't250_ArtifactPreProcessingQualityControl'))
end;
if ~exist('auc','file')
        addpath(fullfile(BIOSIG_PATH, 't490_EvaluationCriteria'));
end;
if ~exist('auc','file')
        addpath(fullfile(BIOSIG_PATH, 't490_EvaluationCriteria'));
end;
if ~exist('mod_optimal_detection_filter','file')
        addpath(fullfile(BIOSIG_PATH, 't340_MachineLearningOptimalDetection'));
end;


%% check prerequisites
if ~exist('sumskipnan','file') || ~exist('roc','file') || isnan(mean([1,3,NaN])) || (sum(isnan(detrend([1,NaN,3])))>1)
        addpath(fullfile(BIOSIG_PATH, '../../NaN/src'));
        addpath(fullfile(BIOSIG_PATH, '../../NaN/inst'));
end
if ~exist('sumskipnan','file') || ~exist('roc','file') || isnan(mean([1,3,NaN])) || (sum(isnan(detrend([1,NaN,3])))>1)
	addpath(fullfile(pwd,'NaN','src'));
        addpath(fullfile(pwd,'NaN','inst'));
	fprintf(2,'Warning: sumskipnan() is missing - installing NaN-toolbox\n');
end
tmp=roc(randn(100,1),[1:100]'<50);
if ~isfield(tmp,'H_kappa')
	error('roc() from NaN-toolbox is missing - install NaN-toolbox');
end

if ~exist('mexSLOAD','file')
        error('mexSLOAD is missing - mexBiosig or Biosig need to be installed');
end
if ~exist('detect_spikes_bursts','file')
        fprintf(2,'Error: detect_spikes_bursts() is missing - Biosig need to be installed\n');
end


% addpath /fs3/home/schloegl/projects/minis/Xiaomin

% addpath('/fs3/home/schloegl/src/minidet-method-mode-release-0.90/core/')
%  addpath('/fs3/home/schloegl/src/minidet-method-mode-release-0.90/analysis/')
%%%  https://github.com/shababo/psc-detection-dist
 

HIGHPASS=0;
LOWPASS=1000;
if ~exist('WINLEN','var')
	WINLEN=4; 	% window size [ms] for scoring trace 
end;
MAXLAG=1000;	% filter length [samples]
Fs = 25000;	% target sampling rate [Hz]


datapath  = 'data';
mkdir(datapath);	% used to write preprocessed data files, for later use - not manadatory
DATAFILES = {
        '21-Feb-2018_001.dat'
        '23-Nov-2017_002.dat'
        'XM160126-02.dat'
        'XM160218-02.dat'
        'XM_2016-12-14_001.dat'
        'XM_2017-03-28.dat'
};

outdir='output';
mkdir(outdir);
fid0=fopen(fullfile(outdir,'summary.txt'),'w');
fid1=fopen(fullfile(outdir,'summary001.txt'),'w');
fid2=fopen(fullfile(outdir,'summary003.txt'),'a');
fprintf(fid2,'dataset#\tDuration(S1,E1)\t#PSP(S1,E1)\tEventRate(S1,E1)\tmean(log(IEI11))\tstd(log(IEI11))\t Duration(S2,E1)\t#PSP(S2,E1)\tEventRate(S2,E1)\tmean(log(IEI21))\tstd(log(IEI21))\t Duration(S1,E2)\t#PSP(S1,E2)\tEventRate(S1,E2)\tmean(log(IEI12))\tstd(log(IEI12))\t Duration(S2,E2)\t#PSP(S2,E2)\tEventRate(S2,E2)\tmean(log(IEI22))\tstd(log(IEI22))\n');

NAN_BLOCK=repmat(NaN,Fs,1);
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
        G1(k1) = ceil(k1/6);
        G2(k1) = mod(k1-1,6)+1;

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

        tix01 = sort(EVT01.EVENT.POS((EVT01.EVENT.TYP==10) | (EVT01.EVENT.TYP==11) )/EVT01.EVENT.SampleRate);
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

	
        fprintf(fid0,'== %s\n\t%f\t%f\t%f\t%f\t%f\n\t%f\t%f\t%f\t%f\t%f\n',datafile,[min(tix11),max(tix11),tt1,min(tix21),max(tix21),min(tix12),max(tix12),tt2,min(tix22),max(tix22)]);
        fprintf(fid1,'== %s\n\t%f\t%f\t%f\t%f\t%f\n\t%f\t%f\t%f\t%f\t%f\n',[f,e],[min(tix11),max(tix11),tt1,min(tix21),max(tix21),min(tix12),max(tix12),tt2,min(tix22),max(tix22)]);

        fprintf(fid0,'E1:\t%s\n\t%s\nE2:\t%s\n\t%s\n',evtfile01,evtfile01b,evtfile2,evtfile3);
        fprintf(fid0,'Duration [s]\t%f\t%f\n\t\t%f\t%f\n',tix11(end)-tix11(1),tix21(end)-tix21(1),tix12(end)-tix12(1),tix22(end)-tix22(1));
        fprintf(fid0,'NoOfPSPs\t%d\t%d\n\t\t%d\t%d\n',length(tix11),length(tix21),length(tix12),length(tix22));
        fprintf(fid0,'EventRate [Hz]\t%f\t%f\n\t\t%f\t%f\n',length(tix11)/(tix11(end)-tix11(1)),length(tix21)/(tix21(end)-tix21(1)),length(tix12)/(tix12(end)-tix12(1)),length(tix22)/(tix22(end)-tix22(1)));

        [tmp,f1,e1] = fileparts(evtfile01);
        [tmp,f1b,e1b] = fileparts(evtfile01);
        [tmp,f2,e2] = fileparts(evtfile2);
        [tmp,f3,e3] = fileparts(evtfile3);
        fprintf(fid1,'E1:\t%s\n\t%s\nE2:\t%s\n\t%s\n',[f1,e1],[f1b,e1b],[f2,e2],[f3,e3]);
        fprintf(fid1,'Duration [s]\t%f\t%f\n\t\t%f\t%f\n',tix11(end)-tix11(1),tix21(end)-tix21(1),tix12(end)-tix12(1),tix22(end)-tix22(1));
        fprintf(fid1,'NoOfPSPs\t%d\t%d\n\t\t%d\t%d\n',length(tix11),length(tix21),length(tix12),length(tix22));
        fprintf(fid1,'EventRate [Hz]\t%f\t%f\n\t\t%f\t%f\n',length(tix11)/(tix11(end)-tix11(1)),length(tix21)/(tix21(end)-tix21(1)),length(tix12)/(tix12(end)-tix12(1)),length(tix22)/(tix22(end)-tix22(1)));


	iei11=diff(tix11);iei11=iei11(isfinite(iei11));
	iei21=diff(tix21);iei21=iei21(isfinite(iei21));
	iei12=diff(tix12);iei12=iei12(isfinite(iei12));
	iei22=diff(tix22);iei22=iei22(isfinite(iei22));
	

        fprintf(fid2,'%i\t%.1f\t%d\t%.2f\t%.3f\t%.3f\t %.1f\t%d\t%.2f\t%.3f\t%.3f\t %.1f\t%d\t%.2f\t%.3f\t%.3f\t %.1f\t%d\t%.2f\t%.3f\t%.3f\n',k1, [ [tix11(end)-tix11(1), tix21(end)-tix21(1), tix12(end)-tix12(1), tix22(end)-tix22(1)]; [length(tix11), length(tix21), length(tix12), length(tix22)]; [length(tix11)/(tix11(end)-tix11(1)), length(tix21)/(tix21(end)-tix21(1)), length(tix12)/(tix12(end)-tix12(1)), length(tix22)/(tix22(end)-tix22(1))]; exp([mean(log(iei11)), mean(log(iei21)), mean(log(iei12)), mean(log(iei22)) ]); [std(log(iei11)), std(log(iei21)), std(log(iei12)), std(log(iei22)) ] ] );

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         load data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [HDR]    = mexSOPEN(datafile);
        chan     = find(strcmp('Vmon-1',HDR.Label));
        if isempty(chan) chan=1; end        
        [data,HDR]  = mexSLOAD(datafile,0);
        HDR.SampleRate = round(HDR.SampleRate);
        S = data(:,chan);        
        [tmp_p,tmp_f,tmp_e]=fileparts(datafile);


        %%%%% apply AP detection  %%%%%
	AP = detect_spikes_bursts(HDR,data,'chan',chan);
        % AP = detect_spikes_bursts(datafile, chan);
        AP.EVENT.CHN(:)=1; % raw data is stored into channel 1 when writing file 

        if 0
                % use original sampling rate
                Fs = HDR.SampleRate;
                HDR.EVENT.POS = round((HDR.EVENT.POS-1)*Fs/HDR.EVENT.SampleRate)+1;
                HDR.EVENT.DUR = round(HDR.EVENT.DUR*Fs/HDR.EVENT.SampleRate);
                HDR.EVENT.SampleRate = Fs;
        else

                % downsample to 25kHz
                DIV = round(HDR.SampleRate/Fs);
                S = rs(S,DIV,1);
                HDR.EVENT.POS = round((HDR.EVENT.POS-1)*Fs/HDR.EVENT.SampleRate)+1;
                if isfield(HDR.EVENT,'DUR');
                        HDR.EVENT.DUR = round(HDR.EVENT.DUR*Fs/HDR.EVENT.SampleRate);
                end
                HDR.SampleRate = Fs; 
                HDR.EVENT.SampleRate = Fs;

                AP.EVENT.POS = round((AP.EVENT.POS-1)*Fs/AP.EVENT.SampleRate)+1;
                if isfield(HDR.EVENT,'DUR');
                        AP.EVENT.DUR = round(AP.EVENT.DUR*Fs/AP.EVENT.SampleRate);
                end
                AP.EVENT.SampleRate = Fs;
        end;

        tix11 = round(tix11*Fs);
        tix21 = round(tix21*Fs);
        tix12 = round(tix12*Fs);
        tix22 = round(tix22*Fs);

        %S=data(:,chan);        
        S0 = S;        % backup for unfiltered data
        selpos = sort([1;HDR.EVENT.POS(HDR.EVENT.TYP==hex2dec('7ffe')); size(S,1)+1]);
        TemplateLength=MAXLAG;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       remove marked artifacts
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%  artifact markers  %%%%%
        EVT = HDR.EVENT;
        if ~isfield(EVT,'DUR')
                EVT.DUR=zeros(size(EVT.POS));
        end
        aix = find(EVT.TYP==20);
        artifact = [EVT.POS(aix),EVT.POS(aix)+EVT.DUR(aix)];
        Stmp = S; 
        for k=1:length(aix)
                Stmp(EVT.POS(aix):EVT.POS(aix)+EVT.DUR(aix),:)=NaN;
        end;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%       remove AP's and spikelets
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	apix = find(AP.EVENT.TYP==hex2dec('201'));
	Stmp(AP.EVENT.POS(apix),:) = NaN;       % trough the template (40ms, 1000samples) the NaN will be expanded to [-10:+30] ms
	[ix,iy] = meshgrid(AP.EVENT.POS(apix),[-round(TemplateLength/4):round(TemplateLength*3/4)]);
	Stmp(ix(:)+iy(:)) = NaN;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%  50 Hertz Notch and Lowpass %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if (HIGHPASS>0)
	        for k2 = 1:length(selpos)-1,
        	        d = Stmp(selpos(k2):selpos(k2+1)-1,:);
                	dnan=isnan(d);
	                if any(dnan)
        	                d(dnan) = mean(d);
                	end
	                D = fft(d); 
        	        ff = [0:length(d)-1]*Fs/length(d);
                	%fix = ((ff > 49.5) & (ff < 50.5)) | (ff > Fs/2);
	                fix = (HIGHPASS < ff) | (ff > LOWPASS);
        	        D( fix ) = 0; 
                	D(1) = D(1)/2;
	                d = real(ifft(D))*2;
        	        if any(dnan)
                	        d(dnan)=NaN;
	                end
        	        Stmp(selpos(k2):selpos(k2+1)-1,:)=d;
	        end;
        else
                B=fir1 (100, 2*LOWPASS/Fs);
                Stmp=filter(B,1,Stmp(end:-1:1));
                Stmp=filter(B,1,Stmp(end:-1:1));
	end

        %%%%% check duration %%%%%
        EVT.TYP(EVT.POS < ((tix2(1)  -0.4)*Fs)) = 0;
        EVT.TYP(EVT.POS > ((tix3(end)+0.4)*Fs)) = 0;
        pos = EVT.POS(EVT.TYP>0);
        L1  = (max(pos)-min(pos))/Fs;
        pos = EVT.TimeStamp(EVT.TYP>0);
        L2  = (max(pos)-min(pos))*24*3600;
        fidres = fopen(sprintf('%s/%s.check.txt',outdir,f),'w');
        fprintf(fidres, '#sanity check of recording length\nrecording_duration: %f s\nwallclock: %f\nratio %f\n',L1,L2,L1/L2);
        fprintf(fidres, 'Number of AP (and spikelets?): %d\n',length(apix));
        fprintf(fidres, 'artifact_duration: %f\n',sum(EVT.DUR(aix))/EVT.SampleRate);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         run EPSP detection, blank AP's 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %assert(isequal(HDR.PhysDim(chan),{'V'}))

        % use scorings from both expoerts and use only common intervals for the analysis
        % t2 indicates the first scoring segment, t3 indicates the 2nd segment 
        dt = median(diff([tix11;tix21;tix12;tix22]));
        tt = max([tix11(1)-dt,tix12(1)-dt,1]) : min([tix21(end)+dt,tix22(end)+dt,length(S)]);
        t2 = max([tix11(1)-dt,tix12(1)-dt,1]) : min([tix11(end)+dt,tix12(end)+dt,length(S)]);
        t3 = max([tix21(1)-dt,tix22(1)-dt,1]) : min([tix21(end)+dt,tix22(end)+dt,length(S)]);

	if fidres>2, fclose(fidres); end

if isempty(t3) || isempty(t2) 
        fprintf(fid0,'# ERROR: t2 or t3 is empty %d/%d\n',length(t2),length(t3));
        continue
end

        %%% Detrend data 
        s = Stmp; 
        s(1:tt(1)-1)=NaN; s(tt(end)+1:end)=NaN;  
        S(1:tt(1)-1)=NaN; S(tt(end)+1:end)=NaN;  
        s = detrend(s);

        tix  = [t2(1):t2(end),t3(1):t3(end)]';
        tix1 = [t2(1):t2(end)]';
        tix2 = [t3(1):t3(end)]';
        % s(tix(isnan(s(tix)))) = mean(s(tix));

        if 0, G1(k1)==1, %~isempty(FLAG_DAVID.BLANK) 
                selpos = sort([1;HDR.EVENT.POS(HDR.EVENT.TYP==hex2dec('7ffe')); size(S,1)+1]);
                sweepLength=diff(selpos)/Fs;
                if all(sweepLength(1)==sweepLength) sweepLength=sweepLength(1); 
                else
                        fprintf(fid0,'# ERROR: sweeplength not equal in file \n#\t',datafile);
                        fprintf(fid0,'%f ',sweepLength);
                        HIS=histo3(sweepLength(:));
                        find(sweepLength~=sweepLength(1))
                        sweepLength=median(sweepLength);
                end 

                FLAG_DAVID.BLANK=[0.7, sweepLength-0.3];
                %FLAG_DAVID.SELPOS=sweepLength;

                % data only from FLAG_DAVID_BLANK(1):FLAG_DAVID.BLANK(2) e.g. 0.5-19.5 s is used

                tmpix = mod(tix11,round(sweepLength*Fs))/Fs;
                tix11 = tix11((FLAG_DAVID.BLANK(1) < tmpix) & (tmpix<FLAG_DAVID.BLANK(2)));

                tmpix = mod(tix12,round(sweepLength*Fs))/Fs;
                tix12 = tix12((FLAG_DAVID.BLANK(1) < tmpix) & (tmpix<FLAG_DAVID.BLANK(2)));

                tmpix = mod(tix21,round(sweepLength*Fs))/Fs;
                tix21 = tix21((FLAG_DAVID.BLANK(1) < tmpix) & (tmpix<FLAG_DAVID.BLANK(2)));

                tmpix = mod(tix22,round(sweepLength*Fs))/Fs;
                tix22 = tix22((FLAG_DAVID.BLANK(1) < tmpix) & (tmpix<FLAG_DAVID.BLANK(2)));

        end


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

        %%% mamud_results018_lda_x1.m
        T1=-.050*Fs;
        T2= .050*Fs;

        if 0
                ; %noop

        elseif 1, %  store clean data to GDF so that we can avoid doing the
		  %  cleanup again when applying the classifier to the whole data set
		[HDR.FILE.Path,HDR.FILE.Name,HDR.FILE.Ext]=fileparts(HDR.FileName);
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
		% copyfile(evtfile01,'data');
		% copyfile(evtfile2,'data');
		% copyfile(evtfile3,'data');


        elseif 1, % trigger and store to GDFmamud_results018_lda_x1.m
                clear H1; 
                ixx=[tix11;tix21; inf];
                tmpix=diff(ixx)>T1;
                ixx=ixx(tmpix(1:end-1) & tmpix(2:end));

                [XX,sz]=trigg(S0, ixx, T1,T2);
                XX=squeeze(reshape(XX,sz));
                tmp=ixx/Fs;
                save('-ascii',sprintf('mamud/stimfit/%s.e1.stimfit.txt',f),'tmp');
                H1.FileName=sprintf('mamud/stimfit/%s.e1.stimfit.gdf',f);
                H1.TYPE='GDF';
                H1.VERSION=3;
                H1.NS=1; 
                H1.SampleRate=Fs; 
                H1.NRec=sz(3);
                H1.SPR=sz(2);
                H1.PhysMax = HDR.PhysMax(chan);
                H1.PhysMin = HDR.PhysMin(chan);
                H1.DigMax  = HDR.DigMax(chan);
                H1.DigMin  = HDR.DigMin(chan);
                H1.PhysDimCode = HDR.PhysDimCode(chan);
                H1.Label = HDR.Label(chan);
                H1.Cal  = HDR.Cal(chan);
                H1.Off  = HDR.Off(chan);
                H1.GDFTYP  = HDR.GDFTYP(chan);
                H1.EVENT.POS=[1:sz(3)-1]'*sz(2)+1;
                H1.EVENT.TYP=repmat(hex2dec('7ffe'),sz(3)-1,1);
                H1.EVENT.CHN=ones(sz(3)-1,1);
                H1.EVENT.DUR=zeros(sz(3)-1,1);
                H1.EVENT.TYP=repmat(hex2dec('7ffe'),sz(3)-1,1);
                H1.FLAG.UCAL=0;
                if isfield(H1,'AS'); H1=rmfield(H1,'AS'); end
                H1=sopen(H1,'w'); H1=swrite(H1,XX(:)); H1=sclose(H1);

                clear H1; 
                ixx=[tix12;tix22;inf];
                tmpix=diff(ixx)>T1;
                ixx=ixx(tmpix(1:end-1) & tmpix(2:end));

                [XX,sz]=trigg(S0, ixx, T1,T2);
                XX=squeeze(reshape(XX,sz));
                tmp=ixx/Fs;
                save('-ascii',sprintf('mamud/stimfit/%s.e2.stimfit.txt',f),'tmp');
                H1.FileName=sprintf('mamud/stimfit/%s.e2.stimfit.gdf',f);
                H1.TYPE='GDF';
                H1.VERSION=3;
                H1.NS=1; 
                H1.SampleRate=Fs; 
                H1.NRec=sz(3);
                H1.SPR=sz(2);
                H1.PhysMax = HDR.PhysMax(chan);
                H1.PhysMin = HDR.PhysMin(chan);
                H1.DigMax  = HDR.DigMax(chan);
                H1.DigMin  = HDR.DigMin(chan);
                H1.PhysDimCode = HDR.PhysDimCode(chan);
                H1.Label = HDR.Label(chan);
                H1.Cal  = HDR.Cal(chan);
                H1.Off  = HDR.Off(chan);
                H1.GDFTYP  = HDR.GDFTYP(chan);
                H1.EVENT.POS=[1:sz(3)-1]'*sz(2)+1;
                H1.EVENT.TYP=repmat(hex2dec('7ffe'),sz(3)-1,1);
                H1.EVENT.CHN=ones(sz(3)-1,1);
                H1.EVENT.DUR=zeros(sz(3)-1,1);
                H1.EVENT.TYP=repmat(hex2dec('7ffe'),sz(3)-1,1);
                H1.FLAG.UCAL=0;
                if isfield(H1,'AS'); H1=rmfield(H1,'AS'); end
                H1=sopen(H1,'w'); H1=swrite(H1,XX(:)); H1=sclose(H1);
        end                        

        % identify time offset "toff" between experts, 
        %   and compute interscorer agreement using Cohen's kappa 
        %    before and after correction for time shift
        [xcf,lag]     = xcorr(c1(tix),c2(tix),0.02*Fs,'coeff');
        [tmp,tmpix]   = max(xcf);
        toff1         = tmpix - find(lag==0);                % that's the time offset between scorers
        toff          = lag(tmpix);                % that's the time offset between scorers
        
        RES.lag(k1,1) = toff;

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
        X.isnan = [X.isnan; isnan(Stmp(tix1));             NAN_BLOCK; isnan(Stmp(tix2));             NAN_BLOCK];
end

for k1=1:length(DATAFILES);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%       within-cell cross-validation (XV)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% XV: S1->S2 cross-validiation
	ixtrain=find(X.G2==1 & X.G0==k1);
	ixtest =find(X.G2==2 & X.G0==k1);
	RES12S_C1=mod_optimal_detection_filter(X.S, X.C1, TemplateLength, ixtrain, ixtest);
	RES12S_C2=mod_optimal_detection_filter(X.S, X.C2, TemplateLength, ixtrain, ixtest);
	% XV: S2->S1 cross-validiation
	ixtrain=find(X.G2==2 & X.G0==k1);
	ixtest =find(X.G2==1 & X.G0==k1);
	RES21S_C1=mod_optimal_detection_filter(X.S, X.C1, TemplateLength, ixtrain, ixtest);
	RES21S_C2=mod_optimal_detection_filter(X.S, X.C2, TemplateLength, ixtrain, ixtest);

	% In [1,2], we combined the output from the two testing sets, and 
	% and applied the ROC analysis to this trace. 
	% for the sake of simplicity, this is omitted here. 

	% XV: A1B2->A2B1 cross-validiation
	ixtrain=find(X.G3==1 & X.G0==k1);
	ixtest =find(X.G3==0 & X.G0==k1);
	RES1A2B_C1=mod_optimal_detection_filter(X.S, X.C1, TemplateLength, ixtrain, ixtest);
	RES1A2B_C2=mod_optimal_detection_filter(X.S, X.C2, TemplateLength, ixtrain, ixtest);
	% XV: A2B1->A1B2 cross-validiation
	ixtrain=find(X.G3==0 & X.G0==k1);
	ixtest =find(X.G3==1 & X.G0==k1);
	RES2A1B_C1=mod_optimal_detection_filter(X.S, X.C1, TemplateLength, ixtrain, ixtest);
	RES2A1B_C2=mod_optimal_detection_filter(X.S, X.C2, TemplateLength, ixtrain, ixtest);

	% In [1,2], we combined the output from the two testing sets, and 
	% and applied the ROC analysis to this trace. 
	% for the sake of simplicity, this is omitted here. 

	fprintf(fid0,'\n#####\nResults from export scoring 1 of %s\n',DATAFILES{k1});
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

fprintf(fid0,'\n\n####################\nResults (XV:LOOM) from export scorings E1 and E2 for each cell [AUC]\n');
for k=1:length(DATAFILES),
	% XV: LOOM
	% this is not implemented here, essentially the idea is to 
	% concatenate the data from multiple cells, and use the data
	% of m-1 cells for training, and the data from the remaining cell for
	% testing.  

	ixtrain=find(X.G0~=k & ~isnan(X.G0));
	ixtest=find(X.G0==k);
	RES_LOOM_C1=mod_optimal_detection_filter(X.S, X.C1, TemplateLength, ixtrain, ixtest);
	RES_LOOM_C2=mod_optimal_detection_filter(X.S, X.C2, TemplateLength, ixtrain, ixtest);

	fprintf(fid0,'cell#%d:\t AUC: %g\t%g\n', k, RES_LOOM_C1.test.AUC, RES_LOOM_C2.test.AUC);
	ResTable3(k,1:2)=[RES_LOOM_C1.test.AUC, RES_LOOM_C2.test.AUC];
end

% general classifier obtained from all scored data, one classifier from each of the export scorings. 
ixtrain=find(~isnan(X.G0));
ixtest=[];
RES_C{1}=mod_optimal_detection_filter(X.S, X.C1, TemplateLength, ixtrain, ixtest);
RES_C{2}=mod_optimal_detection_filter(X.S, X.C2, TemplateLength, ixtrain, ixtest);
RES_C{1}.output=[];
RES_C{2}.output=[];
save('output/GeneralClassifier.mat','RES_C','LOWPASS','HIGHPASS','Fs','WINLEN','Mode', '-mat');

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
	[data,HDR]  = mexSLOAD(fullfile('data',[HDR.FILE.Name,'.gdf']),0);
	HDR.SampleRate = round(HDR.SampleRate);
	AP = detect_spikes_bursts(HDR,data,'chan',chan);
	%% make sure data is clean, no AP, no other artifact that could influence the trend
        S = detrend(data);
                                
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
	H0 = [];
	[tmp_p,tmp_f,tmp_e] = fileparts(HDR.FileName);
	H0.FileName = sprintf('output/%s.%s.w%d.e%d.h%d.gdf',tmp_f, 'mod',round(WINLEN*1000),E,WinHann);
	H2 = minidet2gdf([data, D, D>TH], [], pos, WINLEN, Fs, H0, AP);

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
	H0.FileName = sprintf('output/%s.%s.w%d.e%d.h%d.gdf', tmp_f, 'mod', round(WINLEN*1000), E, WinHann);
	H2 = minidet2gdf([s, D, D>TH], [], pos, WINLEN, Fs, H0, AP);
end
end


