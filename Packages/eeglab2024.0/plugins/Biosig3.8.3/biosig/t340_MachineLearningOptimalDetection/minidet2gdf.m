function HDR = minidet2gdf(S,POS,pos,WINLEN,Fs,fn, AP_Events)
% MINIDET2GDF writes results of miniature detection results to a GDF file. 
%   
%
% Usage: 
%     HDR = minidet2gdf(S,POS,pos,Twin,Fs,fn, AP_Events)
%     HDR = minidet2gdf([S, D, D>TH], [], pos, Twin, Fs, H0, []);
%  
%  Input: 
%    S: raw data
%    D: raw detection trace
%    TH: detection threshold, such that (D>TH) produced binary detection trace
%    POS: (optional) export scoring
%    pos: detected event times, typically these are obtained by
%	  pos = get_local_maxima_above_threshold(D,TH,1,0.01*Fs);         
%    Twin: windows length when translating expert scoring into scoring trace. 
%    Fs:  sampling rate in Hertz
%    HDR: out file structure such that 
%       HDR.FileName indicates the output file name. 
% 	The resulting file will be stored in GDF format [2,3], which can be read
%	by Biosig, Sigviewer and other software. 
%    AP_Events: if AP's have been detected
%	and should be indicated as such in the output file.  
%
%  Output: 
%    HDR: header structure of the output file HDR.FileName. 
%
%
% see also: 
%	AUC, MINIDET2GDF, GET_LOCAL_MAXIMA_ABOVE_THRESHOLD, 
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
% [3] A. Schlögl, GDF - A general dataformat for biosignals, Technical Report.
%     Version 2.x (2006) available at: GDF 2.0
% [4] ON K2204:2015 A general dataformat for biosignals, Austrian Standards, 2015


%    Copyright (C) 2017-2021,2023 by Alois Schloegl <alois.schloegl@ist.ac.at>
%    This is part of the BIOSIG-toolbox http://biosig.sf.net/
%
%    BioSig is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    BioSig is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with BioSig.  If not, see <http://www.gnu.org/licenses/>.




% HDR.FileName = sprintf('/fs3/home/schloegl/projects/minis/Xiaomin/optim_mini_detect.%s_%5.3f.gdf',METHODS{k3},WINLEN(k1));
if isstruct(fn)
	HDR=fn;
else
	HDR.FileName = fn;
end; 

if nargin<7
	AP_Events = [];
end;
if isfield(AP_Events,'EVENT');
	AP_Events = AP_Events.EVENT;
end
if ~isfield(AP_Events,'POS') || ~isfield(AP_Events,'TYP') || ~isfield(AP_Events,'DUR') || ~isfield(AP_Events,'CHN')
	AP.POS=[];
	AP.TYP=[];
	AP.DUR=[];
	AP.CHN=[];
	fprintf('Warning, No AP events provided');
else
	AP = AP_Events;
	% select only events 
	ix = ((hex2dec('200') <= AP.TYP) & (AP.TYP < hex2dec('020F')));
	AP.POS = round((AP.POS(ix)-1)*Fs/AP.SampleRate)+1;
	AP.TYP = AP.TYP(ix);
	AP.DUR = round(AP.DUR(ix)*Fs/AP.SampleRate);
	AP.CHN = AP.CHN(ix);
end

VER   = version;
cname = computer;

% select file format 
HDR.TYPE='GDF';  
% description of recording device 
HDR.Manufacturer.Name = 'BioSig';
HDR.Manufacturer.Model = 'demo_mod';
HDR.Manufacturer.Version = '1.1';
HDR.Manufacturer.SerialNumber = '00000000';

% recording identification, max 80 char.
HDR.RID = 'TestFile 001'; %StudyID/Investigation [consecutive number];
HDR.REC.Hospital   = 'Biosig Laboratory';
HDR.REC.Technician = get_current_username();

HDR.REC.Equipment  = 'biosig';
HDR.REC.IPaddr	   = [127,0,0,1];	% IP address of recording system 	
HDR.Patient.Name   = 'anonymous';  
HDR.Patient.Id     = 'X';	
HDR.Patient.Birthday = [0 0 0 0 0 0];
HDR.Patient.Weight = 0; 	% undefined 
HDR.Patient.Height = 0; 	% undefined 
HDR.Patient.Sex    = 0; 	% 0: undefined,	1: male, 2: female 
HDR.Patient.Birthday = zeros(1,6); %    undefined 
HDR.Patient.Impairment.Heart = 0;  %	0: unknown 1: NO 2: YES 3: pacemaker 
HDR.Patient.Impairment.Visual = 0; %	0: unknown 1: NO 2: YES 3: corrected (with visual aid) 
HDR.Patient.Smoking = 0;           %	0: unknown 1: NO 2: YES 
HDR.Patient.AlcoholAbuse = 0; 	   %	0: unknown 1: NO 2: YES 
HDR.Patient.DrugAbuse = 0; 	   %	0: unknown 1: NO 2: YES 
HDR.Patient.Handedness = 0; 	   % 	unknown, 1:left, 2:right, 3: equal

% recording time [YYYY MM DD hh mm ss.ccc]
if ~isfield(HDR,'T0') HDR.T0 = clock;	end

% number of channels
HDR.NS = size(S,2);

assert(HDR.NS==3)


% Duration of one block in seconds
HDR.SampleRate = Fs;
HDR.SPR = 1;
HDR.Dur = HDR.SPR/HDR.SampleRate;

% Samples within 1 block
HDR.AS.SPR = ones(HDR.NS,1);	% samples per block; 0 indicates a channel with sparse sampling 
%HDR.AS.SampleRate = [1000;100;200;100;20;0];	% samplerate of each channel

% channel identification, max 80 char. per channel
HDR.Label={'raw';'detection trace';'x>TH'};

% Transducer, mx 80 char per channel
HDR.Transducer = {'-'; '-'; '-'};

% define datatypes (GDF only, see GDFDATATYPE.M for more details)
HDR.GDFTYP = 16*ones(1,HDR.NS);

% define scaling factors 
p1=max(S);m1=min(S);
p2=p1 + (p1-m1)/100;
m2=m1 + (m1-p1)/100;

HDR.PhysMax = p2;
HDR.PhysMin = m2;
HDR.DigMax  = p2;
HDR.DigMin  = m2;
HDR.FLAG.UCAL = 1; 	% data S is already converted to internal (usually integer) values (no rescaling within swrite);
HDR.FLAG.UCAL = 0; 	% data S will be converted from physical to digital values within swrite. 
% define filter settings 
HDR.Filter.Lowpass = [NaN,NaN,NaN];
HDR.Filter.Highpass = [0,0,0];
HDR.Filter.Notch = [0,0,0];
% define sampling delay between channels  
HDR.TOffset = [0,0,0];

% define physical dimension
HDR.PhysDim = {'a.u.';'1';'1'};	%% must be encoded in unicode (UTF8)
HDR.Impedance = repmat(NaN,1,HDR.NS);         % electrode impedance (in Ohm) for voltage channels 
HDR.fZ = repmat(NaN,1,HDR.NS);                % probe frequency in Hz for Impedance channel

[p,f,e]=fileparts(HDR.FileName);
mini_epsp_position = pos/Fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  	EPSP Amplitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% baseline time interval
	T01=-0.008;
	T02=-0.004;
	% peak time interval
	T03=-0.001;
	T04=+0.001;

	% save triggered data 
	T=[-0.04*Fs,0.04*Fs];
	t=[T(1):T(2)]'/Fs;
	[x,sx]=trigg(S,pos,T(1),T(2));
	x=reshape(x,sx);

	baseline  = mean(x(1,(T01<t)&(t<T02),:),2);
	amp0      = max(x(1,(T03<t)&(t<T04),:),[],2);
	amplitude = squeeze(amp0-baseline);

PLOTFLAG=0;
	if PLOTFLAG
		subplot(221), plot(t,squeeze(x(1,:,:)));title(HDR.Label{1});
		subplot(223), plot(t,squeeze(x(2,:,:)));title(HDR.Label{2});
		xlabel('time t[s]');	
	end

	T=[-0.010*Fs,0.010*Fs];t=[T(1):T(2)]'/Fs;
	[x,sx]=trigg(S,pos,T(1),T(2));
	x=reshape(x,sx);

	if PLOTFLAG
		subplot(222), plot(t,squeeze(x(1,:,:)));title(HDR.Label{1});
		subplot(224), plot(t,squeeze(x(2,:,:)));title(HDR.Label{2});
		xlabel('time t[s]');	
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  	microStimfit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	TimeConstant = repmat(NaN,size(baseline));

	tt0=0.04*Fs+1;
	baseBeg=-200+tt0;
	baseEnd=-100+tt0;
	peakBeg=-100+tt0;
	peakEnd=100+tt0;
	meanN=25+tt0;
	dir=1; 
	plotFlag=0;
	fitFlag=1;
	%thres, thresFlag

	if 0, %for k=1:length(pos)
		trace = squeeze(x(1,:,k));
		[base, baseSD, peak, tPeak, myf, fExp, ampl, offs, tau, fitResult, threshold, gr1, gr2, t20Int, t20Real, t80Int, t80Real, t50AInt, t50AReal, t50BInt, t50BReal, t0Real, tMaxSlopeRiseReal, tThreshold, yThreshold, maxSlopeRiseReal]= microstimfit4om(trace, baseBeg,  baseEnd, peakBeg, peakEnd, meanN, dir, plotFlag, fitFlag);

		TimeConstant(k)=tau; 
	end


%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  	OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if PLOTFLAG
	print('-dpng',fullfile(p,[f,'.minidet.png']));
end
save('-v4',fullfile(p,[f,'.minidet.mat']),'pos','Fs','baseline','amplitude','amp0');
save('-ascii',fullfile(p,[f,'.minidet.txt']),'mini_epsp_position');



HDR.VERSION   = 2.22;        % experimental  
HDR.EVENT=[];
HDR.EVENT.POS = [AP.POS; POS-WINLEN*Fs/2;pos];
% ### 0x021_      miniature post-synaptic events
% 0x0211  EPSP - Excitatory Post-Synaptic Potentials
% 0x0212  IPSP - Inhibitory Post-Synaptic Potentials
% 0x0213  EPSC - Excitatory Post-Synaptic Currents
% 0x0214  IPSC - Inhibitory Post-Synaptic Currents
HDR.EVENT.TYP = [AP.TYP; repmat(10,length(POS),1);repmat(hex2dec('211'),length(pos),1)];
HDR.EVENT.CHN = [AP.CHN; repmat(1,length(POS),1);repmat(2,length(pos),1);];
HDR.EVENT.DUR = [AP.DUR; repmat(WINLEN*Fs,length(POS),1);zeros(length(pos),1)];

if 0, %try,
	mexSSAVE(HDR,S);
else %catch

	% write results file of mod, useful for debugging
	HDR1 = sopen(HDR,'w');
	HDR1 = swrite(HDR1,S);
	HDR1 = sclose(HDR1);

	% write event file
	HDR2=HDR;
	HDR2.NS=0;
	HDR2.FileName = fullfile(p,[f,'.evt']);
	HDR2.EVENT.CHN(:)=0;

	HDR1 = sopen(HDR2,'w');
	HDR1 = sclose(HDR1);
end;


[s0,HDR0] = sload(HDR.FileName);	% test file 

HDR0      = sopen(HDR0.FileName,'r');
[s0,HDR0] = sread(HDR0);
HDR0      = sclose(HDR0); 

% plot(s0-S)


