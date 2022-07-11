function [HDR, s] = detect_sharp_wave_ripple(fn, chan, varargin)
% DETECT_SHARP_WAVE_RIPPLE detects ripples in neural recordings.
%
% HDR = detect_sharp_wave_ripple(filename, chan)
% ... = detect_sharp_wave_ripple(HDR, data)
% ... = detect_sharp_wave_ripple(... ,'-o', outputFilename)
% ... = detect_sharp_wave_ripple(... ,'-e', eventFilename)
% ... = detect_sharp_wave_ripple(... ,'-b', bandpassFilename)
% ... = detect_sharp_wave_ripple(... ,'-s', segmentFilename)
% ... = detect_sharp_wave_ripple(... ,'Threshold', Threshold)
% ... = detect_sharp_wave_ripple(... ,'winlen', winlen)
% ... = detect_sharp_wave_ripple(... ,'trigchan', trigger_channel)
% ... = detect_sharp_wave_ripple(... ,'method', method)
% [HDR, data] = detect_sharp_wave_ripple(...)
%
% Input:
% 	filename: name of source file
%	chan	list of channels that should be analyzed (default is 0: all channels)
%	HDR	header structure obtained by SOPEN, SLOAD, or meXSLOAD
%	data	signal data that should be analyzed
%	Threshold	[default: 6] number of standard deviations
% 	winlen	windowlength [default 400 ms (+-200ms)]
%               should represent the approximate duration of a single event
%		it determines in the default method (0) at which distance in time,
%		time events are considered one or two separate events.
%	trigger_channel: channel used for triggering/segmenting the data
%	coincidenceWindow: [default: 40 ms]
%		if defined, (TYP=4) indicates a ripple where the current-envelope-peak occured
%               within a time window centered on the LFP-envelope peak. If this criteria is not
%               met or the coincidence window is undefined (empty), TYP=3 is used.
%               Only the first two channels are are used, if data set contains more than two
%               channels, a warning is displayed. method=14 does not support coincidence detection
%               and ignores this parameter.
%
%	outputFilename
%		name of file for storing the resulting data with the
%		detected spikes and bursts in GDF format.
%	eventFilename
%		filename to store event inforamation in GDF format. this is similar to
%		the outputFile, except that the signal data is not included and is, therefore,
%		much smaller than the outputFile
%	bandpassFilename
%		filename for bandpass filtered data
%	segmentFilename
%		filename for segmented data, each segment is twice as large as winlen plus one sample
%		going from -winlen*samplerate : +winlen*samplerate
%		This file can be loaded into stimfit.
%	method: 0 : default
%			filter 100-250,
% 			threshold = 6
%               1 or 'maier2011' : see Maier et al 2011 [1]
%			filter 120-300,
% 			threshold = 6
%               2 or 'tukker2013' : see Tukker et al 2013 [2]
%			filter 90-140,
% 			threshold = 5
%               3 or 'lasz2011' : see Lasztoczi et al 2011 [3]
%			filter 90-200,
% 			threshold = 5
%               4, or 'jiangan2015' : see Jian Gan et al 2017 [1]
%			filter 90-200,
% 			threshold = 4	-> HDR.EVENT.TYP=4
% 			threshold = 6	-> HDR.EVENT.TYP=6
%               5  : see Jian Gan et al 2017 [1]
%			filter 100-250,
% 			threshold = 6
%               14 : see Jian Gan et al 2017 [1]
%			filter 100-250,
% 			threshold = 4	-> HDR.EVENT.TYP=4
% 			threshold = 6	-> HDR.EVENT.TYP=6
%
%	Arguments can appear in any order and multiple times (except for filename, chan, HDR and data),
%	In case of conflicting definitions, the latest definition has highest precedence and is used.
%
%
% Output:
%     HDR	header structure as defined in biosig
%     HDR.EVENT includes the detected spikes and bursts.
%     data	signal data, one channel per column
%		between segments, NaN values for 0.1s are introduced
%
% see also:  DETECT_SPIKES_BURSTS, SPIKE2BURSTS, OPTIMUM_ISI_SPIKE_BURST_SEPARATION
%
% References:
% [0] bandpass 100-250 Hz, - 6 S.D. detection, - search with -200 - +200 ms for maximum positive peak for  alignment
%     extract window -200 - +200 ms window of both channels.
% [1] Nikolaus Maier, Alvaro Tejero-Cantero, Anja L. Dorrn, Jochen Winterer, Prateep S. Beed,1 Genela Morris,
%     Richard Kempter, James F.A. Poulet, Christian Leibold, and Dietmar Schmitz
%     Coherent Phasic Excitation during Hippocampal Ripples
%     Neuron 72, 137–152, October 6, 2011 p149.
%     DOI 10.1016/j.neuron.2011.08.016
% [2] John J. Tukker, Balint Lasztoczi, Linda Katona, J. David B. Roberts, Eleftheria K. Pissadaki, Yannis Dalezios,
%     Laszlo Marton, Limei Zhang, Thomas Klausberger, and Peter Somogyi1
%     Distinct Dendritic Arborization and In Vivo Firing Patterns of Parvalbumin-Expressing Basket Cells in the Hippocampal Area CA3
%     The Journal of Neuroscience, April 17, 2013 • 33(16):6809 – 6825 • p6810
%     J. Neurosci., April 17, 2013 • 33(16):6809 – 6825 • 6811
% [3] Balint Lasztoczi, John J. Tukker, Peter Somogyi, and Thomas Klausberger
%     Terminal Field and Firing Selectivity of Cholecystokinin-
%     Expressing Interneurons in the Hippocampal CA3 Area
%     The Journal of Neuroscience, December 7, 2011 • 31(49):18073–18093 • 18074
% [4] Maria Pangalos, José R. Donoso, Jochen Winterer, Aleksandar R. Zivkovic, Richard Kempter, Nikolaus Maier, and Dietmar Schmitz
%     Recruitment of oriens-lacunosum-moleculare interneurons during hippocampal ripples.
%     4398–4403 | PNAS | March 12, 2013 | vol. 110 | no. 11
%     www.pnas.org/cgi/doi/10.1073/pnas.1215496110
% [5] Jian Gan, Shih-ming Weng, Alejandro J. Pernía-Andrade, Jozsef Csicsvari, Peter Jonas
%     Phase-Locked Inhibition, but Not Excitation, Underlies Hippocampal Ripple Oscillations in Awake Mice In Vivo
%     Neuron, Volume 93, Issue 2, 18 January 2017, Pages 308-314
%     https://doi.org/10.1016/j.neuron.2016.12.018


%    Copyright (C) 2014,2015,2016,2020 by Alois Schloegl <alois.schloegl@ist.ac.at>
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



if nargin < 1,
	help detect_sharp_wave_ripple;
end;

if nargin < 2,
	chan = 0;
end;

%%%%% default settings %%%%%
dT = .2;		%%% window length [s]
Threshold = 6; 		%% S.D.
outFile = [];
evtFile = [];
bandpassFile = [];
segFile = [];
trigChan= 0;
method=0;
coincidenceWindow = [];

%%%%% analyze input arguments %%%%%
k = 1;
K = 0;
while k <= length(varargin)
	if ischar(varargin{k})
		if (strcmp(varargin{k},'-o'))
			k = k + 1;
			outFile = varargin{k};
		elseif (strcmp(varargin{k},'-e'))
			k = k + 1;
			evtFile = varargin{k};
		elseif (strcmp(varargin{k},'-b'))
			k = k + 1;
			bandpassFile = varargin{k};
		elseif (strcmp(varargin{k},'-s'))
			k = k + 1;
			segFile = varargin{k};
		elseif (strcmpi(varargin{k},'Threshold'))
			k = k + 1;
			Threshold = varargin{k};
		elseif (strcmpi(varargin{k},'winlen'))
			k = k + 1;
			dT = varargin{k};
		elseif (strcmpi(varargin{k},'coincidenceWindow'))
			k = k + 1;
			coincidenceWindow = varargin{k};
		elseif (strcmpi(varargin{k},'trigchan'))
			k = k + 1;
			trigChan = varargin{k};
		elseif (strncmpi(varargin{k},'method',6))
			k = k + 1;
			method = varargin{k};
		elseif (strcmpi(varargin{k},'bandpass'))
			k = k + 1;
			if (numel(varargin{k}) == 2) && isnumeric(varargin{k})
				bandpassFilter = varargin{k};
			else
				warning(sprintf('paramater bandpass must contain exactly 2 numeric values'));
			end
		else
			warning(sprintf('unknown input argument <%s>- ignored',varargin{k}));
		end;
	elseif isnumeric(varargin{k})
		K = K + 1
		switch (K)
		case {1}
			Threshold = varargin{k}; %% [S.D.]
		case {2}
			dT = varargin{k};		%%% smoothing window length [s]
		otherwise
			warning(sprintf('unknown input argument <%f> - ignored',varargin{k}));
		end;
	end;
	k = k+1;
end;

if ischar(method)
	method = strmatch(lower(method),{'maier2011', 'tukker2013', 'lasz2011', 'jiangan2015'});
end

switch (method)
	case 1 % Maier et al. p.149, [1]
		bandpassFilter = [120,300];
		smoothingwindow = 0.010;
		Threshold = 6;
	case 2 % Tukker et al. 2013 [2]
		bandpassFilter = [90,140];
		smoothingwindow = 0.010;
		Threshold = 5;
	case 3 % Lasztoczi et al. 2011 [3]
		bandpassFilter = [90,200];
		smoothingwindow = 0.010;
		Threshold = 5;
	case 4 % Jian Gan, March 2015 [5]
		bandpassFilter = [90,200];
		smoothingwindow = 0.010;
		Threshold = 4;
	case 14 % Jian Gan, May 2015 [5]
		bandpassFilter = [100,250];
		smoothingwindow = 0.010;
		Threshold = 4;
	otherwise
		method = 0;
		warning('method not specified, use default method');
		bandpassFilter = [100,250];
		Threshold = 6;
end

Fs = 20000; 	% assumed samplerate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if ischar(fn) && exist(fn,'file') && any(size(chan)==1)
		winlen = Fs*.1;
		[s, HDR] = sload(fn, 0, 'NUMBER_OF_NAN_IN_BREAK', winlen);
		if Fs < HDR.SampleRate,
			winlen   = HDR.SampleRate * .1;
			[s, HDR] = sload(fn, 0, 'NUMBER_OF_NAN_IN_BREAK', winlen);
		end;
		if chan==0, chan = 1:HDR.NS; end;
	elseif isstruct(fn)
		HDR = fn;
		s = chan;
		HDR.NS = size(s,2);
		chan = 1:HDR.NS;
	else
		help(mfilename);
	end


	EVENT = HDR.EVENT;
	if ~isfield(EVENT,'DUR');
		EVENT.DUR = zeros(size(EVENT.POS));
	end;
	if ~isfield(EVENT,'CHN');
		EVENT.CHN = zeros(size(EVENT.TYP));
	end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Set Parameters for ripple Detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% compute bandpower, default windowlength for smoothing is 1s.
	% bp = bandpower(center(s), fs, [    bandpassFilter(:),b2(:)]', [], 5);
	S = fft(center(s,1));
	f = [0:size(S,1)-1]*HDR.SampleRate/size(S,1);

	% frequency band bandpassFilter
	SS1 = zeros(size(S));
	%ix  = ((bandpassFilter(1) <= f) & (f <= 149)) | ((151 <= f ) & (f <= bandpassFilter(2)));
	ix  = (bandpassFilter(1) <= f ) & (f <= bandpassFilter(2));
	SS1(ix,:) = S(ix,:);

	Y = real(ifft(SS1));
	BPSignal=Y;


	D = sparse(size(Y));
	T = repmat(NaN,0,3);
	POS = [];
	CHN = [];
	TYP = [];

	if isempty(trigChan) || (trigChan==0)
		%warning('triggerChannel has not been defined');
		trigChan = 1:size(Y,2);
	else
		if ~isempty(coincidenceWindow)
			error('coincidenceWindow and triggerChannel arguments are mutual exclusive - use only one or the other!')
		end
	end;

	%% compute the envelope
	if (method==1)  % Maier et al. p.149, [1]
			w = smoothingwindow*HDR.SampleRate;
			E = filter(ones(w,1),w,abs(Y));
	elseif any(method==[1:3]),	% Tukker et al. 2013 [2], Lasztoczi et al. 2011 [3]
			w = smoothingwindow*HDR.SampleRate;
			E = sqrt(filter(ones(w,1), w, Y.^2));
	elseif any(method==4),	% Jian Gan, 2015
			w = smoothingwindow*HDR.SampleRate;
			E = sqrt(filtfilt(ones(w,1), w, Y.^2));
	elseif any(method==[14,5]),	% Jian Gan, 2015
			w = smoothingwindow*HDR.SampleRate;
			W = hamming(w);
			E = sqrt(filtfilt(W,sum(W), Y.^2));
	end

	% Detection
	for k = chan(:)',
		if any(method==[1:3])
			% start and end point detection
			% "integrated power (i.e. RMS) trace crossed the mean + 1 SD line"
			% because data is bandpassed, mean of raw data is zero, and mean of RMS is 1 SD.
			% thus threshold of 2 SD = 2 RMS is used.
			% It seems that methods [1-3] use the same way of computing the duration
			d2 = E(:,k) > 2 * rms(Y(:,k));
			d3 = diff([0;d2;0]);
			ix1 = find(d3 > 0);	% onset
			ix2 = find(d3 < 0);	% offset

			d1 = E(:,k) > Threshold * std(Y(:,k));
			d3  = diff([0;d1;0]);
			ix1th = find(d3 > 0);	% onset
			ix2th = find(d3 < 0);	% offset
			for kk = 1:length(ix1th)
				[mx, tix] = max( E( ix1th(kk) : ix2th(kk)-1, k ) );
				tix = tix - 1 + ix1th(kk);
				POS = [POS; tix - w/2];
				CHN = [CHN; k];
				T = [T; max(ix1(ix1<tix)), min(ix2(ix2>tix)), k]; 	% start, and end of
			end

		elseif any(method==[4,14]),	% Jian Gan, 2015
			RMS = rms(Y(:,k));
			d1 = E(:,k) > Threshold * RMS;
			d2 = E(:,k) > 2 * RMS;

			d3 = diff([0;d2;0]);

			%% detection with hysteresis
			ix1 = find(diff([0;d1;0]) > 0);	% onset, exceed threshold
			ix2 = find(d3 < 0);	% offset
			ix3 = find(d3 > 0);	% begin of event - used for determining duration

			tt=0; state=0;
			while (tt <= size(Y,1))
				if (state==0)
					%find next onset
					onset = min(ix1(ix1>tt));
					tt=onset;
					state=1;
				else
					%% find next offset
					offset = min(ix2(ix2>tt));
					tt=offset;
					state=0;

					% enter ripple into event table
					[mx, tix] = max( E( onset : offset-1, k ) );
					tix = tix - 1 + onset;
					POS = [POS; tix];
					CHN = [CHN; k];
					typ = 4 + 2*(mx > 6*RMS);
					TYP = [TYP; typ];
					T   = [T; max(ix3(ix3<=tix)), min(ix2(ix2>tix)), k]; 	% start, and end of

				end
			end
		end;

		if (method==1)  % Maier et al. p.149, [1] excludes events with a duration < 12 ms
		if (size(T,1)>0)
			ix = find((T(:,2)-T(:,1)) < 0.012*HDR.SampleRate);
			POS(ix,:)=[];
			CHN(ix,:)=[];
			TYP(ix,:)=[];
			T(ix,:)=[];
		end;
		end

		if (method==0)	% legacy method
			ix = find(abs(Y(:,k)) > Threshold * std(Y(:,k)));
			m  = 1;
			while m <= length(ix);
				ixx = ix( (ix(m) <= ix) & (ix < (ix(m) + dT*HDR.SampleRate)) );
				T   = [T; ixx(1), ixx(end), k];
				[mx,tix] = max(Y(ixx(1):ixx(end),k));
				POS = [POS; tix - 1 + ixx(1)];
				CHN = [CHN; k];
				TYP = [TYP; 3];
				m   = m + length(ixx);
			end;
		end
	end;

	%% select peaks only if there are peaks in both channels within the coincidence window according to [4] Supporting Information p.1
	if (method < 14)
	TYP=repmat(3,size(POS));
	if ~isempty(coincidenceWindow)
		assert(HDR.NS>=2);
		if (HDR.NS>2)
			fprintf(2,'Warning: (%s) coincidence dection is based on first two channels',fn);
		end;
		for k = 1:min(2,HDR.NS)
			for pp=1:length(POS)
				% extract envelope of other channel around +-50 ms , and find maximum envelope
				tix = max(1,POS(pp) - 0.050*HDR.SampleRate/2):min(size(E,1),POS(pp) + 0.050*HDR.SampleRate/2);
				[tmp,tix2] = max(E(tix,HDR.NS+1-k));
				if abs(tix(tix2)-POS(pp)) < coincidenceWindow*HDR.SampleRate/2,
					TYP(pp)=4;
				end
			end
		end
	end
	end;

	HDR.EVENT.POS = [T(:,1); T(:,2); POS; max(1, POS - dT * HDR.SampleRate/2); min(size(s,1), POS + dT*HDR.SampleRate/2)];
	HDR.EVENT.CHN = [T(:,3); T(:,3); CHN; CHN; CHN];
	HDR.EVENT.TYP = [repmat(hex2dec('0001'),size(T,1),1); repmat(hex2dec('8001'),size(T,1),1); TYP; repmat(hex2dec('0002'),size(POS,1),1); repmat(hex2dec('8002'),size(POS,1),1)];
	HDR.EVENT.DUR = zeros(size(HDR.EVENT.POS));
	%HDR.EVENT.CodeDesc = {'AboveThresholdWindow','selectedSegment','PeakEnvelopeOfSWR','SWRinBothChannelsWithinWindow'};
	if isfield(HDR.EVENT,'TimeStamp')
		HDR.EVENT     = rmfield(HDR.EVENT,'TimeStamp');
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Output ripple detection in GDF file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	HDR.T0 = 0;
	HDR = rmfield(HDR,'AS');
	H = HDR;
	if ~isempty(outFile)
		%%% write data to output
		H.TYPE  = 'GDF';
		H.VERSION = 3;
		%[p,f,e]=fileparts(fn);
		H.FILE = [];
		H.FileName  = outFile;
		H.FILE.Path = '';
		H.PhysMax   = max(s);
		H.PhysMin   = min(s);
		H.DigMax(:) =  2^15-1;
		H.DigMin(:) = -2^15;
		H.GDFTYP(:) = 3;
		H.FLAG.UCAL = 0;
		H.NRec = size(s,1);
		H.SPR = 1;
		H.NS  = size(s,2);
		H.Dur = 1/H.SampleRate;
		H = sopen(H,'w');
		if (H.FILE.FID < 0)
			fprintf(2,'Warning can not open file <%s> - GDF file can not be written\n',H.FileName);
		else
			H = swrite(H,s);
			H = sclose(H);
		end;
	end;

	H = HDR;
	if ~isempty(bandpassFile)
		%%% write data to output
		H.TYPE  = 'GDF';
		H.VERSION = 3;
		%[p,f,e]=fileparts(fn);
		H.FILE = [];
		H.FileName  = bandpassFile;
		H.FILE.Path = '';
		H.PhysMax(:)= 12;
		H.PhysMin(:)= -12;
		H.DigMax(:) = 12;
		H.DigMin(:) = -12;
		H.GDFTYP(:) = 16;
		H.Filter.LowPass(:) = bandpassFilter(1);
		H.Filter.HighPass(:) = bandpassFilter(2);
		H.PhysDimCode(:) = 512;
		H = rmfield(H,'PhysDim');
		H.FLAG.UCAL = 0;
		for ch=1:H.NS,
			H.Label{ch} = sprintf('std(%s)',H.Label{ch});
		end;
		H.NRec = size(BPSignal,1);
		H.SPR = 1;
		H.NS  = size(BPSignal,2);
		H.Dur = 1/H.SampleRate;
		H = sopen(H,'w');
		if (H.FILE.FID < 0)
			fprintf(2,'Warning can not open file <%s> - GDF file can not be written\n',H.FileName);
		else
			H = swrite(H,zscore(BPSignal,1));
			H = sclose(H);
		end;
	end;

	H = [];
	if ~isempty(segFile)
		[data,sz] = trigg([s, BPSignal, E],POS(TYP==4 & CHN==2), -dT*HDR.SampleRate, +dT*HDR.SampleRate, 0);

	if isempty(data)
		fprintf(2,'no matching ripple was identified, can not write %s\n',segFile);
		delete(segFile);
	else
		H.NS = sz(1);
		H.SPR= sz(2);
		H.NRec=sz(3);

		%%% write data to output
		H.TYPE  = 'GDF';
		H.VERSION = 2;
		%[p,f,e]=fileparts(fn);
		H.FILE = [];
		H.FileName  = segFile;
		H.FILE.Path = '';
		PhysMax=max(data');
		PhysMin=min(data');
		H.PhysMax = PhysMax + 0.1*(PhysMax-PhysMin);
		H.PhysMin = PhysMin - 0.1*(PhysMax-PhysMin);
		H.DigMax  = H.PhysMax;
		H.DigMin  = H.PhysMin;

		H.FLAG.UCAL = 1;

		H.SampleRate = HDR.SampleRate;
		H.GDFTYP = repmat(16, H.NS, 1);
		H.AS.SPR = repmat(H.SPR,H.NS,1);
		H.Label=HDR.Label;
		%H=rmfield(H,'Dur');
		%H=rmfield(H,'AS');
		%H=rmfield(H,'ELEC');
		%H=rmfield(H,'Cal');
		%H=rmfield(H,'Off');
		%H3.AS.SPR = repmat(H3.SPR,H3.NS,1);

		H.Filter.LowPass(3:4) = bandpassFilter(1);
		H.Filter.HighPass(3:4) = bandpassFilter(2);
		H.Filter.LowPass(5:6) = 1/smoothingwindow;
		H.Filter.HighPass(5:6) = 0;
		H.Filter.Notch(3:4) = HDR.Filter.Notch(1:2);
		H.Filter.Notch(5:6) = 0;

		H.PhysDimCode = HDR.PhysDimCode([1,2,1,2,1,2]);
		H.Label = HDR.Label;
		for ch=1:HDR.NS,
			H.Label{ch+HDR.NS} = sprintf('BP(%s)',H.Label{ch});
			H.Label{ch+2*HDR.NS} = sprintf('RMS(BP(%s))',H.Label{ch});
		end;

		H.EVENT=[];
		H.EVENT.POS = [1:H.NRec-1]'*H.SPR+1;[1:H.NRec]';
		H.EVENT.TYP = repmat(hex2dec('7ffe'), H.NRec-1,1);
		%H.EVENT.POS = [HDR.EVENT.POS;[1:H.NRec]'*H.SPR-floor(H.SPR/2)];
		%H.EVENT.TYP = [HDR.EVENT.POS;TYP];

%		H.LeadIdCode(1:H.NS)=0;
		H = sopen(H,'w');
		if (H.FILE.FID < 0)
			fprintf(2,'Warning can not open file <%s> - GDF file can not be written\n',H.FileName);
		else
			H = swrite(H,data');
			H = sclose(H);
		end;
	end;
	end;


	H = HDR;
	if ~isempty(evtFile)
		%%% write data to output
		H.TYPE  = 'GDF';
		%[p,f,e]=fileparts(fn);
		H.FILE = [];
		H.FileName  = evtFile;
		H.FILE.Path = '';
		H.NRec = 0;
		H.SPR = 1;
		H.NS  = size(s,2);
		H = sopen(H,'w');
		H.EVENT.SampleRate = H.SampleRate;
		if (H.FILE.FID < 0)
			fprintf(2,'Warning can not open file <%s> - GDF file can not be written\n',H.FileName);
		else
			H = sclose(H);
		end;
	end;


