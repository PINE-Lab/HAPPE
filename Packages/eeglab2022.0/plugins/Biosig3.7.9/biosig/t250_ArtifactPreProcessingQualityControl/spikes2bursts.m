function [H, s] = spikes2bursts(fn, varargin)
% spikes2bursts convert spike trains into bursts. 
%  Spikes with an interspike interval smaller  
%  than 75 ms are considered a burst. The results are stored as an 
% event table. 
%
% HDR = spikes2bursts(filename)
% ... = spikes2bursts(HDR)
% ... = spikes2bursts(... [, dT_Burst ])
% ... = spikes2bursts(... [, dT_Burst [, dT_Exclude] ])
% ... = spikes2bursts(... , 'dT_Burst', dT_Burst)
% ... = spikes2bursts(... , 'dT_Exclude', dT_Exclude)
% ... = spikes2bursts(filename, ... , '-o',outputFilename)
% ... = spikes2bursts(... , '-e',eventFilename)
% ... = spikes2bursts(... , '-b',burstFilename)
%   
% Input: 
%	filename  name of file containing spikes in the event table
%	HDR	header structure obtained by SOPEN, SLOAD, or meXSLOAD
%	dT_Burst	[default: 50e-3 s] am inter-spike-interval (ISI) exceeding this value, 
%		marks the beginning of a new burst
%               'OPTIMUM_ISI' will use the function OPTIMUM_ISI_SPIKE_BURST_SEPARATION
%               for identifying dT_Burst.
%	dT_Exclude an interspike interval smaller than this value, indicates a
%		double detection, and the second detection is deleted.
%		in case of several consecutive ISI's smaller than this value,
%		all except the first spikes are deleted.
%	outputFilename
%		filename to store the original data together with the newly compuated
%		burst information. (Remark: this option is only supported if
% 		filename (not HDR) is given as first argument 
%	eventFilename
%		filename to store event inforamation in GDF format. this is similar to 
%		the outputFile, except that the signal data is not included and is, therefore,
%		much smaller than the outputFile
%	burstFilename
%		filename for the "burst table", containing basic properties of each burst,
%		(it is an ASCII file in <tab>-delimited format)
%
% Output:
%     HDR	header structure as defined in biosig
%     HDR.EVENT includes the detected spikes and bursts. 
%     HDR.BurstTable contains for each burst (each in a row) the following 6 numbers:
%	channel number, sweep number, OnsetTime within sweep [s], 
%	number of spikes within burst, average inter-spike interval (ISI) [ms], 
%	and minimum ISI [ms].
%     data	signal data, one channel per column
%		between segments, NaN values for 0.1s are introduced
%	
% see also:  DETECT_SPIKES_BURSTS, SPIKE2BURSTS, OPTIMUM_ISI_SPIKE_BURST_SEPARATION
%
% References: 
% [1] Janina Kowalski, Jian Gan, Peter Jonas, and Alejandro J. Pernía‐Andrade
%     Intrinsic membrane properties determine hippocampal differential firing pattern in vivo in anesthetized rats
%     Hippocampus. 2016 May; 26(5): 668–682.
%     doi:10.1002/hipo.22550

%	Copyright (C) 2011,2020 by Alois Schloegl <alois.schloegl@ist.ac.at>
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/
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
	help spikes2bursts;
end;

%%%%% default settings %%%%%
dT_Burst   = 50e-3;	%%% smaller ISI is a burst [s]
dT_Exclude = [];	%%% for identifying double detections, spikes with smaller ISI are excluded
outFile    = [];
evtFile    = [];
burstFile  = [];


%%%%% analyze input arguments %%%%%
k = 1;
K = 0;
while k <= length(varargin)
	if ischar(varargin{k})
		if (strcmp(varargin{k},'-e'))
			k = k + 1;
			evtFile = varargin{k};
		elseif (strcmp(varargin{k},'-b'))
			k = k + 1;
			burstFile = varargin{k};
		elseif (strcmp(varargin{k},'-o'))
			k = k + 1;
			outFile = varargin{k};
		elseif (strcmpi(varargin{k},'dT_Burst'))
			k = k + 1;
			dT_Burst = varargin{k};
		elseif (strcmpi(varargin{k},'dT_Exclude'))
			k = k + 1;
			dT_Exclude = varargin{k};
		else
			warning(sprintf('unknown input argument <%s>- ignored',varargin{k}));
		end;
	elseif isnumeric(varargin{k})
		K = K + 1
		switch (K)
		case {1}
			dT_Burst = varargin{k};
		case {2}
			dT_Exclude = varargin{k};
		otherwise
			warning(sprintf('unknown input argument <%f> - ignored',varargin{k}));
		end;
	end;
	k = k+1;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%	load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if ischar(fn) && exist(fn,'file')
		[s, HDR] = sload(fn);
	elseif isstruct(fn)
		HDR = fn;
	else
		help(mfilename); 
	end

	EVENT = HDR.EVENT;
	ix = (EVENT.TYP ~= hex2dec('0202'));
	%if ~all(ix) warning('previously defined bursts are ignored'); end;
	EVENT.POS = EVENT.POS(ix);
	EVENT.TYP = EVENT.TYP(ix);
	if ~isfield(EVENT,'DUR');
		EVENT.DUR = zeros(size(EVENT.POS));
	else
		EVENT.DUR = HDR.EVENT.DUR(ix);
	end; 
	if ~isfield(EVENT,'CHN');
		EVENT.CHN = zeros(size(EVENT.TYP));
	else
		EVENT.CHN = HDR.EVENT.CHN(ix);
	end; 

	Fs = HDR.SampleRate; 
	chan = unique(EVENT.CHN);

	if strcmpi(dT_Burst,'OPTIMUM_ISI');
		dt_Burst = optimum_isi_spike_burst_separation(get_inter_spike_intervals(HDR));
	end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%	Set Parameters for Burst Detection 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	HDR.BurstTable = [];

	for ch = chan(:)';	% look for each channel 	
		OnsetSpike = EVENT.POS((EVENT.CHN==ch) & (EVENT.TYP==hex2dec('0201'))); 	%% spike onset time [samples]
		if isempty(OnsetSpike)
			continue;
		end;

		if ~isempty(dT_Exclude)
			% --- remove double detections < 1 ms
			OnsetSpike = OnsetSpike([1; 1+find(diff(OnsetSpike) > Fs * dT_Exclude)]);
		end; 
		
		%%%%%%% Burst Detection %%%%%%%%%%%%%%%%%%%
		OnsetBurst = OnsetSpike ( [1; 1 + find( diff(OnsetSpike) > Fs * dT_Burst ) ] );

		DUR        = repmat(NaN, length(OnsetBurst), 1);
		BurstTable = repmat(NaN, length(OnsetBurst), 6);
		OnsetBurst(end+1) = inf;

		m2    = 0;
		t0    = [1; EVENT.POS(EVENT.TYP==hex2dec('7ffe'))];
		for m = 1:length(OnsetBurst)-1,	% loop for each burst candidate 
			tmp = OnsetSpike( OnsetBurst(m) <= OnsetSpike & OnsetSpike < OnsetBurst(m+1) );
			d   = diff(tmp);

			ix = sum(t0 <= OnsetBurst(m));  %% number of sweep
			T0 = t0(ix);	% time since beginning of sweep 
			if length(tmp) > 1,
				m2 = m2 + 1;
				DUR(m) = length(tmp)*mean(d);
				BurstTable(m,:) = [ch, ix, (OnsetBurst(m) - T0)/HDR.SampleRate, length(tmp), 1000*mean(d)/HDR.SampleRate, 1000*min(d)/HDR.SampleRate];
			elseif length(tmp)==1
			% 	single spikes are not counted as bursts, DUR(m)==NaN marks them as invalid 
				BurstTable(m,:) = [ch, ix, (OnsetBurst(m) - T0)/HDR.SampleRate, length(tmp), NaN,NaN];
			end;
		end;

		% remove single spike bursts 
		HDR.BurstTable = [HDR.BurstTable; BurstTable];

		EVENT.TYP = [EVENT.TYP; repmat(hex2dec('0202'), size(DUR))];
		EVENT.POS = [EVENT.POS; OnsetBurst(1:end-1)];
		EVENT.DUR = [EVENT.DUR; DUR];
		EVENT.CHN = [EVENT.CHN; repmat(ch, size(DUR)) ];
		if isfield(EVENT,'TimeStamp')
			%%% TODO: these should be properly computed %%%
			EVENT.TimeStamp = [EVENT.TimeStamp; repmat(NaN, size(DUR)) ];
		end;
	end; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%	Output 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	HDR.EVENT = EVENT;
	H = HDR;	%% return value
	if ~isempty(burstFile),
		fid = fopen(burstFile, 'w');
		if (fid<0) 
			fprintf(2,'Warning can not open file <%s> - burst table not written\n',burstFile); 
		else 
			fprintf(fid,'#Burst\t#chan\t#sweep\tOnset [s] \t#spikes\tavgISI [ms] \tminISI [ms] \n');
			for m = 1:size(HDR.BurstTable,1);
				fprintf(fid,'%3d\t%d\t%d\t%9.5f\t%d\t%6.2f\t%6.2f\n', m, HDR.BurstTable(m,:) );
			end;
			if (fid>3) fclose(fid); end;
		end;
	end;

	if ~isempty(outFile) && exist('s','var')
		%%% write data to output
		HDR.TYPE    = 'GDF';
		HDR.VERSION = 3;
		HDR.FILE    = [];
		HDR.FileName  = outFile;
		HDR.FILE.Path = '';
		HDR.Dur     = 1/HDR.SampleRate;
		HDR         = sopen(HDR, 'w');
		if (HDR.FILE.FID < 0) 
			fprintf(2,'Warning can not open file <%s> - GDF file can not be written\n',HDR.FileName); 
		else
			HDR = swrite(HDR,s);
			HDR = sclose(HDR);
 		end;
	end;
	WarnState = warning;
	warning('off'); 
	if ~isempty(evtFile)
		%%% write data to output
		HDR.TYPE  = 'EVENT';
		HDR.VERSION = 3;
		%[p,f,e]=fileparts(fn);
		HDR.FILE = [];
		HDR.FileName  = evtFile;
		HDR.FILE.Path = '';
		HDR.NRec = 0;
		HDR.SPR = 0;
		HDR.Dur = 1/HDR.SampleRate;
		HDR = rmfield(HDR,'AS');
		HDR = sopen(HDR, 'w');
		if (HDR.FILE.FID<0) 
			fprintf(2,'Warning can not open file <%s> - EVT file can not be written\n',HDR.FileName); 
		else
			HDR = sclose(HDR);
 		end;
	end;
	warning(WarnState); 

