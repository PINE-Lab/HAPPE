function [H2,HDR,s] = ecg_wave_analysis(fn,arg2,arg3,varargin)
% ECG_WAVE_ANALYSIS - analysis basic ECG waves, like P-, QRS- and T-wave
%
% Currently, the function is incomplete, only a framework for the event
% handling is implemented, but no actuall detection algorithms.
%
%   HDR = ecg_analysis(fn,chan,Mode)
%   ... = ecg_analysis(fn, 0, Mode, '-o',outputFilename)
%   ... = ecg_analysis(... ,'-e',eventFilename)
%   HDR = qrsdetect(s,Fs,Mode) 
%
% INPUT
%   	fn	filename
%   	chan    channel number of ecg data
%		if chan is empty, all channels which contain ECG, ecg, or EKG in HDR.Label 
%		are used. 
%   	s       ecg signal data 
%   	Fs      sample rate 
%	outputFilename
%		name of file for storing the resulting data with the
%		detected spikes and bursts in GDF format.
%		
%	eventFilename
%		filename to store event inforamation in GDF format. this is similar to 
%		the outputFile, except that the signal data is not included and is, therefore,
%		much smaller than the outputFile
%
% OUTPUT
%   HDR.EVENT  event table containing beginning and end of P, QRS and T wave
%
%
% see also: SLOAD, QRSDETECT
%
% Reference(s):
% [1] The CSE Working Party. Recommendations for measurement standards in
%     quantitvative electrocardiography. European Heart Journal (1985) 6,
%     815-825
%
%


% Copyright (C) 2015 by Alois Schloegl <alois.schloegl@ist.ac.at>
% This is part of the BIOSIG-toolbox http://biosig.sf.net/

% BIOSIG is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% Version 3 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Library General Public License for more details.
%
% You should have received a copy of the GNU Library General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc., 59 Temple Place - Suite 330,
% Boston, MA  02111-1307, USA.


chan = 0; 
MODE = 2;       % default
if (nargin==2) 
	if isnumeric(arg2)
		chan = arg2;
	else
		MODE = arg3;
	end;
elseif (nargin>=3) 
	chan = arg2;
	MODE = arg3;
end;

%%%%% default settings %%%%%
outFile = [];
evtFile = [];

%%%%% analyze input arguments %%%%%
k = 1;
while k <= length(varargin)
	if ischar(varargin{k})
		if (strcmp(varargin{k},'-o'))
			k = k + 1;
			outFile = varargin{k};
		elseif (strcmp(varargin{k},'-e'))
			k = k + 1;
			evtFile = varargin{k};
		else
			warning(sprintf('unknown input argument <%s>- ignored',varargin{k}));
		end;
	end;
	k = k+1;
end;


[s,HDR]=sload(fn,chan);

fiducialpoint = HDR.EVENT.POS(HDR.EVENT.TYP==hex2dec('0501'));
if isempty(fiducialpoint)
	H2 = qrsdetect(s,HDR.SampleRate);
	fiducialpoint = HDR.EVENT.POS(HDR.EVENT.TYP==hex2dec('0501'));
end; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO: implement detection of Onset and Offset of P, QRS, and T waves     %
%
% If you are interested: Contact http://biosig.sf.net
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,'ECG_ANALYSIS is currently not implemented. You can help implementing. As a starting point read [1].\n');


MDC_ECG_TIME_START_P   = [];
MDC_ECG_TIME_END_P     = [];
MDC_ECG_TIME_START_QRS = [];
MDC_ECG_TIME_END_QRS   = [];
MDC_ECG_TIME_START_T   = [];
MDC_ECG_TIME_END_T     = [];


EVT = HDR.EVENT;
EVT = add_event(EVT, hex2dec('0502'), MDC_ECG_TIME_START_P);
EVT = add_event(EVT, hex2dec('8502'), MDC_ECG_TIME_END_P);
EVT = add_event(EVT, hex2dec('0503'), MDC_ECG_TIME_START_QRS);
EVT = add_event(EVT, hex2dec('8505'), MDC_ECG_TIME_END_QRS);
EVT = add_event(EVT, hex2dec('0506'), MDC_ECG_TIME_START_T);
EVT = add_event(EVT, hex2dec('8506'), MDC_ECG_TIME_END_T);


HDR.EVENT=EVENT;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(outFile)
		%%% write data to output
		HDR.TYPE  = 'GDF';
		HDR.VERSION = 3;
		%[p,f,e]=fileparts(fn);
		HDR.FILE = [];
		HDR.FileName  = outFile;
		HDR.FILE.Path = '';
		HDR.PhysMax   = max(s);
		HDR.PhysMin   = min(s);
		HDR.DigMax(:) =  2^15-1;
		HDR.DigMin(:) = -2^15;
		HDR.GDFTYP(:) = 3;
		HDR.FLAG.UCAL = 0;
		HDR.NRec = size(s,1);
		HDR.SPR = 1;
		HDR.NS  = size(s,2);
		HDR.Dur = 1/HDR.SampleRate;
		HDR = rmfield(HDR,'AS');
		HDR.EVENT = H2.EVENT; 
		HDR = sopen(HDR,'w');
		if (HDR.FILE.FID < 0) 
			fprintf(2,'Warning can not open file <%s> - GDF file can not be written\n',HDR.FileName);
		else
			HDR = swrite(HDR,s);
			HDR = sclose(HDR);
		end; 
end;

if ~isempty(evtFile)
		H2 = HDR;
		%%% write data to output
		H2.TYPE  = 'EVENT';
		H2.VERSION = 3;
		%[p,f,e]=fileparts(fn);
		H2.FILE = [];
		H2.FileName  = evtFile;
		H2.FILE.Path = '';
		H2.NRec = 0;
		H2.SPR = 0;
		H2.Dur = 1/HDR.SampleRate;
		%H2 = rmfield(HDR,'AS');
		H2 = sopen(H2, 'w');
		if (H2.FILE.FID<0) 
			fprintf(2,'Warning can not open file <%s> - EVT file can not be written\n',H2.FileName); 
		else
			H2 = sclose(H2);
 		end;
end;


