function [H2,HDR,s] = qrsdetect(fn,arg2,arg3,varargin)
% QRSDETECT - detection of QRS-complexes
%
%   HDR = qrsdetect(fn,chan,Mode)
%   ... = qrsdetect(fn, 0, Mode, '-o',outputFilename)
%   ... = qrsdetect(... ,'-e',eventFilename)
%   HDR = qrsdetect(s,Fs,Mode) 
%
% INPUT
%   	fn	filename
%   	chan    channel number of ecg data
%		if chan is empty, all channels which contain ECG, ecg, or EKG in HDR.Label 
%		are used. 
%   	s       ecg signal data 
%   	Fs      sample rate 
%   	Mode    optional - default is 2
%               1: method [1] is used
%		2: method [2] is used
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
%   HDR.EVENT  fiducial points of qrs complexes	
%
%
% see also: PROCESSING, EVENTCODES.TXT, SLOAD 
%
% Reference(s):
% [1] M.-E. Nygards, L. Sörnmo, Delineation of the QRS complex using the envelope of the e.c.g
%       Med. & Biol. Eng. & Comput., 1983, 21, 538-547.
% [2] V. Afonso, W. Tompkins, T. Nguyen, and S. Luo, "ECG beat detection using filter banks."
% 	IEEE Trans. Biomed. Eng. 46(2):192-202, Feb. 1999.

% $Id$
% Copyright (C) 2000-2003,2006,2009,2011 by Alois Schloegl <alois.schloegl@gmail.com>
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


if isnumeric(fn),
	s = fn;
	HDR.SampleRate = arg2;
	HDR.NS = size(s,2);
	chan = 1:size(s,2);
	CHAN = 1:size(s,2);
else
	if ~isempty(outFile) 
		[s,HDR] = sload(fn);
		if chan==0 || isempty(chan);
			chan = unique([strmatch('ecg',HDR.Label),strmatch('ECG',HDR.Label),strmatch('EKG',HDR.Label)]);
		end;	
		CHAN = chan;
	else 
		HDR = sopen(fn,'r',chan);HDR = sclose(HDR);
		if chan==0 || isempty(chan);
			chan = unique([strmatch('ecg',HDR.Label),strmatch('ECG',HDR.Label),strmatch('EKG',HDR.Label)]);
		end;	
		[s,HDR] = sload(fn,chan);
		CHAN = chan;
		chan = 1:length(CHAN);
	end; 
end;

ET = [];
for ch = 1:length(chan)',
	k = chan(ch);

	if 0,

        elseif MODE==2,     % QRS detection based on Afonso et al. 1999
                POS = nqrsdetect(s(:,k),HDR.SampleRate);

        elseif MODE==1,     % QRS detection based on Hilbert transformation. For details see [1]
                y   = processing({'ECG_envelope',HDR.SampleRate},s(:,k));
                TH  = quantile(y,.90);

                POS = gettrigger(y,TH);	% fiducial point

        else %if Mode== ???     % other QRS detection algorithms
		fprintf(2,'Error QRSDETECT: Mode=%i not supported',Mode);
		return;
        end;

        % difference between R-peak and fiducial point
        [t,sz] = trigg(s(:,k),POS,floor(-HDR.SampleRate),ceil(HDR.SampleRate));
        [tmp,ix] = max(abs(mean(reshape(t,sz(2:3)),2)));
        delay = HDR.SampleRate + 1 - ix;

        ET = [ET; [POS-delay, repmat([hex2dec('0501'),CHAN(ch),0], size(POS,1),1)]];
end;


[tmp,ix] = sort(ET(:,1));
if isfield(HDR,'T0')
	H2.T0 = HDR.T0;
end;
if isfield(HDR,'Patient')
	H2.Patient = HDR.Patient;
end;
H2.EVENT.POS = ET(ix,1);
H2.EVENT.TYP = ET(ix,2);
H2.EVENT.CHN = ET(ix,3);
H2.EVENT.DUR = ET(ix,4);
H2.EVENT.SampleRate = HDR.SampleRate;
H2.TYPE = 'EVENT';

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


