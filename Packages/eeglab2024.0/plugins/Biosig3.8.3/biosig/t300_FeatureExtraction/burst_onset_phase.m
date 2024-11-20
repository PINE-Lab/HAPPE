function [RES, HDR] = burst_onset_phase(fn, ch, band, EVT, varargin)
% BURST_ONSET_PHASE computes the phase of burst onsets 
%
% ... = burst_onset_phase(fn, ch, [f1, f2], EVT)
% ... = burst_onset_phase(..., '-o', outputFilename)
% ... = burst_onset_phase(s, HDR, [f1, f2], EVT)
% [RES, HDR] = burst_onset_phase(...)
% 
% Input: 
%  fn	filename 
%  ch	channel number(s), default=0 (i.e. all)	
%  s	signal data 
%  HDR  header structure
%  [f1,f2]   edge frequencies of bandpass filter
%  EVT  header structure containing events of 
%       burst onset (0x0202) and 
%       start of new segments (0x7ffe)	
%
% Output:
%  RES  a complex matrix containing normalized Amplitude and Phase 
%	for each burst onset (number of rows) and for 
%	each channel (number of columns). The amplitude is 
% 	normalized with the rms of the filtered signal. 
%  angle(RES) returns the phase at burst onset

%	$Id$
%	Copyright (C) 2011 by Alois Schloegl <alois.schloegl@gmail.com>
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

if ~isnumeric(band) || numel(band)~=2 
	error('arg3 (band) must have 2 numeric elements');
else
	f1 = band(1); 
	f2 = band(2); 
end;


if ischar(EVT) && exist(EVT,'file')
	[s,EVT] = sload(EVT);
elseif isstruct(EVT) && isfield(EVT,'POS') && isfield(EVT,'TYP')
	tmp = EVT; 
	EVT = []; 
	EVT.EVENT=tmp;
end; 
if ~isfield(EVT,'EVENT')
	error('arg4 (EVT) does not contain event structure');
end; 

outFile = [];
if length(varargin)>1 && strcmp(varargin{1},'-o'),
	outFile = varargin{2};
end; 


winlen = 20000;

[s, HDR]    = sload(fn, ch, 'NUMBER_OF_NAN_IN_BREAK', winlen);
ixSegStart0 = sort([1; HDR.EVENT.POS(HDR.EVENT.TYP==hex2dec('7ffe')); size(s,1)+1]);
ixSegStart1 = sort([1; EVT.EVENT.POS(EVT.EVENT.TYP==hex2dec('7ffe')); size(s,1)+1]);

if length(ixSegStart0)-length(ixSegStart1)
	error('number of segments differ between raw data and event file');
	return; 
end

if ~isempty(outFile)
	%% keep copy of s for outFile 
	s0 = s; 
end

%%% handle missing values
iv = isnan(s);
s(iv) = 0;


%%% transform to spectral domain, filter, and transform 
s = fft(s,[],1); 
f = (0 : size(s,1) - 1) * HDR.SampleRate / size(s,1);
s( f < f1 | f2 < f, :) = 0; 
if f1==0, s(f==0, : ) = s(f==0, : ) / 2; end; 
s = 2 * ifft(s,[],1); 


s(iv) = NaN; 


%%% rearrange marker to account for NUMBER_OF_NAN_IN_BREAK
tixOn = EVT.EVENT.POS; %(EVT.EVENT.TYP==hex2dec('0202'));
ix    = []; 
for k = 1:length(ixSegStart0)-1,
	ix0 = ixSegStart0(k); 
	ix1 = ixSegStart1(k);
	t   = tixOn( ixSegStart1(k) <= tixOn & tixOn < ixSegStart1(k+1) );
	ix  = [ix; ix0 - ix1 + t(:)];
end; 
EVT.EVENT.POS = ix; 


%%% hack: workaround to nonsense markers 
if any( ix > size(s,1) )
	fprintf(2, '\nWarning: events after end of data found\n');
	for k = find( ix > size(s,1) )
		fprintf(2, '\t%d\t%x\n', EVT.EVENT.POS(k), EVT.EVENT.TYP(k));
	end; 
	ix = ix( ix <= size(s,1) );
end; 


%%% return output: extract phase at burst onset time 
tix = EVT.EVENT.POS( EVT.EVENT.TYP==hex2dec('0202') );
RES = s( tix( tix < size(s,1) ), :) * diag(1/rms(s,1));



if ~isempty(outFile)
	HDR2.EVENT = EVT.EVENT; 
	HDR2.TYPE  = 'GDF'; 
	HDR2.FileName = outFile; 
	HDR2.FILE  = [];
	HDR2.Patient = HDR.Patient; 
	HDR2.T0    = HDR.T0; 
	HDR2.NS    = length(ch)*2; 
	HDR2.SPR   = 1; 
	HDR2.NRec  = size(s0,1);
	HDR2.SampleRate = HDR.SampleRate; 
	HDR2.Label = {HDR.Label{ch}; HDR.Label{ch}};	
	HDR2.DigMax  = [HDR.DigMax(ch);  HDR.DigMax(ch)];
	HDR2.DigMin  = [HDR.DigMin(ch);  HDR.DigMin(ch)];
	HDR2.PhysMax = [HDR.PhysMax(ch); HDR.PhysMax(ch)];
	HDR2.PhysMin = [HDR.PhysMin(ch); HDR.PhysMin(ch)];
	HDR2.GDFTYP  = repmat(16, 1, HDR.NS);
	HDR2.PhysDimCode = [HDR.PhysDimCode(ch);HDR.PhysDimCode(ch)];
	HDR2.Filter.LowPass  = [HDR.Filter.LowPass(ch);	 repmat(f2, length(ch), 1)]; 
	HDR2.Filter.HighPass = [HDR.Filter.HighPass(ch); repmat(f1, length(ch), 1)]; 
	HDR2.Filter.Notch    = [HDR.Filter.Notch(ch);    HDR.Filter.Notch(ch)]; 
	HDR2.FLAG.UCAL = 0; 

	HDR2 = sopen(HDR2,'w');
	HDR2 = swrite(HDR2, [s0, s]);
	HDR2 = sclose(HDR2); 
end; 
