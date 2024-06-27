function [ISI] = get_inter_spike_intervals(F)
% GET_INTER_SPIKE_INTERVALS compute series of interspike 
%  intervals from the eventtable of  a biosignal data set 
%  loaded with SOPEN or SLOAD. 
%
% Usage:
%    [ISI] = get_inter_spike_intervals(F)
% 
% Input:  
% 	F filename
% Output:  
%	ISI series (vector) of interspike intervals in seconds  
%
% See also: SOPEN, SLOAD, DETECT_SPIKES_BURSTS, SPIKE2BURSTS, OPTIMUM_ISI_SPIKE_BURST_SEPARATION

%    Copyright (C) 2014 by Alois Schloegl <alois.schloegl@ist.ac.at>
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



if ischar(F) && exist(F,'file');
	HDR = sopen(F); 
	HDR = sclose(HDR); 

elseif isstruct(HDR) && isfield(HDR,'EVENT') && isfield(HDR.EVENT,'POS') && isfield(HDR.EVENT,'TYP')
	; %% everything fine - continue
else 
	error('invalid input arguments'); 
end; 

%% select only spikes and segment breaks. 
ix = ( ( HDR.EVENT.TYP == hex2dec('0201') ) | ( HDR.EVENT.TYP == hex2dec('7ffe') ) );
EVENT.TYP = HDR.EVENT.TYP(ix);
EVENT.POS = HDR.EVENT.POS(ix);

[EVENT.POS, ixx] = sort(EVENT.POS);
EVENT.TYP = EVENT.TYP(ixx);

%% select spike posititions, segment breaks remain NaN
EVENT.POS(EVENT.TYP == hex2dec('7ffe')) = NaN;

%% compute interspike intervals, segment breaks result in NaN. 
d = diff(EVENT.POS/HDR.EVENT.SampleRate);

%% remove meaningless intervals caused by segment breaks, and normalize to seconds. 
ISI = d(~isnan(d))/HDR.SampleRate;

