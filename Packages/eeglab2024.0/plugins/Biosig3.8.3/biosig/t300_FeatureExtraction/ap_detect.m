function [H2,HDR,s] = ap_detect(fn,arg2,arg3)
% AP_DETECT - detection of action potentionsl 
%
%   HDR = ap_detect(fn,chan,Mode)
%   HDR = ap_detect(s,Fs,Mode)
%
% INPUT
%   fn	        filename
%   chan        channel number of ecg data
%   s           ap signal data 
%   Fs          sample rate 
%   Mode        optional - default is 1
%               1: exceed slope of 20 V/s
%		2: 
% OUTPUT
%   HDR.EVENT  fiducial points of the action potential	
%
%
% see also: QRSDETECT, EVENTCODES.TXT, SLOAD, 
%
% Reference(s):
%
%

% $Id$
% Copyright (C) 2011 by Alois Schloegl <alois.schloegl@gmail.com>
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
MODE = 1;       % default: 
if (nargin==2) 
	if isnumeric(arg2)
		chan = arg2;
	else
		MODE = arg3; 
	end;
elseif (nargin==3) 
	chan = arg2; 
	MODE = arg3;
end;

if isnumeric(fn),
	s = fn; 
	HDR.SampleRate = arg2;
	chan = 1:size(s,2); 
else
	[s,HDR] = sload(fn,chan);
	chan = find(bitand(HDR.PhysDimCode,hex2dec('ffe0'))==4256); 	% select Voltage channels; 4256==physicalunits('V')
end;	

dT = 0.0005;	% window length
B = [1,zeros(1, HDR.SampleRate * dT - 1), -1];

ET = [];
for k = 1:length(chan),
	if 0,

        elseif MODE==1,     % QRS detection based on Hilbert transformation. For details see [1]
		[pu,scale] = physicalunits(HDR.PhysDimCode(k));
		TH = dT*20/(scale);
		POS = find(diff(filter(B,1,s(:,chan(k))) > TH) > 0);	% fiducial point
		

        else %if Mode== ???     % other QRS detection algorithms
		fprintf(2,'Error AP_DETECT: Mode=%i not supported',Mode);
		return; 
        end;

        ET = [ET; [POS, repmat([hex2dec('0010'),chan(k),0], size(POS,1),1)]];
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

