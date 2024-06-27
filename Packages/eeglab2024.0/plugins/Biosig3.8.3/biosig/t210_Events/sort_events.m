function ET = sort_events(EVT)
% SORT_EVENTS event table according to startpostion of each event
%
% Usage:
%    HDR = sort_events(HDR)
%    HDR.EVENT = sort_events(HDR.EVENT)
%
%    Sorts the event table stored in HDR.EVENT according to
%    HDR.EVENT.POS. If several events occur at the same time,
%    the ordering is undefined.

%    Copyright (C) 2015 by Alois Schloegl <alois.schloegl@ist.ac.at>	
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


if nargin<1
	ET = [];
	return
end

%%%%%% CHECK INPUT %%%%%%
if isfield(EVT,'EVENT')
	ET = EVT.EVENT;
else
	ET = EVT;
end;

if ~isfield(ET,'POS') || ~isfield(ET,'TYP')
	error('Event table does not contain POS or TYP')
end; 
assert(numel(ET.TYP)==numel(ET.POS));


% sort event table before extracting HDR.Classlabel and HDR.TRIG
[ET.POS,ix] = sort(ET.POS);
ET.TYP = ET.TYP(ix);
if isfield(ET,'CHN')
	assert(numel(ET.CHN)==numel(ET.POS));
	ET.CHN = ET.CHN(ix);
end;    	    
if isfield(ET,'DUR')
	assert(numel(ET.DUR)==numel(ET.POS));
	ET.DUR = ET.DUR(ix);
end;
if isfield(ET,'TimeStamp')
	assert(numel(ET.TimeStamp)==numel(ET.POS));
	ET.TimeStamp = ET.TimeStamp(ix);
end; 	

% OUTPUT
if isfield(EVT,'EVENT')
	EVT.EVENT = ET;
	ET = EVT;
end;
