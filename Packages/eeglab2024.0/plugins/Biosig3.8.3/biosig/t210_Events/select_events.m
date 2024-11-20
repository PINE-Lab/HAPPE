function EVT = select_events(EVT, TYPLIST, CHNLIST)
% SELECT_EVENT select specific types. 
%
% Usage:
%   EVT = select_events(EVT, TYPLIST, CHNLIST)
%     select all events where TYP matches ANY in TYPLIST OR CHN matches any CHNLIST element
%     unselects all events where TYP is not in TYPLIST AND CHN is not in CHNLIST
%     if you want events that where TYP is in TYPLIST AND CHN is in CHNLIST you need to do
%
%   EVT = select_events(select_event(EVT, TYPLIST, []), [], CHNLIST)
%     select events that match TYPLIST AND CHNLIST
%
%   EVT = unselect_events(EVT, TYPLIST, CHNLIST)
%     unselect all events where TYP matches ANY in TYPLIST and CHN matches any CHNLIST element
%     selects all events where either TYP is not contained in TYPLIST or CHN is not contained in CHNLIST
%
%   EVT = unselect_events(unselect_event(EVT, TYPLIST, []), [], CHNLIST)
%     select events that are neither in TYPLIST nor in CHNLIST 
%

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


%%%%%% CHECK INPUT %%%%%%
if isfield(EVT,'EVENT')
	ET = EVT.EVENT;
else
	ET = EVT;
end;


%%%%%% PROCESSING %%%%%%
flag = repmat(logical(0),size(ET.TYP));

for typ=TYPLIST(:)',
	flag = flag | (ET.TYP==typ);
end;

if isfield(ET,'CHN')
	for chn=CHNLIST(:)',
		flag = flag | (ET.CHN==chn);
	end;
end;

ET.TYP(~flag)=[];
ET.POS(~flag)=[];
if isfield(ET,'CHN')
	ET.CHN(~flag)=[];
end;
if isfield(ET,'DUR')
	ET.DUR(~flag)=[];
end;
if isfield(ET,'TimeStamp')
	ET.TimeStamp(~flag)=[];
end;

%%%%%% OUTPUT %%%%%%
if isfield(EVT,'EVENT')
	EVT.EVENT=ET;
else
	EVT = ET;
end;

