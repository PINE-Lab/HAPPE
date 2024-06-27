function EVT = add_events(EVT, TYP, POS, CHN, DUR)
% ADD_EVENTS adds events to event table 
%
% Usage:
%   EVT = add_events(EVT1, EVT2)
%     merge eventtables EVT1 and EVT2
%   EVT = add_events(EVT, TYP, POS)
%   EVT = add_events(EVT, TYP, POS, CHN, DUR)
%     add events described by TYP and POS to event table
%
%   Event tables with TimeStamp information is currently not supported

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

if (nargin==2) && isstruct(TYP)
	EVT2=TYP;

elseif (nargin==3) && isnumeric(TYP) && isnumeric(POS)
	if (numel(TYP)==numel(POS))
		EVT2.TYP=TYP(:);
		EVT2.POS=POS(:);
	elseif numel(TYP)==1
		EVT2.TYP=repmat(TYP,numel(POS),1);
		EVT2.POS=POS(:);
	else
		error('input arguments not supported');
	end		

elseif (nargin==5) && isnumeric(TYP) && isnumeric(POS) && isnumeric(CHN) && isnumeric(DUR)
	assert(numel(TYP)==numel(POS));
	assert(numel(DUR)==numel(POS) || isempty(DUR));
	assert(numel(CHN)==numel(POS) || isempty(CHN));
	EVT2.TYP=TYP(:);
	EVT2.POS=POS(:);
	EVT2.DUR=DUR(:);
	EVT2.CHN=CHN(:);
else
	error('input arguments not supported')
end

if isfield(EVT,'TimeStamp') || isfield(EVT2,'TimeStamp')
	error('TimeStamps are not supported')
end;

%%%%%% CHECK INPUT %%%%%%
if isfield(EVT,'EVENT')
	ET = EVT.EVENT;
else
	ET = EVT;
end;


if ~isempty(ET.CHN) || ~isempty(ET.DUR) || ~isempty(EVT2.CHN) || ~isempty(EVT2.DUR)
	if isempty(ET.CHN)
		ET.CHN=zeros(size(ET.POS));
	end
	if isempty(EVT2.CHN)
		EVT2.CHN=zeros(size(EVT2.POS));
	end
	if isempty(ET.DUR)
		ET.DUR=zeros(size(ET.POS));
	end
	if isempty(EVT2.DUR)
		EVT2.DUR=zeros(size(EVT2.POS));
	end
	ET.CHN=[ET.CHN; EVT2.CHN];
	ET.DUR=[ET.DUR; EVT2.DUR];
end;
ET.POS=[ET.POS; EVT2.POS];
ET.TYP=[ET.TYP; EVT2.TYP];

