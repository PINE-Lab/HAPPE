% 
% T250: Quality Control and Artifact Processing
% ---------------------------------------------
%	ARTIFACT_SELECTION	converts artifact scorings into trial selections
% 	EEG2HIST	calculates histogram
% 	GETTRIGGER	gets trigger points
% 	TRIGG		extract fixed-length trials around trigger points	
%	DETECT_MUSCLE	detection of muscle artefacts using an inverse filter
%	DETECT_SPIKES_BURSTS detect spikes and bursts in neural recordings
%       REGRESS_EOG     reduce EOG artifacts with regression analysis 
%       REMOVE5060HZ    methods for removing line interference
% 

% 
% Copyright (C) 2003,2004,2005,2006,2011 by Alois Schloegl <alois.schloegl@gmail.com>
% $Id$
% This is part of the BIOSIG-toolbox http://biosig.sf.net/
%
% LICENSE:
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program; if not, write to the Free Software
%     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA 
% 
% 
