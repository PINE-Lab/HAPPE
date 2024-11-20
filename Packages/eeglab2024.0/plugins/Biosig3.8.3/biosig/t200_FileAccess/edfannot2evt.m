function HDR = edfannot2evt(HDR)
% EDFANNOT2EVT converts the EDF+ annotation channel into an event table
%
%  
% see also: SLOAD, SOPEN

%   Copyright (C) 2012-2017,2020 by Alois Schloegl <alois.schloegl@gmail.com>
%    This is part of the BIOSIG-toolbox https://biosig.sourceforge.io/
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

if isfield(HDR,'EDFplus') && isfield(HDR.EDFplus,'ANNONS'),

	warning('EDFANNOT2EVT.M is currently not supported! If there are problems reading EDF files,');
	warning('download and install mexSLOAD.mex* from https://pub.ist.ac.at/~schloegl/src/mexbiosig/')

	sz = size(HDR.EDFplus.ANNONS);
	N  = 0;
	Desc = {};
	TYP  = [];
	POS  = [];
	DUR  = [];
	CHN  = [];
	TimeStamp = [];
	HDR_T0 = NaN;
	try
		HDR_T0 = datenum(HDR.T0);
	end

	for k = 1:sz(2);
		t     = HDR.EDFplus.ANNONS(:,k)';

		tt0   = 0;
		while any(t)
			[t1,t] = strtok(t, 0);

			ix = find(t1==20);
			if isempty(ix)
				t0 = biosig_str2double(t1);
				N        = N + 1;
				tt0      = t0(1);
				POS(N,1) = HDR.SPR * (k-1) + 1;
				DUR(N,1) = 0;
				TYP(N,1) = hex2dec('7ffe');
				TimeStamp(N,1) = HDR_T0 + t0/(24*60);
				continue;
			end
			s1 = t1(1:ix(1)-1);
			s2 = t1(ix(1)+1:ix(2)-1);
			s1(s1==21)=0;
			t0 = biosig_str2double(s1);
			if ((length(ix)>=2) && ((ix(1)+1)==ix(2)) && strcmp(HDR.reserved1(2:5),'DF+D') ) 
				%% time keeping annotation 
				N        = N + 1;
				tt0      = t0(1);
				POS(N,1) = HDR.SPR * (k-1) + 1;
				DUR(N,1) = 0;
				TYP(N,1) = hex2dec('7ffe'); 
				TimeStamp(N,1) = HDR_T0 + t0/(24*60);
			elseif (length(ix)==2 && ix(1)+1 < ix(2))
				N        = N + 1;
				TYP(N,1) = 1; 
				%if all(s2(2:2:end)==0) s2 = s2(1:2:end); end; %% unicode to ascii - FIXME 
				Desc{N}  = s2;
				if length(t0)>1
					DUR(N,1)  = t0(2) * HDR.EVENT.SampleRate;
				else
					DUR(N,1)  = 0;
				end;
				POS(N,1)  = 1 + round((t0(1)-tt0) * HDR.EVENT.SampleRate);
				TimeStamp(N,1) = HDR_T0 + (t0(1)-tt0)/(24*60*60);
			end;
		end; 
	end; 

	HDR.EVENT.POS = POS(:);
	HDR.EVENT.DUR = DUR(:);
	HDR.EVENT.TYP = TYP(:);
	HDR.EVENT.CHN = zeros(N,1);
	HDR.EVENT.TimeStamp = TimeStamp;

	ix = find(TYP < 256);
	if any(ix),
		[HDR.EVENT.CodeDesc, HDR.EVENT.CodeIndex, HDR.EVENT.TYP(ix)] = unique(Desc(ix)');
	end;

        %% TODO: use eventcodes.txt for predefined event types e.g. QRS->0x501

end

