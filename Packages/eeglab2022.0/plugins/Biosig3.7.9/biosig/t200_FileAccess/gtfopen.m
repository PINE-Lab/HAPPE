function [HDR]=gtfopen(HDR,PERMISSION,arg3,arg4,arg5,arg6)
% GTFOPEN reads files in the Galileo Transfer format (GTF) 
%
% HDR = gtfopen(Filename,PERMISSION);
%
% HDR contains the Headerinformation and internal data
%
% see also: SOPEN, SREAD, SSEEK, STELL, SCLOSE, SWRITE, SEOF


% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
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

% Copyright (C) 2005,2020 by Alois Schloegl <alois.schloegl@gmail.com>
%    This is part of the BIOSIG-toolbox https://biosig.sourceforge.io/

        HDR.FILE.FID = fopen(HDR.FileName,'rb','ieee-le');
        
        % read 3 header blocks 
        HDR.H1 = fread(HDR.FILE.FID,[1,512],'uint8');
        HDR.H2 = fread(HDR.FILE.FID,[1,15306],'int8');
        HDR.H3 = fread(HDR.FILE.FID,[1,8146],'uint8');           
        
        HDR.L1 = char(reshape(HDR.H3(1:650),65,10)');                   
        HDR.L2 = char(reshape(HDR.H3(650+(1:20*16)),16,20)');
        HDR.L3 = reshape(HDR.H3(1070+32*3+(1:232*20)),232,20)';
        
        HDR.Label = char(reshape(HDR.H3(1071:1070+32*3),3,32)');        % channel labels
        
        [H.i8,count] = fread(HDR.FILE.FID,inf,'int8');
        fclose(HDR.FILE.FID);
        

	[t,status] = biosig_str2double(char([HDR.H1(35:36),32,HDR.H1(37:39)]));
        if ~any(status) & all(t>0)
                HDR.NS = t(1); 
                HDR.SampleRate = t(2); 
        else
                fprintf(2,'ERROR SOPEN (%s): Invalid GTF header.\n',HDR.FileName);
		HDR.TYPE = 'unknown';
                return; 
        end
        HDR.Dur  = 10; 
        HDR.SPR  = HDR.Dur*HDR.SampleRate; 
        HDR.bits = 8; 
	HDR.GDFTYP = repmat(1,HDR.NS,1);
        HDR.TYPE = 'native'; 
        HDR.THRESHOLD = repmat([-127,127],HDR.NS,1);    % support of overflow detection
        HDR.FILE.POS = 0; 
        HDR.Label = HDR.Label(1:HDR.NS,:);
        
        HDR.AS.bpb = (HDR.SampleRate*240+2048);
	HDR.GTF.Preset = HDR.H3(8134);	% Preset

        t2 = 9248+(0:floor(count/HDR.AS.bpb)-1)*HDR.AS.bpb;
        HDR.NRec = length(t2);
        [s2,sz]  = trigg(H.i8,t2,1,HDR.SampleRate*240);
        HDR.data = reshape(s2,[HDR.NS,sz(2)/HDR.NS*HDR.NRec])';
        
        [s4,sz]  = trigg(H.i8,t2-85,0,1);
        sz(sz==1)= []; 
        x  = reshape(s4,sz)';   
        HDR.GTF.timestamp = (x+(x<0)*256)*[1;256];      % convert from 2*int8 in 1*uint16
        
        HDR.GTF.i8 = H.i8; 
        HDR.GTF.t2 = t2; 

	t2 = HDR.GTF.t2; 
	[s4,sz] = trigg(HDR.GTF.i8,t2,-2047,0);
	sz(sz==1)=[]; if length(sz)<2,sz = [sz,1]; end;
	s4 = reshape(s4,sz);

	tau  = [0.01, 0.03, 0.1, 0.3, 1];
	LowPass = [30, 70];

	%% Scaling 
	Sens = [.5, .7, 1, 1.4, 2, 5, 7, 10, 14, 20, 50, 70, 100, 140, 200]; 
	x    = reshape(s4(13:6:1932,:),32,HDR.NRec*HDR.Dur);
	Cal  = Sens(x(1:HDR.NS,:)+1)'/4;
	HDR.data  = HDR.data.*Cal(ceil((1:HDR.SampleRate*HDR.NRec*HDR.Dur)/HDR.SampleRate),:);
        HDR.Calib = sparse(2:HDR.NS+1,1:HDR.NS,1);
	HDR.PhysDim = 'uV';


