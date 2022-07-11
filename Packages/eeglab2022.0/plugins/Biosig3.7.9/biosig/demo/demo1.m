% DEMO 1 - identifies QRS-complexes and computes HRV parameters 

%	$Id$
%	Copyright (C) 2000-2003,2005,2009 by Alois Schloegl <alois.schloegl@gmail.com>
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

% This library is free software; you can redistribute it and/or
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


% load file

[F,P]=uigetfile('*.*','Pick an ECG file');

if (F==0)
	for ext = {'atr','dat','hea'}
		f = ['100.',ext{1}]
		url=fullfile('https://www.physionet.org/physiobank/database/mitdb/',f);
		fprintf(1,'Example ECG data is automatically downloaded from PhysioBank: http://www.physionet.org/physiobank/database/#ecg\n');
		[S, SUCCESS] = urlread(url);
		if SUCCESS,
			fid=fopen(f,'w');
			fwrite(fid,S,'char');
			fclose(fid);
		else
			error('download of example file failed');
		end;
	end;
	P='';
	F='100.hea';
end; 	

HDR  = sopen(fullfile(P,F),'r');
if HDR.NS > 1,
        CHAN = sort([strmatch('ECG',HDR.Label);strmatch('EKG',HDR.Label);strmatch('ecg',HDR.Label);strmatch('Ecg',HDR.Label)]);
        if length(CHAN)~=1,
                HDR = sclose(HDR);
                fprintf(1,'The selected file contains the following channels: \n');
                for k = 1:HDR.NS,
                        fprintf(1,'%3i: %s\n',k,HDR.Label{k});
                end;
                CHAN = input('Which channel should be used for QRS-detection? ');
        end;
        HDR = sopen(fullfile(P,F),'r',CHAN);
else 
        CHAN = 1;       
end;
[s,HDR] = sread(HDR);
HDR = sclose(HDR);


% QRS-Detection
try
   H2 = qrsdetect(s,HDR.SampleRate,2);
catch
   H2 = qrsdetect(s,HDR.SampleRate,1);
end;
% resampling to 4 Hz using the Berger algorithm 
[HRV,RRI] = berger(H2,4);
% compute HRV parameters 
[X] = heartratevariability(H2);
% Extract QRS-info according to BIOSIG/T200/EVENTCODES.TXT
idx = find(H2.EVENT.TYP == hex2dec('0501'));
qrsindex = H2.EVENT.POS(idx)/H2.EVENT.SampleRate; 

% displays detection
subplot(211)
plot((1:size(s,1))/HDR.SampleRate,s,'-',qrsindex,-ones(size(qrsindex)),'x');
xlabel('time t[s]');
ylabel(sprintf('%s [%s]',HDR.Label{CHAN},HDR.PhysDim{CHAN}));

subplot(212)
semilogy((qrsindex(1:end-1)+qrsindex(2:end))/2,diff(qrsindex));
ylabel('RRI [s]');
xlabel('time t[s]');

%% write output file 
H2.FileName = [HDR.FILE.Name,'.evt']; 
if exist(H2.FileName,'file')
        H2.FileName = input(['File ', H2.FileName, ' does already exist. What filename should be used to store detected QRS events ? '],'s');  
end;
fprintf(1,'File %s contains the detected QRS events\n',H2.FileName) 
H2.NRec = 0;
H2.T0 = HDR.T0; 
H2=sopen(H2,'w'); H2=sclose(H2); 

