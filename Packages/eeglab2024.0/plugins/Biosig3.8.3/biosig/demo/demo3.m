% Demostration for generating EDF/BDF/GDF-files
% DEMO3 is part of the biosig-toolbox
%    and it tests also Matlab/Octave for its correctness. 
% 

%	Copyright (C) 2000-2005,2006,2007,2008,2011,2013 by Alois Schloegl <alois.schloegl@gmail.com>
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

x = randn(10000,6)+(1:1e4)'*ones(1,6); % test data
x = (1:1e4)'*ones(1,6)/1000; % test data
x = reshape(mod(1:6e4,100),6,1e4)'; x(:,6)=NaN;

clear HDR;

VER   = version;
cname = computer;

% select file format 
HDR.TYPE='GDF';  
% HDR.TYPE='EDF';
%HDR.TYPE='BDF'; 
%HDR.TYPE='CFWB';
%HDR.TYPE='CNT';

% set Filename
HDR.FileName = ['TEST_',VER([1,3]),cname(1:3),'_e1.',HDR.TYPE];

% description of recording device 
HDR.Manufacturer.Name = 'BioSig';
HDR.Manufacturer.Model = 'demo3.m';
HDR.Manufacturer.Version = '2.84';
HDR.Manufacturer.SerialNumber = '00000000';

% recording identification, max 80 char.
HDR.RID = 'TestFile 001'; %StudyID/Investigation [consecutive number];
HDR.REC.Hospital   = 'BioSig Test Lab'; 
HDR.REC.Technician = get_current_username();

HDR.REC.Equipment  = 'biosig';
HDR.REC.IPaddr	   = [127,0,0,1];	% IP address of recording system 	
HDR.Patient.Name   = 'anonymous';  
HDR.Patient.Id     = 'P0000';	
HDR.Patient.Birthday = [1951 05 13 0 0 0];
HDR.Patient.Weight = 0; 	% undefined 
HDR.Patient.Height = 0; 	% undefined 
HDR.Patient.Sex    = 'f'; 	% 0: undefined,	1: male, 2: female 
HDR.Patient.Birthday = zeros(1,6); %    undefined 
HDR.Patient.Impairment.Heart = 0;  %	0: unknown 1: NO 2: YES 3: pacemaker 
HDR.Patient.Impairment.Visual = 0; %	0: unknown 1: NO 2: YES 3: corrected (with visual aid) 
HDR.Patient.Smoking = 0;           %	0: unknown 1: NO 2: YES 
HDR.Patient.AlcoholAbuse = 0; 	   %	0: unknown 1: NO 2: YES 
HDR.Patient.DrugAbuse = 0; 	   %	0: unknown 1: NO 2: YES 
HDR.Patient.Handedness = 0; 	   % 	unknown, 1:left, 2:right, 3: equal

% recording time [YYYY MM DD hh mm ss.ccc]
HDR.T0 = clock;	

% number of channels
HDR.NS = size(x,2);

% Duration of one block in seconds
HDR.SampleRate = 1000;
HDR.SPR = 10000;
HDR.Dur = HDR.SPR/HDR.SampleRate;

% Samples within 1 block
HDR.AS.SPR = [1000;100;200;100;20;0];	% samples per block; 0 indicates a channel with sparse sampling 
%HDR.AS.SampleRate = [1000;100;200;100;20;0];	% samplerate of each channel

% channel identification, max 80 char. per channel
HDR.Label={'chan 1  ';'chan 2  ';'chan 3  ';'chan 4  ';'chan 5  ';'NEQS    '};

% Transducer, mx 80 char per channel
HDR.Transducer = {'Ag-AgCl ';'Airflow ';'xyz     ';'        ';'        ';'Thermome'};

	% define datatypes (GDF only, see GDFDATATYPE.M for more details)
HDR.GDFTYP = 3*ones(1,HDR.NS);

% define scaling factors 
HDR.PhysMax = [100;100;100;100;100;100];
HDR.PhysMin = [0;0;0;0;0;0];
HDR.DigMax  = repmat(2^15-1,size(HDR.PhysMax));
HDR.DigMin  = repmat(1-2^15,size(HDR.PhysMax));
HDR.FLAG.UCAL = 1; 	% data x is already converted to internal (usually integer) values (no rescaling within swrite);
HDR.FLAG.UCAL = 0; 	% data x will be converted from physical to digital values within swrite. 
% define filter settings 
HDR.Filter.Lowpass = [0,0,0,NaN,NaN,NaN];
HDR.Filter.Highpass = [100,100,100,NaN,NaN,NaN];
HDR.Filter.Notch = [0,0,0,0,0,0];
% define sampling delay between channels  
HDR.TOffset = [0:5]*1e-6;


% define physical dimension
HDR.PhysDim = {'uV';'mV';'%';'Ohm';'-';'°C'};	%% must be encoded in unicode (UTF8)
HDR.Impedance = [5000,50000,NaN,NaN,NaN,NaN];         % electrode impedance (in Ohm) for voltage channels 
HDR.fZ = [NaN,NaN,NaN,400000,NaN,NaN];                % probe frequency in Hz for Impedance channel

t = [100:100:size(x,1)]';
%HDR.NRec = 100;
HDR.VERSION = 2.22;        % experimental  
HDR.EVENT.POS = t;
HDR.EVENT.TYP = t/100;


if 1, 
HDR.EVENT.CHN = repmat(0,size(t));
HDR.EVENT.DUR = repmat(1,size(t));
HDR.EVENT.VAL = repmat(NaN,size(t));
%% This defines the sparse samples in channel 6
ix = 6:5:60; 
HDR.EVENT.CHN(ix) = 6; 
HDR.EVENT.VAL(ix) = 37.3+ix; %round(100*rand(size(ix))); % HDR.EVENT.TYP(ix) becomes 0x7fff
ix = 8; 
%% The following sparse samples are not valid because channel 5 is not defined as sparse (see HDR.AS.SPR)
HDR.EVENT.CHN(ix) = 5; % not valid because #5 is not sparse sampleing
HDR.EVENT.VAL(ix) = 37.4;
HDR.EVENT.TYP(~isnan(HDR.EVENT.VAL)) = hex2dec('7fff');
end; 

if 0, %try,
	mexSSAVE(HDR,x);
else %catch
	HDR1 = sopen(HDR,'w');
	HDR1 = swrite(HDR1,x);
	HDR1 = sclose(HDR1);
end;

%
[s0,HDR0] = sload(HDR.FileName);	% test file 

HDR0=sopen(HDR0.FileName,'r');
[s0,HDR0]=sread(HDR0);
HDR0=sclose(HDR0); 

%plot(s0-x)


