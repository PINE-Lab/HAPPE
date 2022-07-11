% SLOPE_EVALUATION quantify the error on slope estimation caused by 
%         different noise sources, and sampling rates (see [3] for details). 
%
% evaluate slope analysis, and its influence of windowlength, noise level, sampling rate, 
%
% change sampling rate and noise level, estimated max slope as function of 

%
% References: 
%   discussion on issue 31: http://code.google.com/p/stimfit/issues/detail?id=31
%
% 
% REFERENCE(S):
% [1] Segundo J Guzman , Alois Schlögl and Christoph Schmidt-Hieber
%     Stimfit: quantifying electrophysiological datwith Python
%     Frontiers in Neuroinformatics. (submitted)


% Copyright (C) 2013 by Alois Schloegl <alois.schloegl@ist.ac.at>	
% This is part of the BIOSIG-toolbox http://biosig.sf.net/

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



clear

fn = '/fs3/group/jonasgrp/Jose/publication.Stimfit/AP_interneuron.npy';
fid = fopen(fn,'r')
s = fgetl(fn);
ftell(fid);
fseek(fid,length(s)+1,'bof');
d=fread(fid,inf,'double')/1000;	%% convert from mV to V
fclose(fid)


wintime=.05e-3;	%
H.SampleRate=1e9;
H.NS = 1; 
H.PhysDim={'V'};
t  = [1:length(d)]'/H.SampleRate;	% time axis

figure(10)
g=plot(t(1:1000:end),d(1:1000:end),'g')
set(g,'linewidth',3);


kB=1.3806488e-23; % J/K
T = 273+27;  % K 
B = H.SampleRate/2;

TH = 50; 
M = 100; % number of repetitions

if 1, 

%% List of Impedances used for noise generation
Z = [0,2e6,5e6,10e6,20e6,50e6,100e6];
%Z = [10e6,50e6,100e6];
%Z = [1,10,100]*1e6;

%% List of frequencies 
F = [300,200,150,100,70,50,30,20,15,10,7,5]*1e3;
F = [400,250,200,125,100,80,50,40,25,20,12.5,10,8,5]*1e3;
%F = [200,100,50,20,10,5]*1e3;

%% pre-allocate memory 
maxSlope0 = repmat(NaN,[M,length(F),1]);
maxSlope1 = repmat(NaN,[M,length(F),length(Z)]);
maxSlope2 = repmat(NaN,[M,length(F),length(Z)]);
t1 = repmat(NaN,[M,length(F),length(Z)]);
t2 = repmat(NaN,[M,length(F),length(Z)]);
A1 = repmat(NaN,[M,length(F),length(Z)]);
A2 = repmat(NaN,[M,length(F),length(Z)]);
dQ = 2*(2^-16);	      %% 2 V, 16 bit  

%% intialize random generator
rand('state',1234);
jit = ceil(rand(1,M)*2e5);

for m = 1:M;	% do M permutations 

N  = randn(size(d))*sqrt(4*kB*T*B);	% generate impedance noise
for idxF = 1:length(F)	% for all frequencies 


	fs = F(idxF);	% sampling rate
	div = ceil(H.SampleRate/F(idxF));	% divisor
	
	%% downsampling
	tt = rs(t(jit:end),div,1);
	dd = rs(d(jit:end),div,1);

	maxSlope0(m,idxF,1)=max(diff(dd))*fs; 	% max slope of noise-free data
	dydt0 = diff(dd)*H.SampleRate/div;
	tix0 = find(dydt0 > TH, 1, 'first');
	t0(idxF,1) = sum(tt(tix0+[0,1]))/2;	% detection time
	A0(idxF,1) = sum(dd(tix0+[0,1]))/2;	% Amplitude 

for iz = 1:length(Z)
	%% add impedance noise, downsampling with anti-aliasing lowpasss filter, and add quantization noise		
	s  = round(rs(d+N*sqrt(Z(iz)),div,1)/dQ)*dQ;

	if exist('mexSSAVE','file'),
		%% if biosig toolbox is available, save to gdf file 
		H2.FileName = sprintf('EvalSlopeData_%03i_%03ikHz_%03iMOhm.gdf',m,fs/1000,Z(iz)/1e6); 
		H2.TYPE='GDF';
		[H2.NRec,H2.NS] = size(s); 
		H2.SPR = 1; 
		H2.SampleRate = round(fs); 
		H2.GDFTYP = 3;
		H2.PhysDimCode = physicalunits('mV');
		H2.DigMax = 2^15;
		H2.DigMin = -2^15;
		H2.PhysMax = 1.0;
		H2.PhysMin = -1.0;
		H2.FLAG.UCAL = 	0;
		mexSSAVE(H2,s(:));	
		%H2 = sopen(H2,'w');
		%swrite(H2,int16(round(s(:)/dQ)));
		%sclose(H2);	
		%[s3,H3]=sload(H2.FileName);
		%[s4,H4]=mexSLOAD(H2.FileName);
	end; 

	%% compute maxSlope using single sample and fixedWindow method	
	dydt1 = diff(s)*fs;
	maxSlope1(m,idxF,iz) = max(dydt1);
	winlen = max(1,round(wintime*fs));
	dydt2 =  (s(winlen+1:end)-s(1:end-winlen))*fs/winlen;
	maxSlope2(m,idxF,iz) = max(dydt2);
	
	%% compute detection time and amplitude for both methods 
	tix1 = find(dydt1 > TH, 1, 'first');
	tix2 = find(dydt2 > TH, 1, 'first');	
	t1(m,idxF,iz) = sum(tt(tix1+[0,1]))/2;
	t2(m,idxF,iz) = sum(tt(tix2+[0,winlen]))/2;
	
	A1(m,idxF,iz) = sum(s(tix1+[0,1]))/2;
	A2(m,idxF,iz) = sum(s(tix2+[0,winlen]))/2;
end;
end; 
end; 

%maxSlope0
[se0,m0]=sem(maxSlope0,1);
[se1,m1]=sem(maxSlope1,1);
[se2,m2]=sem(maxSlope2,1);
m1=squeeze(m1);
m2=squeeze(m2);
se1=squeeze(se1)*sqrt(M);
se2=squeeze(se2)*sqrt(M);


q0 = quantile(maxSlope0,[.05,.5,.95],1);
q1 = quantile(maxSlope1,[.05,.5,.95],1);
q2 = quantile(maxSlope2,[.05,.5,.95],1);
c1m=squeeze(q1(2,:,:));
c1l=squeeze(q1(1,:,:))-c1m;
c1u=squeeze(q1(3,:,:))-c1m;
c2m=squeeze(q2(2,:,:));
c2l=squeeze(q2(1,:,:))-c2m;
c2u=squeeze(q2(3,:,:))-c2m;

RES.maxSlope.true              = max(diff(d))*H.SampleRate; 
RES.maxSlope.noisefree.c       = maxSlope0;
RES.maxSlope.noisefree.average = m0;
RES.maxSlope.noisefree.ci     = se0;
RES.maxSlope.noisefree.p95    = squeeze(q0(3,:,:));
RES.maxSlope.noisefree.p50    = squeeze(q0(2,:,:));
RES.maxSlope.noisefree.p05    = squeeze(q0(1,:,:));

RES.maxSlope.singlesample.average = m1;
RES.maxSlope.singlesample.ci     = se1;
RES.maxSlope.singlesample.p95    = squeeze(q1(3,:,:));
RES.maxSlope.singlesample.p50    = squeeze(q1(2,:,:));
RES.maxSlope.singlesample.p05    = squeeze(q1(1,:,:));
RES.amplitude.singlesample.p95   = squeeze(quantile(A1,.95,1));
RES.amplitude.singlesample.p50   = squeeze(quantile(A1,.5,1));
RES.amplitude.singlesample.p05   = squeeze(quantile(A1,.05,1));
RES.detectTime.singlesample.p95  = squeeze(quantile(t1,.95,1));
RES.detectTime.singlesample.p50  = squeeze(quantile(t1,.5,1));
RES.detectTime.singlesample.p05  = squeeze(quantile(t1,.05,1));

RES.maxSlope.fixedinterval.average = m2;
RES.maxSlope.fixedinterval.ci    = se2;
RES.maxSlope.fixedinterval.p95   = squeeze(q2(3,:,:));
RES.maxSlope.fixedinterval.p50   = squeeze(q2(2,:,:));
RES.maxSlope.fixedinterval.p05   = squeeze(q2(1,:,:));
RES.amplitude.fixedinterval.p95  = squeeze(quantile(A2,.95,1));
RES.amplitude.fixedinterval.p50  = squeeze(quantile(A2,.5,1));
RES.amplitude.fixedinterval.p05  = squeeze(quantile(A2,.05,1));
RES.detectTime.fixedinterval.p95 = squeeze(quantile(t2,.95,1));
RES.detectTime.fixedinterval.p50 = squeeze(quantile(t2,.5,1));
RES.detectTime.fixedinterval.p05 = squeeze(quantile(t2,.05,1));
RES.FList = F; 
RES.ZList = Z;
RES.dQ    = dQ; 

outfile = sprintf('slope_estimates_%03i',wintime*1e6);
save('-mat', [outfile,'.mat'], 'RES','F','Z','maxSlope0','A0','t0');
save('-mat', 'mat.mat', 'outfile');


fid = 1; 
fid = fopen([outfile,'.tsv'],'w+'); 
fprintf(fid, '# Temperature: %f K\n',T);
fprintf(fid, '\n# List of Impedances [Ohm]\n');
fprintf(fid, '%.9g\t',Z);
fprintf(fid, '\n\n# List of Sampling Rates\n');
fprintf(fid, '%.9g\t',F);

fprintf(fid, '\n\n# True maximum slope\n');
fprintf(fid, '%.9g\t',RES.maxSlope.true );
#fprintf(fid, '\n\n# Maximum slope of noise free data\n');
#fprintf(fid, '%.9g\t',RES.maxSlope.noisefree);
fprintf(fid, '\n\n# Median of Maximum slope from single-noise free estimation \n');
for k=1:size(RES.maxSlope.noisefree.p50,1)
        fprintf(fid,"%.9g\t",RES.maxSlope.noisefree.p50(k,:));
        fprintf(fid,"\n");
end; 
fprintf(fid, '\n\n# 95 percentile of Maximum slope from single-noise free estimation \n');
for k=1:size(RES.maxSlope.noisefree.p50,1)
        fprintf(fid,"%.9g\t",RES.maxSlope.noisefree.p95(k,:));
        fprintf(fid,"\n");
end; 
fprintf(fid, '\n\n#  5 percentile of Maximum slope from single-noise free estimation \n');
for k=1:size(RES.maxSlope.noisefree.p50,1)
        fprintf(fid,"%.9g\t",RES.maxSlope.noisefree.p05(k,:));
        fprintf(fid,"\n");
end; 

fprintf(fid, '\n\n# Median of Maximum slope from single sample estimation \n');
for k=1:size(RES.maxSlope.fixedinterval.p50,1)
        fprintf(fid,"%.9g\t",RES.maxSlope.singlesample.p50(k,:));
        fprintf(fid,"\n");
end; 
fprintf(fid, '\n# 95 percentile of Maximum slope from single sample estimation \n');
for k=1:size(RES.maxSlope.fixedinterval.p95,1)
        fprintf(fid,"%.9g\t",RES.maxSlope.singlesample.p95(k,:));
        fprintf(fid,"\n");
end; 
fprintf(fid, '\n#  5 percentile of Maximum slope from single sample estimation \n');
for k=1:size(RES.maxSlope.fixedinterval.p05,1)
        fprintf(fid,"%.9g\t",RES.maxSlope.singlesample.p05(k,:));
        fprintf(fid,"\n");
end; 

fprintf(fid, '\n# Median of Maximum slope from fixed window (%f ms) sample estimation \n', wintime*1000 );
for k=1:size(RES.maxSlope.fixedinterval.p50,1)
        fprintf(fid,"%.9g\t",RES.maxSlope.fixedinterval.p50(k,:));
        fprintf(fid,"\n");
end; 
fprintf(fid, '\n# 95 percentile of Maximum slope from fixed window (%f ms)  sample estimation \n', wintime*1000);
for k=1:size(RES.maxSlope.fixedinterval.p95,1)
        fprintf(fid,"%.9g\t",RES.maxSlope.fixedinterval.p95(k,:));
        fprintf(fid,"\n");
end; 
fprintf(fid, '\n#  5 percentile of Maximum slope from fixed window (%f ms)  sample estimation \n', wintime*1000);
for k=1:size(RES.maxSlope.fixedinterval.p05,1)
        fprintf(fid,"%.9g\t",RES.maxSlope.fixedinterval.p05(k,:));
        fprintf(fid,"\n");
end; 
fclose(fid); 



%clear

%else 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  OUTPUT: display results
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load mat.mat
load('-mat', [outfile,'.mat']);


%maxSlope2
figure(1)
hold off
%plot(F'/1000,m2,'b',F'/1000,m1,'r',F',maxSlope0,'k+-')
g = errorbar(repmat(F'/1000,1,length(Z)),RES.maxSlope.singlesample.average,RES.maxSlope.singlesample.ci,'r',repmat(F'/1000,1,length(Z)),RES.maxSlope.fixedinterval.average,RES.maxSlope.fixedinterval.ci,'b');
for k=1:length(Z),set(g(k),  'linewidth',3,'color',[1,[.5-k/14,.5-k/14]*1.5]);end;
for k=1:length(Z),set(g(k+length(Z)),'linewidth',3,'color',[[.5-k/14,.5-k/14]*1.5,1]);end;
hold on
h=plot(F'/1000,maxSlope0,'k-',[F(1),F(end)]/1000,max(diff(d))*H.SampleRate*[1,1],'--g');
xlabel('sampling rate [kHz]')
ylabel('estimated maximum slope [V/s] (mean and c.i. of 1 s.d.)')
set(h(1),'linewidth',3);
set(gca,'box','off')

print('-depsc2', [outfile,'.eps']);
print('-dpdf', [outfile,'.pdf']);
print('-dmeta', [outfile,'.wmf']);
print('-dpng', [outfile,'.png']);


figure(2)
hold off
%plot(F'/1000,m2,'b',F'/1000,m1,'r',F',maxSlope0,'k+-')
g=errorbar(repmat(F'/1000,1,length(Z)),RES.maxSlope.singlesample.p50,RES.maxSlope.singlesample.p50-RES.maxSlope.singlesample.p05,RES.maxSlope.singlesample.p95-RES.maxSlope.singlesample.p50,'b');
for k=1:length(Z),set(g(k),'linewidth',1,'color',[1,[.5-k/14,.5-k/14]*1.5]);end;
set(g([1:length(Z)]),'linewidth',3);
hold on
g=errorbar(repmat(F'/1000,1,length(Z)),RES.maxSlope.fixedinterval.p50,RES.maxSlope.fixedinterval.p50-RES.maxSlope.fixedinterval.p05,RES.maxSlope.fixedinterval.p95-RES.maxSlope.fixedinterval.p50,'b');
for k=1:length(Z),set(g(k),'linewidth',1,'color',[[.5-k/14,.5-k/14]*1.5,1]);end;
set(g([1:length(Z)]),'linewidth',3);
g=errorbar(repmat(F'/1000,1,1),RES.maxSlope.noisefree.p50,RES.maxSlope.fixedinterval.p50-RES.maxSlope.fixedinterval.p05,RES.maxSlope.fixedinterval.p95-RES.maxSlope.fixedinterval.p50,'k');
set(g,'linewidth',3);
%h=plot(F'/1000,maxSlope0,'k-',[F(1),F(end)]/1000,max(diff(d))*H.SampleRate*[1,1],'--g');
h=plot([F(1),F(end)]/1000,max(diff(d))*H.SampleRate*[1,1],'--g');
xlabel('sampling rate [kHz]')
ylabel('maximum slope [V/s] (median and 90% c.i.)')
set(h(1),'linewidth',3);
set(gca,'box','off')
set(gca,'xlim',[0,200],'ylim',[200,600])


if 0,
figure(2)
hold off
%plot(F'/1000,m2,'b',F'/1000,m1,'r',F',maxSlope0,'k+-')
g=errorbar(repmat(F'/1000,1,length(Z)),RES.maxSlope.singlesample.p50-RES.maxSlope.noisefree(:,ones(1,length(Z))),RES.maxSlope.singlesample.p50-RES.maxSlope.noisefree(:,ones(1,length(Z)))-RES.maxSlope.singlesample.p05,RES.maxSlope.singlesample.p95-RES.maxSlope.noisefree(:,ones(1,length(Z)))-RES.maxSlope.singlesample.p50,'b');
for k=1:length(Z),set(g(k),'linewidth',1,'color',[1,[.5-k/14,.5-k/14]*1.5]);end;
set(g([1:length(Z)]),'linewidth',3);
hold on
g=errorbar(repmat(F'/1000,1,length(Z)),RES.maxSlope.fixedinterval.p50-RES.maxSlope.noisefree(:,ones(1,length(Z))),RES.maxSlope.fixedinterval.p50-RES.maxSlope.noisefree(:,ones(1,length(Z)))-RES.maxSlope.fixedinterval.p05,RES.maxSlope.fixedinterval.p95-RES.maxSlope.noisefree(:,ones(1,length(Z)))-RES.maxSlope.fixedinterval.p50,'b');
for k=1:length(Z),set(g(k),'linewidth',1,'color',[[.5-k/14,.5-k/14]*1.5,1]);end;
set(g([1:length(Z)]),'linewidth',3);
hold on
h=plot(F'/1000,maxSlope0,'k-',[F(1),F(end)]/1000,max(diff(d))*H.SampleRate*[1,1],'--g');
xlabel('sampling rate [kHz]')
ylabel('maximum slope [V/s] (median and 90% c.i.)')
set(h(1),'linewidth',3);
set(gca,'box','off')
set(gca,'xlim',[0,200],'ylim',[200,600])
end;



outfile = sprintf('slope_estimates_%03i.med90ci',wintime*1e6);
print('-depsc2', [outfile,'.eps']);
print('-dpdf', [outfile,'.pdf']);
print('-dmeta', [outfile,'.wmf']);
print('-dpng', [outfile,'.png']);


figure(4)
hold off
%plot(F'/1000,m2,'b',F'/1000,m1,'r',F',maxSlope0,'k+-')
g=errorbar(repmat(F'/1000,1,length(Z)),RES.amplitude.singlesample.p50, RES.amplitude.singlesample.p50-RES.amplitude.singlesample.p05, RES.amplitude.singlesample.p95-RES.amplitude.singlesample.p50, 'r'); 
for k=1:length(Z),set(g(k),'linewidth',1,'color',[1,[.5-k/14,.5-k/14]*1.5]);end;
set(g([1:3:length(Z)]),'linewidth',3);
hold on
g=errorbar(repmat(F'/1000,1,length(Z)),RES.amplitude.fixedinterval.p50,RES.amplitude.fixedinterval.p50-RES.amplitude.fixedinterval.p05,RES.amplitude.fixedinterval.p95-RES.amplitude.fixedinterval.p50,'b');
for k=1:length(Z),set(g(k),'linewidth',1,'color',[[.5-k/14,.5-k/14]*1.5,1]);end;
set(g([1:3:length(Z)]),'linewidth',3);
h=plot(F'/1000,A0,'k-');
xlabel('sampling rate [kHz]')
ylabel('Amplitude [V] (median and 90% c.i.)')
set(h(1),'linewidth',3);
set(gca,'box','off')

set(gca,'ylim',[-.05,-.03])


outfile = sprintf('amplitude_estimates_%03i.med90ci',wintime*1e6);
print('-depsc2', [outfile,'.eps']);
print('-dpdf', [outfile,'.pdf']);
print('-dmeta', [outfile,'.wmf']);
print('-dpng', [outfile,'.png']);


figure(6)
hold off
%plot(F'/1000,m2,'b',F'/1000,m1,'r',F',maxSlope0,'k+-')
g=errorbar(repmat(F'/1000,1,length(Z)),RES.detectTime.singlesample.p50, RES.detectTime.singlesample.p50-RES.detectTime.singlesample.p05, RES.detectTime.singlesample.p95-RES.detectTime.singlesample.p50, 'r'); 
for k=1:length(Z),set(g(k),'linewidth',1,'color',[1,[.5-k/14,.5-k/14]*1.5]);end;
set(g([1:3:length(Z)]),'linewidth',3);
hold on
g=errorbar(repmat(F'/1000,1,length(Z)),RES.detectTime.fixedinterval.p50,RES.detectTime.fixedinterval.p50-RES.detectTime.fixedinterval.p05,RES.detectTime.fixedinterval.p95-RES.detectTime.fixedinterval.p50,'b');
for k=1:length(Z),set(g(k),'linewidth',1,'color',[[.5-k/14,.5-k/14]*1.5,1]);end;
set(g([1:3:length(Z)]),'linewidth',3);
h=plot(F'/1000,t0,'k-');
xlabel('sampling rate [kHz]')
ylabel('time [s] (median and 90% c.i.)')
set(h(1),'linewidth',3);
set(gca,'box','off')
set(gca,'ylim',[.007,.008])

outfile = sprintf('onsettime_estimates_%03i.med90ci',wintime*1e6);
print('-depsc2', [outfile,'.eps']);
print('-dpdf', [outfile,'.pdf']);
print('-dmeta', [outfile,'.wmf']);
print('-dpng', [outfile,'.png']);




disp('% List of Impedances')
disp(Z)
disp('% List of Sampling Frequencies')
disp(F)
sprintf('%% Maximum Slope after downsampling, without any Impedance or Quantisation noise. The average of N=%i repeated generation of Z-noise is shown.',M)
disp(maxSlope0)
disp('% Maximum Slope obtained by single sample differences')
disp(RES.maxSlope.singlesample.average)
disp('% S.D. of Maximum Slope estimate obtained by single sample differences.\nThe first column (Z=0), shows the influence of the Q-noise, only')
disp(max(dQ/sqrt(12)*F'*ones(1,length(Z)),RES.maxSlope.singlesample.ci))
sprintf('%% Maximum Slope obtained by a fixed window length of 0.05 ms. The average of N=%i repeated generation of Z-noise is shown.',M);
disp(RES.maxSlope.fixedinterval.average)
disp('% S.D. of Maximum Slope estimate obtained by a fixed window length of 0.05 ms\nThe first column (Z=0), shows the influence of the Q-noise, only')
disp(max(dQ/sqrt(12)*min(F',1/wintime)*ones(1,length(Z)),RES.maxSlope.fixedinterval.ci))


[F',maxSlope0]
[0,Z ;F',RES.maxSlope.singlesample.average]
[0,Z ;F',RES.maxSlope.singlesample.ci]
[0,Z ;F',RES.maxSlope.fixedinterval.average]
[0,Z ;F',RES.maxSlope.fixedinterval.ci]


end;

