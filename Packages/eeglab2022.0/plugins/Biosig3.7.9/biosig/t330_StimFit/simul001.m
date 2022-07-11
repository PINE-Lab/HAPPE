% simulate EPSC with varying parameters 

% Copyright (C) 2013,2019 Alois Schlögl, IST Austria 

% REFERENCES: 
% [1] Jose Guzman, Alois Schlögl, Christoph Schmidt-Hieber
%     Stimfit: quantifying electrophysiological data with Python.
%     Front. Neuroinform. 8:16, 2014
%     available online: doi: http://dx.doi.org/10.3389/fninf.2014.00016
%     https://pub.ist.ac.at/~schloegl/publications/GuzmanEtAl2014.fninf-08-00016.pdf

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

clear 

%%% sampleing rate
Fs = 20000; % kHz
t = -10:1000/Fs:100;    % ms 

%%% Model
M  = @(b,A,t0,t1,t2)(A*(exp((t0-t)/t1)-exp((t0-t)/t2)).*(t>t0)+b);
M1 = @(b,A,t0,t1,t2)(A*(1-exp((t0-t)/t1)).*exp((t0-t)/t2).*(t>t0)+b);

b  = [-20:10:20]'; % pA
A  = [10,30,100,300,1000,2000]'; % pA
t0 = [0:1:4]';  % ms
t1 = [0.1,1.4,0.2,.5,.7,.1,1.4]'; % ms
t2 = [3,5,7,10,15]'; % ms


DIM = [length(b),length(A),length(t0),length(t1),length(t2)]
ix = reshape(1:prod(DIM),DIM);

y = repmat(NaN,prod(DIM),length(t)); 
for k= 1:prod(DIM); 
    [bix, Aix, t0ix, t1ix, t2ix] = ind2sub(DIM, k);
    y(k,:)=M(b(bix),A(Aix),t0(t0ix),t1(t1ix),t2(t2ix));
end;     


% plot(t,y')
% set(gca,'ylim',[-2000,2000]);
 
 
HDR.NS = 1;
HDR.TYPE='GDF'; 
HDR.T0 = now;
HDR.Label = {'Simulated EPSC'};
HDR.Transducer = {'Octave'};
HDR.GDFTYP = 16; 
HDR.PhysDim = {'pA'}; 
HDR.PhysDimCode = physicalunits({'pA'});
[HDR.NRec, HDR.SPR] = size(y); 
HDR.SampleRate = Fs; 
HDR.PhysMax =  2*max(A); 
HDR.PhysMin = -2*max(A); 
HDR.DigMax  = HDR.PhysMax; 
HDR.DigMin  = HDR.PhysMin;


HDR.EVENT.SampleRate = Fs; 
HDR.EVENT.N   = HDR.NRec-1; 
HDR.EVENT.POS = [1:prod(DIM)-1]*length(t); 
HDR.EVENT.TYP = repmat(hex2dec('7ffe'),size(HDR.EVENT.POS)); 

HDR.FileName = 'test01.gdf'; 
HDR = sopen(HDR,'w'); 
fwrite(HDR.FILE.FID,y','float'); 
HDR = sclose(HDR);

if exist('copytest01.out','file')
    [bix, Aix, t0ix, t1ix, t2ix] = ind2sub(DIM, [1:prod(DIM)]');
    X = [b(bix),A(Aix),t0(t0ix),t1(t1ix),t2(t2ix)];
    labelx={'b','A','t0','t1','t2'};
    fid=fopen('origpara_test01.txt','w');
    fprintf(fid,'%s\t',labelx{:});	
    fprintf(fid,'\n');
    fprintf(fid,'%g\t%g\t%g\t%g\t%g\n',X');
    fclose(fid); 	
 
    fid = fopen('copytest01.out','r');
    labely=strsplit(fgetl(fid),"\t");
    fclose(fid);
    Y=dlmread('copytest01.out',' \t',1,0);
 
    fid=fopen('alltest01.txt','w');
    fprintf(fid,'%s\t',labelx{:});	
    fprintf(fid,'%s\t',labely{:});	
    fprintf(fid,'\n');
    fprintf(fid,'%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n',[X,Y]');
    fclose(fid); 	
    fid=fopen('alltest01.noroundingerror.txt','w');
    fprintf(fid,'%s\t',labelx{:});	
    fprintf(fid,'%s\t',labely{:});	
    fprintf(fid,'\n');
    fprintf(fid,'%g\t%g\t%g\t%g\t%g\t%.16g\t%.16g\t%.16g\t%.16g\t%.16g\t%.16g\n',[X,Y]');
    fclose(fid); 	

 
    Z.data = [X,Y]; 
    Z.datatype = 'SCATTER';
    Z.Labels = [labelx,labely];
 
    [r,p]=corrcoef(X,Y);
    [tmp,ix]=max(abs(r));
 

    sel1 = abs(X(:,4)-Y(:,1))<1e-6;
    sel3 = abs(X(:,3)-Y(:,3)+10)<1e-6;
    sel4 = abs(X(:,2)+Y(:,5))<1e-4;
    sel  = all([sel1,sel3,sel4],2);
    %sel = abs(X(:,3)-Y(:,1))>-1;
    figure(1); clf;
    for k=1:length(ix),
        subplot(3,2,k);
    	y = Y(sel,k)*sign(r(ix(k),k));
    	x = X(sel,ix(k));
        plot(x,y,'x');
        xlabel(labelx(ix(k)));	
        ylabel(labely(k));	
    end;
    print -dpng -f1 f1.png

    fprintf(1,'Label\tcorrelation\tp-value\tbias\trelative error\n');
    for k=1:length(ix),
    	y = Y(sel,k)*sign(r(ix(k),k));
    	x = X(sel,ix(k));
	e = x-y;
	figure(k+1); 
	plot([x+y]/2,e,'x');	
	ylabel('error');
	%legend([labelx{ix(k)},'+',labely{k}]);
	xlabel([labelx{ix(k)},'+',labely{k}]);
	bias   = mean(e); 	
	relerr = std(e)/sqrt(std(x)*std(y));
        fprintf(1,'%s-%s:\t%.3f %g\t%f\t%f\n',labelx{ix(k)},labely{k},r(ix(k),k),p(ix(k),k),bias,relerr);
    end
end;
 
