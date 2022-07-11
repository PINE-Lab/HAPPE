% SIMULATE_EPSP  generates a large number sweeps of EPSP data using
%         different models, and parameters for validation of Stimfit model
%         fitting algorithms [3] 
% For each model, a separate GDF file with sweeps of varying model 
%    paramters is generated. 

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

%%% sampleing rate
Fs = 20000; % kHz
t = -10:1000/Fs:100;    % ms 

for kk=1:8

clear P; 

switch (kk) 
case 1  %%% monoexponential free fit 
	M  = @(t0,A,t1,b)(A*exp((t0-t)/t1).*(t>t0)+b);
	b  = [-20:10:20]'; % pA
	A  = [10,30,100,300,1000,2000]'; % pA
	t0 = [0]';  % ms
	t1 = [0.1,1.4,0.2,.5,.7,.1,1.4,3,5,7,10,15]'; % ms

	DIM = [length(b),length(A),length(t0),length(t1)]
	ix = reshape(1:prod(DIM),DIM);

	y = repmat(NaN,prod(DIM),length(t)); 
	for k= 1:prod(DIM); 
    		[bix, Aix, t0ix, t1ix] = ind2sub(DIM, k);
		y(k,:)=M(b(bix),A(Aix),t0(t0ix),t1(t1ix));
	end;     

case 2  %%% monoexponential - fixed baseline 
	M  = @(t0,A,t1,b)(A*exp((t0-t)/t1).*(t>t0)+b);
	b  = [0]'; % pA
	A  = [10,30,100,300,1000,2000]'; % pA
	t0 = [0]';  % ms
	t1 = [0.1,1.4,0.2,.5,.7,.1,1.4,3,5,7,10,15]'; % ms

	DIM = [length(b),length(A),length(t0),length(t1)]
	ix = reshape(1:prod(DIM),DIM);

	y = repmat(NaN,prod(DIM),length(t)); 
	for k= 1:prod(DIM); 
    		[bix, Aix, t0ix, t1ix] = ind2sub(DIM, k);
		y(k,:)=M(b(bix),A(Aix),t0(t0ix),t1(t1ix));
	end;     

case 3  %%% monoexponential + delay - fixed baseline 
	M  = @(b,A,t0,t1)(A*exp((t0-t)/t1).*(t>t0)+b);
	b  = [0]'; % pA
	A  = [10,30,100,300,1000,2000]'; % pA
	t0 = [0:1:4]';  % ms
	t1 = [0.1,1.4,0.2,.5,.7,.1,1.4,3,5,7,10,15]'; % ms

	DIM = [length(b),length(A),length(t0),length(t1)]
	ix = reshape(1:prod(DIM),DIM);

	y = repmat(NaN,prod(DIM),length(t)); 
	for k= 1:prod(DIM); 
    		[bix, Aix, t0ix, t1ix] = ind2sub(DIM, k);
		y(k,:)=M(b(bix),A(Aix),t0(t0ix),t1(t1ix));
	end;     

case 4  %%% biexponential + freefit 
	M  = @(t0,A1,t1,A2,t2,b)(A1*exp((t0-t)/t1 + A2*exp((t0-t)/t2 )).*(t>t0)+b);

	b  = [-20:10:20]'; % pA
	A1 = [10,30,100,300,1000,2000]'; % pA
	A2 = [10,30,100,300,1000,2000]'; % pA
	t0 = [0]';  % ms
	t1 = [0.1,1.4,0.2,.5,.7,.1,1.4,3,5,7,10,15]'; % ms
	t2 = NaN;
	A1 = [-A1;A1]; 
	A2 = [-A2;A2]; 

	DIM = [length(b),length(t0),length(A1),length(t1),length(A2),length(t2)]
	ix = reshape(1:prod(DIM),DIM);

	y = repmat(NaN,prod(DIM),length(t)); 
	for k= 1:prod(DIM); 
    		[bix, t0ix, A1ix, t1ix, A2ix, t2ix] = ind2sub(DIM, k);
		y(k,:)=M(b(bix),t0(t0ix),A1(A1ix),t1(t1ix),A1(A1ix),t2(t2ix));
	end;     

case 5  %%% biexponential, offset fixed 
	M  = @(t0,A1,t1,A2,t2,b)(A1*exp((t0-t)/t1 + A2*exp((t0-t)/t2 )).*(t>t0)+b);
	b  = [0]'; % pA
	A1 = [10,30,100,300,1000,2000]'; % pA
	A2 = [10,30,100,300,1000,2000]'; % pA
	t0 = [0]';  % ms
	t1 = [0.1,1.4,0.2,.5,.7,.1,1.4,3,5,7,10,15]'; % ms
	t2 = NaN;
	A1 = [-A1;A1]; 
	A2 = [-A2;A2]; 

	DIM = [length(b),length(t0),length(A1),length(t1),length(A2),length(t2)]
	ix = reshape(1:prod(DIM),DIM);

	y = repmat(NaN,prod(DIM),length(t)); 
	for k= 1:prod(DIM); 
    		[bix, t0ix, A1ix, t1ix, A2ix, t2ix] = ind2sub(DIM, k);
		y(k,:)=M(b(bix),t0(t0ix),A1(A1ix),t1(t1ix),A1(A1ix),t2(t2ix));
	end;     

case 6  %%% Model: biexponential with delay 
	M  = @(b,A,t0,t1,t2)(A*(exp((t0-t)/t1)-exp((t0-t)/t2)).*(t>t0)+b);
	M1 = @(b,A,t0,t1,t2)(A*(1-exp((t0-t)/t1)).*exp((t0-t)/t2).*(t>t0)+b);

	b  = [-20:10:20]'; % pA
	A  = [10,30,100,300,1000,2000]'; % pA
	t0 = [0:1:4]';  % ms
	t2 = [0.1,1.4,0.2,.5,.7,.1,1.4]'; % ms
	t1 = [3,5,7,10,15]'; % ms

	DIM = [length(b),length(A),length(t0),length(t1),length(t2)]
	ix = reshape(1:prod(DIM),DIM);

	y = repmat(NaN,prod(DIM),length(t)); 
	for k= 1:prod(DIM); 
    		[bix, Aix, t0ix, t1ix, t2ix] = ind2sub(DIM, k);
		y(k,:)=M(b(bix),A(Aix),t0(t0ix),t1(t1ix),t2(t2ix));
	end;     


case 7  %%% triexponential + freefit 
	M  = @(t0,A1,t1,A2,t2,A3,t3,b)( (A1*exp((t0-t)/t1) + A2*exp((t0-t)/t2)  + A3*exp((t0-t)/t3) ).*(t>t0) + b);
	b  = [-20:10:20]'; % pA
	A1 = [10,30,100,300,1000,2000]'; % pA
	A2 = [10,30,100,300,1000,2000]'; % pA
	A3 = [10,30,100,300,1000,2000]'; % pA
	t0 = [0]';  % ms
	t1 = [3,5,7,10,15]'; % ms
	t2 = [0.1,1.4,0.2,.5,.7,.1,1.4]'; % ms
	t3 = [t1;t2]; t3=t3(1:2:end);
	A1 = [-A1;A1](1:2:end); 
	A2 = [-A2;A2](2:2:end); 
	A3 = [-A3;A3](1:2:end); 

	DIM = [length(t0),length(A1),length(t1),length(A2),length(t2),length(A3),length(t3),length(b)];
	ix = reshape(1:prod(DIM),DIM);

	y = repmat(NaN,prod(DIM),length(t)); 
	for k= 1:prod(DIM); 
	    [t0ix, A1ix, t1ix, A2ix, t2ix, A3ix, t3ix, bix] = ind2sub(DIM, k);
	    y(k,:)=M(t0(t0ix),A1(A1ix),t1(t1ix),A2(A2ix),t2(t2ix),A3(A3ix),t3(t3ix),b(bix));
	end;     

case 8  %%% triexponential + freefit, fix offset 
	M  = @(t0,A1,t1,A2,t2,A3,t3,b)( (A1*exp((t0-t)/t1) + A2*exp((t0-t)/t2)  + A3*exp((t0-t)/t3) ).*(t>t0) + b);
	b  = [0]'; % pA
	A1 = [10,30,100,300,1000,2000]'; % pA
	A2 = [10,30,100,300,1000,2000]'; % pA
	A3 = [10,30,100,300,1000,2000]'; % pA
	t0 = [0]';  % ms
	t1 = [3,5,7,10,15]'; % ms
	t2 = [0.1,1.4,0.2,.5,.7,.1,1.4]'; % ms
	t3 = [t1;t2]; t3=t3(1:2:end);
	A1 = [-A1;A1](1:2:end); 
	A2 = [-A2;A2](2:2:end); 
	A3 = [-A3;A3](1:2:end); 

	DIM = [length(t0),length(A1),length(t1),length(A2),length(t2),length(A3),length(t3),length(b)];
	ix = reshape(1:prod(DIM),DIM);

	y = repmat(NaN,prod(DIM),length(t)); 
	for k= 1:prod(DIM); 
	    [t0ix, A1ix, t1ix, A2ix, t2ix, A3ix, t3ix, bix] = ind2sub(DIM, k);
	    y(k,:)=M(t0(t0ix),A1(A1ix),t1(t1ix),A2(A2ix),t2(t2ix),A3(A3ix),t3(t3ix),b(bix));
	end;     

end; 

clear HDR
HDR.NS = 1;
HDR.TYPE='GDF'; 
HDR.T0 = now;
HDR.Label = {'Simulated EPSC'};
HDR.Transducer = {'Octave'};
HDR.GDFTYP = 16; 
HDR.PhysDim = {'pA'}; 
HDR.PhysDimCode = physicalunits('pA');
[HDR.NRec, HDR.SPR] = size(y); 
HDR.SampleRate = Fs; 
HDR.PhysMax =  7000
HDR.PhysMin = -7000; 
HDR.DigMax  = HDR.PhysMax; 
HDR.DigMin  = HDR.PhysMin;


HDR.EVENT.SampleRate = Fs; 
HDR.EVENT.N   = HDR.NRec-1; 
HDR.EVENT.POS = [1:prod(DIM)-1]*length(t)+1;
HDR.EVENT.TYP = repmat(hex2dec('7ffe'),size(HDR.EVENT.POS)); 

HDR.FileName = sprintf('test.%i.%i.gdf',prod(DIM),kk); 
HDR = sopen(HDR,'w'); 
fwrite(HDR.FILE.FID,y','float'); 
HDR = sclose(HDR);


end; 
