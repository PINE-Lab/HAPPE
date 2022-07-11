p=pwd
cd ..
install
cd(p)

if exist('OCTAVE_VERSION','builtin')
	pkg load mexbiosig
	pkg load optim
else
	addpath /fs3/home/schloegl/src/biosig-code/biosig4matlab/t330_StimFit
	addpath /fs3/home/schloegl/src/biosig-code/biosig4c++/mex
end

simul001; 
b  = [-20:10:20]'; % pA
A  = [10,30,100,300,1000,2000]'; % pA
t0 = [0:1:4]';  % ms
t1 = [0.1,1.4,0.2,.5,.7,.1,1.4]'; % ms
t2 = [3,5,7,10,15]'; % ms


DIM = [length(b),length(A),length(t0),length(t1),length(t2)]
ix = reshape(1:prod(DIM),DIM);

params = repmat(NaN,prod(DIM),5);
for k= 1:prod(DIM);
    [bix, Aix, t0ix, t1ix, t2ix] = ind2sub(DIM, k);
    params(k,:)=[A(Aix),b(bix),t0(t0ix),t1(t1ix),t2(t2ix)];
end;


% load data 
[data,HDR]=mexSLOAD('test01.gdf', 1);
Fs=round(HDR.SampleRate);

% run microstimfit
default.t1=round(-Fs*10e-3);
default.t2=round(+Fs*100e-3);
default.baseBegin=round(-Fs*7e-3);
default.baseEnd=round(-Fs*2e-3);
default.peakBegin=round(-Fs*1e-3);
default.peakEnd=round(+Fs*10e-3);
default.fitEnd=round(+Fs*50e-3);
default.meanN=1;
default.dir=0;
default.plotFlag=0;
default.baseFlag=0;
default.fitFlag=1;
default.thres=0.03;
default.fitFlag=1;
default.thresFlag=0;

[results, opt] = microstimfit(data, HDR.SampleRate, [0:5249]*2201+201, default);


subplot(3,1,1)
plot([results.data(:,17+1),-params(:,1)])
title('comarison between true and estimated model parameters')
legend({'A estimated','A'},'box','off')
subplot(3,1,2)
plot([results.data(:,17+2),params(:,2)])
legend({'b estimated','b'},'box','off')
subplot(3,1,3)
plot([results.data(:,17+3),params(:,5)*Fs/1000])
legend({'tau estimated','tau2'},'box','off')


