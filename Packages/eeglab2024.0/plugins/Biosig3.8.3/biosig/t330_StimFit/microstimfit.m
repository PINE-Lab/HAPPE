function [results, option] = microstimfit(data, Fs, evtpos, varargin)
% microStimfit for analysis of PSC's in Matlab. 
%  It is a re-implementention of Peter Jonas' implementation for Wolfram Mathematica, 
%  but it supports also (i) a continues data stream, with (ii) multiple 
%  events, and (3) default parameter settings, to allow an easier start. 
%
% Prerequisite: 
% [1] Biosig for Octave and Matlab
%    http://biosig.sourceforge.net/documentation.html
% [2a] Matlab, and Optimization toolbox 
% [2b] Octave, and Optim package
%
% Usage: 
%  [...] = microstimfit(data, samplerate, event_pos, ...)
%  [...] = microstimfit(filename, channel, event_pos, ...)
%  [results, opt] = microstimfit(filename, channel, event_pos, option)
%  [results, opt] = microstimfit(data, samplerate, event_pos, option)
%
%  Input arguments:
%      filename 	name of biosig file
%      channel 		channel number 
%      data 		vector of data samples with 	 
%      samplerate	sampling rate of data 
%      event_pos	position (in samples) of the event in the data 
%
%      option.t1	start of analysing window [in samples relative to event_pos]
%      option.t2	end of analysing window [in samples relative to event_pos]
%      option.baseBegin	[in samples relative to event_pos]
%      option.baseEnd	[in samples relative to event_pos]
%      option.peakBegin	[in samples relative to event_pos]
%      option.peakEnd	[in samples relative to event_pos]
%      option.meanN	smoothing window, in order to get the peak
%      option.dir	direction of peak( 1: up, -1: down, 0: up or down)
%      option.plotFlag	1: draw figure and wait for keystrock;	0: no plotting; -1: draw figure and continue
%      optioon.baseFlag: 0 (default): mean (average) of baseline
%                        1 : median of baseline
%      option.fitFlag	0 [default]: no fitting
%			1: fit decay with monoexponential model 'a*exp(-x/tau)+offset' (from peakTime to t2)
%			2: fit decay with biexponential model 'a1*exp(-x/tau1)+a2*exp(-x/tau2) ' (from peakTime to t2)
%			3: fit decay with biexponential model 'a1*exp(-x/tau1)+a2*exp(-x/tau2)+offset' (from peakTime to t2)
%      option.thres	threshold value for AP (in relative units or unit voltage or current per sample interval), and 
%      option.thresFlag	controls the threshold criterion (0 = relative to max slope, 1 = absolute). 
%      option.fitEnd	data length for fitting, starting from the peak time 
%
%  Output arguments:
%	opt		same as option, adding the default values for unspecified option fields
%	results.label	cell array, containing the label for each column in results.data
%	results.data	table containing the resulting values, NaN indicates the result is not available. 
%
%  Figure plots:
%	Crosshair indicates peak, points: 
%	black = 0%, red = 20%, green = 80%, blue = 50%, cyan = max dV / dt, magenta = threshold. 
%
% REFERENCES: 
% [1] Jose Guzman, Alois Schlögl, Christoph Schmidt-Hieber
%     Stimfit: quantifying electrophysiological data with Python.
%     Front. Neuroinform. 8:16, 2014
%     available online: doi: http://dx.doi.org/10.3389/fninf.2014.00016
%     https://pub.ist.ac.at/~schloegl/publications/GuzmanEtAl2014.fninf-08-00016.pdf


% Copyright (C) 2017-2019,2022 Alois Schlögl, IST Austria <alois.schloegl@ist.ac.at>

% TODO: 
% smoothing window length 
% 	outer and inner rise time
% 	adding more fitting models
% 
% DONE:
%	peak smoothing
%	BASELINE: median, IQR
%	smoothing is applied only for peak identification
%	validations succesful with 
%		octave 4.4.1, optim 1.6.0, 
%		matlab R2019b with curve fitting toolbox 

opt.out.base=[];
opt.out.baseSD=[];
opt.out.peak=[];
opt.out.tPeak=[];
opt.out.myf=[];
opt.out.fExp=[];
opt.out.ampl=[];
opt.out.offs=[];
opt.out.tau=[];
opt.out.fitResult=[];
opt.out.threshold=[];
opt.out.gr1=[];
opt.out.gr2=[];
opt.out.t20Int=[];
opt.out.t20Real=[];
opt.out.t80Int=[];
opt.out.t80Real=[];
opt.out.t50AInt=[];
opt.out.t50AReal=[];
opt.out.t50BInt=[];
opt.out.t50BReal=[];
opt.out.t0Real=[];
opt.out.tMaxSlopeRiseReal=[];
opt.out.tThreshold=[];
opt.out.yThreshold=[];
opt.out.maxSlopeRiseReal=[];

if nargin<3,
	error('number of arguments must larger than 2');
end

if ischar(data) && exist(data,'file') && isscalar(Fs)
	filename=data; 
	chan=Fs
	[data,HDR]=sload(filename,chan);
	Fs = HDR.SampleRate; 
elseif isnumeric(data) && (sum(size(data)>1)==1) && isscalar(Fs) && (Fs>0)
	; %noop
else
	error('invalid input arguments - use either %s(file,channel,..) or %s(data,Fs,..) ', mfilename(), mfilename());
end

% default values
default.t1=round(-Fs*100e-3);
default.t2=round(+Fs*100e-3);
default.baseEnd=round(-Fs*10e-3);
default.peakBeg=round(-Fs*10e-3);
default.peakEnd=round(+Fs*20e-3);
default.fitEnd=round(+Fs*50e-3);
default.meanN=round(50e-6*Fs);	% see [1], Fig. 3-C
default.dir=0;
default.plotFlag=0;
default.baseFlag=0;
default.fitFlag=0;
default.thres=0.03;
default.thresFlag=0;

if (nargin>3) && isstruct(varargin{1})
	option = varargin{1};
	F      = fieldnames(default);
	for k  = 1:length(F)
		f = F{k};
		if ~isfield(option,f)
			option = setfield(option, f, getfield(default,f));
		end
	end
else
    option = default;	
    flag_numeric_only = 1;
    for k = 1:length(varargin)
	if ischar(varargin{k})
		flag_numeric_only = 0;
	end
	if flag_numeric_only,
		switch (k)
		case {1}
			option.baseBegin=varargin{k};
		case {2}
			option.baseEnd=varargin{k};
		case {3}
			option.peakBegin=varargin{k};
		case {4}
			option.peakEnd=varargin{k};
		case {5}
			option.meanN=varargin{k};
		case {6}
			option.dir=varargin{k};
		case {6}
			option.plotFlag=varargin{k};
		case {7}
			option.fitFlag=varargin{k};
		case {8}
			option.thres=varargin{k};
		case {9}
			option.thresFlag=varargin{k};
		end
	elseif ischar(varargin{k}) && isnumeric(varargin{k+1})
		setfield(option,varargin{k},varargin{k+1});
	else
		error(sprintf('%s-th input argument invalid/unsupported',k))
	end
    end	
end

if (evtpos~=round(evtpos))
	warning('event_pos is not integer')
end

N = length(evtpos);
results.label = {'tix','baseline','baseSD','-','peak','tPeak', ...
		 't20Real', 't80Real', 't50AInt', 't50AReal', 't50BInt', 't50BReal', 't0Real', ...
		 'tMaxSlopeRiseReal', 'maxSlopeRiseReal', 'tThreshold', 'yThreshold', ...
		 'tMaxSlopeDecayReal', 'maxSlopeDecayReal', ...
		};

numFixLabels = length(results.label);

if (option.fitFlag>0)
if exist('OCTAVE_VERSION','builtin') && ~exist('lsqcurvefit','file')
	try
		pkg load optim
	end
end
end

if (option.fitFlag==1)
	results.label(end+[1:3]) = {'a','offset','tau [s]'};
	option.fitfun = @(p, x) p(1) * exp (-x*p(3)) + p(2);	% A * exp(-t * invTau) + B
	% default.fitfunJ = @(p, x) ((exp (-x*p(3)) - x*p(1)*exp (-x*p(3))) + 1);	% A * exp(-t*invTau) + B
	results.fitResults=repmat(NaN,N,3);
elseif (option.fitFlag==2)
	results.label(end+[1:4]) = {'A1','A2', 'tau1 [s]', 'tau2 [s]'};
	option.fitfun = @(p, x) (p(1) * exp (-x*p(3)) + p(2) * exp (-x*p(4)));	% A1 * exp(-t*invTau1) + A2 * exp(-t*invTau2)
	results.fitResults=repmat(NaN,N,4);
	results.weightedTau=repmat(NaN,N,1);
elseif (option.fitFlag==3)
	results.label(end+[1:5]) = {'A1','A2', 'tau1 [s]', 'tau2 [s]','offset'};
	option.fitfun = @(p, x) (p(1) * exp (-x*p(3)) + p(2) * exp (-x*p(4)) + p(5));	% A1 * exp(-t*invTau1) + A2 * exp(-t*invTau2) + offset
	results.fitResults=repmat(NaN,N,5);
	results.weightedTau=repmat(NaN,N,1);
elseif (option.fitFlag==0)
	;
else
	warning('fitting function not defined')
end

results.data = repmat(NaN,N,length(results.label)); 


if 0, (option.meanN~=1)
	% triangular window
	n=option.meanN;
	sdata=filter(ones(n,1),n,data(end:-1:1));
	sdata=filter(ones(n,1),n,sdata(end:-1:1));
elseif (option.meanN~=1)
	% this is closer to what stimfit is doing 
	n=option.meanN;
	n2=floor(n/2);
	sdata=filter(ones(n,1),n,data);
	sdata=[sdata(n2+1:end);repmat(NaN,n2,size(sdata,2))];
else
	sdata=data;
end

markerSize = 16;
lineWidth  = 2;

[base,baseSZ]=trigg(data, evtpos, option.t1, option.t2);
results.datatype='STIMFIT';
results.samples=reshape(base,baseSZ);
results.time = [option.t1:option.t2]'/Fs;
results.Fs = Fs;
results.options=option; 

[base,baseSZ]=trigg(data,evtpos,option.baseBegin, option.baseEnd);
[peakRegion,peakSZ]=trigg(sdata,evtpos,option.peakBegin, option.peakEnd);

results.baselineMean   = mean(reshape(base,baseSZ),2);
results.baselineMedian = median(reshape(base,baseSZ),2);
results.baselineStdDev = std(reshape(base,baseSZ),[],2);
results.baselineIQR    = iqr(reshape(base,baseSZ),2);

if option.baseFlag==0
	results.data(:,2)  = results.baselineMedian;
	results.data(:,3)  = results.baselineStdDev;
	Baseline = results.baselineMean;
	BaseSD   = results.baselineStdDev;
else
	results.data(:,2)  = results.baselineMedian;
	results.data(:,3)  = results.baselineIQR;
	Baseline = results.baselineMedian;
	BaseSD   = results.baselineIQR;
end

peakRegion = reshape(peakRegion,peakSZ);
if option.dir==1,
	[peak, tPeak2] = max(peakRegion,[],2);
	peak = peak - Baseline;
elseif option.dir==-1,
	[peak, tPeak2] = min(peakRegion,[],2);
	peak = peak - Baseline;
elseif option.dir==0,
	tmp = peakRegion - repmat(Baseline,[1,peakSZ(2),1]);
	[peak, tPeak2] = max(abs(tmp), [], 2);
	for k1=1:size(tmp,1)
	for k3=1:size(tmp,3)
		peak(k1,1,k3) = tmp(k1,tPeak2(k1,k3),k3);
	end; end; 
end

results.data(:,1) = evtpos;
results.data(:,5) = peak;
% FIXME: works only correctly only for single channel, 
%    add additional checks, or fix this. 
results.data(:,6) = tPeak2(:) + option.peakBegin + evtpos(:) - 1;

results.PeakAmplitude=peak;
results.peakTime = (tPeak2(:) + option.peakBegin - 1)/Fs;

for k=1:N;
	ix = evtpos(k);

	if ( ((option.t1+ix) < 1) || ((option.t2+ix)>size(data,1)) ) continue; end; 

	base = Baseline(k);
	results.data(k,2) = Baseline(k);
	results.data(k,3) = BaseSD(k);

	peakRegion = sdata( max(option.peakBegin+ix,1) : min(ix+option.peakEnd,size(data,1)) ) - base;

	if option.dir && any(xor(peak < 0, option.dir < 0))
		warning('microStimfit: peak has wrong direction  - this is strange'); 
	end

	peakCor     = abs(peak(k));
	peakRegion2 = sign(peak(k))*peakRegion;
	peakRegion3 = sign(peak(k))*diff(peakRegion);

	t20Int  = find( diff(peakRegion2 >  0.2 * peakCor) > 0); t20Int (t20Int  > tPeak2(k))=[];
	t80Int  = find( diff(peakRegion2 >  0.8 * peakCor) > 0); t80Int (t80Int  > tPeak2(k))=[];
	t50AInt = find( diff(peakRegion2 >  0.5 * peakCor) > 0); t50AInt(t50AInt > tPeak2(k))=[];
	t50BInt = find( diff(peakRegion2 >= 0.5 * peakCor) < 0); t50BInt(t50BInt < tPeak2(k))=[];

	t20Real = ix + option.peakBegin - 1; 
	t80Real = ix + option.peakBegin - 1;
	t50AReal= ix + option.peakBegin - 1;
	t50BReal= ix + option.peakEnd   - 1;

	if (length(t20Int)>0)
		t20Real  = t20Int  + (0.2*peakCor-peakRegion2(t20Int))./peakRegion3(t20Int) + ix + option.peakBegin-1; 
	end;
	if (length(t80Int)>0)
		t80Real  = t80Int  + (0.8*peakCor-peakRegion2(t80Int))./peakRegion3(t80Int) + ix + option.peakBegin-1;
	end;
	if (length(t50AInt)>0)
		t50AReal = t50AInt + (0.5*peakCor-peakRegion2(t50AInt))./peakRegion3(t50AInt) + ix + option.peakBegin-1;
	end;
	if (length(t50BInt)>0)
		t50BReal = t50BInt + (0.5*peakCor-peakRegion2(t50BInt))./peakRegion3(t50BInt) + ix + option.peakBegin-1;
	end;
	t0Real = mean(t20Real) - (mean(t80Real)-mean(t20Real))/3;

	maxSlopeRiseReal  = max(peakRegion3);
	tMaxSlopeRiseReal = find(peakRegion3==maxSlopeRiseReal) + ix + option.peakBegin - 0.5;
	maxSlopeDecayReal = min(peakRegion3);
	tMaxSlopeDecayReal= find(peakRegion3==maxSlopeDecayReal) + ix + option.peakBegin - 0.5;

	t50AReal = mean(t50AReal);
	tMaxSlopeRiseReal(2:end) = [];	% First only
	t0Real(2:end) = [];		% First only

	results.data(k,7)  = mean(t20Real);
	results.data(k,8)  = mean(t80Real);
	results.data(k,10) = mean(t50AReal);
	results.data(k,12) = mean(t50BReal);
	results.data(k,13) = mean(t0Real);
	results.data(k,14) = mean(tMaxSlopeRiseReal);
	results.data(k,15) = mean(maxSlopeRiseReal);
	results.data(k,18) = mean(tMaxSlopeDecayReal);
	results.data(k,19) = mean(maxSlopeDecayReal);

	% peakRegion4
	if option.thresFlag>0
		threshold = option.thres;
	else
		threshold = option.thres * maxSlopeRiseReal;
	end; 

	myList = find( diff(peakRegion3 > threshold) > 0 );
	if numel(myList), 
		myList=myList(1);
		tThreshold = myList + ix + option.peakBegin + .5;
		yThreshold = sum(peakRegion(myList + [0;1]))/2;
		results.data(k,16:17) = [tThreshold,yThreshold];
	else 
		tThreshold = NaN;
		yThreshold = NaN;
	end

	%traceCor   = data - base; %CHECK_MMA%
	tix  = [option.t1:option.t2] + ix;
	d    = data(tix);
	myf  = @(x)interp1([option.t1:option.t2], d, x, 'linear');
	YLIM = [min(d);max(d)];

	t = [0:option.fitEnd]';
	if (option.fitFlag==1) % && (abs(peak)>0))
		% results.label(end+[1:3]) = {'a','offset','tau [s]'};
		% default.fitfun = @(p, x) p(1) * exp (-x*p(3)) + p(2);	% A * exp(-t * invTau) + B

		fitResult = [];
		fExp      = option.fitfun;

		pInit = [peak(k), base, Fs/(mean(t50BReal) - mean(t50AReal))];
		LB = [-inf, -inf, 0];
		UB = [+inf, +inf, Fs];

		%% settings = optimset ('lbound',[-inf,-inf,0]','ubound',[inf,inf,Fs*10]');
		%opts = optimset('Algorithm','levenberg-marquardt');
		%opts = optimset('Jacobian',default.fitfunJ);
		%opts = optimset('Algorithm','levenberg-marquardt','Jacobian',default.fitfunJ);
		try
			t = [0:option.fitEnd]';
			decay = data(results.peakTime(k)*Fs+evtpos(k) + t);
			[fitResult, RESNORM, RESIDUAL, EXITFLAG, OUTPUT, LAMBDA, JACOBIAN] = lsqcurvefit (option.fitfun, pInit', t/Fs, decay, LB, UB);
			results.fitResults(k,1:3) = [fitResult(1), fitResult(2), 1/fitResult(3)];
			results.data(k, numFixLabels+[1:3])    = [fitResult(1), fitResult(2), 1/fitResult(3)];
		catch
			warning('fitting failed; optimization toolbox or optim package missing')
			fitResult=[];
		end

	elseif (option.fitFlag==2) % && (abs(peak)>0))
		% results.label(end+[1:3]) = {'A1','A2', 'tau1 [s]', 'tau2 [s]'};
		% default.fitfun = @(p, x) p(1) * exp (-x*p(3)) + p(2) * exp (-x*p(4));	% A1 * exp(-t*invTau1) + A2 * exp(-t*invTau2)

		fitResult = [];
		fExp      = option.fitfun;
		invtau1   = Fs/(mean(t50BReal) - mean(t50AReal));
		pInit = [peak(k), peak(k)/2, invtau1, invtau1/2];
		LB = [-inf, -inf, 0, 0];
		UB = [+inf, +inf, Fs, Fs];

		try
			t = [0:option.fitEnd-results.peakTime(k)*Fs]';
			decay = data(results.peakTime(k)*Fs+evtpos(k) + t) - base;
			[fitResult, RESNORM, RESIDUAL, EXITFLAG, OUTPUT, LAMBDA, JACOBIAN] = lsqcurvefit (option.fitfun, pInit', t/Fs, decay, LB, UB);
			results.fitResults(k,1:4) = [fitResult(1), fitResult(2), 1/fitResult(3), 1/fitResult(4)];
			% (A1/invTau1 + A2/invTau2)/(A1+A2)
			results.weightedTau(k,1)  = sum(fitResult(1:2)./fitResult(3:4)) / sum(fitResult(1:2));
		catch
			warning('fitting failed; optimization toolbox or optim package missing')
			fitResult=[];
		end

	elseif (option.fitFlag==3) % && (abs(peak)>0))
		% results.label(end+[1:3]) = {'A1','A2', 'tau1 [s]', 'tau2 [s]'};
		% default.fitfun = @(p, x) p(1) * exp (-x*p(3)) + p(2) * exp (-x*p(4)) + offset;	% A1 * exp(-t*invTau1) + A2 * exp(-t*invTau2) + offset

		fitResult = [];
		fExp      = option.fitfun;
		invtau1   = Fs/(mean(t50BReal) - mean(t50AReal));
		pInit = [peak(k), peak(k)/2, invtau1, invtau1/2, base];
		LB = [-inf, -inf, 0, 0, -inf];
		UB = [+inf, +inf, Fs, Fs, +inf];

		try
			t = [0:option.fitEnd-results.peakTime(k)*Fs]';
			decay = data(results.peakTime(k)*Fs+evtpos(k) + t);
			[fitResult, RESNORM, RESIDUAL, EXITFLAG, OUTPUT, LAMBDA, JACOBIAN] = lsqcurvefit (option.fitfun, pInit', t/Fs, decay, LB, UB);
			results.fitResults(k,1:5) = [fitResult(1), fitResult(2), 1/fitResult(3), 1/fitResult(4), fitResult(5)];
			% (A1/invTau1 + A2/invTau2)/(A1+A2)
			results.weightedTau(k,1)  = sum(fitResult(1:2)./fitResult(3:4)) / sum(fitResult(1:2));
		catch
			warning('fitting failed; optimization toolbox or optim package missing')
			fitResult=[];
		end
	end

	if (option.plotFlag),
		h = plot([option.t1:option.t2]/Fs, data([option.t1:option.t2]+ix,:),'-', ...
			(t0Real-ix)/Fs, base,'kx',...
			(t20Real-ix)/Fs, myf(t20Real-ix),'rx', ...
			(t80Real-ix)/Fs, myf(t80Real-ix),'gx', ...
			(t50AReal-ix)/Fs, myf(t50AReal-ix),'b.', ...
			(t50BReal-ix)/Fs, myf(t50BReal-ix),'b.', ...
			(tMaxSlopeRiseReal-ix)/Fs, myf(tMaxSlopeRiseReal-ix),'c.', ...
			(tMaxSlopeDecayReal-ix)/Fs, myf(tMaxSlopeDecayReal-ix),'c.', ...
			(tThreshold-ix)/Fs, myf(tThreshold-ix), 'm.',  ...
			[option.t1;option.t2]'/Fs, base*[1;1],'g-',  ...
			[option.baseBegin;option.baseEnd]'/Fs, base*[1;1],'g-', ...
			[option.t1;option.t2]/Fs, [1;1]*peak(k)+base,'r-',  ...
			[option.peakBegin;option.peakEnd]/Fs, [1;1]*peak(k)+base,'r-', ...
			results.peakTime(k)*[1;1], YLIM,'r-',  ...
			results.peakTime(k)+t/Fs, myf(t+results.peakTime(k)*Fs),'b-');

		set(h(2),'markersize',markerSize)
		set(h(3),'markersize',markerSize)
		set(h(4),'markersize',markerSize)
		set(h(5),'markersize',markerSize)
		set(h(6),'markersize',markerSize)
		set(h(7),'markersize',markerSize)
		set(h(8),'markersize',markerSize)

		set(h(end-4),'linewidth',lineWidth);
		set(h(end-2),'linewidth',lineWidth);
		set(h(end-1),'linewidth',lineWidth);
		set(h(end),'linewidth',lineWidth);
		title(sprintf('%i/%i',k,N))
		
		if any(option.fitFlag==[1:3]) && ~isempty(fitResult)
			hold on;
			y=option.fitfun(fitResult,t/Fs);
			if (option.fitFlag==2) y=y+base; end
			plot(t/Fs+results.peakTime(k),y,'c','linewidth',4);
			%plot(x/Fs,fitResult(x),'c','linewidth',4);
			hold off
			drawnow
			if option.plotFlag > 0,
				fprintf(1,'Press any key');
				pause;
			end
		end
	end
end

results.MaxSlope = results.data(:,15);		% for backwards compatibility
results.MaxRiseSlope = results.data(:,15);	% same as .MaxSlope
results.MaxDecaySlope = results.data(:,19);
results.BaseLine = Baseline;
results.BaseLineDev = BaseSD;

results.OnsetTime = (results.data(:,13)-evtpos)/Fs;
results.RiseTime = diff(results.data(:,[7,8]),[],2)/Fs;
results.HalfWidth = diff(results.data(:,[10,12]),[],2)/Fs;
results.PeakTime = diff(results.data(:,[1,6]),[],2)/Fs;
results.MaxSlopeTime   = diff(results.data(:,[1,14]),[],2)/Fs;	% for backwards compatibility
results.tMaxRiseSlope  = diff(results.data(:,[1,14]),[],2)/Fs;	% same  .MaxSlopeTime
results.tMaxDecaySlope = diff(results.data(:,[1,18]),[],2)/Fs;
results.thresholdTime = diff(results.data(:,[1,16]),[],2)/Fs;


