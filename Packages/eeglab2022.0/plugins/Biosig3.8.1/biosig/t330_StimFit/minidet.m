function [results, traceTab, average] = minidet(trace,template,rCrit,aCrit,refract) 
% Miniature EPSC detection by template matching 
%
% Usage: 
%    [dtrace] = minidet(trace,template)
%    [dtrace,a] = minidet(trace,template)
%    [results, traceTab, average] = minidet(trace,template,rCrit,aCrit,refract)
%
%  Input: 
%    dtrace :   raw sampling data
%    template: default is a bi-exponential template with tauRise=0.5ms, 
%              tauDecay=5ms, a 10ms baseline, for data sampled with 20kHz
%          template = F([0:0.05:49.9]', 10, .5, 5);
% 	        with F defined as            
% 	   F = @(t,delta,tauRise,tauDecay) (-double(t>=delta) .* ...
%	       (exp((delta-t)/tauRise) - exp((delta-t)/tauDecay) ) / ...
%              ( (tauRise/tauDecay)^(tauDecay/(tauDecay-tauRise)) - ...
%                  (tauRise/tauDecay)^(tauRise/(tauDecay-tauRise)) ) );
%          alternative templates can be described as 
%	       template = F([0:1000/SampleRate:49.9]', delay, tauRise, tauDecay);
%    rCrit: 	default: 0.5
%    aCrit: 	default: 10
%    refract: refractory period
%  Output:
%    dtrace: raw detection trace
%    result.tEventList: event times
%    traceTab: traces from window [-200,600] samples of each event
%    average:  average trace 

%
% Copyright (C) 2020,2022 Alois Schlögl, ISTA
%
% Re-Implemented in Matlab/Octave after the Mathematica version of Peter Jonas. 
% File : Microminidet.nb
% PJ, 04.11.2020
% Version: 2
%

% References: 
% Jonas P, Major G, Sakmann B. Quantal components of unitary EPSCs at the mossy
%   fibre synapse on CA3 pyramidal cells of rat hippocampus. J Physiol. 1993
%   Dec;472:615-63.
% 
% von Kitzing E, Jonas P, Sakmann B. Quantal analysis of excitatory
%   postsynaptic currents at the hippocampal mossy fiber-CA3 pyramidal cell synapse.
%   Adv Second Messenger Phosphoprotein Res. 1994;29:235-60.
% 
% Clements JD, Bekkers JM.Detection of spontaneous synaptic events with an
%   optimally scaled template.Biophys J.1997 Jul; 73 (1) : 220 - 9. doi : 10.1016/S0006 - 3495 (97) 78062 - 7.
% 
% Pernía - Andrade AJ, Goswami SP, Stickler Y, Fröbe U, Schlögl A, Jonas P.A
%   deconvolution - based method with high sensitivity and temporal resolution for
%   detection of spontaneous synaptic currents in vitro and in vivo.Biophys J.2012
%   Oct 3; 103 (7) : 1429 - 39.

if nargin<2,
	template = [];
end
if (nargin<3) || isempty(rCrit),
	rCrit = .5;
end
if (nargin<4) || isempty(aCrit),
	aCrit = 10;
end
if nargin<5,
	refract = [];
end
if isempty(template),
	F = @(t,delta,tauRise,tauDecay) (-double(t>=delta) .* (exp((delta-t)/tauRise) - exp((delta-t)/tauDecay) ) / ( (tauRise/tauDecay)^(tauDecay/(tauDecay-tauRise)) - (tauRise/tauDecay)^(tauRise/(tauDecay-tauRise)) ) );
	sampleInt = 0.05; %1000/round(HDR.SampleRate);
	%%% PARAM 1: delta
	shift     = 10; % round(10 / sampleInt); 
	%%% PARAM 2 (length), 3 (tauRise), 4 (tauDecay)
	template  = F([0:sampleInt:49.9	]', 10, .5, 5);
	%%% PARAM 5 (rCrit), 6 (aCrit)
end;

%MMA% tracePart =  Partition[myTrace, Length[myTemplate],  1]; (* partitions the trace *) 
%OCT% tracePart=trace(hankel(1:n2,n2:n1));
n1 = length(trace);
n2 = length(template);

%MMA% detectionTrace =  fit3[myTemplate, #] & /@ tracePart; (* computes detection trace *) 
%OCT% SPY = var(tracePart,[],1);
m1   = filter(ones(n2,1), n2, trace);
m2   = filter(ones(n2,1), n2, trace.^2);
SPY  = (m2-m1.*m1)*n2/(n2-1);
%OCT% SP   = center(template)' * center(tracePart,1)/n2;
SP   = filter(template(end:-1:1)-mean(template), n2, trace);
%
SPX  = var(template);
%
a    = SP(n2:end)./SPX;
r    = SP(n2:end)./sqrt(SPX.*SPY(n2:end));

if nargin<4,
	results=r;
	traceTab=a;
	return
end
delay=find(template~=0,1)-2;

%MMA%  tEventList = (* generates event time list computing peaks above threshold *) 
% Flatten@Position[Partition[detectionTrace, 5, 1], 
%    x_ /; (x[[1, 3]] <=  x[[2, 3]] <=   x[[3, 3]] > x[[4, 3]] >  
%        x[[5, 3]]) && (x[[3, 1]] > aCrit) && (x[[3, 3]] > 
%        rCrit), {1}, Heads -> False]  + shift ; 
x = (r(1:end-4)<=r(2:end-3) & r(2:end-3)<=r(3:end-2) & r(3:end-2)>r(4:end-1) & r(4:end-1)>r(5:end) & a(3:end-2)>aCrit & r(3:end-2)>rCrit);
tEventList = find(x)'+delay;

results.tEventList = tEventList; % stimfit(...)
if ~isempty(refract)
	%MMA%  tEventListRev = Mean /@ Split[   tEventList, (Abs[#2 - #1] <=  refract) &] ;  (* split list into sublists containing close elements *)  
	ix = cumsum([1,diff(tEventList)>refract]);
	tEventListRev=zeros(ix(end),1);
	for k=1:ix(end)
		tEventListRev(k) = mean(tEventList(ix==k));
	end
	results.tEventList = tEventListRev; % stimfit(...)
end

%MMA% traceTab =  Take[myTrace, {# - 200 , # + 600 }] & /@ 
%        tEventList; (* generates table of traces around PSC onset points *) \
[x,sx]   = trigg(trace, round(results.tEventList), -delay, 3*delay);
traceTab = reshape(x,sx);

%MMA% average =  Mean /@ Transpose[ traceTab]; (* computes average PSC trace *)
average  = reshape(mean(traceTab,3),sx(1:2))';

% TODO:
%MMA% results =  stimfit[#, 1, 200, 180, 400, 5, -1, 1, 1] & /@  traceTab ; (* runs stimfit analysis on table of traces *)

