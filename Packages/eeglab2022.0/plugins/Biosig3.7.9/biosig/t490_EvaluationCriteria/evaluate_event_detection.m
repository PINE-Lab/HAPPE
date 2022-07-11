function RES = evaluate_event_detection(D, trig)
% EVALUATE_EVENT_DETECTION quantifizes the performance of an 
%  event detection algorithm producing a detection trace D. 
%  The detection trace is compared to the reference trigger values, 
%  
%  A sample-based evaluation is used for measuring the performance
%  and for optimizing unknown parameters, like delay factor (TAU),
%  window length (WINLEN), and the detection threshold (TH). 
%
%  Based on the found detection threshold, parameter of an event-based
%  evaluation are reported. 

%       $Id$
%       Copyright (c) 2011 by  Alois Schloegl <alois.schloegl@gmail.com>
%       This is part of the BIOSIG-toolbox http://biosig.sf.net/
%
% Biosig is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% version 3 of the License, or (at your option) any later version.
%
% Biosig is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Library General Public License for more details.
%
% You should have received a copy of the GNU Library General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc., 59 Temple Place - Suite 330,
% Boston, MA  02111-1307, USA.
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rank sample data, 
% a simple algorithm not resolving ties is used.  
% d = ranks(D);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = zeros(size(D));
[TH, I] = sort(D,1);
for k=1:size(D,2)
	[tmp,d(:,k)]  = sort(I(:,k),1);     % d yields the rank of each element         
%	x(:,k) = r(I(:,k));
end;
d(isnan(D)) = nan;
N = sum(~isnan(D));    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find onset and offset of optimal detection window. 
% This is a 2-dimensional optimization problem 
% with window size larger than 0.
% For optimization, the correlation on the ranked samples is
% used because it is faster than AUC or optimum Kappa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  initialize
IX1 = 0; 
IX2 = 0;
%% reference
d = zscore(d); 

r = zeros(size(D));
r(trig) = 1;
R0 = corrcoef(d,r);
flag_improved = 1; 
while (flag_improved) 
	flag_improved = 0; 

	IX1 = IX1+1;
	t = trig(trig <= size(r,1) - IX1) + IX1;
	r(t) = 1;
	R1 = mean(d.*zscore(r)); %corrcoef(d,r);
	d11 = R1-R0;
	if R1 > R0,
		R0 = R1; 
		flag_improved = 1; 
		d12 = 0;  
	else 
		% no improvement, revert. 
		r(t) = 0;
		IX1 = IX1-1;
		% go opposite direction 
		t = trig(trig <= size(r,1) - IX1) + IX1;
		r(t) = 0;
		IX1 = IX1-1;
		R1  = mean(d.*zscore(r)); %corrcoef(d,r);
		d12 = -(R1 - R0);
		if R1 > R0,
			R0 = R1; 
			flag_improved = 1; 
		else
			% no improvement, revert. 
			IX1 = IX1+1;
			r(t) = 1;
		end 
	end		

	IX2 = IX2-1;
	t = trig(trig > -IX2) + IX2;
	r(t) = 1;
	R2 = mean(d.*zscore(r)); %corrcoef(d,r);
	d21 = R2 - R0;
	if R2 > R0,
		R0 = R2;
		flag_improved = 1; 
		d22 = 0;  
	else 
		% no improvement, revert. 
		r(t) = 0;
		IX2 = IX2+1;
		% go opposite direction 
		t = trig(trig > -IX2) + IX2;
		r(t) = 0;
		IX2 = IX2+1;
		R2  = mean(d.*zscore(r)); %corrcoef(d,r);
		d22 = -(R2 - R0);
		if R2 > R0,
			R0  = R2; 
			flag_improved = 1; 
		else
			% no improvement, revert. 
			IX2 = IX2-1;
			r(t) = 1;
		end 
	end	
%[IX2,IX1],R0,[d11,d21,d22,d12],
	if all([d11,d21,d22,d12] <= 0) break; end; 
	if IX2 > IX1, break end; 	% no convergence
end;
RES.WIN = [IX2,IX1];
RES.corrcoef = R0; 
if RES.WIN(2)<RES.WIN(1), return; end; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ROC analysis: get AUC, and Threshold with maximum kappa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[TH,I] = sort(D); % already done above 
for k=1:size(D,2)
	[tmp,d(:,k)]  = sort(I(:,k),1);     % d yields the rank of each element         
	x(:,k) = r(I(:,k));
end;
d(isnan(D)) = nan;


FNR = cumsum(x==1)/sum(x==1);
TPR = 1-FNR;

TNR = cumsum(x==0)/sum(x==0);
FPR = 1-TNR;

FN = cumsum(x==1);
TP = sum(x==1)-FN;

TN = cumsum(x==0);
FP = sum(x==0)-TN;

ACC = (TP+TN)./N;

%%% compute Cohen's kappa coefficient
N = size(d,1);
% H = [TP,FN;FP,TN];
p_i = [TP+FP,FN+TN];%sum(H,1);
pi_ = [FP+FN,FP+TN];%sum(H,2)';
pe  = sum(p_i.*pi_,2)/(N*N);  % estimate of change agreement
kap = (ACC-pe)./(1-pe);

% area under the ROC curve
RES.AUC = -diff(FPR)' * (TPR(1:end-1)+TPR(2:end))/2;

% Youden index
YI = (TP+TN)/N;

RES.YI = YI; 
RES.ACC = ACC; 
RES.KAPPA = kap; 
RES.TH = TH;

[tmp,ix] = max(kap);
RES.THRESHOLD.maxKAPPA = TH(ix); 

% find optimal threshold 
[tmp,ix] = max(YI);
RES.THRESHOLD.maxYI = TH(ix);
[tmp,ix] = max(kap);
RES.THRESHOLD.maxKAPPA = TH(ix); 
[tmp,ix] = max(ACC);
RES.THRESHOLD.maxACC = TH(ix); 
RES.H = [TP(ix),FN(ix);FP(ix),TN(ix)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample-based analysis 
%    compute the confusion matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1 = D > RES.THRESHOLD.maxKAPPA; 
v  = ~isnan(D);

[ix,iy] = meshgrid(trig, IX2:IX1);
ix = ix(:) + iy(:);
ix = ix(0 < ix & ix <= N);
x0 = full(sparse(ix, 1, 1, N, 1)>0);

RES.d.sum = [sum(x0),sum(x1),size(x0),size(x1),sum(v)];

D = histo4([x0(v),x1(v)]);
RES.D = D; 
RES.K = kappa(full(sparse(D.X(:,1)+1, D.X(:,2)+1, D.H, 2, 2)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Event-based analysis 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% see spikerror.m from analysis for JK

%% TODO: support for window with defined onset and offset


TP = 0; 
%TN = 0;
FP = 0;
FN = 0;

t = 0; 
k1 = 1; k0=1;

t0 = trig + IX2;
t1 = find(diff([0;x1])>0);
winlen = IX1-IX2;
t0(end+1)=inf;
t1(end+1)=inf;
while (t1(k1)<inf && t0(k0)<inf)
        t = min(t0(k0), t1(k1)); 
        
        if t0(k0) + winlen < t1(k1), 
                k0=k0+1;
                FN=FN+1;
        elseif t1(k1) + winlen < t0(k0), 
                k1=k1+1;
                FP=FP+1;
        else 
                k0=k0+1;
                k1=k1+1;
                TP=TP+1;        
        end; 
end; 

RES.EVT.H = [TP,FP;FN,NaN];

