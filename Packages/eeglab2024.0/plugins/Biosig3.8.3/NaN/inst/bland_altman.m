function RES = bland_altman(data,group,arg3)
% BLAND_ALTMANN shows the Bland-Altman plot of two columns of measurements
%   and computes several summary results. 
%
%   bland_altman(m1, m2 [,group])
%   bland_altman(data [, group])
%   R = bland_altman(...)
% 
%   m1,m2 are two colums with the same number of elements
%	containing the measurements. m1,m2 can be also combined 
%       in a single two column data matrix. 
%   group [optional] indicates which measurements belong to the same group
%	This is useful to account for repeated measurements.  
%
%
% References:
% [1] JM Bland and DG Altman, Measuring agreement in method comparison studies. 
%       Statistical Methods in Medical Research, 1999; 8; 135. 
%       doi:10.1177/09622802990080204
% [2] P.S. Myles, Using the Bland– Altman method to measure agreement with repeated measures
%	British Journal of Anaesthesia 99(3):309–11 (2007)
%	doi:10.1093/bja/aem214

%	$Id$
%	Copyright (C) 2010,2011 by Alois Schloegl <alois.schloegl@gmail.com>	
%       This function is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
% of the  License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street - Fifth Floor, Boston, MA 02110-1301, USA.

if nargin<2, group = []; end;
if nargin<3, arg3 = []; end;

if (size(data,2)==1)
	data  = [data, group];
	group = arg3;
end; 		


D = data * [1;-1];
M = data * [1;1]/2;

RES.corrcoef = corrcoef(data(:,1),data(:,2),'spearman');
[REs.cc,RES.p] = corrcoef(M,D,'spearman');
if (RES.p<0.05)
	warning('A regression model according to section 3.2 [1] should be used');
	%% TODO: implement support for this type of data. 	
	RES.a = [ones(size(data,1),1),D]\M;
	RES.b = [ones(size(data,1),1),M]\D;
end; 	

if isempty(group)
	G = [1:size(data,1)]';
	m = ones(size(data,1),1);
	d = D;
	RES.Bias = mean(d,1);
	RES.Var  = var(d);
	
elseif ~isempty(group)
	%% TODO: this is not finished  	
	warning('analysis of data with repetitions is experimental - it might yield incorrect results - you are warned.!')
	[G,I,J] = unique (group);
	R       = zeros(size(data));
	m       = repmat(NaN,length(G),1);
	n       = repmat(NaN,length(G),1);
	d       = repmat(NaN,length(G),1);
	d2      = repmat(NaN,length(G),1);
	data2   = repmat(NaN,length(G),size(data,2));
	SW2     = repmat(NaN,length(G),size(data,2));
	for i   = 1:length(G),
		ix         = find(group==G(i));
		n(i)       = length(ix);
%		IX((i-1)*N+1:i*N) = ix(ceil(rand(N,1)*n(i)));		
		
		[R(ix,:), data2(i,:)] = center(data(ix,:),1);
		d(i)       = mean(D(ix,:),1);
		m(i)       = mean(M(ix,:),1);
		d2(i)      = mean(D(ix,:).^2,1);
		RES.SW2(i,:)   = var(data(ix,:),[],1);
		RES.avg(i,:)   = mean(data(ix,:),1);
	end;

	W = 1./n(J);
	RES.SSW = sumskipnan(R.^2,1,W);
	RES.SSB = var(data,[],1,W)*sum(W)*(sum(W)-1);
	RES.sigma2_w = RES.SSW/(sum(W)*(length(G)-1));
	RES.sigma2_u = RES.SSB/(sum(W)*(length(G)-1)) - RES.sigma2_w/(length(G));
	RES.group  = bland_altman(data2); 	% FIXME: this plot shows incorrect interval, it does not account for the group/repeated samples.
	RES.repeatability_coefficient1 = 2.77*sqrt(var(R,1,1)); 	% variance with factor group removed
	RES.repeatability_coefficient = 2.77*sqrt(mean(SW2,1)); 	% variance with factor group removed
	RES.std_d_ = std(d);
	RES.std_D_ = std(D);
	RES.std_m_ = std(m);	
	
	RES.n = n;
	return; 

	D = d; 
	M = m;
%	RES.sigma2_dw = 

	RES.Bias = mean(d,1,[],n);
end; 


plot(M,D,'o', [min(M),max(M)]', [0,0]','k--', [min(M),max(M)]', [1,1,1; 0,1.96,-1.96]'*[RES.Bias;std(D)]*[1,1], 'k-');
xlabel('mean');
ylabel('difference (m1-m2)');

