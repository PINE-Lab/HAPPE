function [SIGMA,PHI,OMEGA,m0,m1,m2,cout] = wackermann(S,UC,A,cini)
% Wackermann calculates the global field strength SIGMA,  the global
%    frequency PHI and a measure of spatial complexity OMEGA [1]. 
%    A stationary and a time-varying (adaptive) estimator is implemented 
%   
%  [SIGMA,PHI,OMEGA] = wackermann(...)
%  
%  [...] = wackermann(S,0)
%       calculates stationary Wackermann parameter 
%  [...] = wackermann(S,UC) with 0<UC<1,
%       calculates time-varying Wackermann parameter using 
%       exponential window 
%  [...] = wackermann(S,N) with N>1,
%       calculates time-varying Wackermann parameter using 
%       rectangulare window of length N
%  [...] = wackermann(S,B,A) with B>=1 oder length(B)>1,
%       calulates time-varying Wackermann parameters using 
%       transfer function B(z)/A(z) for windowing 
%
%       S       data (each channel is a column)
%       UC      update coefficient 
%       B,A     filter coefficients (window function) 
%
% Remark: estimating of Omega requires the eigenvalues, the adaptive estimator  
%   utilized the adaptive eigenanalysis method [2].   
%
% see also: TDP, BARLOW, HJORTH
%
% REFERENCE(S):
% [1] Jiri Wackermann, Towards a quantitative characterization of
%     functional states of the brain: from the non-linear methodology to the
%     global linear descriptor. International Journal of Psychophysiology, 34 (1999) 65-80.
% [2] Bin Yang, Projection approximation subspace tracking. 
%     IEEE Trans. on Signal processing, vol. 43, no. 1 jan. 2005, pp. 95-107.
% [3] Saito N., Kuginuki T., Yagyu T., Kinoshita T., Koenig T., Pascual-Marqui
%	R. D., Kochi K., Wackermann J., and Lemann D., 1998. 
%	Global, regional and local measures of complexity of multichannel EEG in acute, neurolepticnaive,
%	first-break schizophrenics. Society of Biological Psychiatry. 43:794–802.

%	$Id$
%	Copyright (C) 2004,2005,2008,2009 by Alois Schloegl <alois.schloegl@gmail.com>
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

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


[N,K] = size(S); 	% number of electrodes K, number of samples N
S=center(S,1);

if nargin<2, 
        UC = 0; 
end;
FLAG_ReplaceNaN = 0;
[N,K] = size(S); 	% number of electrodes K, number of samples N
if nargin<3;
        if UC==0,
                                
        elseif UC>=1,
                B = ones(1,UC);
                A = UC;
        elseif UC<1,
                FLAG_ReplaceNaN = 1;
                B = UC; 
                A = [1, UC-1];
        end;
else
        B = UC;    
end;

d0 = sum(S.*S,2);
d1 = sum(diff([zeros(1,K);S],[],1).^2,2);

[k1,k2] = meshgrid(1:size(S,2));
k1 = k1(:); k2 = k2(:); 
d2 = S(:,k1).*S(:,k2);  

if UC==0,
        m0 = mean(d0);
        m1 = mean(d1);
        m2 = mean(d2);
elseif 0, UC>=1
        m0 = rs(d0,UC,1);
        m1 = rs(d1,UC,1);
        m2 = rs(d2,UC,1);
else
        if FLAG_ReplaceNaN;
                d0(isnan(d0)) = 0;
                d1(isnan(d1)) = 0;
                d2(isnan(d2)) = 0;
        end;
        %m0 = mean(sumsq(S,2));
        m0 = filter(B,A,d0)./filter(B,A,real(~isnan(d0)));
        %m1 = mean(sumsq(diff(S,[],1),2));
        m1 = filter(B,A,d1)./filter(B,A,real(~isnan(d1)));
        m2 = filter(B,A,d2)./filter(B,A,real(~isnan(d2)));
end;

SIGMA = sqrt(m0/K);

PHI   = sqrt(m1./m0)/(2*pi);

OMEGA = repmat(NaN,size(m0));
if ~exist('OCTAVE_VERSION','builtin'),

	% this branch is equivalent to the next one, but faster on Matlab
	for k = 1:size(m2,1),
	        if all(isfinite(m2(k,:))), 
	                L = eig(reshape(m2(k,:), [K,K]));
	                L = L./sum(L);
	                if all(L>0)
	                        OMEGA(k) = -sum(L.*log(L));
        	        end;
	        end;
	end;        

elseif exist('OCTAVE_VERSION','builtin'),

	% this branch is equivalent to the previous one, but faster on octave
	rows_m2 = size(m2, 1);
	m3 = permute (reshape (m2, [rows_m2, K, K]), [2, 3, 1]);
	idx = all (isfinite (m2), 2);
	t = cellfun (@eig, mat2cell (m3 (:, :, idx), K, K, ones(1, sum(idx))),'UniformOutput', false);
	t = [t{:}];
	idx2 = all(t>0);
	t = t(:,idx2) ./ [ones(K,1) * sum(t(:,idx2))];
	t = sum (t .* log (t));
	idx = find(idx); 
	OMEGA(idx(idx2)) = t;

else 
%% THIS CODE IS CURRENTLY USED FOR TESTING ONLY 


	% alternative adaptive estimation of eigenvalues based on [2] 
	% implemented by Carmen Vidaurre
	

	if (nargin<4) %%,isempty(fieldnames(cini)))
		Cini=eye(size(S,2));
		W=Cini;
		P=inv(W);
		cxyt=W*Cini;
		cxxt=Cini;
		icxxt=inv(Cini);
		icyyt=P;
	else
		cxyt=cini.cxyt;
		cxxt=cini.cxxt;
		icxxt=inv(cini.cxxt);
		icyyt=cini.icyyt;
	end;
	wt=cxyt*icyyt;
	[wt]=orth(wt);
	d=zeros(size(S));
	cont=0;
	for j=1:size(S,1)
    		xt=S(j,:);
		if all(isfinite(isnan(xt)))
			yt=wt'*xt';
			cxyt=(1-UC)*cxyt+UC*xt'*yt';
			vt=icxxt*xt';
			icxxt=1/(1-UC)*(icxxt-(vt*vt'/(((1-UC)/UC)+xt*vt)));
			icxxt=(icxxt+icxxt')/2;
			vt=icyyt*yt;
			icyyt=1/(1-UC)*(icyyt-(vt*vt'/(((1-UC)/UC)+yt'*vt)));
			icyyt=(icyyt+icyyt')/2;
			wt=cxyt*icyyt;
			[wt]=orth(wt);
			d(j,:)=diag(wt'*icxxt*wt)';
			[d(j,:) ii]=sort(d(j,:),2,'descend');
			wt=wt(:,ii);
			L=d(j,:);
			L=L./sum(L);
			cxxt=(1-UC)*cxxt+UC*xt'*xt;
    			if all(L>0)
			    	OMEGA(j) = -sum(L.*log(L));
			else
				cont=cont+1;
			end;
		elseif (j>1) 	
			OMEGA(j)=OMEGA(j-1);	
		end; 		
	end;
	cout.cxyt=cxyt;
	cout.cxxt=cxxt;
	cout.icyyt=icyyt;
end; 
