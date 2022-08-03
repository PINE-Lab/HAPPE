function [kap,se,H,z,p0,SA,R]=kappa(d,c,arg3,w)
% KAPPA estimates Cohen's kappa coefficient
%   and related statistics 
%
% [...] = kappa(d1,d2);
%	NaN's are handled as missing values and are ignored
% [...] = kappa(d1,d2,'notIgnoreNAN');
%	NaN's are handled as just another Label.
% [kap,sd,H,z,ACC,sACC,MI] = kappa(...);
% X = kappa(...);
%
% d1    data of scorer 1 
% d2    data of scorer 2 
%
% kap	Cohen's kappa coefficient point
% se	standard error of the kappa estimate
% H	Concordance matrix, i.e. confusion matrix
% z	z-score
% ACC	overall agreement (accuracy) 
% sACC	specific accuracy 
% MI 	Mutual information or transfer information (in [bits])
% X 	is a struct containing all the fields above
%       For two classes, a number of additional summary statistics including 
%         TPR, FPR, FDR, PPV, NPF, F1, dprime, Matthews Correlation coefficient (MCC) or
%	Phi coefficient (PHI=MCC), Specificity and Sensitivity, Youden index (YI)
%       are provided. Note, the positive category must the larger label (in d and c), otherwise 
%       the confusion matrix becomes transposed and the summary statistics are messed up. 
%
%
% Reference(s):
% [1] Cohen, J. (1960). A coefficient of agreement for nominal scales. Educational and Psychological Measurement, 20, 37-46.
% [2] J Bortz, GA Lienert (1998) Kurzgefasste Statistik f|r die klassische Forschung, Springer Berlin - Heidelberg. 
%        Kapitel 6: Uebereinstimmungsmasze fuer subjektive Merkmalsurteile. p. 265-270.
% [3] http://www.cmis.csiro.au/Fiona.Evans/personal/msc/html/chapter3.html
% [4] Kraemer, H. C. (1982). Kappa coefficient. In S. Kotz and N. L. Johnson (Eds.), 
%        Encyclopedia of Statistical Sciences. New York: John Wiley & Sons.
% [5] http://ourworld.compuserve.com/homepages/jsuebersax/kappa.htm
% [6] http://en.wikipedia.org/wiki/Receiver_operating_characteristic

%    Copyright (c) 1997-2021 by Alois Schloegl <alois.schloegl@gmail.com>
%    This function is part of the NaN-toolbox
%    http://pub.ist.ac.at/~schloegl/matlab/NaN/
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


mode.ignoreNAN = 1; 
kk = [];
if nargin>2
	if ischar(arg3)
		if strcmpi(arg3,'notIgnoreNAN')
			mode.ignoreNAN = 0; 
		 end
	else 
		kk = arg3; 
	end
end; 		 
if nargin<4
	w = [];
end; 	

if nargin>1,
	d = d(:);
	c = c(:);
	
	tmp    = [d;c];
	maxtmp = max(tmp);
	tmp(isnan(tmp)) = maxtmp+1;
	[X.Label,i,j]   = unique(tmp);
	c = j(1+numel(d):end);
	d = j(1:numel(d));
	kk = max(j);
	maxCLASS = kk - any(tmp>maxtmp);

	if mode.ignoreNAN,
		if any(j > maxCLASS)
%			fprintf(2,'Warning KAPPA: some elements are NaN. These are handled as missing values and are ignored.\n');
%			fprintf(2,'If NaN should be handled as just another label, use kappa(..,''notIgnoreNaN'').\n');
			ix = find((c<=maxCLASS) & (d<=maxCLASS));
			d = d(ix); c=c(ix);
			if ~isempty(w), w = w(ix); end; 
			kk = kk - 1;
		end;
		X.Label(X.Label>maxtmp) = []; 
	else 
		X.Label(X.Label>maxtmp) = NaN; 
	end;

	if isempty(w)
		H = full( sparse (d, c, 1, kk, kk) );
	elseif ~isempty(w),
		H = full( sparse (d, c, w, kk, kk) );
	end;

else
	X.Label = 1:min(size(d));
    	H = d(X.Label,X.Label);

end;

s = warning; 
warning('off');

N = sum(H(:)); 
p0  = sum(diag(H))/N;  %accuracy of observed agreement, overall agreement 
%OA = sum(diag(H))/N);

p_i = sum(H,1);
pi_ = sum(H,2)';

SA  = 2*diag(H)'./(p_i+pi_); % specific agreement 

pe  = (p_i*pi_')/(N*N);  % estimate of change agreement

px  = sum(p_i.*pi_.*(p_i+pi_))/(N*N*N);

%standard error 
kap = (p0-pe)/(1-pe);
sd  = sqrt((pe+pe*pe-px)/(N*(1-pe*pe)));

%standard error 
se  = sqrt((p0+pe*pe-px)/N)/(1-pe);
if ~isreal(se),
	z = NaN;
else
        z = kap/se;
end

if ((1 < nargout) && (nargout<7)) 
	warning(s);	
	return; 
end; 

% Nykopp's entropy
pwi = sum(H,2)/N;                       % p(x_i)
pwj = sum(H,1)/N;                       % p(y_j)
pji = H./repmat(sum(H,2),1,size(H,2));  % p(y_j | x_i) 
R   = - sumskipnan(pwj.*log2(pwj)) + sumskipnan(pwi'*(pji.*log2(pji)));

if (nargout>1), return; end; 

X.kappa = kap; 
X.kappa_se = se; 
X.data = H;
X.H    = X.data;
X.z    = z; 
X.ACC  = p0; 
X.sACC = SA;
X.MI   = R;
X.SSE  = sum(X.data(:))-trace(X.data);
X.datatype = 'confusion';

if length(H)==2,
	% see http://en.wikipedia.org/wiki/Receiver_operating_characteristic
  	% Note that the confusion matrix used here uses more positive values in
	% the 2nd row and column, moreover the true values are indicated by
	% rows (transposed). Thus, in summary H(1,1) and H(2,2) are exchanged 
	% as compared to the wikipedia article.  
	X.TP  = H(2,2);
	X.TN  = H(1,1);
	X.FP  = H(1,2);
	X.FN  = H(2,1);
	X.FNR = H(2,1) / sum(H(2,:));
	X.FPR = H(1,2) / sum(H(1,:));
	X.TPR = H(2,2) / sum(H(2,:));
	X.PPV = H(2,2) / sum(H(:,2));
	X.NPV = H(1,1) / sum(H(:,1));
	X.FDR = H(1,2) / sum(H(:,2));
	X.MCC = det(H) / sqrt(prod([sum(H), sum(H')]));
	X.PHI = X.MCC; 
	X.F1  = 2 * X.TP / (sum(H(2,:)) + sum(H(:,2)));
	X.Sensitivity = X.TPR;	%% hit rate, recall
	X.Specificity = 1 - X.FPR;
	X.Precision   = X.PPV;
	X.dprime = norminv(X.TPR) - norminv(X.FDR);
	X.YI = X.Sensitivity + X.Specificity - 1;  % Youden Index

	% statistical significance test of Matthews' correlation coefficient
	NN = sum(H(:));
	R = X.MCC;
	R(isnan(R)) = 0;
	tmp = 1 - R.*R;
	tmp(tmp<0) = 0;		% prevent tmp<0 i.e. imag(t)~=0
	t   = R.*sqrt(max(NN-2,0)./tmp)
	sig = tcdf(t,NN-2);
	sig = 2 * min(sig,1 - sig);
	X.MCC_p_value = sig;
end;

kap = X;  
warning(s);
