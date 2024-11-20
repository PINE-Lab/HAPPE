function [V,D] = csp(ECM,arg2,arg3)
% CSP computes common spatial patterns
% 	this version supports multiple classes using a One-vs-Rest scheme
%
% [V] = csp(ECM)
% [V] = csp(X,Y)
% [V] = csp(...,Mode)
%
% ECM(k,:,:) is the extended covariance matrices (see COVM) for class k  
% X,Y	are matrices of the two classes (one channel per column)
%	the number of columns must be the same for X and Y
% Mode  = 'CSP0' uses common diagonalization 
%	= 'CSP3' solves generalized eigenvalue problem
% V	each column represents one CSP component. 
%
% REFERENCES: 
% [1] Koles ZJ, Soong AC.
% 	EEG source localization: implementing the spatio-temporal decomposition approach.
% 	Electroencephalogr Clin Neurophysiol. 1998 Nov;107(5):343-52
% [2] Ramoser, H.; Muller-Gerking, J.; Pfurtscheller, G.;
% 	Optimal spatial filtering of single trial EEG during imagined hand movement
%	Rehabilitation Engineering, IEEE Transactions on [see also IEEE Trans. on Neural Systems and Rehabilitation]
%	Volume 8,  Issue 4,  Dec. 2000 Page(s):441 - 446 
% [3] Dornhege, G.; Blankertz, B.; Curio, G.; Muller, K.-R.;
%    	Boosting bit rates in noninvasive EEG single-trial classifications by feature combination and multiclass paradigms
% 	Biomedical Engineering, IEEE Transactions on
% 	Volume 51,  Issue 6,  June 2004 Page(s):993 - 1002
%	Digital Object Identifier 10.1109/TBME.2004.827088 
% [4] Lemm, S.; Blankertz, B.; Curio, G.; Muller, K.-R.;
%    	Spatio-Spectral Filters for Improving the Classification of Single Trial EEG
% 	Biomedical Engineering, IEEE Transactions on
% 	Volume 52,  Issue 9,  Sept. 2005 Page(s):1541 - 1548
%	Digital Object Identifier 10.1109/TBME.2005.851521 

%	Copyright (C) 2007,2008,2009,2019 by Alois Schloegl <alois.schloegl@gmail.com>
%	This is part of the BIOSIG-toolbox http://biosig.sf.net/

p = 2; 

Mode =[];
sz = size(ECM); 
if (length(sz)<3) && isnumeric(arg2),
	X = ECM;
	ECM = permute(cat(3,covm(X,'E'),covm(arg2,'E')),[3,1,2]);
end; 
sz = size(ECM); 
if (nargin>1) && ischar(arg2)
	Mode = arg2; 
end; 
if (nargin>1) && ischar(arg2)
	Mode = arg3; 
end; 
if isempty(Mode),
	Mode = 'CSP3';
end; 	

sz=size(ECM);
COV = zeros(sz([2,3,1])-[1,1,0]);
for k=1:sz(1),
	[mu,sd,COV(:,:,k),xc,N,R2]=decovm(squeeze(ECM(k,:,:)));
end; 


V = repmat(NaN,sz(2)-1,2*sz(1));
if 0,

elseif strcmpi(Mode,'CSP0');  
	% common diagonalization 
	[P,D] = eig(squeeze(sum(COV,3))); 
	P = diag(sqrt(1./diag(D)))*P';

	% class specific components
	V = repmat(NaN,sz(2)-1,2*sz(1));
	d = V(1,:);
	for k = 1:sz(1), 
		C = P * squeeze(COV(:,:,k)) * P';
  		[R,d1]  = eig((C+C')/2);
		[d1,ix] = sort(diag(d1)); 
	
	  	V(:,2*k+[-1:0]) = P'*R(:,ix([1,end]));
	  	d(1,2*k+[-1:0]) = d1(ix([1,end]))';
	end; 

elseif strcmpi(Mode,'CSP3');  
	%% do actual CSP calculation as generalized eigenvalues
	%% R = permute(COV,[2,3,1]);
	for k = 1:sz(1),
		[W,D] = eig(COV(:,:,k),sum(COV,3));
		V(:,2*k+[1-p:0]) = W(:,[1,end]);
	end;
end; 

%!test
%! c1=[ones(30,1), randn(30, 50)]; c2=[ones(30,1),randn(30, 50)];
%! [v]=csp(c1'*c1, c2'*c2);
%! [v]=csp(c1'*c1, c2'*c2);
%! assert(all(size(v)==[51,4]))
%! assert(all(size(v)==[51,4]),'ECM0')
%! assert(all(size(v)==[51,4]),'ECM3')
%! c1=randn(30, 50); c2=randn(30, 50);
%! [v]=csp(c1, c2);
%! assert(all(size(v)==[50,4]))
%! assert(all(size(v)==[50,4]),'ECM0')
%! assert(all(size(v)==[50,4]),'ECM3')

