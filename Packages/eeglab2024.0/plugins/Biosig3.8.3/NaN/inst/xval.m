function [R,CC]=xval(D,classlabel,MODE,arg4)
% XVAL is used for crossvalidation 
%
%  [R,CC] = xval(D,classlabel)
%  .. = xval(D,classlabel,CLASSIFIER)
%  .. = xval(D,classlabel,CLASSIFIER,type)
%  .. = xval(D,{classlabel,W},CLASSIFIER)
%  .. = xval(D,{classlabel,W,NG},CLASSIFIER)
% 
%  example: 
%      load_fisheriris;    %builtin iris dataset      
%      C = species;
%      K = 5; NG = [1:length(C)]'*K/length(C);
%      [R,CC] = xval(meas,{C,[],NG},'NBC');            
%
% Input:
%    D:	data features (one feature per column, one sample per row)
%    classlabel	labels of each sample, must have the same number of rows as D. 
% 		Two different encodings are supported: 
%		{-1,1}-encoding (multiple classes with separate columns for each class) or
%		1..M encoding. 
% 		So [1;2;3;1;4] is equivalent to 
%			[+1,-1,-1,-1;
%			[-1,+1,-1,-1;
%			[-1,-1,+1,-1;
%			[+1,-1,-1,-1]
%			[-1,-1,-1,+1]
%		Note, samples with classlabel=0 are ignored. 
%
%    CLASSIFIER can be any classifier supported by train_sc (default='LDA')
%       {'REG','MDA','MD2','QDA','QDA2','LD2','LD3','LD4','LD5','LD6','NBC','aNBC','WienerHopf', 'RDA','GDBC',
%	 'SVM','RBF','PSVM','SVM11','SVM:LIN4','SVM:LIN0','SVM:LIN1','SVM:LIN2','SVM:LIN3','WINNOW'}
%       these can be modified by ###/GSVD, ###/sparse and ###/DELETION. 
%	   /DELETION removes in case of NaN's either the rows or the columns (which removes less data values) with any NaN
%	   /sparse and /GSVD preprocess the data an reduce it to some lower-dimensional space. 
%       Hyperparameters (like alpha for PLA, gamma/lambda for RDA, c_value for SVM, etc) can be defined as 
% 	CLASSIFIER.hyperparameter.alpha, etc. and 
% 	CLASSIFIER.TYPE = 'PLA' (as listed above). 
%       See train_sc for details.
%    W:	weights for each sample (row) in D. 
%	default: [] (i.e. all weights are 1)
%	number of elements in W must match the number of rows of D 
%    NG: used to define the type of cross-valdiation
% 	Leave-One-Out-Method (LOOM): NG = [1:length(classlabel)]' (default)
% 	Leave-K-Out-Method: NG = ceil([1:length(classlabel)]'/K)
%	K-fold XV:  NG = ceil([1:length(classlabel)]'*K/length(classlabel))
%	group-wise XV (if samples are not indepentent) can be also defined here
%	samples from the same group (dependent samples) get the same identifier
%	samples from different groups get different classifiers
%    TYPE:  defines the type of cross-validation procedure if NG is not specified 
%	'LOOM'  leave-one-out-method
%       k	k-fold crossvalidation
%
% OUTPUT: 
%    R contains the resulting performance metric
%    R.OUTPUT classifier output
%    R.CL	class labels
%    R.H	confusion matrix
%    R.kappa	Cohen's kappa coefficient
%    R.ACC	accuracy 
%    R.ERR	error rate
%        ... and a number of addition parameters obtained from the confusionm matrix [2]
%    R.noXV.{...} contains the result w/o crossvalidation
%    CC contains the classifier  
%
%    plota(R) shows the confusion matrix of the results
%
% see also: TRAIN_SC, TEST_SC, CLASSIFY, PLOTA
%
% References: 
% [1] R. Duda, P. Hart, and D. Stork, Pattern Classification, second ed. 
%       John Wiley & Sons, 2001. 
% [2] A. Schlögl, J. Kronegg, J.E. Huggins, S. G. Mason;
%       Evaluation criteria in BCI research.
%       (Eds.) G. Dornhege, J.R. Millan, T. Hinterberger, D.J. McFarland, K.-R.Müller;
%       Towards Brain-Computer Interfacing, MIT Press, 2007, p.327-342

%	Copyright (C) 2008,2009,2010,2016 by Alois Schloegl
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

if (nargin<3) || isempty(MODE),
	MODE = 'LDA';
end;
if ischar(MODE)
        tmp = MODE;
        clear MODE;
        MODE.TYPE = tmp;
elseif ~isfield(MODE,'TYPE')
        MODE.TYPE='';
end;

sz = size(D);
NG = [];
W = [];

if iscell(classlabel)
        % hack to handle NaN's in unique(...)
        c  = classlabel{:,1};
        ix = find(~isnan(c));
        C  = c;
        [b,i,C(ix)] = unique(c(ix));
        if size(classlabel,2)>1,
                W = [classlabel{:,2}]; 
        end; 
	if size(classlabel,2)>2,
		[Label,tmp1,NG] = unique(classlabel{:,3});
	end;
elseif size(classlabel,2)>1,
	%% group-wise classvalidation
	C = classlabel(:,1);
	W = classlabel(:,2);
	if size(classlabel,2)==2,
	        warning('This option defines W and NG in an ambigous way - use instead xval(D,{C,[],NG},...) or xval(D,{C,W},...)'); 
	else
		[Label,tmp1,NG] = unique(classlabel(:,3));
	end;
else
	C = classlabel;	
end; 
if all(W==1), W = []; end;
if sz(1)~=size(C,1),
        error('length of data and classlabel does not fit');
end;


if isempty(NG)
if (nargin<4) || strcmpi(arg4,'LOOM')
	%% LOOM 
	NG = (1:sz(1))';

elseif isnumeric(arg4)
	if isscalar(arg4)  
	% K-fold XV
		NG = ceil((1:length(C))'*arg4/length(C));
	elseif length(arg4)==2,
		NG = ceil((1:length(C))'*arg4(1)/length(C));
	end;

end;
end;

sz = size(D);
if sz(1)~=length(C),
        error('length of data and classlabel does not fit');
end;
if ~isfield(MODE,'hyperparameter')
        MODE.hyperparameter = [];
end

cl     = repmat(NaN,size(classlabel,1),1);
output = repmat(NaN,size(classlabel,1),max(C));
for k = 1:max(NG),
 	ix = find(~any(isnan(C),2) & (NG~=k));
	if isempty(W),	
		CC = train_sc(D(ix,:), C(ix), MODE);
	else
		CC = train_sc(D(ix,:), C(ix), MODE, W(ix));
	end;
 	ix = find(NG==k);
	r  = test_sc(CC, D(ix,:));
	cl(ix,1)     = r.classlabel;
	output(ix,:) = r.output;
end; 

%R = kappa(C,cl,'notIgnoreNAN',W);
R = kappa(C,cl,[],W);
%R2 = kappa(R.H);

R.OUTPUT=output;
R.CL=cl;
R.ERR = 1-R.ACC; 
if isnumeric(R.Label)
	R.Label = cellstr(int2str(R.Label)); 
end;

if nargout>1,
	% final classifier 
 	ix = find(~any(isnan(C),2));
	if isempty(W), 
		CC = train_sc(D(ix,:), C(ix), MODE);
	else	
		CC = train_sc(D(ix,:), C(ix), MODE,W);
	end; 	
	CC.Labels = 1:max(C);
	%CC.Labels = unique(C);

	% test w/o cross-validation
	r  = test_sc(CC, D(ix,:));
	R.noXV = kappa(C(ix), r.classlabel,[],W);
	R.noXV.OUTPUT = r.output;
	R.noXV.CL     = r.classlabel;
	R.noXV.CC     = CC; 
end;

