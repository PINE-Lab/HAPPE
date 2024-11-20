function [R]=test_sc(CC,D,mode,classlabel)
% TEST_SC: apply statistical and SVM classifier to test data 
%
%  R = test_sc(CC,D,TYPE [,target_Classlabel]) 
%       R.output     	output: "signed" distance for each class. 
%		This represents the distances between sample D and the separating hyperplane
%		The "signed distance" is possitive if it matches the target class, and 
%		and negative if it lays on the opposite side of the separating hyperplane. 
%       R.classlabel 	class for output data
%  The target class is optional. If it is provided, the following values are returned. 
%       R.kappa 	Cohen's kappa coefficient
%       R.ACC   	Classification accuracy 
%       R.H     	Confusion matrix 
%
% The classifier CC is typically obtained by TRAIN_SC. If a statistical 
% classifier is used, TYPE can be used to modify the classifier. 
%    TYPE = 'MDA'    mahalanobis distance based classifier
%    TYPE = 'MD2'    mahalanobis distance based classifier
%    TYPE = 'MD3'    mahalanobis distance based classifier
%    TYPE = 'GRB'    Gaussian radial basis function 
%    TYPE = 'QDA'    quadratic discriminant analysis
%    TYPE = 'LD2'    linear discriminant analysis
%    TYPE = 'LD3', 'LDA', 'FDA, 'FLDA'   (Fisher's) linear discriminant analysis
%    TYPE = 'LD4'    linear discriminant analysis
%    TYPE = 'GDBC'   general distance based classifier
% 
% see also: TRAIN_SC
%
% References: 
% [1] R. Duda, P. Hart, and D. Stork, Pattern Classification, second ed. 
%       John Wiley & Sons, 2001.

%	Copyright (C) 2005,2006,2008,2009,2010,2016 by Alois Schloegl <alois.schloegl@gmail.com>
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

if nargin<3,
        mode = [];
end;
[t1,t] = strtok(CC.datatype,':');
[t2,t] = strtok(t,':');
[t3] = strtok(t,':');
if ~strcmp(t1,'classifier'), return; end; 

if isfield(CC,'prewhite')
        D = D*CC.prewhite(2:end,:) + CC.prewhite(ones(size(D,1),1),:);
        CC = rmfield(CC,'prewhite');
end;

POS1 = [strfind(CC.datatype,'/gsvd'),strfind(CC.datatype,'/sparse'),strfind(CC.datatype,'/delet')];

if 0,


elseif strcmp(CC.datatype,'classifier:nbpw')
	error('NBPW not implemented yet')
	%%%% Naive Bayesian Parzen Window Classifier %%%%
	d = repmat(NaN,size(D,1),size(CC.MEAN,1));
	for k = 1:size(CC.MEAN,1)
		z = (D - CC.MEAN(repmat(k,size(D,1),1),:)).^2 ./ (CC.VAR(repmat(k,size(D,1),1),:));
		z = z + log(CC.VAR(repmat(k,size(D,1),1),:)); % + log(2*pi);
		d(:,k) = sum(-z/2, 2) + log(mean(CC.N(k,:)));
	end;
	d = exp(d-log(mean(sum(CC.N,1)))-log(2*pi)/2);


elseif strcmp(CC.datatype,'classifier:nbc')
	%%%% Naive Bayesian Classifier %%%%
	d = repmat(NaN,size(D,1),size(CC.MEAN,1));
	for k = 1:size(CC.MEAN,1)
		z = (D - CC.MEAN(repmat(k,size(D,1),1),:)).^2 ./ (CC.VAR(repmat(k,size(D,1),1),:));
		z = z + log(CC.VAR(repmat(k,size(D,1),1),:)); % + log(2*pi);
		d(:,k) = sum(-z/2, 2) + log(mean(CC.N(k,:)));
	end; 
	d = exp(d-log(mean(sum(CC.N,1)))-log(2*pi)/2);


elseif strcmp(CC.datatype,'classifier:anbc')
	%%%% Augmented Naive Bayesian Classifier %%%%
	d = repmat(NaN,size(D,1),size(CC.MEAN,1));
	for k = 1:size(CC.MEAN,1)
		z = (D*CC.V - CC.MEAN(repmat(k,size(D,1),1),:)).^2 ./ (CC.VAR(repmat(k,size(D,1),1),:));
		z = z + log(CC.VAR(repmat(k,size(D,1),1),:)); % + log(2*pi);
		d(:,k) = sum(-z/2, 2) + log(mean(CC.N(k,:)));
	end; 
	d = exp(d-log(mean(sum(CC.N,1)))-log(2*pi)/2);


elseif strcmp(CC.datatype,'classifier:statistical:rda')
	% Friedman (1989) Regularized Discriminant analysis
	if isfield(CC,'hyperparameter') && isfield(CC.hyperparameter,'lambda')  && isfield(CC.hyperparameter,'gamma')
		D = [ones(size(D,1),1),D];  % add 1-column
		lambda = CC.hyperparameter.lambda;
		gamma  = CC.hyperparameter.gamma;
		d = repmat(NaN,size(D,1),size(CC.MD,3));
                ECM = CC.MD./CC.NN; 
                NC = size(ECM); 
                ECM0 = squeeze(sum(ECM,3));  %decompose ECM
                [M0,sd,COV0] = decovm(ECM0);
                for k = 1:NC(3);
			[M,sd,s,xc,N] = decovm(squeeze(ECM(:,:,k)));
                	s = ((1-lambda)*N*s+lambda*COV0)/((1-lambda)*N+lambda);
                	s = (1-gamma)*s+gamma*(trace(s))/(NC(2)-1)*eye(NC(2)-1);
                        ir  = [-M;eye(NC(2)-1)]*inv(s)*[-M',eye(NC(2)-1)];  % inverse correlation matrix extended by mean
                        d(:,k) = -sum((D*ir).*D,2); % calculate distance of each data point to each class
                end;
	else 
		error('QDA: hyperparamters lambda and/or gamma not defined')
	end;


elseif strcmp(CC.datatype,'classifier:csp')
	d = filtfilt(CC.FiltB,CC.FiltA,(D*CC.csp_w).^2);
	R = test_sc(CC.CSP,log(d));	% LDA classifier of 
	d = R.output; 


elseif strcmp(CC.datatype,'classifier:svm:lib:1vs1') || strcmp(CC.datatype,'classifier:svm:lib:rbf');
        nr = size(D,1);
	[cl] = svmpredict_mex(ones(nr,1), D, CC.model, '-q');   %Use the classifier
        %Create a pseudo tsd matrix for bci4eval
	d = full(sparse(1:nr,cl,1,nr,CC.model.nr_class));

elseif isfield(CC,'weights'); %strcmpi(t2,'svm') || (strcmpi(t2,'statistical') & strncmpi(t3,'ld',2)) ;
        % linear classifiers like: LDA, SVM, LPM 
        %d = [ones(size(D,1),1), D] * CC.weights;
        d = repmat(NaN,size(D,1),size(CC.weights,2));
        for k = 1:size(CC.weights,2),
                d(:,k) = D * CC.weights(2:end,k) + CC.weights(1,k);
        end;


elseif ~isempty(POS1)	% GSVD, sparse & DELETION
        CC.datatype = CC.datatype(1:POS1(1)-1);
        r = test_sc(CC, D*sparse(CC.G));
        d = r.output; 


elseif strcmp(t2,'statistical');
	if isempty(mode)
		mode.TYPE = upper(t3); 
	else
		tmp=mode;
		clear mode;
		mode.TYPE=tmp;                
	end;
	D = [ones(size(D,1),1),D];  % add 1-column
	W = repmat(NaN, size(D,2), size(CC.MD,3));

        if 0,
        elseif strcmpi(mode.TYPE,'LD2'),
                %d = ldbc2(CC,D);
                ECM = CC.MD./CC.NN; 
                NC = size(ECM); 
                ECM0 = squeeze(sum(ECM,3));  %decompose ECM
                [M0] = decovm(ECM0);
                for k = 1:NC(3);
                        ecm = squeeze(ECM(:,:,k));
                        [M1,sd,COV1] = decovm(ECM0-ecm);
                        [M2,sd,COV2] = decovm(ecm);
                        w     = (COV1+COV2)\(M2'-M1')*2;
                        w0    = -M0*w;
                        W(:,k) = [w0; w];
                end;
                d = D*W;
        elseif strcmpi(mode.TYPE,'LD3') || strcmpi(mode.TYPE,'FLDA');
                %d = ldbc3(CC,D);
                ECM = CC.MD./CC.NN; 
                NC = size(ECM); 
                ECM0 = squeeze(sum(ECM,3));  %decompose ECM
                [M0,sd,COV0] = decovm(ECM0);
                for k = 1:NC(3);
                        ecm = squeeze(ECM(:,:,k));
                        [M1] = decovm(ECM0-ecm);
                        [M2] = decovm(ecm);
                        w     = COV0\(M2'-M1')*2;
                        w0    = -M0*w;
                        W(:,k) = [w0; w];
                end;
                d = D*W;
        elseif strcmpi(mode.TYPE,'LD4');
                %d = ldbc4(CC,D);
                ECM = CC.MD./CC.NN; 
                NC = size(ECM); 
                ECM0 = squeeze(sum(ECM,3));  %decompose ECM
                M0 = decovm(ECM0);
                for k = 1:NC(3);
                        ecm = squeeze(ECM(:,:,k));
                        [M1,sd,COV1,xc,N1] = decovm(ECM0-ecm);
                        [M2,sd,COV2,xc,N2] = decovm(ecm);
                        w     = (COV1*N1+COV2*N2)\((M2'-M1')*(N1+N2));
                        w0    = -M0*w;
                        W(:,k) = [w0; w];
                end;
                d = D*W;
        elseif strcmpi(mode.TYPE,'MDA');
                d = repmat(NaN,size(D,1),length(CC.IR)); 
                for k = 1:length(CC.IR);
                        d(:,k) = -sum((D*CC.IR{k}).*D,2); % calculate distance of each data point to each class
                end;
        elseif strcmpi(mode.TYPE,'MD2');
                d = repmat(NaN,size(D,1),length(CC.IR)); 
                for k = 1:length(CC.IR);
                        d(:,k) = sum((D*CC.IR{k}).*D,2); % calculate distance of each data point to each class
                end;
                d = -sqrt(d);
        elseif strcmpi(mode.TYPE,'GDBC');
                d = repmat(NaN,size(D,1),length(CC.IR)); 
                for k = 1:length(CC.IR);
                        d(:,k) = sum((D*CC.IR{k}).*D,2) + CC.logSF7(k); % calculate distance of each data point to each class
                end;
                d = exp(-d/2);
        elseif strcmpi(mode.TYPE,'MD3');
                d = repmat(NaN,size(D,1),length(CC.IR)); 
                for k = 1:length(CC.IR);
                        d(:,k) = sum((D*CC.IR{k}).*D,2) + CC.logSF7(k); % calculate distance of each data point to each class
                end;
                d = exp(-d/2);
                d = d./repmat(sum(d,2),1,size(d,2));  % Zuordungswahrscheinlichkeit [1], p.601, equ (18.39)
        elseif strcmpi(mode.TYPE,'QDA');
                d = repmat(NaN,size(D,1),length(CC.IR)); 
                for k = 1:length(CC.IR);
                        % [1] (18.33) QCF - quadratic classification function  
                        d(:,k) = -(sum((D*CC.IR{k}).*D,2) - CC.logSF5(k)); 
                end;
        elseif strcmpi(mode.TYPE,'QDA2');
                d = repmat(NaN,size(D,1),length(CC.IR)); 
                for k = 1:length(CC.IR);
                        % [1] (18.33) QCF - quadratic classification function  
                        d(:,k) = -(sum((D*(CC.IR{k})).*D,2) + CC.logSF4(k)); 
                end;
        elseif strcmpi(mode.TYPE,'GRB');     % Gaussian RBF
                d = repmat(NaN,size(D,1),length(CC.IR)); 
                for k = 1:length(CC.IR);
                        d(:,k) = sum((D*CC.IR{k}).*D,2); % calculate distance of each data point to each class
                end;
                d = exp(-sqrt(d)/2);
        elseif strcmpi(mode.TYPE,'GRB2');     % Gaussian RBF
                d = repmat(NaN,size(D,1),length(CC.IR)); 
                for k = 1:length(CC.IR);
                        d(:,k) = sum((D*CC.IR{k}).*D,2); % calculate distance of each data point to each class
                end;
                d = exp(-d);
        elseif strcmpi(mode.TYPE,'MQU');     % Multiquadratic 
                d = repmat(NaN,size(D,1),length(CC.IR)); 
                for k = 1:length(CC.IR);
                        d(:,k) = sum((D*CC.IR{k}).*D,2); % calculate distance of each data point to each class
                end;
                d = -sqrt(1+d);
        elseif strcmpi(mode.TYPE,'IMQ');     % Inverse Multiquadratic 
                d = repmat(NaN,size(D,1),length(CC.IR)); 
                for k = 1:length(CC.IR);
                        d(:,k) = sum((D*CC.IR{k}).*D,2); % calculate distance of each data point to each class
                end;
                d = (1+d).^(-1/2);
        elseif strcmpi(mode.TYPE,'Cauchy');     % Cauchy RBF
                d = repmat(NaN,size(D,1),length(CC.IR)); 
                for k = 1:length(CC.IR);
                        d(:,k) = sum((D*CC.IR{k}).*D,2); % calculate distance of each data point to each class
                end;
                d = 1./(1+d);
        else 
                error('Classifier %s not supported. see HELP TRAIN_SC for supported classifiers.',mode.TYPE);
        end;
else
        fprintf(2,'Error TEST_SC: unknown classifier\n');
        return;
end;

if size(d,2)>1,
	[tmp,cl] = max(d,[],2);
	cl = CC.Labels(cl); 
	cl(isnan(tmp)) = NaN; 
elseif size(d,2)==1,
	cl = (d<0) + 2*(d>0);
	cl(isnan(d)) = NaN; 
end; 	

R.output = d; 
R.classlabel = cl; 

if nargin>3,
        [R.kappa,R.sd,R.H,z,R.ACC] = kappa(classlabel(:),cl(:));
end;
