function [CC]=train_sc(D,classlabel,MODE,W)
% Train a (statistical) classifier
% 
%  CC = train_sc(D,classlabel)
%  CC = train_sc(D,classlabel,MODE)
%  CC = train_sc(D,classlabel,MODE, W)
%	weighting D(k,:) with weight W(k) (not all classifiers supported weighting)
%
% CC contains the model parameters of a classifier which can be applied 
%   to test data using test_sc. 
%   R = test_sc(CC,D,...) 
%
%   D		training samples (each row is a sample, each column is a feature)	
%   classlabel	labels of each sample, must have the same number of rows as D. 
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
%  The following classifier types are supported MODE.TYPE
%    'MDA'      mahalanobis distance based classifier [1]
%    'MD2'      mahalanobis distance based classifier [1]
%    'MD3'      mahalanobis distance based classifier [1]
%    'GRB'      Gaussian radial basis function     [1]
%    'QDA'      quadratic discriminant analysis    [1]
%    'LD2'      linear discriminant analysis (see LDBC2) [1]
%		MODE.hyperparameter.gamma: regularization parameter [default 0] 
%    'LD3', 'FDA', 'LDA', 'FLDA'
%               linear discriminant analysis (see LDBC3) [1]
%		MODE.hyperparameter.gamma: regularization parameter [default 0] 
%    'LD4'      linear discriminant analysis (see LDBC4) [1]
%		MODE.hyperparameter.gamma: regularization parameter [default 0] 
%    'LD5'      another LDA (motivated by CSP)
%		MODE.hyperparameter.gamma: regularization parameter [default 0] 
%    'RDA'      regularized discriminant analysis [7]
%		MODE.hyperparameter.gamma: regularization parameter 
%		MODE.hyperparameter.lambda =
%		gamma = 0, lambda = 0 : MDA
%		gamma = 0, lambda = 1 : LDA [default]
%		Hint: hyperparameter are used only in test_sc.m, testing different 
%		the hyperparameters do not need repetitive calls to train_sc, 
%		it is sufficient to modify CC.hyperparameter before calling test_sc. 	
%    'GDBC'     general distance based classifier  [1]
%    ''         statistical classifier, requires Mode argument in TEST_SC	
%    '###/DELETION'  if the data contains missing values (encoded as NaNs), 
%		a row-wise or column-wise deletion (depending on which method 
%		removes less data values) is applied;  
%    '###/GSVD'	GSVD and statistical classifier [2,3], 
%    '###/sparse'  sparse  [5] 
%		'###' must be 'LDA' or any other classifier 
%    'PLS'	(linear) partial least squares regression 
%    'REG'      regression analysis;
%    'WienerHopf'	Wiener-Hopf equation  
%    'NBC'	Naive Bayesian Classifier [6]     
%    'aNBC'	Augmented Naive Bayesian Classifier [6]
%    'NBPW'	Naive Bayesian Parzen Window [9]     
%
%    'PLA'	Perceptron Learning Algorithm [11]
%		MODE.hyperparameter.alpha = alpha [default: 1]
%		 w = w + alpha * e'*x
%    'LMS', 'AdaLine'  Least mean squares, adaptive line element, Widrow-Hoff, delta rule 
%		MODE.hyperparameter.alpha = alpha [default: 1]
%    'Winnow2'  Winnow2 algorithm [12]
%
%    'PSVM'	Proximal SVM [8] 
%		MODE.hyperparameter.nu  (default: 1.0)
%    'LPM'      Linear Programming Machine
%                 uses and requires train_LPM of the iLog CPLEX optimizer 
%		MODE.hyperparameter.c_value = 
%    'CSP'	CommonSpatialPattern is very experimental and just a hack
%		uses a smoothing window of 50 samples.
%    'SVM','SVM1r'  support vector machines, one-vs-rest
%		MODE.hyperparameter.c_value = 
%    'SVM11'    support vector machines, one-vs-one + voting
%		MODE.hyperparameter.c_value = 
%    'RBF'      Support Vector Machines with RBF Kernel
%		MODE.hyperparameter.c_value = 
%		MODE.hyperparameter.gamma = 
%    'SVM:LIB'    libSVM [default SVM algorithm)
%    'SVM:bioinfo' uses and requires svmtrain from the bioinfo toolbox        
%    'SVM:OSU'   uses and requires mexSVMTrain from the OSU-SVM toolbox 
%    'SVM:LOO'   uses and requires svcm_train from the LOO-SVM toolbox 
%    'SVM:Gunn'  uses and requires svc-functios from the Gunn-SVM toolbox 
%    'SVM:KM'    uses and requires svmclass-function from the KM-SVM toolbox 
%    'SVM:LINz'  LibLinear [10] (requires train.mex from LibLinear somewhere in the path)
%            z=0 (default) LibLinear with -- L2-regularized logistic regression
%            z=1 LibLinear with -- L2-loss support vector machines (dual)
%            z=2 LibLinear with -- L2-loss support vector machines (primal)
%            z=3 LibLinear with -- L1-loss support vector machines (dual)
%    'SVM:LIN4'  LibLinear with -- multi-class support vector machines by Crammer and Singer
%    'DT'	decision tree - not implemented yet.  
%
% {'REG','MDA','MD2','QDA','QDA2','LD2','LD3','LD4','LD5','LD6','NBC','aNBC','WienerHopf','LDA/GSVD','MDA/GSVD', 'LDA/sparse','MDA/sparse', 'PLA', 'LMS','LDA/DELETION','MDA/DELETION','NBC/DELETION','RDA/DELETION','REG/DELETION','RDA','GDBC','SVM','RBF','PSVM','SVM11','SVM:LIN4','SVM:LIN0','SVM:LIN1','SVM:LIN2','SVM:LIN3','WINNOW', 'DT'};
%
% CC contains the model parameters of a classifier. Some time ago,     
% CC was a statistical classifier containing the mean 
% and the covariance of the data of each class (encoded in the 
%  so-called "extended covariance matrices". Nowadays, also other 
% classifiers are supported. 
%
% see also: TEST_SC, COVM, ROW_COL_DELETION
%
% References: 
% [1] R. Duda, P. Hart, and D. Stork, Pattern Classification, second ed. 
%       John Wiley & Sons, 2001. 
% [2] Peg Howland and Haesun Park,
%       Generalizing Discriminant Analysis Using the Generalized Singular Value Decomposition
%       IEEE Transactions on Pattern Analysis and Machine Intelligence, 26(8), 2004.
%       dx.doi.org/10.1109/TPAMI.2004.46
% [3] http://www-static.cc.gatech.edu/~kihwan23/face_recog_gsvd.htm
% [4] Jieping Ye, Ravi Janardan, Cheong Hee Park, Haesun Park
%       A new optimization criterion for generalized discriminant analysis on undersampled problems.
%       The Third IEEE International Conference on Data Mining, Melbourne, Florida, USA
%       November 19 - 22, 2003
% [5] J.D. Tebbens and P. Schlesinger (2006), 
%       Improving Implementation of Linear Discriminant Analysis for the Small Sample Size Problem
%	Computational Statistics & Data Analysis, vol 52(1): 423-437, 2007
%       http://www.cs.cas.cz/mweb/download/publi/JdtSchl2006.pdf
% [6] H. Zhang, The optimality of Naive Bayes, 
%	 http://www.cs.unb.ca/profs/hzhang/publications/FLAIRS04ZhangH.pdf
% [7] J.H. Friedman. Regularized discriminant analysis. 
%	Journal of the American Statistical Association, 84:165–175, 1989.
% [8] G. Fung and O.L. Mangasarian, Proximal Support Vector Machine Classifiers, KDD 2001.
%        Eds. F. Provost and R. Srikant, Proc. KDD-2001: Knowledge Discovery and Data Mining, August 26-29, 2001, San Francisco, CA.
% 	p. 77-86.
% [9] Kai Keng Ang, Zhang Yang Chin, Haihong Zhang, Cuntai Guan.
%	Filter Bank Common Spatial Pattern (FBCSP) in Brain-Computer Interface.
%	IEEE International Joint Conference on Neural Networks, 2008. IJCNN 2008. (IEEE World Congress on Computational Intelligence). 
%	1-8 June 2008 Page(s):2390 - 2397
% [10] R.-E. Fan, K.-W. Chang, C.-J. Hsieh, X.-R. Wang, and C.-J. Lin. 
%       LIBLINEAR: A Library for Large Linear Classification, Journal of Machine Learning Research 9(2008), 1871-1874. 
%       Software available at http://www.csie.ntu.edu.tw/~cjlin/liblinear 
% [11] http://en.wikipedia.org/wiki/Perceptron#Learning_algorithm
% [12] Littlestone, N. (1988) 
%       "Learning Quickly When Irrelevant Attributes Abound: A New Linear-threshold Algorithm" 
%       Machine Learning 285-318(2)
% 	http://en.wikipedia.org/wiki/Winnow_(algorithm)
 
%	Copyright (C) 2005-2010,2019 by Alois Schloegl <alois.schloegl@gmail.com>
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

        if nargin<2,
                error('insufficient input arguments\n\tusage: train_sc(D,C,...)\n');
        end
        if nargin<3, MODE = 'LDA'; end
        if nargin<4, W = []; end
        if ischar(MODE)
                tmp = MODE;
                clear MODE;
                MODE.TYPE = tmp;
        elseif ~isfield(MODE,'TYPE')
                MODE.TYPE='';
        end

        if isfield(MODE,'hyperparameters') && ~isfield(MODE,'hyperparameter'),
                %% for backwards compatibility, this might become obsolete
                warning('MODE.hyperparameters are used, You should use MODE.hyperparameter instead!!!');
                MODE.hyperparameter = MODE.hyperparameters;
        end

        sz = size(D);
        if sz(1)~=size(classlabel,1),
                error('length of data and classlabel does not fit');
        end

        % remove all NaN's
        if 1,
                % several classifier can deal with NaN's, there is no need to remove them.
        elseif isempty(W)
                %% TODO: some classifiers can deal with NaN's in D. Test whether this can be relaxed.
                %ix = any(isnan([classlabel]),2);
                ix = any(isnan([D,classlabel]),2);
                D(ix,:) = [];
                classlabel(ix,:)=[];
                W = [];
        else
                %ix = any(isnan([classlabel]),2);
                ix = any(isnan([D,classlabel]),2);
                D(ix,:)=[];
                classlabel(ix,:)=[];
                W(ix,:)=[];
                warning('support for weighting of samples is still experimental');
        end

        sz = size(D);
        if sz(1)~=length(classlabel),
                error('length of data and classlabel does not fit');
        end
        if ~isfield(MODE,'hyperparameter')
                MODE.hyperparameter = [];
        end

	if ~isfield(MODE,'RowDeletion')
		MODE.RowDeletion = 0;
	end
	if MODE.RowDeletion,
		D = D(~any(isnan(D),2),:);
	end

        if 0,
                ;
        elseif ~isempty(strfind(lower(MODE.TYPE),'/delet'))
                POS1 = find(MODE.TYPE=='/');
                [rix,cix] = row_col_deletion(D);
                if ~isempty(W), W=W(rix); end
                CC   = train_sc(D(rix,cix),classlabel(rix,:),MODE.TYPE(1:POS1(1)-1),W);
                CC.G = sparse(cix, 1:length(cix), 1, size(D,2), length(cix));
                if isfield(CC,'weights')
                        W = [CC.weights(1,:); CC.weights(2:end,:)];
                        CC.weights = sparse(size(D,2)+1, size(W,2));
                        CC.weights([1,cix+1],:) = W;
                        CC.datatype = ['classifier:statistical:',lower(MODE.TYPE)];
                else
                        CC.datatype = [CC.datatype,'/delet'];
                end

        elseif ~isempty(strfind(lower(MODE.TYPE),'nbpw'))
                error('NBPW not implemented yet')
                %%%% Naive Bayesian Parzen Window Classifier.
                [classlabel,CC.Labels] = CL1M(classlabel);
                for k = 1:length(CC.Labels),
                        [d,CC.MEAN(k,:)] = center(D(classlabel==CC.Labels(k),:),1);
                        [CC.VAR(k,:),CC.N(k,:)] = sumskipnan(d.^2,1);
                        h2_opt = (4./(3*CC.N(k,:))).^(2/5).*CC.VAR(k,:);
                        %%% TODO
                end


        elseif ~isempty(strfind(lower(MODE.TYPE),'nbc'))
                %%%% Naive Bayesian Classifier
                if ~isempty(strfind(lower(MODE.TYPE),'anbc'))
                        %%%% Augmented Naive Bayesian classifier.
                        [CC.V,L] = eig(covm(D,'M',W));
                        D = D*CC.V;
                else
                        CC.V = eye(size(D,2));
                end
                [classlabel,CC.Labels] = CL1M(classlabel);
                for k = 1:length(CC.Labels),
                        ix = classlabel==CC.Labels(k);
                        %% [d,CC.MEAN(k,:)] = center(D(ix,:),1);
                        if ~isempty(W)
                                [s,n] = sumskipnan(D(ix,:),1,W(ix));
                                CC.MEAN(k,:) = s./n;
                                d = D(ix,:) - CC.MEAN(repmat(k,sum(ix),1),:);
                                [CC.VAR(k,:),CC.N(k,:)] = sumskipnan(d.^2,1,W(ix));
                        else
                                [s,n] = sumskipnan(D(ix,:),1);
                                CC.MEAN(k,:) = s./n;
                                d = D(ix,:) - CC.MEAN(repmat(k,sum(ix),1),:);
                                [CC.VAR(k,:),CC.N(k,:)] = sumskipnan(d.^2,1);
                        end
                end
                CC.VAR = CC.VAR./max(CC.N-1,0);
                CC.datatype = ['classifier:',lower(MODE.TYPE)];


        elseif ~isempty(strfind(lower(MODE.TYPE),'lpm'))
                if ~isempty(W)
                        error('Error TRAIN_SC: Classifier (%s) does not support weighted samples.',MODE.TYPE);
                end
                % linear programming machine
                % CPLEX optimizer: ILOG solver, ilog cplex 6.5 reference manual http://www.ilog.com
                MODE.TYPE = 'LPM';
                if ~isfield(MODE.hyperparameter,'c_value')
                        MODE.hyperparameter.c_value = 1;
                end
                [classlabel,CC.Labels] = CL1M(classlabel);

                M = length(CC.Labels);
                if M==2, M=1; end   % For a 2-class problem, only 1 Discriminant is needed
                for k = 1:M,
                        %LPM = train_LPM(D,(classlabel==CC.Labels(k)),'C',MODE.hyperparameter.c_value);
                        LPM = train_LPM(D',(classlabel'==CC.Labels(k)));
                        CC.weights(:,k) = [-LPM.b; LPM.w(:)];
                end
                CC.hyperparameter.c_value = MODE.hyperparameter.c_value;
                CC.datatype = ['classifier:',lower(MODE.TYPE)];


        elseif ~isempty(strfind(lower(MODE.TYPE),'pla')),
                % Perceptron Learning Algorithm

                [rix,cix] = row_col_deletion(D);
                [CL101,CC.Labels] = cl101(classlabel);
                M = size(CL101,2);
                weights   = sparse(length(cix)+1,M);

                %ix = randperm(size(D,1));      %% randomize samples ???
                if ~isfield(MODE.hyperparameter,'alpha')
                        if isfield(MODE.hyperparameter,'alpha')
                                alpha = MODE.hyperparameter.alpha;
                        else
                                alpha = 1;
                        end
                        for k = rix(:)',
                                %e = ((classlabel(k)==(1:M))-.5) - sign([1, D(k,cix)] * weights)/2;
                                e = CL101(k,:) - sign([1, D(k,cix)] * weights);
                                weights = weights + alpha * [1,D(k,cix)]' * e ;
                        end

                else %if ~isempty(W)
                        if isfield(MODE.hyperparameter,'alpha')
                                W = W*MODE.hyperparameter.alpha;
                        end
                        for k = rix(:)',
                                %e = ((classlabel(k)==(1:M))-.5) - sign([1, D(k,cix)] * weights)/2;
                                e = CL101(k,:) - sign([1, D(k,cix)] * weights);
                                weights = weights + W(k) * [1,D(k,cix)]' * e ;
                        end
                end
                CC.weights  = sparse(size(D,2)+1,M);
                CC.weights([1,cix+1],:) = weights;
                CC.datatype = ['classifier:',lower(MODE.TYPE)];


        elseif  ~isempty(strfind(lower(MODE.TYPE),'adaline')) || ~isempty(strfind(lower(MODE.TYPE),'lms')),
                % adaptive linear elemente, least mean squares, delta rule, Widrow-Hoff,

                [rix,cix] = row_col_deletion(D);
                [CL101,CC.Labels] = cl101(classlabel);
                M = size(CL101,2);
                weights  = sparse(length(cix)+1,M);

                %ix = randperm(size(D,1));      %% randomize samples ???
                if isempty(W)
                        if isfield(MODE.hyperparameter,'alpha')
                                alpha = MODE.hyperparameter.alpha;
                        else
                                alpha = 1;
                        end
                        for k = rix(:)',
                                %e = (classlabel(k)==(1:M)) - [1, D(k,cix)] * weights;
                                e = CL101(k,:) - sign([1, D(k,cix)] * weights);
                                weights = weights + alpha * [1,D(k,cix)]' * e ;
                        end

                else %if ~isempty(W)
                        if isfield(MODE.hyperparameter,'alpha')
                                W = W*MODE.hyperparameter.alpha;
                        end
                        for k = rix(:)',
                                %e = (classlabel(k)==(1:M)) - [1, D(k,cix)] * weights;
                                e = CL101(k,:) - sign([1, D(k,cix)] * weights);
                                weights = weights + W(k) * [1,D(k,cix)]' * e ;
                        end
                end
                CC.weights  = sparse(size(D,2)+1,M);
                CC.weights([1,cix+1],:) = weights;
                CC.datatype = ['classifier:',lower(MODE.TYPE)];


        elseif ~isempty(strfind(lower(MODE.TYPE),'winnow'))
                % winnow algorithm
                if ~isempty(W)
                        error('Classifier (%s) does not support weighted samples.',MODE.TYPE);
                end

                [rix,cix] = row_col_deletion(D);
                [CL101,CC.Labels] = cl101(classlabel);
                M = size(CL101,2);
                weights  = ones(length(cix),M);
                theta = size(D,2)/2;

                for k = rix(:)',
                        e = CL101(k,:) - sign(D(k,cix) * weights - theta);
                        weights = weights.* 2.^(D(k,cix)' * e);
                end

                CC.weights = sparse(size(D,2)+1,M);
                CC.weights(cix+1,:) = weights;
                CC.datatype = ['classifier:',lower(MODE.TYPE)];

        elseif ~isempty(strfind(lower(MODE.TYPE),'pls')) || ~isempty(strfind(lower(MODE.TYPE),'reg'))
                % 4th version: support for weighted samples - work well with unequally distributed data:
                % regression analysis, can handle sparse data, too.

                if nargin<4,
                        W = [];
                end
                [rix, cix] = row_col_deletion(D);
                wD = [ones(length(rix),1),D(rix,cix)];

                if ~isempty(W)
                        %% wD = diag(W)*wD
                        W = W(:);
                        for k=1:size(wD,2)
                                wD(:,k) = W(rix).*wD(:,k);
                        end
                end
                [CL101, CC.Labels] = cl101(classlabel(rix,:));
                M = size(CL101,2);
                CC.weights = sparse(sz(2)+1,M);

                %[rix, cix] = row_col_deletion(wD);
                [q,r] = qr(wD,0);

                if isempty(W)
                        CC.weights([1,cix+1],:) = r\(q'*CL101);
                else
                        CC.weights([1,cix+1],:) = r\(q'*(W(rix,ones(1,M)).*CL101));
                end
                %for k = 1:M,
                %       CC.weights(cix,k) = r\(q'*(W.*CL101(rix,k)));
                %end
                CC.datatype = ['classifier:statistical:',lower(MODE.TYPE)];


        elseif ~isempty(strfind(MODE.TYPE,'WienerHopf'))
                % Q: equivalent to LDA
                % equivalent to Regression, except regression can not deal with NaN's
                [CL101,CC.Labels] = cl101(classlabel);
                M = size(CL101,2);
                CC.weights = sparse(size(D,2)+1,M);
                cc = covm(D,'E',W);
                %c1 = classlabel(~isnan(classlabel));
                %c2 = ones(sum(~isnan(classlabel)),M);
                %for k = 1:M,
                %       c2(:,k) = c1==CC.Labels(k);
                %end
                %CC.weights  = cc\covm([ones(size(c2,1),1),D(~isnan(classlabel),:)],2*real(c2)-1,'M',W);
                CC.weights  = cc\covm([ones(size(D,1),1),D],CL101,'M',W);
                CC.datatype = ['classifier:statistical:',lower(MODE.TYPE)];


        elseif ~isempty(strfind(lower(MODE.TYPE),'/gsvd'))
                if ~isempty(W)
                        error('Classifier (%s) does not support weighted samples.',MODE.TYPE);
                end
                % [2] Peg Howland and Haesun Park, 2004
                %       Generalizing Discriminant Analysis Using the Generalized Singular Value Decomposition
                %       IEEE Transactions on Pattern Analysis and Machine Intelligence, 26(8), 2004.
                %       dx.doi.org/10.1109/TPAMI.2004.46
                % [3] http://www-static.cc.gatech.edu/~kihwan23/face_recog_gsvd.htm

                [classlabel,CC.Labels] = CL1M(classlabel);
                [rix,cix] = row_col_deletion(D);

                Hw = zeros(length(rix)+length(CC.Labels), length(cix));
                Hb = [];
                m0 = mean(D(rix,cix));
                K = length(CC.Labels);
                N = zeros(1,K);
                for k = 1:K,
                        ix   = find(classlabel(rix)==CC.Labels(k));
                        N(k) = length(ix);
                        [Hw(ix,:), mu] = center(D(rix(ix),cix));
                        %Hb(k,:) = sqrt(N(k))*(mu(k,:)-m0);
                        Hw(length(rix)+k,:) = sqrt(N(k))*(mu-m0);  % Hb(k,:)
                end
                try
                        [P,R,Q] = svd(Hw,'econ');
                catch   % needed because SVD(..,'econ') not supported in Matlab 6.x
                        [P,R,Q] = svd(Hw,0);
                end
                t = rank(R);

                clear Hw Hb mu;
                %[size(D);size(P);size(Q);size(R)]
                R = R(1:t,1:t);
                %P = P(1:size(D,1),1:t);
                %Q = Q(1:t,:);
                [U,E,W] = svd(P(1:length(rix),1:t),0);
                %[size(U);size(E);size(W)]
                clear U E P;
                %[size(Q);size(R);size(W)]

                %G = Q(1:t,:)'*[R\W'];
                G = Q(:,1:t)*(R\W');   % this works as well and needs only 'econ'-SVD
                %G = G(:,1:t);  % not needed

                % do not use this, gives very bad results for Medline database
                %G = G(:,1:K); this seems to be a typo in [2] and [3].
                CC = train_sc(D(:,cix)*G,classlabel,MODE.TYPE(1:find(MODE.TYPE=='/')-1));
                CC.G = sparse(size(D,2),size(G,2));
                CC.G(cix,:) = G;
                if isfield(CC,'weights')
                        CC.weights  = sparse([CC.weights(1,:); CC.G*CC.weights(2:end,:)]);
                        CC.datatype = ['classifier:statistical:', lower(MODE.TYPE)];
                else
                        CC.datatype = [CC.datatype,'/gsvd'];
                end


        elseif ~isempty(strfind(lower(MODE.TYPE),'sparse'))
                if ~isempty(W)
                        error('Classifier (%s) does not support weighted samples.',MODE.TYPE);
                end
                % [5] J.D. Tebbens and P.Schlesinger (2006),
                %       Improving Implementation of Linear Discriminant Analysis for the Small Sample Size Problem
                %       http://www.cs.cas.cz/mweb/download/publi/JdtSchl2006.pdf

                [classlabel,CC.Labels] = CL1M(classlabel);
                [rix,cix] = row_col_deletion(D);

                warning('sparse LDA is sensitive to linear transformations')
                M = length(CC.Labels);
                G  = sparse([],[],[],length(rix),M,length(rix));
                for k = 1:M,
                        G(classlabel(rix)==CC.Labels(k),k) = 1;
                end
                tol  = 1e-10;

                G    = train_lda_sparse(D(rix,cix),G,1,tol);
                CC.datatype = 'classifier:slda';
                POS1 = find(MODE.TYPE=='/');
                %G = v(:,1:size(G.trafo,2)).*G.trafo;
                %CC.weights = s * CC.weights(2:end,:) + sparse(1,1:M,CC.weights(1,:),sz(2)+1,M);

                CC = train_sc(D(rix,cix)*G.trafo,classlabel(rix),MODE.TYPE(1:POS1(1)-1));
                CC.G = sparse(size(D,2),size(G.trafo,2));
                CC.G(cix,:) = G.trafo;
                if isfield(CC,'weights')
                        CC.weights = sparse([CC.weights(1,:); CC.G*CC.weights(2:end,:)]);
                        CC.datatype = ['classifier:statistical:',lower(MODE.TYPE)];
                else
                        CC.datatype = [CC.datatype,'/sparse'];
                end

        elseif ~isempty(strfind(lower(MODE.TYPE),'rbf'))
                if ~isempty(W)
                        error('Classifier (%s) does not support weighted samples.',MODE.TYPE);
                end

                % Martin Hieden's RBF-SVM
                if exist('svmpredict_mex','file')==3,
                        MODE.TYPE = 'SVM:LIB:RBF';
                else
                        error('No SVM training algorithm available. Install LibSVM for Matlab.\n');
                end
                CC.options = '-t 2 -q';   %use RBF kernel, set C, set gamma
                if isfield(MODE.hyperparameter,'gamma')
                        CC.options = sprintf('%s -c %g', CC.options, MODE.hyperparameter.c_value);  % set C
                end
                if isfield(MODE.hyperparameter,'c_value')
                        CC.options = sprintf('%s -g %g', CC.options, MODE.hyperparameter.gamma);  % set C
                end

                % pre-whitening
                [D,r,m]=zscore(D,1);
                CC.prewhite = sparse(2:sz(2)+1,1:sz(2),r,sz(2)+1,sz(2),2*sz(2));
                CC.prewhite(1,:) = -m.*r;

                [classlabel,CC.Labels] = CL1M(classlabel);
                CC.model = svmtrain_mex(classlabel, sparse(D), CC.options);    % Call the training mex File
                CC.datatype = ['classifier:',lower(MODE.TYPE)];

        elseif ~isempty(strfind(lower(MODE.TYPE),'svm11'))
                if ~isempty(W)
                        error('Classifier (%s) does not support weighted samples.',MODE.TYPE);
                end
                % 1-versus-1 scheme
                if ~isfield(MODE.hyperparameter,'c_value')
                        MODE.hyperparameter.c_value = 1;
                end

                CC.options=sprintf('-c %g -t 0 -q',MODE.hyperparameter.c_value);  %use linear kernel, set C
                CC.hyperparameter.c_value = MODE.hyperparameter.c_value;

                % pre-whitening
                [D,r,m]=zscore(D,1);
                CC.prewhite = sparse(2:sz(2)+1,1:sz(2),r,sz(2)+1,sz(2),2*sz(2));
                CC.prewhite(1,:) = -m.*r;

                [classlabel,CC.Labels] = CL1M(classlabel);
                CC.model = svmtrain_mex(classlabel, D, CC.options);    % Call the training mex File

                FUN = 'SVM:LIB:1vs1';
                CC.datatype = ['classifier:',lower(FUN)];


        elseif ~isempty(strfind(lower(MODE.TYPE),'psvm'))
                if ~isempty(W)
                        %%% error('Classifier (%s) does not support weighted samples.',MODE.TYPE);
                        warning('Classifier (%s) in combination with weighted samples is not tested.',MODE.TYPE);
                end
                if ~isfield(MODE,'hyperparameter')
                        nu = 1;
                elseif isfield(MODE.hyperparameter,'nu')
                        nu = MODE.hyperparameter.nu;
                else
                        nu = 1;
                end
                [m,n] = size(D);
                [CL101,CC.Labels] = cl101(classlabel);
                CC.weights = sparse(n+1,size(CL101,2));
                M = size(CL101,2);
                for k = 1:M,
                        d = sparse(1:m,1:m,CL101(:,k));
                        H = d * [ones(m,1),D];
                        %%% r = sum(H,1)';
                        r = sumskipnan(H,1,W)';
                        %%% r = (speye(n+1)/nu + H' * H)\r; %solve (I/nu+H’*H)r=H’*e
                        [HTH, nn] = covm(H,H,'M',W);
                        r = (speye(n+1)/nu + HTH)\r; %solve (I/nu+H’*H)r=H’*e
                        u = nu*(1-(H*r));
                        %%% CC.weights(:,k) = u'*H;
                        [c,nn] = covm(u,H,'M',W);
                        CC.weights(:,k) = c';
                end
                CC.hyperparameter.nu = nu;
                CC.datatype = ['classifier:',lower(MODE.TYPE)];

        elseif ~isempty(strfind(lower(MODE.TYPE),'svm:lin4'))
                if ~isfield(MODE.hyperparameter,'c_value')
                        MODE.hyperparameter.c_value = 1;
                end

                [classlabel,CC.Labels] = CL1M(classlabel);
                M = length(CC.Labels);
                CC.weights = sparse(size(D,2)+1,M);

                [rix,cix] = row_col_deletion(D);

                % pre-whitening
                [D,r,m]=zscore(D(rix,cix),1);
                sz2 = length(cix);
                s = sparse(2:sz2+1,1:sz2,r,sz2+1,sz2,2*sz2);
                s(1,:) = -m.*r;

                CC.options = sprintf('-s 4 -B 1 -c %f -q', MODE.hyperparameter.c_value);      % C-SVC, C=1, linear kernel, degree = 1,
		% TODO: support for weight vector got lost
                %model = train(W,classlabel, sparse(D), CC.options);    % C-SVC, C=1, linear kernel, degree = 1,
                model = train(classlabel, sparse(D), CC.options);    % C-SVC, C=1, linear kernel, degree = 1,
                weights = model.w([end,1:end-1],:)';

                CC.weights([1,cix+1],:) = s * weights(2:end,:) + sparse(1,1:M,weights(1,:),sz2+1,M); % include pre-whitening transformation
                CC.weights([1,cix+1],:) = s * CC.weights(cix+1,:) + sparse(1,1:M,CC.weights(1,:),sz2+1,M); % include pre-whitening transformation
                CC.hyperparameter.c_value = MODE.hyperparameter.c_value;
                CC.datatype = ['classifier:',lower(MODE.TYPE)];


        elseif ~isempty(strfind(lower(MODE.TYPE),'svm'))

                if ~isfield(MODE.hyperparameter,'c_value')
                        MODE.hyperparameter.c_value = 1;
                end
                if any(MODE.TYPE==':'),
                        % nothing to be done
                elseif exist('train','file')==3,
                        MODE.TYPE = 'SVM:LIN';        %% liblinear
                elseif exist('svmtrain_mex','file')==3,
                        MODE.TYPE = 'SVM:LIB';
                elseif (exist('svmtrain','file')==3),
                        MODE.TYPE = 'SVM:LIB';
                        fprintf(1,'You need to rename %s to svmtrain_mex.mex !! \n Press any key to continue !!!\n',which('svmtrain.mex'));
                elseif exist('svmtrain','file')==2,
                        MODE.TYPE = 'SVM:bioinfo';
                elseif exist('mexSVMTrain','file')==3,
                        MODE.TYPE = 'SVM:OSU';
                elseif exist('svcm_train','file')==2,
                        MODE.TYPE = 'SVM:LOO';
                elseif exist('svmclass','file')==2,
                        MODE.TYPE = 'SVM:KM';
                elseif exist('svc','file')==2,
                        MODE.TYPE = 'SVM:Gunn';
                else
                        error('No SVM training algorithm available. Install OSV-SVM, or LOO-SVM, or libSVM for Matlab.\n');
                end

                %%CC = train_svm(D,classlabel,MODE);
                [CL101,CC.Labels] = cl101(classlabel);
                M = size(CL101,2);
                [rix,cix] = row_col_deletion(D);
                CC.weights = sparse(sz(2)+1, M);

                % pre-whitening
                [D,r,m]=zscore(D(rix,cix),1);
                sz2 = length(cix);
                s = sparse(2:sz2+1,1:sz2,r,sz2+1,sz2,2*sz2);
                s(1,:) = -m.*r;

                for k = 1:M,
                        cl = CL101(rix,k);
                        if strncmp(MODE.TYPE, 'SVM:LIN',7);
                                if isfield(MODE,'options')
                                        CC.options = MODE.options;
                                else
                                        t = 1;
                                        if length(MODE.TYPE)>7, t=str2num(MODE.TYPE(8:end)); end
                                        if ((t<0) || (t>13)) t=0; end
                                        CC.options = sprintf('-s %i -B 1 -c %f -q',t, MODE.hyperparameter.c_value);      % C-SVC, C=1, linear kernel, degree = 1,
                                end
				% TODO: support for weight vector got lost
                                %model = train(W, cl, sparse(D), CC.options);    % C-SVC, C=1, linear kernel, degree = 1,
                                model = train(cl, sparse(D), CC.options);    % C-SVC, C=1, linear kernel, degree = 1,
                                w = model.w(:,1:end-1)';
                                Bias = model.w(:,end)';

                        elseif strcmp(MODE.TYPE, 'SVM:LIB');    %% tested with libsvm-mat-3.2.3
                                if isfield(MODE,'options')
                                        CC.options = MODE.options;
                                else
                                        CC.options = sprintf('-s 0 -c %f -t 0 -d 1 -q', MODE.hyperparameter.c_value);      % C-SVC, C=1, linear kernel, degree = 1,
                                end
                                model = svmtrain_mex(cl, sparse(D), CC.options);    % C-SVC, C=1, linear kernel, degree = 1,
                                w = -cl(1) * model.SVs' * model.sv_coef;  %Calculate decision hyperplane weight vector
                                % ensure correct sign of weight vector and Bias according to class label
                                Bias = -model.rho * cl(1);

                        elseif strcmp(MODE.TYPE, 'SVM:bioinfo');
                                % SVM classifier from bioinformatics toolbox.
                                % Settings suggested by Ian Daly, 2011-06-06
                                options = optimset('Display','iter','maxiter',20000, 'largescale','off');
                                CC.SVMstruct = svmtrain(D, cl, 'AUTOSCALE', 0, 'quadprog_opts', options, 'Method', 'LS', 'kernel_function', 'polynomial');
                                Bias = -CC.SVMstruct.Bias;
                                w = -CC.SVMstruct.Alpha'*CC.SVMstruct.SupportVectors;

                        elseif strcmp(MODE.TYPE, 'SVM:OSU');
                                [AlphaY, SVs, Bias] = mexSVMTrain(D', cl', [0 1 1 1 MODE.hyperparameter.c_value]);    % Linear Kernel, C=1; degree=1, c-SVM
                                w = -SVs * AlphaY'*cl(1);  %Calculate decision hyperplane weight vector
                                % ensure correct sign of weight vector and Bias according to class label
                                Bias = -Bias * cl(1);

                        elseif strcmp(MODE.TYPE, 'SVM:LOO');
                                [a, Bias, g, inds]  = svcm_train(D, cl, MODE.hyperparameter.c_value); % C = 1;
                                w = D(inds,:)' * (a(inds).*cl(inds)) ;

                        elseif strcmp(MODE.TYPE, 'SVM:Gunn');
                                [nsv, alpha, Bias,svi]  = svc(D, cl, 1, MODE.hyperparameter.c_value); % linear kernel, C = 1;
                                w = D(svi,:)' * alpha(svi) * cl(1);
                                Bias = mean(D*w);

                        elseif strcmp(MODE.TYPE, 'SVM:KM');
                                [xsup,w1,Bias,inds] = svmclass(D, cl, MODE.hyperparameter.c_value, 1, 'poly', 1); % C = 1;
                                w = -D(inds,:)' * w1;

                        else
                                fprintf(2,'Error TRAIN_SVM: no SVM training algorithm available\n');
                                return;
                        end

                        CC.weights(1,k) = -Bias;
                        CC.weights(cix+1,k) = w;
                end
                CC.weights([1,cix+1],:) = s * CC.weights(cix+1,:) + sparse(1,1:M,CC.weights(1,:),sz2+1,M); % include pre-whitening transformation
                CC.hyperparameter.c_value = MODE.hyperparameter.c_value;
                CC.datatype = ['classifier:',lower(MODE.TYPE)];


        elseif ~isempty(strfind(lower(MODE.TYPE),'csp'))
                CC.datatype = ['classifier:',lower(MODE.TYPE)];
                [classlabel,CC.Labels] = CL1M(classlabel);
                CC.MD = repmat(NaN,[sz(2)+[1,1],length(CC.Labels)]);
                CC.NN = CC.MD;
                for k = 1:length(CC.Labels),
                        %% [CC.MD(k,:,:),CC.NN(k,:,:)] = covm(D(classlabel==CC.Labels(k),:),'E');
                        ix = classlabel==CC.Labels(k);
                        if isempty(W)
                                [CC.MD(:,:,k),CC.NN(:,:,k)] = covm(D(ix,:), 'E');
                        else
                                [CC.MD(:,:,k),CC.NN(:,:,k)] = covm(D(ix,:), 'E', W(ix));
                        end
                end
                ECM = CC.MD./CC.NN;
                W   = csp(ECM,'CSP3');
                %%% ### This is a hack ###
                CC.FiltA = 50;
                CC.FiltB = ones(CC.FiltA,1);
                d   = filtfilt(CC.FiltB,CC.FiltA,(D*W).^2);
                CC.csp_w = W;
                CC.CSP = train_sc(log(d),classlabel);


        else          % Linear and Quadratic statistical classifiers
                CC.datatype = ['classifier:statistical:',lower(MODE.TYPE)];
                [classlabel,CC.Labels] = CL1M(classlabel);
                CC.MD = repmat(NaN,[sz(2)+[1,1],length(CC.Labels)]);
                CC.NN = CC.MD;
                for k = 1:length(CC.Labels),
                        ix = classlabel==CC.Labels(k);
                        if isempty(W)
                                [CC.MD(:,:,k),CC.NN(:,:,k)] = covm(D(ix,:), 'E');
                        else
                                [CC.MD(:,:,k),CC.NN(:,:,k)] = covm(D(ix,:), 'E', W(ix));
                        end
                end

                ECM = CC.MD./CC.NN;
                NC  = size(CC.MD);
                if strncmpi(MODE.TYPE,'LD',2) || strncmpi(MODE.TYPE,'FDA',3) || strncmpi(MODE.TYPE,'FLDA',3),

                        %if NC(1)==2, NC(1)=1; end                % linear two class problem needs only one discriminant
                        CC.weights = repmat(NaN,NC(2),NC(3));     % memory allocation
                        type = MODE.TYPE(3)-'0';

                        ECM0 = squeeze(sum(ECM,3));  %decompose ECM
                        for k = 1:NC(3);
                                ix = [1:k-1,k+1:NC(3)];
                                dM = CC.MD(:,1,k)./CC.NN(:,1,k) - sum(CC.MD(:,1,ix),3)./sum(CC.NN(:,1,ix),3);
                                switch (type)
                                        case 2          % LD2
                                                ecm0 = (sum(ECM(:,:,ix),3)/(NC(3)-1) + ECM(:,:,k));
                                        case 4          % LD4
                                                ecm0 = 2*(sum(ECM(:,:,ix),3) + ECM(:,:,k))/NC(3);
                                                % ecm0 = sum(CC.MD,3)./sum(CC.NN,3);
                                        case 5          % LD5
                                                ecm0 = ECM(:,:,k);
                                        case 6          % LD6
                                                ecm0 = sum(CC.MD(:,:,ix),3)./sum(CC.NN(:,:,ix),3);
                                        otherwise       % LD3, LDA, FDA
                                                ecm0 = ECM0;
                                end
                                if isfield(MODE.hyperparameter,'gamma')
                                        ecm0 = ecm0 + mean(diag(ecm0))*eye(size(ecm0))*MODE.hyperparameter.gamma;
                                end

                                CC.weights(:,k) = ecm0\dM;

                        end
                        %CC.weights = sparse(CC.weights);

                elseif strcmpi(MODE.TYPE,'RDA');
                        if isfield(MODE,'hyperparameter')
                                CC.hyperparameter = MODE.hyperparameter;
                        end
                        % default values
                        if ~isfield(CC.hyperparameter,'gamma')
                                CC.hyperparameter.gamma = 0;
                        end
                        if ~isfield(CC.hyperparameter,'lambda')
                                CC.hyperparameter.lambda = 1;
                        end
                else
                        ECM0 = sum(ECM,3);
                        nn = ECM0(1,1,1);       % number of samples in training set for class k
                        XC = squeeze(ECM0(:,:,1))/nn;           % normalize correlation matrix
                        M  = XC(1,2:NC(2));             % mean
                        S  = XC(2:NC(2),2:NC(2)) - M'*M;% covariance matrix

                        try
                                [v,d]=eig(S);
                                U0 = v(diag(d)==0,:);
                                CC.iS2 = U0*U0';
                        end

                        %M  = M/nn; S=S/(nn-1);
			v=eig(S);
			if all(v > eps)
	                        ICOV0 = inv(S);
			else
				ICOV0 = repmat(NaN,size(S));
			end
                        CC.iS0 = ICOV0;
                        % ICOV1 = zeros(size(S));
                        for k = 1:NC(3),
                                %[M,sd,S,xc,N] = decovm(ECM{k});  %decompose ECM
                                %c  = size(ECM,2);
                                nn = ECM(1,1,k);% number of samples in training set for class k
                                XC = squeeze(ECM(:,:,k))/nn;% normalize correlation matrix
                                M  = XC(1,2:NC(2));% mean
                                S  = XC(2:NC(2),2:NC(2)) - M'*M;% covariance matrix
                                %M  = M/nn; S=S/(nn-1);

                                %ICOV(1) = ICOV(1) + (XC(2:NC(2),2:NC(2)) - )/nn
				%if all(eig(S) > eps)	% only when S is positive definite
				v=eig(S);
				if all(v > eps)
					iS = inv(S);
					detS = det(S);
				else
					iS = repmat(NaN,size(S));
					detS = NaN;
				end

                                CC.M{k}   = M;
                                CC.IR{k}  = [-M;eye(NC(2)-1)]*iS*[-M',eye(NC(2)-1)];  % inverse correlation matrix extended by mean
                                CC.IR0{k} = [-M;eye(NC(2)-1)]*ICOV0*[-M',eye(NC(2)-1)];  % inverse correlation matrix extended by mean
                                d = NC(2)-1;
                                if exist('OCTAVE_VERSION','builtin')
                                        S = full(S);
                                end
                                CC.logSF(k)  = log(nn) - d/2*log(2*pi) - detS/2;
                                CC.logSF2(k) = -2*log(nn/sum(ECM(:,1,1)));
                                CC.logSF3(k) = d*log(2*pi) + log(det(S));
                                CC.logSF4(k) = log(detS) + 2*log(nn);
                                CC.logSF5(k) = log(detS);
                                CC.logSF6(k) = log(detS) - 2*log(nn/sum(ECM(:,1,1)));
                                CC.logSF7(k) = log(detS) + d*log(2*pi) - 2*log(nn/sum(ECM(:,1,1)));
                                CC.logSF8(k) = sum(log(svd(S))) + log(nn) - log(sum(ECM(:,1,1)));
                                CC.SF(k) = nn/sqrt((2*pi)^d * detS);
                                %CC.datatype='LLBC';
                        end
                end
        end
end

function [CL101,Labels] = cl101(classlabel)
        %% convert classlabels to {-1,1} encoding

        if (all(classlabel>=0) && all(classlabel==fix(classlabel)) && (size(classlabel,2)==1))
                M = max(classlabel);
                if 0, % M==2;
			% this approach is disabled for now because it currently 
  			% is not compatible with all supported classifiers (e.g. NBC). 
                        CL101 = (classlabel==2)-(classlabel==1);
                else
                        CL101 = zeros(size(classlabel,1),M);
                        for k=1:M,
                                %% One-versus-Rest scheme
                                CL101(:,k) = 2*real(classlabel==k) - 1;
                        end
                end
                CL101(isnan(classlabel),:) = NaN; %% or zero ???

        elseif all((classlabel==1) | (classlabel==-1)  | (classlabel==0) )
                CL101 = classlabel;
                M = size(CL101,2);
        else
                classlabel,
                error('format of classlabel unsupported');
        end
        Labels = 1:M;
        return;
end


function [cl1m, Labels] = CL1M(classlabel)
        %% convert classlabels to 1..M encoding
        if (all(classlabel>=0) && all(classlabel==fix(classlabel)) && (size(classlabel,2)==1))
                cl1m = classlabel;

        elseif all((classlabel==1) | (classlabel==-1) | (classlabel==0) )
                CL101 = classlabel;
                M = size(classlabel,2);
                if any(sum(classlabel==1,2)>1)
                        warning('invalid format of classlabel - at most one category may have +1');
                end
                if (M==1),
                        cl1m = (classlabel==-1) + 2*(classlabel==+1);
                else
                        [tmp, cl1m] = max(classlabel,[],2);
                        if any(tmp ~= 1)
                                warning('some class might not be properly represented - you might what to add another column to classlabel = [max(classlabel,[],2)<1,classlabel]');
                        end
                        cl1m(tmp<1)= 0; %% or NaN ???
                end
        else
                classlabel
                error('format of classlabel unsupported');
        end
        Labels = 1:max(cl1m);
        return;
end
