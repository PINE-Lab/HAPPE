function RES=mod_optimal_detection_filter(data, C, maxlag, ixtrain, ixtest, method)
% MOD_OPTIMAL_DETECTION_FILTER for detections miniature events [1,2]. The parameters
%    of the detection method (e.g. filter coefficents, detection thresheld, and the 
%    time shift) are estimated from events defined as "training set", and 
%    and was tested on events declared as "test set". Several performances 
%    measures (including AUC and Cohen's kappa) are reported. In case the training
%    and test set did not overlap, crossvalidation is implemented and
%    it can be expeted that the the obtained performance can be reached 
%    also on unseen data. Several methods are implemented, the default and 
%    recommended method ('wopt2') is based on an a-causal Wiener filter 
%    as described in [1,2]. The method tries different time shift values for 
%    training the wiener filter, and uses the filter that leads to the 
%    largest AUC. 
%
% Usage:
%    RES = get_minidet_classification(data, C, maxlag, ixtrain, ixtest, method): 
%    Input:  
%	data: raw data trace (single column vector) 
%  		the data may contain NaN's indicating missing values.
%       C:    scoring trace 
%       maxlag: maximum lag numbers for computing the correlation functions 
% 		it should be the some of the filter length (typically 1000 samples)
% 		and the maximum timeshift.
%	ixtrain: list of sample numbers used for training the detector
%
%	ixtest:  list of sample numbers for testing the performance of the detector
%
% 	method: 'wopt2' (a.k.a. 'MOD') is the default method, and the recommended method
%		because it provides the best performance [1] so far.  
%		other approaches were used for comparison and are not recommended. 
%   Output:
%	RES.CLASSIFIER.delay
%		time shift for making the a-causal filter causal.  
%	RES.CLASSIFIER.A
%		coefficients of the Wiener filter  used for 
%		u = filter(RES.CLASSIFIER.A, 1, data)
%
% 	RES.WOPT.DLIST =
%		list of time shift values used to train the optimum filter
%		currently this list is [-maxlag/4:25:maxlag], 
%	RES.WOPT.AUC   = AUC;
%		AUC result of the corresponding time shift values. 
%	plot(RES.WOPT.DLIST/Fs, RES.WOPT.AUC) can be used to produce 
%       Suppl Figure 3 [1].  
%
% see also: 
%	AUC, MINIDET2GDF, GET_LOCAL_MAXIMA_ABOVE_THRESHOLD, 
%
% Requirements: 
% - Octave or Matlab
%	http://www.octave.org
% - NaN toolbox
%       https://pub.ist.ac.at/~schloegl/matlab/NaN/	
%	https://octave.sourceforge.io/nan/index.html
% - Biosig for Octave and Matlab
%	http://biosig.sourceforge.net/download.html
%
% References: 
% [1] Zhang X, Schlögl A, Vandael D, Jonas P, 
%     MOD: A novel machine-learning optimal-filtering method for accurate and efficient
%        detection of subthreshold synaptic events in vivo.
%     Journal of Neuroscience Methods, 2021.
%     doi:10.1016/j.jneumeth.2021.109125
% [2] Zhang X., Schloegl A., Vandael D., Jonas P. (2020)
%      A novel machine learning-based method for accurate and efficient detection of
%      subthreshold synaptic events in vivo and in vitro (in revision).
%      Journal of Neuroscience Methods. 

%    Copyright (C) 2016-2022 by Alois Schloegl <alois.schloegl@ist.ac.at>
%    This is part of the BIOSIG-toolbox http://biosig.sf.net/
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


if nargin<6,
	method.preproc='wiener2';
elseif ischar(method)
	tmp=method; 
	clear method; 
	method.preproc=tmp;
end;


%NAN_BLOCK=repmat(NaN,maxlag,1);
%data = [NAN_BLOCK; data; NAN_BLOCK];
%C    = [NAN_BLOCK; C;    NAN_BLOCK];
%ixtrain = ixtrain+maxlag; 
%ixtest  = ixtest+maxlag;

if isempty(ixtest)
	ixtest  = ixtrain; 
	ixtest0 = [2*maxlag:size(data,1)]';
else
	ixtest0 = ixtest;
end


% if isfield(method,'preproc') && strcmp(method.preproc,'wopt2') && (size(data,2)==1)
if 1
	ID0=tic(); cput0=cputime(); 

	Rdc=repmat(NaN,2*maxlag+1,1);
	Rdd=repmat(NaN,2*maxlag+1,1);

	dl=length(data); 
	if (any(ixtest > dl-maxlag) || any(ixtest <= maxlag))
		warning('ixtest should be limited to the range MAXLAG+1:length(data)-MAXLEN')
	end
	if (any(ixtrain > dl-maxlag) || any(ixtrain <= maxlag))
		warning('ixtrain should be limited to the range MAXLAG+1:length(data)-MAXLEN')
	end

	c = C-mean(C(ixtrain));

	TTLabel{1} = 'length-of-data';
	TTLabel{2} = 'number-of-lags';
	TTLabel{3} = 'length-ixtrain';
	TTLabel{4} = 'sum-ixtrain';
	TT=[size(data,1), 3*maxlag,length(ixtrain),sum(ixtrain)];
 
	ID1=tic(); cput1=cputime(); 

	if exist('accovf_mex','file')
		[Sxx,Nxx,Sxy,Nxy]=accovf_mex(data, c, 2*maxlag, ixtrain);
		Rdc = (Sxy./Nxy);
		Rdd = (Sxx./Nxx);
	else
		data((end+1) : (ixtest0(end)+maxlag*1.5)) = NaN;
		for k = -maxlag*2:maxlag*2,
			Rdc(k+1+maxlag*2) = mean(data(ixtrain-k).*c(ixtrain));
			Rdd(k+1+maxlag*2) = mean(data(ixtrain-k).*data(ixtrain));
		end;
	end

	TT=[TT, toc(ID1), cputime()-cput1]; ID1=tic(); cput1=cputime();

	TTLabel{end+1} = 'wallclock-acf+ccf';
	TTLabel{end+1} = 'cputime-acf+ccf';

	TT2Label={'wallclock-toeplitz'};
	TT2Label{end+1} = 'cputime-toeplitz';
	TT2Label{end+1} = 'wallclock-filter(*2)';
	TT2Label{end+1} = 'cputime-filter(*2)';
	TT2Label{end+1} = 'wallclock-auc(*2)';
	TT2Label{end+1} = 'cputime-auc(*2)';

	A=[];
	B=[];
	%DLIST=0:5:1000; 
        if isfield(method,'DLIST')
                DLIST=method.DLIST; % 10 ms
        else
        	DLIST=-maxlag/4:25:maxlag; 
        end
	AUC=repmat(NaN,length(DLIST),4);

	for k=1:length(DLIST)
		ID2=tic(); cput2=cputime(); 

		delay = DLIST(k);
		A = toeplitz(Rdd(maxlag*2+(0:maxlag))) \ Rdc(maxlag*2+1-delay+[0:maxlag]);
		B = toeplitz(Rdd(maxlag*2+(0:maxlag))) \ Rdc(maxlag*2+1+delay-[0:maxlag]);

		TT2(k,1:2)=[toc(ID2), cputime()-cput2]; ID2=tic(); cput2=cputime(); 
                 
                out1a = filter(A,1,data);
                %out1b = filter(B,1,data);

		TT2(k,3:4)=[toc(ID2), cputime()-cput2]; ID2=tic(); cput2=cputime(); 

                %AUC(k,1:2)=[auc(out1a(ixtest+delay),C(ixtest)), auc(out1b(ixtest+delay),C(ixtest))];

		ixtmp      = ixtest( 0<(ixtest+delay) & (ixtest+delay)<=length(C) );
		tmp        = roc(out1a(ixtmp+delay),C(ixtmp));
		AUC(k,1)   = tmp.AUC;

		TT2(k,5:6) = [toc(ID2), cputime()-cput2]; ID2=tic(); cput2=cputime();

		if 0, %try
	                out2a = filter(A,1,data);
        	        out2b = filter(B,1,data);
        	        AUC(k,3:4)=[auc(out2a(ixtest+delay),C(ixtest)), auc(out2b(ixtest+delay),C(ixtest))];
		end

		RES.WOPT.A(:,k) = A;
		RES.WOPT.B(:,k) = B;
	end

	TT=[TT, k, toc(ID1), cputime()-cput1]; ID1=tic(); cput1=cputime(); 
	TTLabel{end+1}='number-of-delays';
	TTLabel{end+1}='wallclock-find-delay(wiener+auc)';
	TTLabel{end+1}='cputime-find-delay(wiener+auc)';

	RES.WOPT.DLIST = DLIST(:);
	RES.WOPT.AUC   = AUC;
	% AUC,
	[tmpAUC,tmpix]=max(AUC,[],1);
	[tmpm,mix] = max(tmpAUC);
	if 1, %mix==1
		RES.CLASSIFIER.A = RES.WOPT.A(:,tmpix(mix));
		RES.CLASSIFIER.delay = DLIST(tmpix(mix));
	elseif mix==2
		RES.CLASSIFIER.A = RES.WOPT.B(:,tmpix(mix));
		RES.CLASSIFIER.delay = DLIST(tmpix(mix));
	elseif mix==3
		RES.CLASSIFIER.A = RES.WOPT.A(:,tmpix(mix));
		RES.CLASSIFIER.delay = -DLIST(tmpix(mix));
	elseif mix==4
		RES.CLASSIFIER.A = RES.WOPT.B(:,tmpix(mix));
		RES.CLASSIFIER.delay = -DLIST(tmpix(mix));
	end
	RES.WOPT.delay=RES.CLASSIFIER.delay;

	ID1=tic(); cput1=cputime(); 

	out   = filter(RES.CLASSIFIER.A, 1, data);

	TT=[TT, toc(ID1), cputime()-cput1]; ID1=tic(); cput1=cputime(); 
	TTLabel{end+1}='wallclock-filtering-all-data';
	TTLabel{end+1}='cputime-filtering-all-data';

	ixtrain   = ixtrain(0 < (ixtrain+RES.CLASSIFIER.delay));
	ixtrain   = ixtrain((ixtrain+RES.CLASSIFIER.delay)<=length(out));
	ixtest0   = ixtest0(0 < (ixtest0+RES.CLASSIFIER.delay));
	ixtest0   = ixtest0((ixtest0+RES.CLASSIFIER.delay)<=length(out));
	ROC.train = roc(out(ixtrain+RES.CLASSIFIER.delay), C(ixtrain));
	ROC.test  = roc(out(ixtest0+RES.CLASSIFIER.delay), C(ixtest0));
	RES.CLASSIFIER.CM        = ROC.test.H_kappa;
	RES.CLASSIFIER.THRESHOLD = ROC.test.THRESHOLD.maxKAPPA;
	RES.CLASSIFIER.AUC       = ROC.test.AUC;

        RES.train.AUC  = ROC.train.AUC;
        RES.test.AUC   = ROC.test.AUC;
        RES.THRESHOLD  = ROC.train.THRESHOLD.maxKAPPA;
        RES.train.CM   = ROC.train.H_kappa;
        RES.test.CM    = ROC.test.H_kappa;
        [RES.train.Kappa,sd,H] = kappa(ROC.train.H_kappa);
        [RES.test.Kappa,sd,H]  = kappa(ROC.test.H_kappa);

	TT=[TT toc(ID1), cputime()-cput1]; ID1=tic(); cput1=cputime(); 
	TTLabel{end+1}='wallclock-auc';
	TTLabel{end+1}='cputime-auc';

	RES.output = out(ixtest0+RES.CLASSIFIER.delay);

	TT=[TT, toc(ID0), cputime()-cput0, mean(TT2)/2]; 
	TTLabel{end+1}='wallclock-total';
	TTLabel{end+1}='cputime-total';

	TTLabel{end+1}='wallclock-toeplitz-average(x2)';
	TTLabel{end+1}='cputime-toeplitz-average(x2)';
	TTLabel{end+1}='wallclock-filter-training-data-average(x2)';
	TTLabel{end+1}='cputime-filter-training-data-average(x2)';
	TTLabel{end+1}='wallclock-auc-average(x2)';
	TTLabel{end+1}='cputime-auc-average(x2)';

	RES.TT=TT;
	RES.TT2=TT2;
	RES.TTLabel=TTLabel;
	RES.TT2Label=TT2Label;

        fprintf(1,'TimeLabels:\t');
	fprintf(1,'%s\t',TTLabel{:});
	fprintf(1,'%s\nTiming:\t');
	fprintf(1,'%g\t',TT);
        fprintf(1,'\n');

else 
	whos
	method
	error('invalid input argumments')

end


