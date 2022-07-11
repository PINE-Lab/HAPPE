function [CLASS,ERR,POSTERIOR,LOGP,COEF]=classify(sample,training,classlabel,TYPE)
% CLASSIFY classifies sample data into categories 
% defined by the training data and its group information 
%
%  CLASS = classify(sample, training, group) 
%  CLASS = classify(sample, training, group, TYPE) 
%  [CLASS,ERR,POSTERIOR,LOGP,COEF] = CLASSIFY(...) 
%
%  CLASS contains the assigned group. 
%  ERR is the classification error on the training set weighted by the 
%	prior propability of each group. 
%
%  The same classifier as in TRAIN_SC are supported. 
%
% ATTENTION: no cross-validation is applied, therefore the 
%    classification error is too optimistic (overfitting). 
%    Use XVAL instead to obtain cross-validated performance. 
% 
% see also: TRAIN_SC, TEST_SC, XVAL
%
% References: 
% [1] R. Duda, P. Hart, and D. Stork, Pattern Classification, second ed. 
%       John Wiley & Sons, 2001. 

%	$Id$
%	Copyright (C) 2008,2009 by Alois Schloegl <alois.schloegl@gmail.com>
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

if nargin<4
	TYPE = 'linear';
end; 

if strcmp(TYPE,'linear')
	TYPE = 'LDA';
elseif strcmp(TYPE,'quadratic')
	TYPE = 'QDA2'; % result is closer to Matlab 
elseif strcmp(TYPE,'diagLinear')
	TYPE = 'NBC';
elseif strcmp(TYPE,'diagQuadratic')
	TYPE = 'NBC';
elseif strcmp(TYPE,'mahalanobis')
	TYPE = 'MDA';
end; 	

[group,I,classlabel] = unique(classlabel);

CC = train_sc(training,classlabel,TYPE); 
R  = test_sc(CC,sample);
CLASS = group(R.classlabel); 

if nargout>1,
	R  = test_sc(CC,training,[],classlabel);
	ERR = 1-R.ACC; 
end; 

if nargout>2,
	warning('output arguments POSTERIOR,LOGP and COEF not supported')
	POSTERIOR = [];
	LOGP = [];
	COEF = [];
end; 

