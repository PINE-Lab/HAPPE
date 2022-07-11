function [W] = bss(data,Mode,M,maxlag)
%  BSS is a wrapper for various Blind Source Separation algorithms
%
%  W = bss(data, Mode)
%  W = bss(data, Mode, M)
%  W = bss(data, Mode, [], maxlag)
%  W = bss(data, Mode, M, maxlag)
%
%  W: 	   unmixing matrix
%  data:   signal data (each column is a channel)
%  M: 	   [optional] number of components.
%  maxlag: maximum lag for time delayed separation 	
%  Mode:   algorithm used. Currently are supported: 
%	PCA
%	FastICA [6]
%	JADE [1-3]
%	NGCA
% 	FFDIAG [4-5]
%	TDSEP (old not recommended)
%	TDSEP1 [4-5]
%	TDSEP3 [4-5]
%
% References:
% [1] Cardoso, Jean-François; Souloumiac, Antoine (1993).
%   Blind beamforming for non-Gaussian signals.
%   IEE Proceedings F (Radar and Signal Processing). 140 (6): 362–370.
% [2] http://perso.telecom-paristech.fr/~cardoso/guidesepsou.html
% [3] http://perso.telecom-paristech.fr/~cardoso/Algo/Jade/jade.m
% [4] A. Ziehe, G.Nolte, K-R. Mueller
%   A Fast Algorithm for Joint Diagonalization with Non-orthogonal
%   Transformations and its Application to Blind Source Separation.
%   Journal of Machine Learning Research 5 (2004) 777–800
%   http://www.jmlr.org/papers/volume5/ziehe04a/ziehe04a.pdf
% [5] http://www.user.tu-berlin.de/aziehe/code/
% [6] https://research.ics.aalto.fi/ica/fastica/

%	Copyright (C) 2007,2018 by Alois Schloegl
%	This is part of the BIOSIG-toolbox http://biosig.sf.net/

% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% Version 3 of the License, or (at your option) any later version.


if ~exist('jade','file')
	addpath([getenv('HOME'),'/matlab/other/jade']);
end;
if ~exist('NGCA','file')
	addpath('/home/neuro/schalo/cvs/neuro_cvs/Tex/codename_spok/code/JMLR_paper'); 
end; 
if ~exist('tdsep','file')
%	addpath('/home/neuro/schalo/cvs/neuro_cvs/matlab/meg_lab/'); 
end; 
if ~exist('tdsep0','file')
	addpath('/home/neuro/schalo/cvs/neuro_cvs/matlab/bci/ica/'); 
end; 
if ~exist('tdsep3','file')
	addpath('/home/neuro/schalo/matlab/other/motoaki/'); 
end; 
if ~exist('jadeR','file')
	addpath('/home/data/share/ica/'); 
end; 
if ~exist('ffdiag2','file')
	addpath('/home/data/share/ica/ffdiag/'); 
end; 

r = rank(data); 
m = size(data,2);
if nargin<3,
	M = r;
end; 	
if isempty(M),
	M = r;
end; 	

Mode = upper(Mode); 

[u,s,v]=svd(data,0); 
V = v*inv(s);

data= center(data);
V0 = inv(sqrtm(covm(data,'D')));   %% whitening or sphering
data = data*V0';

if 0, 
elseif strcmp(Mode,'PCA')
	[u,s,W] = svd(data,M); % get k (i.e. FLAG.PCA) largests PCA's

elseif strcmp(Mode,'JADE')
	W = jade(data', M)'; 

elseif strcmp(Mode,'FastICA')
	%[A,W] = fastica(data'); 
	[A,W] = fastica(data','numOfIC',M); 
	W=W'; 

elseif strcmp(Mode,'NGCA')
	[W, projdata, signalmatrix] = NGCA(data',M);	

elseif strcmp(Mode,'FFDIAG')
	C   = zeros(m,m,maxlag+1);
	for k = 1:maxlag+1,
		C(:,:,k) = tdCorr(data', k-1);
	end
	[W,d] = ffdiag2(C);
elseif strcmp(Mode,'TDSEP')
	[W,d] = tdsep(data',0:maxlag);

elseif strcmp(Mode,'TDSEP0')
	[W] = tdsep0(data',0:maxlag);

elseif strcmp(Mode,'TDSEP3')
	[W,d] = tdsep3(data',0:maxlag);
end;

W = V0*W; 


