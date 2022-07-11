function [o] = bci3eval(x1,x2,DIM)
% BCI3eval evaluations a BCI-result as suggested in [1,2].
%   - It returns the classification error, the signal to noise ratio,
%   the mutual information, as well as mean, standard error and
%   standard deviation for both classes.
%   - time course of these resulting parameters are supported
%   - Missing values can be encoded as NaN.
%
%
% [o] = bci3eval(x1, x2 [,DIM])
%
% x1 is the bci output for class 1
% x2 is the bci output for class 2
% o is a struct with various results
%
%
% see also: SUMSKIPNAN, PLOTA
%
% REFERENCES:
%  [1] Schlögl A., Neuper C. Pfurtscheller G.
%	Estimating the mutual information of an EEG-based Brain-Computer-Interface
%	Biomedizinische Technik 47(1-2): 3-8, 2002.
%  [2] A. Schlögl, C. Keinrath, R. Scherer, G. Pfurtscheller,
%	Information transfer of an EEG-based Bran-computer interface.
%	Proceedings of the 1st International IEEE EMBS Conference on Neural Engineering, pp.641-644, Mar 20-22, 2003.
%  [3] Evaluation of the dataset III of the BCI-competition 2003.
%	http://ida.first.gmd.de/~blanker/competition/results/TR_BCI2003_III.pdf
%  [4] BIOSIG.SF.NET


%    Copyright (C) 2003 by Alois Schloegl <alois.schloegl@gmail.com>

%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

%    Copyright (C) 2000-2003 by  Alois Schloegl <alois.schloegl@gmail.com>


if nargin<3,
        DIM=min(find(size(x1)>1));
        if isempty(DIM), DIM=1; end;
end;


% classification error
if DIM==1,
	o.ERR = (1-mean(sign([-x1;x2]),DIM))/2;
elseif DIM==2,
	o.ERR = (1-mean(sign([-x1,x2]),DIM))/2;
end;


%%%%% 2nd order statistics
[i1.SUM,i1.N,i1.SSQ] = sumskipnan(x1,DIM);
[i2.SUM,i2.N,i2.SSQ] = sumskipnan(x2,DIM);

o.MEAN1 = i1.SUM./i1.N;	% mean
v1    = i1.SSQ-i1.SUM.*o.MEAN1;	% n*var
o.SD1 = sqrt(v1./i1.N); % standard deviation
%o.SE1 = sqrt(v1)./i1.N; % standard error of the mean

o.MEAN2 = i2.SUM./i2.N;
v2    = i2.SSQ-i2.SUM.*o.MEAN2;
o.SD2 = sqrt(v2./i2.N);
%o.SE2 = sqrt(v2)./i2.N;


%%%%% Signal-to-Noise Ratio

	% intra-class variability
if DIM==1,
	vd = var([-x1;x2],[],DIM);
elseif DIM==2,
	vd = var([-x1,x2],[],DIM);
end;

o.SNR = 1/4*(o.MEAN2-o.MEAN1).^2./vd;


%%%%% Mutual Information
o.I   = 1/2*log2(o.SNR+1);


o.datatype = 'TSD_BCI7';  % useful for PLOTA
