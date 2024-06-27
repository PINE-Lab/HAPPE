function [R,sig,ci1,ci2] = partcorrcoef(X,Y,Z,Mode)
% PARTCORRCOEF calculates the partial correlation between X and Y
% after removing the influence of Z.
% X, Y and Z can contain missing values encoded with NaN.
% NaN's are skipped, NaN do not result in a NaN output. 
% (Its assumed that the occurence of NaN's is uncorrelated) 
% The output gives NaN, only if there are insufficient input data.
%
%  The partial correlation  is defined as 
%  pcc(xy|z)=(cc(x,y)-cc(x,z)*cc(y,z))/sqrt((1-cc(x,y)�)*((1-cc(x,z)�)))
%
%
% PARTCORRCOEF(X [,Mode]);
%      calculates the (auto-)correlation matrix of X
% PARTCORRCOEF(X,Y,Z);
% PARTCORRCOEF(X,Y,Z,[]);
% PARTCORRCOEF(X,Y,Z,'Pearson');
% PARTCORRCOEF(X,Y,Z,'Rank');
% PARTCORRCOEF(X,Y,Z,'Spearman');
%
% Mode=[] [default]
%	removes from X and Y the part that can be explained by Z
%	and computes the correlation of the remaining part. 
% 	Ideally, this is equivalent to Mode='Pearson', however, in practice
%	this is more accurate.
% Mode='Pearson' or 'parametric'
% Mode='Spearman'
% Mode='Rank'
%	computes the partial correlation based on cc(x,y),cc(x,z) and cc(y,z) 
%	with the respective mode. 
%
% [R,p,ci1,ci2] = PARTCORRCOEF(...);
%  r is the partialcorrelation matrix
%	r(i,j) is the partial correlation coefficient r between X(:,i) and Y(:,j) 
%	when influence of Z is removed. 
%  p    gives the significance of PCC
%	It tests the null hypothesis that the product moment correlation coefficient is zero 
%       using Student's t-test on the statistic t = r sqrt(N-Nz-2)/sqrt(1-r^2) 
%       where N is the number of samples (Statistics, M. Spiegel, Schaum series).
%  p > alpha: do not reject the Null hypothesis: "R is zero".
%  p < alpha: The alternative hypothesis "R2 is larger than zero" is true with probability (1-alpha).
%  ci1	lower 0.95 confidence interval 
%  ci2	upper 0.95 confidence interval 
%
% see also: SUMSKIPNAN, COVM, COV, COR, SPEARMAN, RANKCORR, RANKS, CORRCOEF
%
% REFERENCES:
% on the partial correlation coefficient 
% [1] http://www.tufts.edu/~gdallal/partial.htm
% [2] http://www.nag.co.uk/numeric/fl/manual/pdf/G02/g02byf.pdf

%    $Id$
%    Copyright (C) 2000-2002,2009 by Alois Schloegl <alois.schloegl@gmail.com>	
%    This function is part of the NaN-toolbox
%    http://pub.ist.ac.at/~schloegl/matlab/NaN/

%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA

% Features:
% + interprets NaN's as missing value
% + Pearson's correlation
% + Spearman's rank correlation
% + Rank correlation (non-parametric, non-Spearman)
% + is fast, using an efficient algorithm O(n.log(n)) for calculating the ranks
% + significance test for null-hypthesis: r=0 
% + confidence interval (0.99) included
% - rank correlation works for cell arrays, too (no check for missing values).
% + compatible with Octave and Matlab


if nargin==3
        Mode=[];
elseif nargin==4,
else
        error('Error PARTCORRCOEF: Missing argument(s)\n');
end;        

if isempty(Z)
	R = corrcoef(X,Y,Mode);

elseif isempty(Mode) 
	if ~isempty(Z)
	for j=1:size(X,2)
		ix = ~any(isnan(Z),2) & ~isnan(X(:,j));
		X(:,j) = X(:,j) - Z*(Z(ix,:)\X(ix,j));
	end; 	
	for j=1:size(Y,2)
		ix = ~any(isnan(Z),2) & ~isnan(Y(:,j));
		Y(:,j) = Y(:,j) - Z*(Z(ix,:)\Y(ix,j));
	end;
	end;
	R = corrcoef(X,Y,Mode);

else 
	rxy = corrcoef(X,Y,Mode);
	rxz = corrcoef(X,Z,Mode);
	if isempty(Y),
        	ryz = rxz;
	else
        	ryz = corrcoef(Y,Z,Mode);
	end;

	%rxy,rxz,ryz 
	R = (rxy-rxz*ryz')./sqrt((1-rxz.^2)*(1-ryz.^2)');
	
end;

if nargout<2, 
        return, 
end;

% SIGNIFICANCE TEST
%warning off; 	% prevent division-by-zero warnings in Matlab.
NN=size(X,1)-size(Z,2);

tmp = 1 - R.*R;
tmp(tmp<0) = 0;		% prevent tmp<0 i.e. imag(t)~=0 
t   = R.*sqrt(max(NN-2,0)./tmp);

if exist('t_cdf','file')
        sig = t_cdf(t,NN-2);
elseif exist('tcdf','file')
        sig = tcdf(t,NN-2);
else
        fprintf('Warning CORRCOEF: significance test not completed because of missing TCDF-function\n')
        sig = repmat(nan,size(R));
end;
sig  = 2 * min(sig,1 - sig);

if nargout<3, 
        return, 
end;


% CONFIDENCE INTERVAL
if exist('flag_implicit_significance','file'),
        alpha = flag_implicit_significance;
else
	alpha = 0.01;        
end;

fprintf(1,'CORRCOEF: confidence interval is based on alpha=%f\n',alpha);

tmp = R;
%tmp(ix1 | ix2) = nan;		% avoid division-by-zero warning
z   = log((1+tmp)./(1-tmp))/2; 	% Fisher's z-transform; 
%sz  = 1./sqrt(NN-3);		% standard error of z
sz  = sqrt(2)*erfinv(1-2*alpha)./sqrt(NN-3);		% confidence interval for alpha of z

ci1 = tanh(z-sz);
ci2 = tanh(z+sz);



