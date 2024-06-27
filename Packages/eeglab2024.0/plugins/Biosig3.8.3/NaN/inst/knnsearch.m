function [idx, dist]=knnsearch(X,Y,varargin)
% KNNSEARCH search for K nearest neighbors
%   and related statistics 
%
%  Usage: 
%     IDX = knnsearch(X,Y);
%	  finds for each element (row) in Y, the nearest 
% 	  of all elements in X, such that 
%         IDX(k) points to X(IDX(k),:) that is nearest to Y(k,:)
% 	  IDX has as many elements as Y has rows
%     [IDX,DIST] = knnsearch(X,Y);
%     ... = knnsearch(...,'k',k);
% 		search for k nearest neighbors (default: k=2)
%     ... = knnsearch(...,'Scale',Scale);
% 		Scaling vector of 'seuclidian' metric
%		default value is std(X)
%     ... = knnsearch(...,'Cov',Cov);
% 	    Cov is the covariance matrix used for Mahalanobis distance
%	    default value is cov(X)
%     ... = knnsearch(...,'Distance',Distance);
% 	the following distance metrics are currently supported:
% 	   'euclidean' [1], 
% 	   'seuclidean', (scaled euclidian)
%	   'minkowski' [3], 
%          'cityblock' or 'manhattan' [4], 
%          'hamming' [5],
%          'mahalanobis' [6],
%          'cosine' [7]
%               (one minus the cosine of the angle between the two samples),
%          'correlation'
%               (one minus the linear correlation between each pair f data vectors),
%          'spearman'
%               (one minus the rank correlation between each pair of data vectors),
%
% SEE ALSO: corrcoef, spearman, rankcorr, cov, std
%
% Reference(s):
%  [1] https://en.wikipedia.org/wiki/K-nearest_neighbors_algorithm
%  [2] https://en.wikipedia.org/wiki/Euclidean_distance
%  [3] https://en.wikipedia.org/wiki/Minkowski_distance
%  [4] https://en.wikipedia.org/wiki/Taxicab_geometry
%  [5] https://en.wikipedia.org/wiki/Hamming_distance
%  [6] https://en.wikipedia.org/wiki/Mahalanobis_distance
%  [7] https://en.wikipedia.org/wiki/Cosine_similarity

%    Copyright (C) 2021 by Alois Schloegl <alois.schloegl@gmail.com>
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


if nargin<2
	error('missing input arguments')
end

if size(X,2)~=size(Y,2)
	error('number of rows in X and Y must match')
end

% default values 
K=1; # number of NN
P=2; # exponent for minkowski distance
Distance='euclidean';
NSMethod='exhaustive';
Scale = []; 
Cov = [];

k=1;
while (k<length(varargin))
	if strcmpi(varargin{k},'k')
		K = varargin{k+1};
		k=k+2;
		continue;
	elseif strcmp(varargin{k},'P')
		P = varargin{k+1};
		k=k+2;
		continue;
	elseif strcmp(varargin{k},'Scale')
		Scale = varargin{k+1};
		if (length(Scale)~=size(X,2))
			error('size of Cov does not match input data');
		end
		k=k+2;
		continue;
	elseif strcmp(varargin{k},'Cov')
		Cov = varargin{k+1};
		if ~all(size(Cov)==size(X,2))
			error('size of Cov does not match input data');
		end
		k=k+2;
		continue;
	elseif strcmpi(varargin{k},'Distance')
		Distance = varargin{k+1};
		k=k+2;
		continue;
	elseif strcmpi(varargin{k},'NSMethod')
		NSMethod = varargin{k+1};
		k=k+2;
		continue;
	else 
		disp(varargin{k});
		fprintf(1,'Warning: input argument is ignored');
	end
	k=k+1;
end

[ix,iy]=meshgrid(1:size(X,1),1:size(Y,1));
if strcmp(Distance,'euclidean')
	D = sqrt(sum((X(ix(:),:)-Y(iy(:),:)).^2,2));

elseif strcmp(Distance,'seuclidean')
	if isempty(Scale), Scale=std(X,[],1); end; 
	IS  = Scale(:).^(-2);
	dxy = X(ix(:),:) - Y(iy(:),:);
	D   = sqrt((dxy.^2)*IS);

elseif strcmp(Distance,'mahalanobis')
	if isempty(Cov), Cov=cov(X(~any(isnan(X),2),:)); end; 
	dxy = X(ix(:),:)-Y(iy(:),:);
	D   = sqrt( sum( (dxy*inv(Cov)).*dxy, 2) );

elseif strcmp(Distance,'minkowski')
	D = sum(abs(X(ix(:),:)-Y(iy(:),:)).^P,2).^(1/P);

elseif strcmp(Distance,'cityblock') || strcmp(Distance,'manhattan')
	D = sum(abs(X(ix(:),:)-Y(iy(:),:)),2);

elseif strcmp(Distance,'cosine')
	sx = sum(X.^2, 2).^(-1/2);
	sy = sum(Y.^2, 2).^(-1/2);
	D  = 1 - sum(X(ix(:),:).*Y(iy(:),:), 2).*sx(ix(:)).*sy(iy(:));

elseif strcmp(Distance,'correlation')
	D = 1 - corrcoef(Y', X');

elseif strcmp(Distance,'spearman')
	D = 1 - corrcoef(Y', X', 'Rank');

elseif strcmp(Distance,'hamming')
	D = mean(abs(X(ix(:),:)~=Y(iy(:),:)),2);

elseif 0, % add more
else
	error(sprintf('distance metric "%s" not supported yet',Distance));
end

D = reshape(D,size(Y,1),size(X,1));
if K==1,
	[dist,idx]=min(D,[],2);
else 
	[dist,idx]=sort(D,2);
	dist=dist(:,1:K);
	idx=idx(:,1:K);
end
