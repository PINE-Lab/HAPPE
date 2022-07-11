function [idx,score] = fss(D,cl,N,MODE)
% FSS - feature subset selection and feature ranking 
%   the method is motivated by the max-relevance-min-redundancy (mRMR) 
%   approach [1]. However, the default method uses partial correlation,
%   which has been developed from scratch. PCCM [3] describes
%   a similar idea, but is more complicated. 
%   An alternative method based on FSDD is implemented, too. 
%    
%  [idx,score] = fss(D,cl) 
%  [idx,score] = fss(D,cl,MODE) 
%  [idx,score] = fss(D,cl,MODE) 
%    
% D 	data - each column represents a feature 
% cl	classlabel   
% Mode 	'Pearson' [default] correlation
%	'rank' correlation 
%       'FSDD' feature selection algorithm based on a distance discriminant [2]
%       %%% 'MRMR','MID','MIQ' max-relevance, min redundancy [1] - not supported yet. 
%
% score score of the feature
% idx	ranking of the feature    
%       [tmp,idx]=sort(-score)
%
% see also: TRAIN_SC, XVAL, ROW_COL_DELETION
%
% REFERENCES:
% [1] Peng, H.C., Long, F., and Ding, C., 
%   Feature selection based on mutual information: criteria of max-dependency, max-relevance, and min-redundancy, 
%   IEEE Transactions on Pattern Analysis and Machine Intelligence, 
%   Vol. 27, No. 8, pp.1226-1238, 2005.
% [2] Jianning Liang, Su Yang, Adam Winstanley, 
%   Invariant optimal feature selection: A distance discriminant and feature ranking based solution, 
%   Pattern Recognition, Volume 41, Issue 5, May 2008, Pages 1429-1439.
%   ISSN 0031-3203, DOI: 10.1016/j.patcog.2007.10.018.
% [3] K. Raghuraj Rao and S. Lakshminarayanan
%   Partial correlation based variable selection approach for multivariate data classification methods
%   Chemometrics and Intelligent Laboratory Systems
%   Volume 86, Issue 1, 15 March 2007, Pages 68-81 
%   http://dx.doi.org/10.1016/j.chemolab.2006.08.007

%	$Id$
%	Copyright (C) 2009,2010 by Alois Schloegl <alois.schloegl@gmail.com>
%       This function is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/

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
%    along with this program; If not, see <http://www.gnu.org/licenses/>.

        if nargin<3
                MODE = [];
                N = [];
        elseif ischar(N)
                MODE = N;
                N = [];
        elseif nargin<4,
                MODE = [];
        end

        if isempty(N), N = size(D,2); end
        score = repmat(NaN,1,size(D,2));

        if 0, %strcmpi(MODE,'MRMR') || strcmpi(MODE,'MID') || strcmpi(MODE,'MIQ');
                %% RMRM/MID/MIQ is not supported
                %% TODO: FIXME 
                 
                [tmp,t] = sort([cl,D]);
                cl = t(:,1:size(cl,2)); 
                D  = t(:,1:size(D,2)); 
                for k = 1:N,
                        V(k) = mi(cl, D(:,k));

                        for m = 1:N,
                                W(k,m) = mi(D(:,m), D(:,k));
                        end
                        MID(k) = V(k) - mean(W(k,:));
                        MIQ(k) = V(k) / mean(W(k,:));
                end

                if  strcmpi(MODE,'MIQ')
                        [score,idx] = sort(MIQ,[],'descend');
                else
                        [score,idx] = sort(MID,[],'descend');
                end

        elseif strcmpi(MODE,'FSDD');
                [b,i,j]=unique(cl);
                for k=1:length(b)
                        n(k,1) = sum(j==k);
                        m(k,:) = mean(D(j==k,:),1);
                        v(k,:) = var(D(j==k,:),1);
                end
                m0 = mean(m,1,n);
                v0 = var(D,[],1);
                s2 = mean(m.^2,1,n) - m0.^2;
                score = (s2 - 2*mean(v,1,n)) ./ v0;
                [t,idx] = sort(-score);

        elseif isempty(MODE) || strcmpi(MODE,'rank') || strcmpi(MODE,'Pearson')
                cl = cat2bin(cl);
                if strcmpi(MODE,'rank'),
                        [tmp,D] = sort(D,1);
                end
                idx = repmat(NaN,1,N);
                for k = 1:N,
                        f = isnan(score);

                        %%%%% compute partial correlation (X,Y|Z)
                        % r = partcorrcoef(cl, D(:,f), D(:,~f)); % obsolete, not very robust

                        %% this is a more robust version 
                        X = cl; Y = D(:,f); Z = D(:,~f);
                        if (k>1)
                                X = X-Z*(Z\X);
                                Y = Y-Z*(Z\Y);
                        end
                        r = corrcoef(X,Y);

                        [s,ix] = max(sumsq(r,1));
                        f      = find(f);
                        idx(k) = f(ix);
                        score(idx(k)) = s;
                end

        end
end

function I = mi(x,y)
        ix  = ~any(isnan([x,y]),2);
        H   = sparse(x(ix),y(ix));
        pij = H./sum(ix); 
        Iij = pij.*log2(pij./(sum(pij,2)*sum(pij,1)));
        Iij(isnan(Iij)) = 0;
        I = sum(Iij(:));
end
