function [rix,cix] = row_col_deletion(d,c,w)
% ROW_COL_DELETION selects the rows and columns for removing any missing values. 
%    A heuristic based on maximizing the number of remaining sample values
%    is used. In other words, if there are more rows than columns, it is 
%    more likely that a row-wise deletion will be applied and vice versa. 
% 
%    [rix,cix] = row_col_deletion(d)
%    [rix,cix] = row_col_deletion(d,c,w)
% 
% Input: 
%    d        data (each row is a sample, each column a feature)
%    c        classlabels (not really used) [OPTIONAL]
%    w        weight for each sample vector [OPTIONAL]   
% Output:
%    rix      selected samples
%    cix      selected columns    
% 
%   d(rix,cix) does not contain any NaN's i.e. missing values      
% 
% see also: TRAIN_SC, TEST_SC
 
%	$Id$
%	Copyright (C) 2009,2010 by Alois Schloegl <alois.schloegl@gmail.com>
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


if nargin > 2,
	if isempty(w) || all(w==w(1)), 
		ix = ~isnan(c);
	else
		ix = ~any(isnan(c) | isnan(w));
	end; 
	d = d(ix,:);	%% ignore samples with invalid c or w
	w = w(ix,:);

elseif nargin > 1,
	d = d(~isnan(c),:);	%% ignore samples with invalid c or w
        w = [];
else
        w = [];        
end; 	


if 0, 
	% decides whether row-wise or column-wise deletion removes less data. 
	% rix and cix are the resulting index vectors
	% either row-wise or column-wise deletion, but not a combination of both, is used. 
	% this is obsolete
	
        	n   = numel(d); 
        	cix = find(~any(isnan(d),1)); 
        	rix = find(~any(isnan(d),2)); 
        	nr  = length(rix)*size(d,2); % number of elements after row-wise deletion
        	nc  = length(cix)*size(d,1); % number of elements after column-wise deletion
        	
        	if (nr>nc)
        		cix = 1:size(d,2);  % select all columns
        		%fprintf(1,'row-wise deletion (%i,%i,%i)\n',n,nr,nc);		
        	else
        		rix = 1:size(d,1);  % select all rows 
        		%fprintf(1,'column-wise deletion (%i,%i,%i)\n',n,nr,nc);		
        	end; 

else

	        %% a mix of row- and column-wise deletion is possible 
		if ~isempty(w) && (abs(sum(w)-1) > log2(N)*eps || any(w<0) || any(~isfinite(w)))
		        error('weight vector must contain only non-negative and finite values');
		end; 
		[N,M] = size(d); 
		rix = ones(N,1); cix = ones(1,M);
		while 1;
                        e = ~isnan(d(rix>0,cix>0));
                        if  ~isempty(w),
                                colCost  = mean(e, 1, w(rix>0)/sum(w(rix>0)))';        % cost of deleting columns
                        else
                                colCost  = mean(e, 1)';        % cost of deleting columns
                        end;         
                        rowCost  = mean(e, 2);           % cost of deleting rows     
                        [tmp,ix] = sort([colCost; rowCost]);
                        
                if abs(tmp(1)-1) < log2(N)*eps, break; end;        % stopping criterion

                        if diff(tmp(1:2))==0, warning('row/col deletion: arbitrary selection [%i,%i]',ix(1:2)); end;  
                        ix = ix(1);
                        if (ix<=sum(cix))
                                tmp = find(cix>0);
                                cix(tmp(ix)) = 0;
                        else
                                tmp = find(rix>0);
                                rix(tmp(ix-sum(cix))) = 0;
                        end;
                end; 
                rix = find(rix);
                cix = find(cix);

end

