function Q=quantile(Y,q,DIM,method)
% QUANTILE calculates the quantiles of histograms and sample arrays.  
%
%  Q = quantile(Y,q)
%  Q = quantile(Y,q,DIM)
%     returns the q-th quantile along dimension DIM of sample array Y.
%     size(Q) is equal size(Y) except for dimension DIM which is size(Q,DIM)=length(Q)
%
%  Q = quantile(HIS,q)
%     returns the q-th quantile from the histogram HIS. 
%     HIS must be a HISTOGRAM struct as defined in HISTO2 or HISTO3.
%     If q is a vector, the each row of Q returns the q(i)-th quantile 
%
% see also: HISTO2, HISTO3, PERCENTILE


%	$Id$
%	Copyright (C) 1996-2003,2005,2006,2007,2009,2011 by Alois Schloegl <alois.schloegl@gmail.com>	
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

if nargin<3,
        DIM = [];
end;
if isempty(DIM),
        DIM = find(size(Y)>1,1);
        if isempty(DIM), DIM = 1; end;
end;


if nargin<2,
	help quantile
        
else
	[q, rix]  = sort(q(:)'); 	% sort quantile values 
	[tmp,rix] = sort(rix);	% generate reverse index 

        SW = isstruct(Y);
        if SW, SW = isfield(Y,'datatype'); end;
        if SW, SW = strcmp(Y.datatype,'HISTOGRAM'); end;
        if SW,                 
                [yr, yc] = size(Y.H);
                Q = repmat(nan,length(q),yc);
                if ~isfield(Y,'N');
                        Y.N = sum(Y.H,1);
                end;
                
                for k1 = 1:yc,
                        tmp = Y.H(:,k1)>0;
                        h = full(Y.H(tmp,k1));
                        t = Y.X(tmp,min(size(Y.X,2),k1));

			N = Y.N(k1);  
			t2(1:2:2*length(t)) = t;
			t2(2:2:2*length(t)) = t;
			x2 = cumsum(h); 
			x(1)=0; 
			x(2:2:2*length(t)) = x2;
			x(3:2:2*length(t)) = x2(1:end-1);

			% Q(q < 0 | 1 < q,:) = NaN;  % already done at initialization
			Q(q==0,k1) = t2(1);
			Q(q==1,k1) = t2(end);
			n = 1; 
		    	for k2 = find( (0 < q) & (q < 1) )
				while (q(k2)*N > x(n)),
					n=n+1; 
				end;

				if q(k2)*N==x(n)
					% mean of upper and lower bound 
					Q(k2,k1) = (t2(n)+t2(n+1))/2;
				else
					Q(k2,k1) = t2(n);
				end;
                        end;
			Q = Q(rix,:);	% order resulting quantiles according to original input q 
                end;


        elseif isnumeric(Y),
		sz = size(Y);
		if DIM>length(sz),
		        sz = [sz,ones(1,DIM-length(sz))];
		end;

		f  = zeros(1,length(q));        
		f( (q < 0) | (1 < q) ) = NaN;
		D1 = prod(sz(1:DIM-1));
		D3 = prod(sz(DIM+1:length(sz)));
		Q  = repmat(nan,[sz(1:DIM-1),length(q),sz(DIM+1:length(sz))]);
		for k = 0:D1-1,
		for l = 0:D3-1,
		        xi = k + l * D1*sz(DIM) + 1 ;
			xo = k + l * D1*length(q) + 1;
		        t  = Y(xi:D1:xi+D1*sz(DIM)-1);
		        t  = t(~isnan(t));
		        N  = length(t); 
			
			if (N==0)
		            	f(:) = NaN;        
			else
				t  = sort(t);
				t2(1:2:2*length(t)) = t; 
				t2(2:2:2*length(t)) = t;
				x = floor((1:2*length(t))/2);
				%f(q < 0 | 1 < q) = NaN;  % for efficiency its defined outside loop 
				f(q==0) = t2(1);
				f(q==1) = t2(end);

				n = 1; 
			    	for k2 = find( (0 < q) & (q < 1) )
					while (q(k2)*N > x(n)),
						n = n+1;
					end;

					if q(k2)*N==x(n)
						% mean of upper and lower bound 
						f(k2) = (t2(n) + t2(n+1))/2;   
					else
						f(k2) = t2(n);
					end; 
				end; 		
			end; 
			Q(xo:D1:xo + D1*length(q) - 1) = f(rix);
		end;
		end;

        else
                fprintf(2,'Error QUANTILES: invalid input argument\n');
                return;
        end;
        
end;

%!assert(quantile(1:10,[.2,.5]),[2.5, 5.5])
%!assert(quantile([1:3,NaN],[.10,.50,.90])==[1,2,3])


