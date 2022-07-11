function [varargout]=statistic(i,DIM,fun)
% STATISTIC estimates various statistics at once.
% 
% R = STATISTIC(x,DIM)
%   calculates all statistic (see list of fun) in dimension DIM
%   R is a struct with all statistics 
%
% y = STATISTIC(x,fun)
%   estimate of fun on dimension DIM
%   y gives the statistic of fun	
%
% DIM	dimension
%	1: STATS of columns
%	2: STATS of rows
% 	N: STATS of  N-th dimension 
%	default or []: first DIMENSION, with more than 1 element
%
% fun	'mean'	mean
%	'std'	standard deviation
%	'var'	variance
%	'sem'	standard error of the mean
%	'rms'	root mean square
%	'meansq' mean of squares
%	'sum'	sum
%	'sumsq'	sum of squares
%	'CM#'	central moment of order #
%	'skewness' skewness 
%	'kurtosis' excess coefficient (Fisher kurtosis)
%	'mad'	mean absolute deviation
%
% features:
% - can deal with NaN's (missing values)
% - dimension argument 
% - compatible to Matlab and Octave
%
% see also: SUMSKIPNAN
%
% REFERENCE(S):
% [1] http://www.itl.nist.gov/
% [2] http://mathworld.wolfram.com/

%	$Id$
%	Copyright (C) 2000-2003,2010 by Alois Schloegl <alois.schloegl@gmail.com>	
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



if nargin==1,
        DIM=[];
	fun=[];        
elseif nargin==2,
        if ~isnumeric(DIM),
	        fun=DIM;
                DIM=[];
        else
                fun=[];
        end
end
if isempty(DIM), 
        DIM = find(size(i)>1,1);
        if isempty(DIM), DIM=1; end;
end;

%R.N   	= sumskipnan(~isnan(i),DIM); 	% number of elements
[R.SUM,R.N,R.SSQ] = sumskipnan(i,DIM);	% sum
%R.S3P  = sumskipnan(i.^3,DIM);		% sum of 3rd power
R.S4P  = sumskipnan(i.^4,DIM);		% sum of 4th power
%R.S5P  = sumskipnan(i.^5,DIM);		% sum of 5th power

R.MEAN 	= R.SUM./R.N;			% mean 
R.MSQ  	= R.SSQ./R.N;   		% mean square
R.RMS  	= sqrt(R.MSQ);			% root mean square
%R.SSQ0	= R.SSQ-R.SUM.*R.MEAN;		% sum square of mean removed
R.SSQ0	= R.SSQ - real(R.SUM).*real(R.MEAN) - imag(R.SUM).*imag(R.MEAN);	% sum square of mean removed

%if flag_implicit_unbiased_estim;    %% ------- unbiased estimates ----------- 
    n1 	= max(R.N-1,0);			% in case of n=0 and n=1, the (biased) variance, STD and SEM are INF
%else
%    n1	= R.N;
%end;

R.VAR  	= R.SSQ0./n1;	     		% variance (unbiased) 
R.STD  	= sqrt(R.VAR);		     	% standard deviation
R.SEM  	= sqrt(R.SSQ0./(R.N.*n1)); 	% standard error of the mean
R.SEV	= sqrt(n1.*(n1.*R.S4P./R.N+(R.N.^2-2*R.N+3).*(R.SSQ./R.N).^2)./(R.N.^3)); % standard error of the variance
R.COEFFICIENT_OF_VARIATION = R.STD./R.MEAN;

q = quantile(i, (1:3)/4, DIM);

%sz=size(i);sz(DIM)=1;
%Q0500=repmat(nan,sz);
%Q0250=Q0500;
%Q0750=Q0500;
%MODE=Q0500;
%for k=1:size(i,2),
%        tmp = sort(i(:,k));
        %ix = find(~~diff([-inf;tmp;inf]))
        %ix2=diff(ix)
        %MODE(k)= tmp(max(ix2)==ix2)
%        Q0500(k) = flix(tmp,R.N(k)/2   + 0.5);
%        Q0250(k) = flix(tmp,R.N(k)/4   + 0.5);
%        Q0750(k) = flix(tmp,R.N(k)*3/4 + 0.5);
%end;
%R.MEDIAN	= Q0500;
%R.Quartiles   	= [Q0250; Q0750];

%R.Skewness.Fisher = (R.CM3)./(R.STD.^3);	%%% same as R.SKEWNESS

%R.Skewness.Pearson_Mode   = (R.MEAN-R.MODE)./R.STD;
%R.Skewness.Pearson_coeff1 = (3*R.MEAN-R.MODE)./R.STD;
%R.Skewness.Pearson_coeff2 = (3*R.MEAN-R.MEDIAN)./R.STD;
%R.Skewness.Bowley = (Q0750+Q0250 - 2*Q0500)./(Q0750-Q0250); % quartile skewness coefficient

R.CM2	= R.SSQ0./n1;
szi = size(i); szm = [size(R.MEAN),1];
i       = i - repmat(R.MEAN,szi./szm(1:length(szi)));
R.CM3 	= sumskipnan(i.^3,DIM)./n1;
R.CM4 	= sumskipnan(i.^4,DIM)./n1;
%R.CM5 	= sumskipnan(i.^5,DIM)./n1;

R.SKEWNESS = R.CM3./(R.STD.^3);
R.KURTOSIS = R.CM4./(R.VAR.^2)-3;
[R.MAD,N]  = sumskipnan(abs(i),DIM);	% mean absolute deviation
R.MAD = R.MAD./n1;

R.datatype = 'STAT Level 3';

tmp = version;
if 0, %str2num(tmp(1))*1000+str2num(tmp(3))*100+str2num(tmp(5:6))<2136,
	% ###obsolete: was needed for Octave version < 2.1.36
        if strcmp(fun(1:2),'CM') 
                oo = str2double(fun(3:length(fun)));
                varargout  = sumskipnan(i.^oo,DIM)./n1;
        elseif isempty(fun)
	        varargout  = R;
        else	            
                varargout  = getfield(R,upper(fun));
        end;
else
	if iscell(fun),  
                for k=1:length(fun),
	                if strcmp(fun{k}(1:2),'CM') 
            	                oo = str2double(fun{k}(3:length(fun{k})));
                    	        varargout{k}  = sumskipnan(i.^oo,DIM)./n1;
                    	else	            
                    		varargout{k}  = getfield(R,upper(fun{k}));
	                end;
                end;
	elseif ischar(fun),
            	if strcmp(fun(1:2),'CM') 
                    	oo = str2double(fun(3:length(fun)));
                	varargout{1}  = sumskipnan(i.^oo,DIM)./n1;
        	else	            
    		        varargout{1}  = getfield(R,upper(fun));
            	end;
	else
		varargout{1} = R;
	end;
end;
