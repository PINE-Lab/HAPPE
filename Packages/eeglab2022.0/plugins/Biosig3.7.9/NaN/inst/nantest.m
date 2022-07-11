% NANTEST checks several mathematical operations and a few 
% statistical functions for their correctness related to NaN's.
% e.g. it checks norminv, normcdf, normpdf, sort, matrix division and multiplication.
%
%
% see also: NANINSTTEST
%
% REFERENCE(S): 
% [1] W. Kahan (1996) Lecture notes on the Status of "IEEE Standard 754 for 
%     Binary Floating-point Arithmetic. 
%


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
%    along with this program; If not, see <http://www.gnu.org/licenses/>.

%	$Id$
%	Copyright (C) 2000-2004,2009 by Alois Schloegl <alois.schloegl@gmail.com>
%       This script is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/

%FLAG_WARNING = warning;
%warning('off');

try
	x = randn([3,4,5]); 
	x(~isnan(x)) = 0;
catch
	fprintf(1,'WARNING: NANTEST fails for 3-DIM matrices. \n');
end;
try
	[s,n] = sumskipnan([nan,1,4,5]);
catch
	fprintf(1,'WARNING: SUMSKIPNAN is not avaible. \n');
end;

% check NORMPDF, NORMCDF, NORMINV
x = [-inf,-2,-1,-.5,0,.5,1,2,3,inf,nan]';
if exist('normpdf','file')==2,
        q(1) = sum(isnan(normpdf(x,2,0)))>sum(isnan(x));
        if q(1),
                fprintf(1,'NORMPDF cannot handle v=0.\n');
                fprintf(1,'-> NORMPDF should be replaced\n');
        end;
end;

if exist('normcdf','file')==2,
        q(2) = sum(isnan(normcdf(x,2,0)))>sum(isnan(x));
        if q(2),
                fprintf(1,'NORMCDF cannot handle v=0.\n');
                fprintf(1,'-> NORMCDF should be replaced\n');
        end;
end;

if ~(any(exist('erfinv') == [2,5]))
       fprintf(1,'ERFINV is not available\n');

elseif exist('norminv','file')==2,
        p = [-inf,-.2,0,.2,.5,1,2,inf,nan];
        q(3) = sum(~isnan(norminv(p,2,0)))<4;
        if q(3),
                fprintf(1,'NORMINV cannot handle  correctly v=0.\n');
                fprintf(1,'-> NORMINV should be replaced\n');
        end;
        q(4) = ~isnan(norminv(0,NaN,0)); 
        q(5) = any(norminv(0.5,[1 2 3],0)~=(1:3));
end;

if exist('tpdf','file')==2,
        q(6) = ~isnan(tpdf(nan,4));
        if q(6),
                fprintf(1,'TPDF(NaN,4) does not return NaN\n');
                fprintf(1,'-> TPDF should be replaced\n');
        end;
end;

if exist('tcdf','file')==2,
        try
                q(7) = ~isnan(tcdf(nan,4));
        catch
                q(7) = 1;
        end;
        if q(7),
                fprintf(1,'TCDF(NaN,4) does not return NaN\n');
                fprintf(1,'-> TCDF should be replaced\n');
        end;
end;

if exist('tinv','file')==2,
        try
                q(8) = ~isnan(tinv(nan,4));
        catch
                q(8) = 1;
        end;
        if q(8),
                fprintf(1,'TINV(NaN,4) does not return NaN\n');
                fprintf(1,'-> TINV should be replaced\n');
        end;
end;

q(9) = isreal(double(2+3i));
if q(9)
	printf('DOUBLE rejects imaginary part\n-> this can affect SUMSKIPNAN\n');
end; 

try 
        x = reshape(1:6,3,2); 
        [cc,nn] = covm(x+i*x,'e');
        q(10) = 0; 
catch
        q(10) = 1; 
end; 

if 0,
%%%%% MOD 
if exist('mod')>1,
        if (mod(5,0))~=0,
                fprintf(1,'WARNING: MOD(x,0) does not return 0.\n');
        end;
        if isnan(mod(5,0)),
                fprintf(1,'WARNING: MOD(x,0) returns NaN.\n');
        end;
        if isnan(mod(5,inf)),
                fprintf(1,'WARNING: MOD(x,INF) returns NaN.\n');
        end;
end;
%%%%% REM 
if exist('rem')>1,
        if (rem(5,0))~=0,
                fprintf(1,'WARNING: REM(x,0) does not return 0.\n');
        end;
        if isnan(rem(5,0)),
                fprintf(1,'WARNING: REM(x,0) returns NaN.\n');
        end;
        if isnan(mod(5,inf)),
                fprintf(1,'WARNING: REM(x,INF) returns NaN.\n');
        end;
end;
end; 

%%%%% NANSUM(NAN) - this test addresses a problem in Matlab 5.3, 6.1 & 6.5
if exist('nansum','file'),
        if isnan(nansum(nan)),
                fprintf(1,'Warning: NANSUM(NaN) returns NaN instead of 0\n');
                fprintf(1,'-> NANSUM should be replaced\n');
        end;
end;
%%%%% NANSUM(NAN) - this test addresses a problem in Matlab 5.3, 6.1 & 6.5
if exist('nanstd','file'),
        if ~isnan(nanstd(0)),
                fprintf(1,'Warning: NANSTD(x) with isscalar(x) returns 0 instead of NaN\n');
                fprintf(1,'-> NANSTD should be replaced\n');
        end;
end;
%%%%% GEOMEAN - this test addresses a problem in Octave
if exist('geomean','file'),
        if isnan(geomean((0:3)')),
                fprintf(1,'Warning: GEOMEAN([0,1,2,3]) NaN instead of 0\n');
                fprintf(1,'-> GEOMEAN should be replaced\n');
        end;
end;
%%%%% HARMMEAN - this test addresses a problem in Octave
if exist('harmmean','file'),
        if isnan(harmmean(0:3)),
                fprintf(1,'Warning: HARMMEAN([0,1,2,3]) NaN instead of 0\n');
                fprintf(1,'-> HARMMEAN should be replaced\n');
        end;
end;
%%%%% BITAND - this test addresses a problem in Octave
if exist('bitand')>1,
        if isnan(bitand(2^33-1,13)),
                fprintf(1,'BITAND can return NaN. \n');
        end;
end;
%%%%% BITSHIFT - this test addresses a problem in Octave
if exist('bitshift','file'),
        if isnan(bitshift(5,30,32)),
                fprintf(1,'BITSHIFT can return NaN.\n');
        end;
end;
%%%%% ALL - this test addresses a problem in some old Octave and FreeMat v3.5
if any(NaN)==1,
	fprintf(1,'WARNING: ANY(NaN) returns 1 instead of 0\n');
end;
if any([])==1,
	fprintf(1,'WARNING: ANY([]) returns 1 instead of 0\n');
end;
%%%%% ALL - this test addresses a problem in some old Octave and FreeMat v3.5
if all(NaN)==0,
	fprintf(1,'WARNING: ALL(NaN) returns 0 instead of 1\n');
end;
if all([])==0,
	fprintf(1,'WARNING: ALL([]) returns 0 instead of 1\n');
end;
	
%%%%% SORT - this was once a problem in Octave Version < 2.1.36, and still is in FreeMat 4.0 %%%%
if ~all(isnan(sort([3,4,NaN,3,4,NaN]))==[0,0,0,0,1,1]), 
        warning('Warning: SORT does not properly handle NaN.');
end;

%%%%% commutativity of 0*NaN	%%% This test adresses a problem in Octave
x=[-2:2;4:8]';
y=x;y(2,1)=nan;y(4,2)=nan;
B=[1,0,2;0,3,1];
if ~all(all(isnan(y*B)==isnan(B'*y')')),
        fprintf(2,'WARNING: 0*NaN within matrix multiplication is not commutative\n');
end;

% from Kahan (1996)
tmp = (0-3*i)/inf;
if isnan(tmp)
        fprintf(2,'WARNING: (0-3*i)/inf results in NaN instead of 0.\n');
end;

%(roots([5,0,0])-[0;0])
%(roots([2,-10,12])-[3;2])
%(roots([2e-37,-2,2])-[1e37;1])
%%%%% check nan/nan   %% this test addresses a problem in Matlab 5.3, 6.1 & 6.5
p    = 4;
tmp1 = repmat(nan, 4);
tmp2 = repmat(nan, 4);
if ispc
        % Octave 4.0.0 on Windows crashes, therefore the test is disabled
	warning('mrdivide (repmat(nan,4), repmat(nan,4)) and mldivide (repmat(nan,4), repmat(nan,4)) not tested because it might crash Octave on Windows.\n');	
else        
try
	tmp1 = repmat(nan,p) / repmat(nan,p);
catch   % exception error in Octave 3.8.2 and later of debian wheezy
	fprintf(2,'mrdivide (repmat(nan,4), repmat(nan,4)) fails with an exception\n');	
end;
try
	tmp2 = repmat(nan,p) \ repmat(nan,p);
catch   % exception error in Octave 3.8.2 and later of debian wheezy
	fprintf(2,'mldivide (repmat(nan,4), repmat(nan,4)) fails with an exception\n');
end
end;
tmp3 = repmat(0,p)/repmat(0,p);
tmp4 = repmat(0,p)\repmat(0,p);
tmp5 = repmat(0,p)*repmat(inf,p);
tmp6 = repmat(inf,p)*repmat(0,p);
x    = randn(100,1)*ones(1,p); y=x'*x; 
tmp7 = y/y;
tmp8 = y\y;

if ~all(isnan(tmp1(:))),
        fprintf(1,'WARNING: matrix division NaN/NaN does not result in NaN\n');
end;
if ~all(isnan(tmp2(:))),
        fprintf(1,'WARNING: matrix division NaN\\NaN does not result in NaN\n');
end;
if ~all(isnan(tmp3(:))),
        fprintf(2,'WARNING: matrix division 0/0 does not result in NaN\n');
end;
if ~all(isnan(tmp4(:))),
        fprintf(2,'WARNING: matrix division 0\\0 does not result in NaN\n');
end;
if ~all(isnan(tmp5(:))),
        fprintf(2,'WARNING: matrix multiplication 0*inf does not result in NaN\n');
end;
if ~all(isnan(tmp6(:))),
        fprintf(2,'WARNING: matrix multiplication inf*0 does not result in NaN\n');
end;
if any(any(tmp7==inf));
        fprintf(2,'WARNING: right division of two singulare matrices return INF\n');
end;
if any(any(tmp8==inf));
        fprintf(2,'WARNING: left division of two singulare matrices return INF\n');
end;

tmp  = [tmp1;tmp2;tmp3;tmp4;tmp5;tmp6;tmp7;tmp8];



%warning(FLAG_WARNING); 


%%%%% QUANTILE TEST 
d = [1 1 2 2 4 4 10 700]'; 
q = [-1,0,.05,.1,.25,.49,.5,.51,.75,.8, .999999,1,2];
r = [ NaN, 1, 1, 1, 1.5, 2, 3, 4, 7, 10, 700, 700, NaN]; 
if any( quantile(d, q)' -  r>0)
	fprintf(1,'Quantile(1): failed\n');
else
	fprintf(1,'Quantile(1): OK\n'); 
end; 
if exist('histo3','file')
	H = histo3(d);
else
	H.X = [1;2;4;10;700];
	H.H = [2;2;2;1;1];
	H.datatype = 'HISTOGRAM'; 
end; 	 
if any( quantile(H, q)' -  r>0)
	fprintf(1,'Quantile(2): failed\n');
else
	fprintf(1,'Quantile(2): OK\n'); 
end; 

