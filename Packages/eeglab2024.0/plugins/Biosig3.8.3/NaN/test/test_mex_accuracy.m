% TEST_MEX_ACCURACY evaluates the accuracy and speed of 
%   different accuracy levels in SUMSKIPNAN_MEX and COVM_MEX
%
% see also: FLAG_ACCURACY_LEVEL, SUMSKIPNAN_MEX, COVM_MEX
%
% Reference:
% [1] David Goldberg, 
%       What Every Computer Scientist Should Know About Floating-Point Arithmetic
%       ACM Computing Surveys, Vol 23, No 1, March 1991. 

% 	Copyright (C) 2009,2010,2019 by Alois Schloegl <alois.schloegl@gmail.com>
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
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA

clear
flag=0;
N = 1e7;
x=randn(N,10)+1e6;

level = flag_accuracy_level;         %% backup original level 
flag_accuracy_level(0);
tic,t=cputime();[cc0,nn0]=covm_mex(x,[],flag);t0=[cputime-t,toc];

flag_accuracy_level(1);
tic,t=cputime();[cc1,nn1]=covm_mex(x,[],flag);t1=[cputime-t,toc];

flag_accuracy_level(2);
tic,t=cputime();[cc2,nn2]=covm_mex(x,[],flag);t2=[cputime-t,toc];

flag_accuracy_level(3);
tic,t=cputime();[cc3,nn3]=covm_mex(x,[],flag);t3=[cputime-t,toc];

tic,t=cputime();cc4=x'*x;nn4=size(x,1);t4=[cputime-t,toc];

flag_accuracy_level(0);
tic,t=cputime();[c0,n0]=sumskipnan_mex(x,1,flag);t0s=[cputime-t,toc];

flag_accuracy_level(1);
tic,t=cputime();[c1,n1]=sumskipnan_mex(x,1,flag);t1s=[cputime-t,toc];

flag_accuracy_level(2);
tic,t=cputime();[c2,n2]=sumskipnan_mex(x,1,flag);t2s=[cputime-t,toc];

flag_accuracy_level(3);
tic,t=cputime();[c3,n3]=sumskipnan_mex(x,1,flag);t3s=[cputime-t,toc];

tic,t=cputime();c4=sum(x,1);n4=size(x,1);t4s=[cputime-t,toc];

flag_accuracy_level(level);        %% restore original level 

cc = {cc0,cc1,cc2,cc3};
c  = {c0,c1,c2,c3};
tt = [t0;t1;t2;t3;t4];
t  = [t0s;t1s;t2s;t3s;t4s];

[status, host] = system ('hostname');
fprintf('hostname: %s\n',host);
fprintf('Sum squared differences between accuracy levels:\n');
fprintf('Level:\t|(0) naive-dou\t|(1) naive-ext\t|(2) kahan-dou \t| (3) kahan-ext\n')
fprintf('error:\t|N*2^-52\t|N*2^-64\t| 2^-52  \t| 2^-64\n')
fprintf('COVM_MEX:\ntime:\t|%f\t|%f\t| %f  \t| %f \t| %f',tt(:,2))
fprintf('\nntime:\t|%f\t|%f\t| %f  \t| %f ',tt(1:4,2)./tt(5,2))
for K1=1:4,
fprintf('\n(%i)\t',K1-1);
for K2=1:4,
	EE(K1,K2)=sum(sum((cc{K1}-cc{K2}).^2));
	E(K1,K2) =sum(sum((c{K1}-c{K2}).^2));
        fprintf('|%8g\t',EE(K1,K2)/nn1(1));
end;
end;
fprintf('\nSUMSKIPNAN_MEX:\n')
fprintf('time:\t|%f\t|%f\t| %f  \t| %f \t| %f',t(:,2))
fprintf('\nntime:\t|%f\t|%f\t| %f  \t| %f ',t(1:4,2)./t(5,2))
for K1=1:4,
fprintf('\n(%i)\t',K1-1);
for K2=1:4,
        fprintf('|%8g\t',E(K1,K2)/n1(1));
end;
end;
fprintf('\n');

