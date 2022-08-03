% Test train_sc and test_sc, weighted samples  



%	$Id$
%	Copyright (C) 2010 by Alois Schloegl <alois.schloegl@gmail.com>	
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


clear 
classifier= {'REG','REG2','MDA','MD2','QDA','QDA2','LD2','LD3','LD4','LD5','LD6','NBC','aNBC','WienerHopf','PLA', 'LMS','LDA/DELETION','MDA/DELETION','NBC/DELETION','RDA/DELETION','RDA','GDBC','SVM','RBF'};% 'LDA/GSVD','MDA/GSVD', 'LDA/GSVD','MDA/GSVD', 'LDA/sparse','MDA/sparse', 

N=1e2;
c=[1:N]'*2>N;

W3 = [ones(1,N/2)/5,ones(1,N/10)];
for l=1:length(classifier),
	fprintf(1,'%s\n',classifier{l});
for k=1:10,

x=randn(N,2);
x=x+[c,c];

ix = 1:0.6*N;

try,
CC = train_sc(x(ix,:),c(ix)+1,classifier{l});
R1 = test_sc(CC,x,[],c+1);

CC = train_sc(x,c+1,classifier{l});
R2 = test_sc(CC,x,[],c+1);

CC = train_sc(x(ix,:),c(ix)+1,classifier{l},W3);
R3 = test_sc(CC,x,[],c+1);

acc1(k,l)=[R1.ACC];
kap1(k,l)=[R1.kappa];
acc2(k,l)=[R2.ACC];
kap2(k,l)=[R2.kappa];
acc3(k,l)=[R3.ACC];
kap3(k,l)=[R3.kappa];
end; 

end;
end; 
 
[se,m]=sem(acc1);m
[se,m]=sem(acc2);m
[se,m]=sem(acc3);m

%[diff(m),diff(m)/sqrt(sum(se.^2))]
%[se,m]=sem(kap);[diff(m),diff(m)/sqrt(sum(se.^2))]

%These are tests to compare varios classiers

return 


N=1e2;
c=[1:N]'*2>N;

for k=1:1000,k

x=randn(N,2);
x=x+[c,c];

ix = 1:0.6*N;
[R1,CC]=xval(x(ix,:),c(ix)+1,'REG');
[R2,CC]=xval(x,c+1,'REG');
[R3,CC]=xval(x(ix,:),c(ix)+1,'LDA');
[R4,CC]=xval(x,c+1,'LDA');

acc(k,1:4)=[R1.ACC,R2.ACC,R3.ACC,R4.ACC];
kap(k,1:4)=[R1.kappa,R2.kappa,R3.kappa,R4.kappa];

end;
 
[se,m]=sem(acc),%[diff(m),diff(m)/sqrt(sum(se.^2))]
%[se,m]=sem(kap);[diff(m),diff(m)/sqrt(sum(se.^2))]

