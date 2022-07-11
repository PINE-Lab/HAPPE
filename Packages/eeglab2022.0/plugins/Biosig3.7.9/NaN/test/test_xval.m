% test_classifier;

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


if 1,
clear
N=100;  % number of samples 
M=10;   % number of features 
classifier = {'SVM:LIB', 'REG', 'MDA', 'MD2', 'QDA', 'QDA2', 'LD2', 'LD3', 'LD4','LD5', 'LD6', 'NBC', 'aNBC', 'WienerHopf', 'LDA/GSVD', 'MDA/GSVD', 'LDA/sparse', 'MDA/sparse', 'PLA', 'LMS', 'LDA/DELETION', 'MDA/DELETION', 'NBC/DELETION', 'RDA/DELETION', 'REG/DELETION', 'RDA', 'GDBC', 'SVM', 'RBF', 'PSVM', 'SVM11', 'SVM:LIN4', 'SVM:LIN0', 'SVM:LIN1', 'SVM:LIN2', 'SVM:LIN3', 'SVM:LIN5', 'SVM:LIN6', 'SVM:LIN7', 'SVM:LIN11', 'SVM:LIN12', 'SVM:LIN13', 'WINNOW' };

x = randn(N,M);         % data
c = ([1:N]'>(N/2))+1;   % classlabel 
%w = [ones(1,N/2)/5,ones(1,N/10),zeros(1,4*N/10)];
w = [];                 % no weightening

x = randn(N,M);
x = x+c*ones(1,M);

if 1,
%x(2:2:N/2,2) = NaN; 
x(2:2:N,2) = NaN; 
x(3,2:2:end) = NaN;
end; 
end; 

for k = 1:length(classifier);
	try,
		[R{k},CC{k}] = xval(x, {c,w}, classifier{k}); 
		fprintf(1,'%8s\t%i\t%5.2f\t%5.2f+-%5.2f\n',classifier{k},sum(R{k}.data(:)),R{k}.ACC*100,R{k}.kappa,R{k}.kappa_se);
		save -v6 debug.mat
	catch,
		R{k} = [];
	end; 
end;

for k = 1:length(R)
	if isempty(R{k})
		fprintf(1,'%8s \t failed\n',classifier{k});
	else
		fprintf(1,'%8s\t%i\t%5.2f\t%5.2f+-%5.2f\n',classifier{k},sum(R{k}.data(:)),R{k}.ACC*100,R{k}.kappa,R{k}.kappa_se);
	end; 	
end

