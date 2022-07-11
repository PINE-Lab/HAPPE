% TEST_FSS test of fss.m 

% Copyright (C) 2009,2010,2016 by Alois Schloegl <alois.schloegl@gmail.com>	
%    This function is part of the NaN-toolbox
%    http://pub.ist.ac.at/~schloegl/matlab/NaN/

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
if ~exist('ue6.mat','file')	
	if strncmp(computer,'PCWIN',5)
		fprintf(1,'Download http://pub.ist.ac.at/~schloegl/LV/SMBS/UE6/ue6.mat and save in local directory %s\nPress any key to continue ...\n',pwd);
		pause;
	else 	
		system('wget http://pub.ist.ac.at/~schloegl/LV/SMBS/UE6/ue6.mat');
	end; 	
end
load ue6; 

N = 10;   % select N highest ranked features
[ix,score] = fss(data, C, N);
end; 

classifier= {'REG','MDA','MD2','QDA','QDA2','LD2','LD3','LD4','LD5','LD6','NBC','aNBC','WienerHopf','LDA/GSVD','MDA/GSVD', 'LDA/sparse','MDA/sparse', 'PLA', 'LMS','LDA/DELETION','MDA/DELETION','NBC/DELETION','RDA/DELETION','REG/DELETION','RDA','GDBC','SVM','PSVM','SVM11','SVM:LIN4','SVM:LIN0','SVM:LIN1','SVM:LIN2','SVM:LIN3'};%,'RBF'

%% compute cross-validated result; 
for k=1:N
        [R1{k},CC1{k}]=xval(data(:,ix(1:k)),C);
end; 
for k=1:length(classifier),
	fprintf(1,'%i:\t%s\n',k,classifier{k});
        [R2{k},CC2{k}]=xval(data(:,ix(1:5)),C,classifier{k});
end; 

fprintf(1,'#\tFeature\tN\tACC [%%]\tKappa+-se\t I [bit]\n'); 
R=R1; 
for k=1:length(R);
        n(k)=sum(R{k}.data(:));
        ACC(k)=R{k}.ACC;
        KAP(k)=R{k}.kappa;
        KAP_Se(k)=R{k}.kappa_se;
        MI(k)=R{k}.MI;
        
        fprintf(1,'%3i:\t%4i\t%i\t%5.2f\t%5.2f+-%5.2f\t%4.2f\n',k,ix(k),n(k),ACC(k),KAP(k),KAP_Se(k),MI(k));
end
R=R2; 
for k=1:length(R);
        n(k)=sum(R{k}.data(:));
        ACC(k)=R{k}.ACC;
        KAP(k)=R{k}.kappa;
        KAP_Se(k)=R{k}.kappa_se;
        MI(k)=R{k}.MI;
        
        fprintf(1,'%3i:\t%8s\t%i\t%5.2f\t%5.2f+-%5.2f\t%4.2f\n',k,classifier{k},n(k),ACC(k),KAP(k),KAP_Se(k),MI(k));
end


%% display 
plot(ACC*100,'x'); 
set(gca,'YLim',[0,100])
ylabel('Accuracy [%]')
title('selection of N out of 2540 features')


