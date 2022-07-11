% TEST_CLASSIFY tests and compares NaN/CLASSIFY.M with the matlab version of CLASSIFY

%	Copyright (C) 2009,2010,2019 by Alois Schloegl <alois.schloegl@gmail.com>	
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
load_fisheriris
SL = meas(51:end,1);
SW = meas(51:end,2);
group = species(51:end);
h1 = gscatter(SL,SW,group,'rb','v^',[],'off');
set(h1,'LineWidth',2)
legend('Fisher versicolor','Fisher virginica','Location','northwest')
       
[X,Y] = meshgrid(linspace(4.5,8),linspace(2,4));
X = X(:); Y = Y(:);

classifiers={'linear','quadratic','diagLinear','diagQuadratic','mahalanobis'};

p = which('train_sc.m'); 
p = fileparts(p); 
rmpath(p); 	
for k=1:length(classifiers)
  [C1,err(1,k),P1,logp1,coeff1] = classify([X Y],[SL SW],group,classifiers{k});
end; 

addpath(p);
for k=1:length(classifiers)
  [C2,err(2,k),P2,logp2,coeff2] = classify([X Y],[SL SW],group,classifiers{k});
end; 


