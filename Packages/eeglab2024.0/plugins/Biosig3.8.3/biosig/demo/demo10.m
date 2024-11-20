%% Implement deconvolution method as used in [1]
%%
%% References: 
%%  [1] Alejandro Javier Pernía-Andrade, Sarit Pati Goswami, Yvonne Stickler, 
%%      Ulrich Fröbe, Alois Schlögl, and Peter Jonas (submitted) 
%%  A deconvolution-based method with high sensitivity and temporal
%%  resolution for detection of spontaneous synaptic currents in vitro and
%%  in vivo. Biophysical Journal Volume 103 October 2012 1–11.

%  $Id$
%  Copyright (C) 2012 by Alois Schloegl  <alois.schloegl@ist.ac.at>  
%  This is part of the BIOSIG-toolbox http://biosig.sf.net/

% BioSig is free software; you can redistribute it and/or
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
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
% USA.


f1 = 'http://pub.ist.ac.at/~schloegl/download/rawdata.gdf';
f2 = 'http://pub.ist.ac.at/~schloegl/download/template.gdf';

while ~exist('rawdata.gdf','file') || ~exist('template.gdf','file')
    fprintf(1,'Download example and template file from here:\n');
    fprintf(1,'\t%s\n\t%s\n',f1,f2);
    fprintf(1,'and save them in local working directory "%s".\n',pwd);
    pause;
end; 


clear; 
B = [0.1,150]; %% bandwidth 

% load raw data
[r,rHDR]=sload('rawdata.gdf');

% load template
[t,tHDR]=sload('template.gdf');

fs = rHDR.SampleRate;

% deconvolution and filteriong
d = signal_deconvolution(r, t, fs, B); 

% get local maxima above threshold
TH = quantile(d,1-1e-3);
pos = get_local_maxima_above_threshold(d,TH);
	
% display brief segment 
xlim = [15,20];
t = [1:length(r)]'/fs;
subplot(2,1,1); plot(t,r); set(gca,'xlim',xlim);
legend({'Raw data'});
subplot(2,1,2); plot(t,d,'g-',t(pos),TH,'ro'); set(gca,'xlim',xlim);
legend({'Deconvolved data','local max above threshold'});

