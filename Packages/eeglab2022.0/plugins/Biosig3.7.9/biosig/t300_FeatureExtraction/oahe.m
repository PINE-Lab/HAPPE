function [Y0,EVENT] = OAHE(S,Fs,varargin)
% OAHE detectes obstructive Apnea/Hypopnea event
%
%    [Y,EVENT] = OAHE(X,Fs)
%    [Y,EVENT] = OAHE(filename,CHAN)
%   ... = OAHE(... ,'-o',outputFilename)
%   ... = OAHE(... ,'-e',eventFilename)
%
% INPUT:
%       X   respiratory channel
%       Fs  sampleing rate
%       filename        source filename 
%       CHAN            respiratory channels for calculating OAHE
%	outputFilename
%		name of file for storing the resulting data with the
%		detected spikes and bursts in GDF format.
%	eventFilename
%		filename to store event inforamation in GDF format. this is similar to 
%		the outputFile, except that the signal data is not included and is, therefore,
%		much smaller than the outputFile
%
% OUTPUT: 
%       Y       detection trace
%       EVENT   event structure as used in BIOSIG  
%
% see also: SVIEWER, SLOAD 
%
%
% REFERENCES:
% [1] AASM Task Force. Sleep-Related Breathing Disorders in Adults: 
%       Recommendations for Syndrome Definition, and Measurement Techniques in Clinical Research. Sleep, 22(5), 1999. 
% [2] Meoli AL, Casey KR, Clark RW, Coleman JA Jr, Fayle RW, Troell RJ, Iber C; Clinical Practice Review Committee. 
%       Hypopnea in Sleep-Disordered Breathing in Adults. Sleep. 2001 Jun 15;24(4):469-70.
% [3] Penzel T.,  Brandenburg U., Fischer J., Jobert M., Kurella B., Mayer G., Nioewerth H.J., Peter J.H., Pollmächer T., Schäfer T., Steinberg R., Trowitzsch E., Warmuth R., Weeß H.-G., Wölk C., Zulley J., 
%       Empfehlungen zur computerunterstützen Aufzeichnung und Auswertung von Polygraphien. Somnologie, 2, 42-48, 1998. 
% [4] Ross SD, Sheinhait IA, Harrison KJ, Kvasz M, Connelly JE, Shea SA, Allen IE. 
%       Systematic review and meta-analysis of the literature regarding the diagnosis of sleep apnea. Sleep. 2000 Jun 15;23(4):519-32.
% [5] Schafer T. Methodik der Atmungsmessung im Schlaf: Kapneographie zur Beurteilung der Ventilation. 
%       Biomed Tech (Berl). 2003 Jun; 48(6):170-5.
% [6] Thurnheer R, Xie X, Bloch KE. Accuracy of nasal cannula pressure recordings for assessment of ventilation during sleep. 
%       Am J Respir Crit Care Med. 2001 Nov 15;164(10 Pt 1):1914-9.
% [7] Whitney CW, Gottlieb DJ, Redline S, Norman RG, Dodge RR, Shahar E, Surovec S, Nieto FJ. 
%       Reliability of scoring respiratory disturbance indices and sleep staging. Sleep. 1998 Nov 1;21(7):749-57.
% [8] Schlögl A, Kemp B, Penzel T, Kunz D, Himanen SL, Varri A, Dorffner G, Pfurtscheller 
%       G. Quality control of polysomnographic sleep data by histogram and entropy analysis.
%       Clin Neurophysiol. 1999 Dec;110(12):2165-70.

%	$Revision: 1.3 $
%	$Id$
%	Copyright (C) 2003,2005 by Alois Schloegl <alois.schloegl@gmail.com>

% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% Version 3 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Library General Public License for more details.
%
% You should have received a copy of the GNU Library General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc., 59 Temple Place - Suite 330,
% Boston, MA  02111-1307, USA.

%%%%% default settings %%%%%
outFile = [];
evtFile = [];

%%%%% analyze input arguments %%%%%
k = 1;
while k <= length(varargin)
	if ischar(varargin{k})
		if (strcmp(varargin{k},'-o'))
			k = k + 1;
			outFile = varargin{k};
		elseif (strcmp(varargin{k},'-e'))
			k = k + 1;
			evtFile = varargin{k};
		else
			warning(sprintf('unknown input argument <%s>- ignored',varargin{k}));
		end;
	end;
	k = k+1;
end;

if ischar(S)
        [S,HDR]=sload(S,Fs);
        Fs = HDR.SampleRate;
else
        HDR.InChanSelect = 1:size(X,1);
end;

[nr,nc]=size(S);
 
nseg  = ceil (nr/(120*Fs));
nseg2 = floor(nr/(  5*Fs));

S = [repmat(nan,125*Fs,nc); S; repmat(nan, nseg*120*Fs-nr, nc)];
for k = 1:nc,
        X = reshape(S(:,k),5*Fs,24*nseg+25);
        
        hi_baseline = max(X);
        lo_baseline = min(X);
        
        ix = hankel(1:nseg2,nseg2:nseg2+25)';
        
        blh = median(hi_baseline(ix(1:24,:)));
        bll = median(lo_baseline(ix(1:24,:)));
        
        mh  = max(hi_baseline(ix(25:26,:)));
        ml  = min(lo_baseline(ix(25:26,:)));
        
        Y   = real((blh-bll)' > 2*(mh-ml)');
        Y(any(isnan([mh;ml;blh;bll]))) = NaN;
        
        Y0(:,k) = Y; 
end;

%%%%% generate event encoding 
Y = Y0; 
ix = isnan(Y(1,:));
Y(1,ix) = 0; 
for k=2:length(Y), 
	ix = isnan(Y(k,:));
	Y(k,ix) = Y(k-1,ix); 
end;
Y(end+1,:)=0; 

EVENT.SampleRate = Fs; 
N = 0;
for k=1:nc;
	tmp1 = diff(Y(:,k));
        tmp2 = find(tmp1>0); 
	L    = length(tmp2);
	EVENT.TYP(N+1:N+L,1) = hex2dec('0401');
	EVENT.POS(N+1:N+L,1) = tmp2*Fs;
	EVENT.DUR(N+1:N+L,1) = (find(tmp1<0)-tmp2)*Fs;
	EVENT.CHN(N+1:N+L,1) = HDR.InChanSelect(k);
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(outFile)
		%%% write data to output
		HDR.TYPE  = 'GDF';
		HDR.VERSION = 3;
		%[p,f,e]=fileparts(fn);
		HDR.FILE = [];
		HDR.FileName  = outFile;
		HDR.FILE.Path = '';
		HDR.PhysMax   = max(s);
		HDR.PhysMin   = min(s);
		HDR.DigMax(:) =  2^15-1;
		HDR.DigMin(:) = -2^15;
		HDR.GDFTYP(:) = 3;
		HDR.FLAG.UCAL = 0;
		HDR.NRec = size(s,1);
		HDR.SPR = 1;
		HDR.NS  = size(s,2);
		HDR.Dur = 1/HDR.SampleRate;
		HDR = rmfield(HDR,'AS');
		HDR.EVENT = EVENT; 
		HDR = sopen(HDR,'w');
		if (HDR.FILE.FID < 0) 
			fprintf(2,'Warning can not open file <%s> - GDF file can not be written\n',HDR.FileName);
		else
			HDR = swrite(HDR,s);
			HDR = sclose(HDR);
		end; 
end;

if ~isempty(evtFile)
		%%% write data to output
		HDR.TYPE  = 'EVENT';
		HDR.VERSION = 3;
		%[p,f,e]=fileparts(fn);
		HDR.FILE = [];
		HDR.FileName  = evtFile;
		HDR.FILE.Path = '';
		HDR.NRec = 0;
		HDR.SPR = 0;
		HDR.Dur = 1/HDR.SampleRate;
		HDR = rmfield(HDR,'AS');
		HDR = sopen(HDR, 'w');
		if (HDR.FILE.FID<0) 
			fprintf(2,'Warning can not open file <%s> - EVT file can not be written\n', HDR.FileName); 
		else
			HDR = sclose(HDR);
 		end;
end;


