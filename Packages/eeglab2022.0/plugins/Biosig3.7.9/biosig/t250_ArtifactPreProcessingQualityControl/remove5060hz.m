function [s,HDR] = remove5060hz(varargin)
%  REMOVE5060HZ removes line interference artifacts
%    from biomedical signal data. Segmented data supported only 
%    with an FFT-based method. 
%
%   [s,HDR] = remove5060hz(fn, Mode)
%   [s,HDR] = remove5060hz(fn, CHAN, Mode)
%   [s,HDR] = remove5060hz(s, HDR, Mode)
%             remove5060hz(..., '-o',outfile)
%  
%
% INPUT:
% 	fn	biosignal file
%	CHAN	channel selection [default: 0 (all)]
%	s	signal data
%	HDR	header structure (contains labels, sampling rate etc)
%	Mode    'PCA', 'PCA 50'	50 Hz PCA/SVD filter
%		'NOTCH'  	50 Hz FIR Notch filter, order=3
%               'FIT'           fit and remove 50 Hz sine/cosine wave
%		'FFT'		fft filter - cancels 50+-0.5 Hz in Frequency domain (default)
%		'FFT3'		fft filter - cancels also 3rd harmonic 50 and 150 Hz 
%		'PCA 60'	60 Hz PCA/SVD filter
%		'NOTCH 60' 	60 Hz FIR Notch filter, order=3
%               'FIT 60'        fit and remove 60 Hz sine/cosine wave
%		'FFT 60'	fft filter - cancels 60+-0.5 Hz in frequency domain
%		'FFT3 60'	fft filter - cancels also 3rd harmonic 60 and 180 Hz 
%       outfile	Name of file for storing corrected data. 
%
% OUTPUT:
%	s	corrected signal data
%	HDR	header structure
%
%
% see also: REGRESS_EOG, EVENTCODES
%
% Reference(s):
%


% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

%	$Id$
%	(C)2006,2009,2012 by Alois Schloegl <alois.schloegl@ist.ac.at>
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

MODE = [];

arg2 = []; arg3=[];
arg1 = varargin{1};
if length(varargin)>1, arg2 = varargin{2}; end;
if length(varargin)>2, arg3 = varargin{3}; end; 

outfile = []; 
k = strmatch('-o',varargin); 
if k, outfile = varargin{k+1}; end; 

if ischar(arg1),
        CHAN = 0;
        if (nargin>1)
                if ischar(arg2)
                        MODE = arg2;
                        CHAN = 0;
                else
                        MODE = '50';
                        CHAN = arg2;
                end;
        end;
        [s,HDR] = sload(arg1,CHAN);

elseif isnumeric(arg1) && isstruct(arg2)
        s = arg1;
        HDR = arg2;
        if (nargin>2)
                MODE = arg3;
        end;

end

if isempty(MODE)
	MODE = 'FFT';
end;

NotchFreq = 50; 	% default 50 Hz Notch
if strfind(upper(MODE),'60')
        NotchFreq = 60;
end;


if isfield(HDR,'EVENT') 
	segs = [1;HDR.EVENT.POS(HDR.EVENT.TYP==hex2dec('7ffe'));size(s,1)];
else
	segs = [1;size(s,1)];
end; 
flagNoSegmentSupport = length(segs)>2; 


if strfind(upper(MODE),'PCA') % strfind(upper(MODE),'SVD')
	if flagNoSegmentSupport, 	
		error('non-continous data is not supported by MODE=%s',MODE);
	end; 

        % fit and remove sine/cosine wave 
        % Assumption: stationarity of line interference

        bw = 1;  % bandwidth of 1 Hz
        r = exp(-bw/HDR.SampleRate);	% convert bandwidth into magnitude of poles
        %B = real(r*poly(exp(i*2*pi/HDR.SampleRate*[NotchFreq:2*NotchFreq:HDR.SampleRate-NotchFreq])));
        B = real(poly(r*exp(i*2*pi*NotchFreq/HDR.SampleRate*[1,-1])));

        %%%%% identify 50/60 Hz component
        [u,d,v]  = svd(filtfilt(1,B,s),0); % s = u*d*v';
        [tmp,ix] = max(diag(d));

        %%%%% regress 50/60 Hz component
        [R,S] = regress_eog([s,u(:,ix)],1:size(s,2),size(s,2)+1);	% regression on filtered component
        s = S(:,1:size(s,2));
        %[R2,s2] = regress_eog(s,1:size(s,2),v(:,ix));			% regression on raw component


elseif ~isempty(strfind(upper(MODE),'FIT'))
	if flagNoSegmentSupport, 	
		error('non-continous data is not supported by MODE=%s',MODE);
	end; 

        % fit and remove sine/cosine wave
        % Assumption: stationarity and locking of line and sampling frequency 
        t = [1:size(s,1)]'/HDR.SampleRate;
        u = exp(j*2*pi*NotchFreq/HDR.SampleRate*[1:size(s,1)]');
        A = covm(u,s);
        s = s - 2*real(u*A);


elseif strfind(upper(MODE),'ICA')

        %%%%% not implemented %%%%%
        error('ICA method not implemented');


elseif ~isempty(strfind(upper(MODE),'NOTCH'))
	if flagNoSegmentSupport, 	
		error('non-continous data is not supported by MODE=%s',MODE);
	end; 

        if ~isempty(strfind(upper(MODE),'+'))
                B = real(poly(exp(i*2*pi/HDR.SampleRate*[NotchFreq,HDR.SampleRate-NotchFreq])));
        else
                warning('NOTCH+ - mode is not recommended');        
                B = real(poly(exp(i*2*pi/HDR.SampleRate*[NotchFreq:2*NotchFreq:HDR.SampleRate-NotchFreq/2])));
        end;
        s = filtfilt(B/sum(B),1,s);


elseif ~isempty(strfind(upper(MODE),'FFT'))
    	% cancel 50/60 Hz
	b = 1/2;   % half bandwidth 

	for k=1:length(segs)-1,	
		%% support segmented data 
		u = fft(s(segs(k):segs(k+1)-1,:),[],1); 
		f = [0,segs(k+1)-segs(k)-1];
		u((NotchFreq-b) < f & f < (NotchFreq+b),:) = 0; 
		u((NotchFreq-b) < (HDR.SampleRate-f) & (HDR.SampleRate-f) < (NotchFreq+b),:) = 0; 
		if ~isempty(strfind(upper(MODE),'FFT3'))
			% cancel 3rd harmonic 150/180 Hz
			u(3*(NotchFreq-b) < f & f < 3*(NotchFreq+b),:) = 0; 
			u(3*(NotchFreq-b) < (HDR.SampleRate-f) & (HDR.SampleRate-f) < 3*(NotchFreq+b),:) = 0; 
		end; 
		s(segs(k):segs(k+1)-1,:) = real(ifft(u,[],1)); 
	end; 


elseif ~isempty(strfind(upper(MODE),'USER-SPECIFIC')) && (abs(HDR.SampleRate-1000)<1e3),
        % Sebastians filter 
        B = load('49__51bs_2751pts_1000Hz.asc');
        s = filtfilt(B,1,s);

else
        warning('unknown method: no correction applied.');

end;


if ~isempty(outfile)
	HDR.FileName = outfile; 
        HDR.TYPE = 'GDF'; 
	HDR.VERSION = 3; 
	HDR.GDFTYP(:) = 16;
	try, 
		mexSSAVE(HDR,s); 
	catch 
		HDR = sopen(HDR,'w'); 
		HDR = swrite(HDR,s); 
		HDR = sclose(HDR); 	
	end; 
end; 	


