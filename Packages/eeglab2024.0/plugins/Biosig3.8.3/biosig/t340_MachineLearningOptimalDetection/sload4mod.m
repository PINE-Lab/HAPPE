function [data,HDR,rawdata] = sload4mod(datafile,chan,Fs,MODE,TemplateLength,BlankingPeriod)
% sload_mod is a function that loads biosignal data and 
%   and does all preprocessing steps for the MOD method. 
% Preprocessing steps consist of 
% - AP detection and blanking of AP's  
% - Resampling to a common target sampling rate (default: 25kHz)
% - detrending and/or highpass filtering
% - lowpass filtering
%
% Usage: 
%  [data,HDR]=sload4mod(filename,chan,Fs,Mode,TemplateLength);
%  [data,HDR]=sload4mod(filename,chan,Fs,Mode,TemplateLength,BlankingPeriod);
%
%
%
% Fs:  target sampling rate (default: 25000 Hz)
% BlankingPeriod: 
%       []: default, no blaninking
%       +t   if t>0, from the beginning of each sweep until t is blanked
%            if t<0  from end-t to the end is blanked	
%       [t1,t2]  interval [t1,t2] is blanked	
%  
% Mode: 'fft:10:1000'
% 	'gauss:10:1000'  Gaussian filter with 10-1000 Hz bandpass
%	'fir1:0:1000'	fir filter
%	'D1a:0:1000'	fir1-forward-backward filter
%	'D1a:-1:1000'	1st order difference filter + fir1-forward-backward filter
%	'D1b:0:1000'	fir1-forward-backward filter
%	'D1b:-1:1000'	1st order difference filter + fir1-forward-backward filter

% TODO: 
% - identify RS pulses
% + downsampling using fft (for Gauss and FFT filter)
% + check whether AP detection returns all events or only AP events 
% + merge AP events with all other events, this is done in minidet2gdf
% - remove marked artifacts - 

if (nargin<2) || isempty(chan)
	chan = 1; %find(strcmp('Adc0',HDR.Label));
end
if (nargin<3) || isempty(Fs)
	Fs = 25000; % kHz
end
if (nargin<6)
	BlankingPeriod = [];
end

if nargin<4,
	Mode.FTYPE = 'fir1';
	Mode.LOWPASS  = 1000; 
	Mode.HIGHPASS = 0; 
	Mode.WINLEN   = 4; 	% window size [ms] for scoring trace 
elseif ischar(MODE)
	C=strsplit(MODE,'_');
	Mode.FTYPE = C{1};
	Mode.HIGHPASS = str2num(C{2}); 
	Mode.LOWPASS = str2num(C{3}); 
	Mode.WINLEN   = 4; 	% window size [ms] for scoring trace 
else
	error('invalid argument Mode ')	
end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       read in data 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	[data,HDR]     = mexSLOAD(datafile,0);
	HDR.SampleRate = round(HDR.SampleRate);
	S = data(:,chan);        
	[HDR.FILE.Path,HDR.FILE.Name,HDR.FILE.Ext]=fileparts(datafile);
        S0 = S;        % backup for unfiltered data

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       remove marked artifacts
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%  artifact markers  %%%%%
        EVT = HDR.EVENT;
        if ~isfield(EVT,'DUR')
                EVT.DUR=zeros(size(EVT.POS));
        end
        aix = find(EVT.TYP==20);
        artifact = [EVT.POS(aix),EVT.POS(aix)+EVT.DUR(aix)];
        Stmp = S; 
        for k=1:length(aix)
                Stmp(EVT.POS(aix):EVT.POS(aix)+EVT.DUR(aix),:)=NaN;
        end;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%       Blank RS pulse 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        selpos = sort([1;HDR.EVENT.POS(HDR.EVENT.TYP==hex2dec('7ffe')); size(S,1)+1]);
        bi = [0,0]; 
        if length(BlankingPeriod)==1, 
        	if BlankingPeriod>0,
        		bi=[0,BlankingPeriod];
        	elseif BlankingPeriod<0,
        		bi=[BlankingPeriod,0];
        	end
        elseif 	length(BlankingPeriod)==2,
        	bi = BlankingPeriod;
        end
	[ix,iy] = meshgrid(selpos-1,round(bi(1)*HDR.SampleRate):round(bi(2)*HDR.SampleRate)-1);
        ix = ix(:)+iy(:); 
        ix = ix( (0 < ix) & (ix <= size(S,1)) );

        Stmp(ix) = NaN;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%       Detect and remove AP's and spikelets
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	AP = detect_spikes_bursts(HDR,data,'chan',chan);
	% AP = detect_spikes_bursts(datafile, chan);
	AP.EVENT.CHN(:)=1; % raw data is stored into channel 1 when writing file 

	DIV = round(HDR.SampleRate/Fs);
	if DIV~=1, 
		% downsample to 25kHz
		% S   = rs(S,DIV,1);
		HDR.EVENT.POS = round((HDR.EVENT.POS-1)*Fs/HDR.EVENT.SampleRate)+1;
		  
		if isfield(HDR.EVENT,'DUR');
			HDR.EVENT.DUR = round(HDR.EVENT.DUR*Fs/HDR.EVENT.SampleRate);
		end
		HDR.EVENT.SampleRate = Fs;

		if ~isempty(AP) && isfield(AP,'EVENT')
		        AP.EVENT.POS = round((AP.EVENT.POS-1)*Fs/AP.EVENT.SampleRate)+1;
		        if isfield(HDR.EVENT,'DUR');
		                AP.EVENT.DUR = round(AP.EVENT.DUR*Fs/AP.EVENT.SampleRate);
		        end
		        AP.EVENT.SampleRate = Fs;
		end
		selpos2=round((selpos-1)/DIV)+1;
	else
		selpos2=selpos;
	end;

	sweepLen= diff(selpos); 
	FLAG_ALL_SWEEPS_HAVE_EQUAL_LENGTH=0; 
	if all(sweepLen==sweepLen(1)), 
		sweepLen=sweepLen(1);
		FLAG_ALL_SWEEPS_HAVE_EQUAL_LENGTH=1; 
       	        ff = [0:sweepLen-1]'*HDR.SampleRate/sweepLen;

                ffix = (Mode.HIGHPASS < ff) | (ff > Mode.LOWPASS);	% used in fft filter

       	        gfix = min(ff, HDR.SampleRate-ff);			% used in gauss filter
		if (Mode.HIGHPASS > 0)
	       	        G0 = exp(-0.5 * (gfix / Mode.HIGHPASS).^2);
	       	else
	       		G0 = 0; 
	       	end	        
		if (Mode.LOWPASS > 0)
       	        	G1 = exp(-0.5 * (gfix / Mode.LOWPASS).^2);
       	        else 
       	        	G1 = 1;	
       	        end
       	        G2 = G1.*(1-G0);
       	        G1 = G1 - G0; 
	end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%  50 Hertz Notch and Lowpass %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(Mode.FTYPE,'fft')
        	%  does the resampling through fft, not rs()
        	Stmp2 = repmat(NaN,size(Stmp,1)/DIV,1);
	        for k2 = 1:length(selpos)-1,
        	        d    = Stmp(selpos(k2):selpos(k2+1)-1,:);
                	dnan = isnan(d);
	                if any(dnan),
				d(dnan) = mean(d);
                	end
	                D  = fft(d); 
	                if ~FLAG_ALL_SWEEPS_HAVE_EQUAL_LENGTH,
	        	        ff = [0:length(d)-1]'*HDR.SampleRate/length(d);
		                ffix = (Mode.HIGHPASS < ff) | (ff > Mode.LOWPASS);	
	        	end
        	        D( 2*ff >= Fs) = 0; 	% Nyquist for downsampling 
                	%fix = ((ff > 49.5) & (ff < 50.5)) | (ff > Fs/2);
        	        D( ffix ) = 0; 
                	D(1) = D(1) / 2;

	                d2   = real(ifft(D(ff<Fs)))*2; 	%% fft-based downsampling to Fs
	                dnan = any(reshape(dnan,DIV,[]),1)';
        	        if any(dnan)
				d2(dnan) = NaN;
	                end
        	        Stmp2(selpos2(k2):selpos2(k2+1)-1,:) = d2;
	        end;
	        Stmp = Stmp2;	
		HDR.SampleRate = Fs; 

        elseif strcmp(Mode.FTYPE,'gauss2')
        	%  does the resampling through fft, not rs()
        	Stmp2 = repmat(NaN,size(Stmp,1)/DIV,1);
	        for k2 = 1:length(selpos)-1,
        	        d    = Stmp(selpos(k2):selpos(k2+1)-1,:);
                	dnan = isnan(d);
	                if any(dnan)
				d(dnan) = mean(d);
                	end
	                D  = fft(d); 
	                if ~FLAG_ALL_SWEEPS_HAVE_EQUAL_LENGTH,
	        	        ff = [0:length(d)-1]'*HDR.SampleRate/length(d);
	        	        gfix = min(ff, HDR.SampleRate-ff);
	        	        if (Mode.HIGHPASS > 0)
		        	        G0 = exp(-0.5 * (gfix / Mode.HIGHPASS).^2);
				else
					G0 = 0;
				end	
	        	        if (Mode.LOWPASS > 0)
					G1 = exp(-0.5 * (gfix / Mode.LOWPASS).^2);
				else
					G1 = 1;
				end
				G2 = G1.*(1-G0);
	        	end
        	        D( 2*ff >= Fs) = 0; 	% Nyquist for downsampling 
                	%fix = ((ff > 49.5) & (ff < 50.5)) | (ff > Fs/2);
        	        D = D.*G2;	% gauss filter in F-domain

	                d2 = real(ifft(D(ff<Fs)))*2; 	%% fft-based downsampling to Fs
	                dnan=any(reshape(dnan,DIV,[]),1)';
        	        if any(dnan)
				d2(dnan)=NaN;
	                end
        	        Stmp2(selpos2(k2):selpos2(k2+1)-1,:) = d2;
	        end;
	        Stmp = Stmp2;	
		HDR.SampleRate = Fs; 

        elseif strcmp(Mode.FTYPE,'gauss')
        	%  does the resampling through fft, not rs()
        	Stmp2 = repmat(NaN,size(Stmp,1)/DIV,1);
	        for k2 = 1:length(selpos)-1,
        	        d    = Stmp(selpos(k2):selpos(k2+1)-1,:);
                	dnan = isnan(d);
	                if any(dnan)
				d(dnan) = mean(d);
                	end
	                D  = fft(d); 
	                if ~FLAG_ALL_SWEEPS_HAVE_EQUAL_LENGTH,
	        	        ff = [0:length(d)-1]'*HDR.SampleRate/length(d);
	        	        gfix = min(ff, HDR.SampleRate-ff);
	        	        if (Mode.HIGHPASS > 0)
		        	        G0 = exp(-0.5 * (gfix / Mode.HIGHPASS).^2);
				else
					G0 = 0;
				end	
	        	        if (Mode.LOWPASS > 0)
					G1 = exp(-0.5 * (gfix / Mode.LOWPASS).^2) - G0;
				else
					G1 = 1 - G0;
				end
	        	end
        	        D( 2*ff >= Fs) = 0; 	% Nyquist for downsampling 
                	%fix = ((ff > 49.5) & (ff < 50.5)) | (ff > Fs/2);
        	        D = D.*G1;	% gauss filter in F-domain

	                d2 = real(ifft(D(ff<Fs)))*2; 	%% fft-based downsampling to Fs
	                dnan=any(reshape(dnan,DIV,[]),1)';
        	        if any(dnan)
				d2(dnan)=NaN;
	                end
        	        Stmp2(selpos2(k2):selpos2(k2+1)-1,:) = d2;
	        end;
	        Stmp = Stmp2;	
		HDR.SampleRate = Fs; 

        elseif strcmp(Mode.FTYPE,'fir1')
		HDR.SampleRate = Fs; 
        	Stmp = rs(Stmp,DIV,1);
	        if (Mode.HIGHPASS<0)
	        	Stmp=[NaN;diff(Stmp)];
	        end
	        if ((0 < Mode.LOWPASS) && (Mode.LOWPASS*2 < Fs))
	                B=fir1 (100, 2*Mode.LOWPASS/Fs);
	                Stmp=filter(B,1,Stmp(end:-1:1));
	                Stmp=filter(B,1,Stmp(end:-1:1));
		end

        elseif strcmp(Mode.FTYPE,'D1a')
		HDR.SampleRate = Fs; 
        	Stmp = rs(Stmp,DIV,1);
	        if (Mode.HIGHPASS<0)
	        	Stmp=[NaN;diff(Stmp)];
	        end
	        if ((0 < Mode.LOWPASS) && (Mode.LOWPASS*2 < Fs))
	                B=fir1 (100, 2*Mode.LOWPASS/Fs);
        	        Stmp=filter(B,1,Stmp(end:-1:1));
        	        Stmp=filter(B,1,Stmp(end:-1:1));
		end
        elseif strcmp(Mode.FTYPE,'D1b')
		HDR.SampleRate = Fs; 
        	Stmp = rs(Stmp,DIV,1);
	        if ((0 < Mode.LOWPASS) && (Mode.LOWPASS*2 < Fs))
	                B=fir1 (100, 2*Mode.LOWPASS/Fs);
	                Stmp=filter(B,1,Stmp(end:-1:1));
	                Stmp=filter(B,1,Stmp(end:-1:1));
		end
	        if (Mode.HIGHPASS<0)
	        	Stmp=[NaN;diff(Stmp)];
	        end

 	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%       Detect and remove AP's and spikelets
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	apix = find(AP.EVENT.TYP==hex2dec('201'));
	Stmp(AP.EVENT.POS(apix),:) = NaN;       % trough the template (40ms, 1000samples) the NaN will be expanded to [-10:+30] ms
	[ix,iy] = meshgrid(AP.EVENT.POS(apix),[-round(TemplateLength/4):round(TemplateLength*3/4)]);
	Stmp(ix(:)+iy(:)) = NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawdata = data;
data = Stmp; 

