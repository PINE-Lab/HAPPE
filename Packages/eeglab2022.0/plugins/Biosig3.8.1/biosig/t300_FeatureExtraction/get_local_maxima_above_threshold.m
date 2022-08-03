function [pos] = get_local_maxima_above_threshold(data, TH, Mode, winlen)
% GET_LOCAL_MAXIMA_ABOVE_THRESHOLD is used to identify the events
% derived from the trigger trace. It has been used to identify
% miniature EPSP and IPSP events. The functionality of Mode=0 resembles a
% function implemented in FBrain [1-3].
% This functions is provided for backwards compatibility, for new projects you
% should consider using findpeaks() instead.
%
% pos = get_local_maxima_above_threshold(data,TH)
% pos = get_local_maxima_above_threshold(data,TH,Mode)
% pos = get_local_maxima_above_threshold(data,TH,Mode,winlen)
%
% Input:
%   data: sample vector of detection trace
%   TH: threshold
%   Mode==0: [default], detects all (local) maxima above threshold [1]
%         1: only one maximum above threshold within the
%               window of size winlen is considered [2].
%         2: only single detections are considered.
%               if two detections with a distance smaller than
%               winlen occur, both are omitted. That might be
%               useful for optain clean candiate templates.
%         3: single detection for period between positive
%               and negative threshold crossing (winlen is ignored)
%         4: single detection for period between positive
%               and negative threshold, however the distance
%               winlen must exceed the distance between the
%               positive and negative threshold crossing
%
%   winlen (Mode=1 only): window length (in number of samples)
%               in which all detections collapse to one event
%
% Output:
%   pos: time points of local maxima above threshold
%
% see also: signal_deconvolution, findpeaks
%
% Reference(s): 
% [1] A. Pernía-Andrade, S.P. Goswami, Y. Stickler, U. Fröbe, A. Schlögl, and P. Jonas (2012)
%     A deconvolution-based method with high sensitivity and temporal resolution for 
%     detection of spontaneous synaptic currents in vitro and in vivo.
%     Biophysical Journal Volume 103 October 2012 1–11.
% [2] Zhang X., Schloegl A., Vandael D., Jonas P. (2020)
%      A novel machine learning-based method for accurate and efficient detection of
%      subthreshold synaptic events in vivo and in vitro (in revision).
%      Journal of Neuroscience Methods.
% [3] Zhang X, Schlögl A, Vandael D, Jonas P (2021),
%     MOD: A novel machine-learning optimal-filtering method for accurate and efficient
%        detection of subthreshold synaptic events in vivo.
%     Journal of Neuroscience Methods, 2021.
%     doi:10.1016/j.jneumeth.2021.109125


% Copyright 2011,2012,2016,2019,2020,2021 Alois Schloegl, IST Austria <alois.schloegl@ist.ac.at>
% This is part of the BIOSIG-toolbox http://biosig.sf.net/ 


if nargin<3,
	Mode=0;
end;
if nargin<4,
	winlen=1;
end;

if any(Mode==[0,2]) && (numel(TH)==1),
	%%% InVitro Data from Sarit
	data = [data;+inf];
	pos = (data(1:end-1) >= TH) & (data(1:end-1) >  data(2:end) ) & (data(1:end-1) > [+inf;data(1:end-2)]);		%% local maxima above threshold
	ix  = (data(1:end-1) >= TH) & (data(1:end-1) == data(2:end) ) & (data(1:end-1) > [+inf;data(1:end-2)]);		%% local maxima above threshold

elseif (Mode==0) && (numel(TH)==2),
	%%% InVivo data from Alejo
	th = repmat(NaN,size(data));
	th(1:end/2) = TH(1);
	th(end/2+1:end) = TH(2);

	data = [data;+inf];
	pos = (data(1:end-1) >= th) & (data(1:end-1) >  data(2:end) ) & (data(1:end-1) > [+inf;data(1:end-2)]);		%% local maxima above threshold
	ix  = (data(1:end-1) >= th) & (data(1:end-1) == data(2:end) ) & (data(1:end-1) > [+inf;data(1:end-2)]);		%% local maxima above threshold

elseif (Mode==1)
	% approach Mode==1
	% only one maximum above threshold within the window of size winlen is considered.
	dix    = diff([0; data>TH; 0]);
	onset  = find(dix>0);
	offset = find(dix<0);

	pos  = [];
	stop = 0;
	while 1,
		k     = find(onset > stop, 1, 'first');
		if isempty(k) break; end
		start = onset(k);

		ix = find((start < offset) & (offset < (start+winlen)));
		if ~isempty(ix)
			stop = offset(ix(end));
		else
			stop = offset(k);
		end;
		[tmp,ix]=max(data(start:stop-1));

		pos(end+1,1)=start+ix-1;
	end;
	return;

elseif (Mode==3)
	d = diff([-inf;data;-inf] > TH);
	onset  = find(d>0);
	offset = find(d<0)-1;
	A      = repmat(NaN,size(onset));
	pos    = repmat(NaN,size(onset));
	for k  = 1:length(onset),
		[A,tix] = max(data(onset(k):offset(k)));
		pos(k)  = onset(k)+tix-1;
	end
	return;

elseif (Mode==4)
	d      = diff([-inf;data;-inf] > TH);
	onset  = find(d > 0);
	offset = find(d < 0);
	A      = repmat(NaN, size(onset));
	pos    = repmat(NaN, size(onset));
	for k  = find((offset - onset)' >= winlen),
		[A,tix] = max(data(onset(k):offset(k)-1));
		pos(k)  = onset(k)+tix-1;
	end
	fprintf(1,'%d spurious events removed\n',sum(isnan(pos)));
	pos = pos(~isnan(pos));
	return;

else
	error('invalid input argument')
end 

k   = 0; 
ix  = find(ix);
while (1) 
	k   = k+1;
	ix2 = ix + k;
	ix  = ix(data(ix2) <= data(ix));
	ix2 = ix + k;
	if (isempty(ix) || ~any(data(ix)==data(ix2))) break; end; 
end
pos = sort([find(pos); ix]);

if (Mode==2)
	% remove those events that are to close together
        % !!! NOTE: In case of too many detection, all events are removed !!!
	tmp = find(diff(pos)<winlen);
	pos(tmp)=NaN;
	pos(tmp+1)=NaN;
	pos = pos(~isnan(pos));
end

return

%!assert(get_local_maxima_above_threshold([0,0,0,1,0]',0),4)
%!assert(get_local_maxima_above_threshold([0,0,0,1,1,0]',0),4)
%!assert(get_local_maxima_above_threshold([0,0,0,1,1,2,0]',0),6)
%!assert(get_local_maxima_above_threshold([0,0,0,3,1,1,0]',0),4)
%!assert(get_local_maxima_above_threshold([0,0,0,3,1,1,2,0]',0),[4;7])
%!assert(get_local_maxima_above_threshold([0,0,0,1,1,1,0]',0),4)
%!assert(get_local_maxima_above_threshold([0,1,2,3,4,3,2,1,0]',3,3),5)
%!assert(get_local_maxima_above_threshold([0,1,2,3,4,3,2,1,0]',2,4),5)
%!assert(get_local_maxima_above_threshold([0,4,1,3,4,3,2,1,0]',2,4,2),5)
%!assert(length(get_local_maxima_above_threshold([0,4,1,3,4,3,2,1,0]',2,4,1)),2)
%!assert(length(get_local_maxima_above_threshold([0,4,1,3,4,3,2,1,0]',2,4,2)),1)

