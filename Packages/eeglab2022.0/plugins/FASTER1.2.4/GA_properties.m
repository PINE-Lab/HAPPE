function list_properties = GA_properties(EEG,eeg_chans,EOG_chans)

% Copyright (C) 2010 Hugh Nolan, Robert Whelan and Richard Reilly, Trinity College Dublin,
% Ireland
% nolanhu@tcd.ie, robert.whelan@tcd.ie
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

list_properties = [];

if length(size(EEG.data)) < 3
	fprintf('Not epoched.\n');
	return;
end

measure = 1;

means = mean(EEG.data(eeg_chans,:),2);

% 1 Epoch's mean deviation from channel means.
for u = 1:size(EEG.data,3)
	list_properties(u,measure) = mean(abs(squeeze(mean(EEG.data(eeg_chans,:,u),2)) - means));
end
measure = measure + 1;

% 2 Epoch variance
list_properties(:,measure) = mean(squeeze(var(EEG.data(eeg_chans,:,:),0,2)));
measure = measure + 1;

% 3 Max amplitude difference
for t = eeg_chans
	for u = 1:size(EEG.data,3)
		ampdiffs(t,u) = max(EEG.data(t,:,u)) - min(EEG.data(t,:,u));
	end
end
list_properties(:,measure) = mean(ampdiffs,1);
measure = measure + 1;

% 4 EOG channel max value
if ~isempty(EOG_chans)
list_properties(:,measure) = max(max(abs(EEG.data(EOG_chans,:,:)),[],2),[],1);
measure = measure + 1;
end

for v = 1:size(list_properties,2)
	list_properties(:,v) = list_properties(:,v) - median(list_properties(:,v));
end