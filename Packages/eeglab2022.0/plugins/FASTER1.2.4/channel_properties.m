function list_properties = channel_properties(EEG,eeg_chans,ref_chan)

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


if ~isstruct(EEG)
	newdata=EEG;
	clear EEG;
	EEG.data=newdata;
	clear newdata;
end

measure = 1;

if ~isempty(ref_chan) && length(ref_chan)==1
    pol_dist=distancematrix(EEG,eeg_chans);
    [s_pol_dist dist_inds] = sort(pol_dist(ref_chan,eeg_chans));
    [s_inds idist_inds] = sort(dist_inds);
end

% TEMPORAL PROPERTIES

% 1 Mean correlation between each channel and all other channels

% Ignore zeroed channels (ie reference channels) to avoid NaN problems
ignore = [];
datacorr = EEG.data;
for u = eeg_chans
	if max(EEG.data(u,:))==0 && min(EEG.data(u,:))==0
		ignore=[ignore u];
	end
end

% Calculate correlations
calc_indices=setdiff(eeg_chans,ignore);
ignore_indices=intersect(eeg_chans,ignore);
corrs = abs(corrcoef(EEG.data(setdiff(eeg_chans,ignore),:)'));
mcorrs=zeros(size(eeg_chans));
for u=1:length(calc_indices)
    mcorrs(calc_indices(u))=mean(corrs(u,:));
end
mcorrs(ignore_indices)=mean(mcorrs(calc_indices));

% Quadratic correction for distance from reference electrode

if (~isempty(ref_chan) && length(ref_chan)==1)
    p = polyfit(s_pol_dist,mcorrs(dist_inds),2);
    fitcurve = polyval(p,s_pol_dist);
    corrected = mcorrs(dist_inds) - fitcurve(idist_inds);

    list_properties(:,measure) = corrected;
else
    list_properties(:,measure) = mcorrs(dist_inds);
end
measure = measure + 1;

% 3 Variance of the channels
vars = var(EEG.data(eeg_chans,:)');
vars(~isfinite(vars))=mean(vars(isfinite(vars)));
% Quadratic correction for distance from reference electrode

if (~isempty(ref_chan) && length(ref_chan)==1)
    p = polyfit(s_pol_dist,vars(dist_inds),2);
    fitcurve = polyval(p,s_pol_dist);
    corrected = vars - fitcurve(idist_inds);

    list_properties(:,measure) = corrected;
else
    list_properties(:,measure) = vars;
end
measure = measure + 1;

% 4 Hurst exponent
for u=1:length(eeg_chans)
    list_properties(u,measure) = hurst_exponent(EEG.data(eeg_chans(u),:));
end

for u = 1:size(list_properties,2)
    list_properties(isnan(list_properties(:,u)),u)=nanmean(list_properties(:,u));
	list_properties(:,u) = list_properties(:,u) - median(list_properties(:,u));
end