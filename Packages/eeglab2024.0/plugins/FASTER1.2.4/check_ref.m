function EEG=check_ref(EEG,ref_chan,num_chans,num_exts)

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

if all(round10(EEG.data(ref_chan,:),10)==0)
	EEG=h_pop_reref(EEG,[],'keepref','on');
end
mediffs=zeros(1,size(EEG.data,3));vars=mediffs;
for v=1:size(EEG.data,3)
	mediffs(v)=median(diff(EEG.data(ref_chan,:,v),1,2),2);
	vars(v)=var(EEG.data(ref_chan,:,v),[],2);
end
bad_subjs=union(find(min_z(mediffs)),find(min_z(vars)));
for v=1:length(bad_subjs)
	tempEEG=h_eeg_interp_spl(pop_rejepoch(EEG,setdiff(1:size(EEG.data,3),bad_subjs(v)),0),ref_chan,'spherical',num_chans+1:num_chans+num_exts);
	EEG.data(:,:,bad_subjs(v))=tempEEG.data;
end