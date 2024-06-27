function [distmatrixpol distmatrixxyz distmatrixproj] = distancematrix(EEG,eeg_chans)

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

num_chans = size(EEG.data,1);
distmatrix = zeros(length(eeg_chans),length(eeg_chans));
distmatrixpol = [];
for chan2tst = eeg_chans;
	for q=eeg_chans
		distmatrixpol(chan2tst,q)=sqrt(((EEG.chanlocs(chan2tst).radius^2)+(EEG.chanlocs(q).radius^2))-(2*((EEG.chanlocs(chan2tst).radius)*...
			(EEG.chanlocs(q).radius)*cosd(EEG.chanlocs(chan2tst).theta - EEG.chanlocs(q).theta))));%calculates the distance between electrodes using polar format
	end
end

locs = EEG.chanlocs;
for u = eeg_chans
	if ~isempty(locs(u).X)
		Xs(u) = locs(u).X;
	else
		Xs(u) = 0;
	end
	if ~isempty(locs(u).Y)
		Ys(u) = locs(u).Y;
	else
		Ys(u) = 0;
		end
	if ~isempty(locs(u).Z)
		Zs(u) = locs(u).Z;
	else
		Zs(u) = 0;
	end
end
Xs = round2(Xs,6);
Ys = round2(Ys,6);
Zs = round2(Zs,6);

for u = eeg_chans
	for v=eeg_chans
		distmatrixxyz(u,v) = dist(Xs(u),Xs(v))+dist(Ys(u),Ys(v))+dist(Zs(u),Zs(v));
	end
end
D = max(max(distmatrixxyz));
distmatrixproj = (pi-2*(acos(distmatrixxyz./D))).*(D./2);
	function d = dist(in1,in2)
		d = sqrt(abs(in1.^2 - in2.^2));
	end

	function num = round2(num,decimal)
		num = num .* 10^decimal;
		num = round(num);
		num = num ./ 10^decimal;
	end
end

