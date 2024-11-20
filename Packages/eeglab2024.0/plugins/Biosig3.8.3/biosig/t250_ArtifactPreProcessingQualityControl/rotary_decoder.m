function [e,v] = rotary_decoder(s,chan)
% ROTARY_DECODER decodes the spatial position from the
%   two TTL signals from the spokes of a wheel
%   This functions has been developed and used in [1,2] 
%   for obtaining the spatial position. An earlier version of 
%   of this functions was called GET_SPATIAL_POSITION.
%
% Usage:
%   [e,v] = rotary_decoder(filename)
%   [e,v] = rotary_decoder(filename,chan)
%   e = rotary_decoder(data)
%
%  Input:
% 	filename: channel 7 and 8 need to contain the 
%		TTL-signals from the spokes
%	data: must a matrix with two columns, that 
%		contain the TTL-signals of the two spokes
%       chan: channels containing the TTL-signals
% 		default: channels with labels 'Adc-9' 
%               and 'Adc-10', or channel [7,8].
%
%  Output:
%       e: spatial position at each sample in time.
%       v: velocity (estimated based on distance differenc within on 25ms) 
%
% REFERENCES:
% [1] Zhang X, Schlögl A, Vandael D, Jonas P, 
%     MOD: A novel machine-learning optimal-filtering method for accurate and efficient detection of subthreshold synaptic events in vivo. 
%     Journal of Neuroscience Methods, 2021.  
%     https://doi.org/10.1016/j.jneumeth.2021.109125
% [2] Xiaomin Zhang, Alois Schlögl, Peter Jonas
%     Selective Routing of Spatial Information Flow from Input to Output in Hippocampal Granule Cells, Neuron, 2020. 
%     https://doi.org/10.1016/j.neuron.2020.07.006 
%     http://www.sciencedirect.com/science/article/pii/S0896627320305237 
% 

%    Copyright (C) 2017-2021 by Alois Schloegl, IST Austria
%    This is part of the BIOSIG-toolbox http://biosig.sf.net/
%
%    BioSig is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    BioSig is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with BioSig.  If not, see <http://www.gnu.org/licenses/>.


	TTL_THRESHOLD = 1.5;
	if nargin<2
		chan=[];
	end

	SampleRate = [];
	if ischar(s) && exist(s,'file'),
		fn    = s;
		[d,H] = mexSLOAD(fn, 0, 'OVERFLOWDETECTION:OFF');
		SampleRate = H.SampleRate;
		
		TTL2 = find(strcmp(H.Label,'Adc-9'));
		TTL3 = find(strcmp(H.Label,'Adc-10'));
		chan = [TTL2,TTL3]; 
		if numel(chan)~=2,
			chan=7:8;
		end

		d  = d(:,chan) > TTL_THRESHOLD;

	elseif isnumeric(s) && size(s,2)==2,
		d  = s > TTL_THRESHOLD;

	else
		error('invalid input arguments');
	end;

	ds = diff(d);

	tix1u = find(ds(:,1) > 0);
	tix2u = find(ds(:,2) > 0);
	tix1d = find(ds(:,1) < 0);
	tix2d = find(ds(:,2) < 0);


if  isempty(tix1d) ||  isempty(tix2d) || isempty(tix1u) ||  isempty(tix2u)
	e = []; 
	return
end

	diameter = 10; 		% diameter of wheel in [cm]
	nspokes  = 500;		% number of spokes
	% 4 state changes per spoke

	e        = zeros(size(d,1),1);
	e(tix1u) = e(tix1u) - 0.5 + d(tix1u,2);	% #1 up   #2 low
	e(tix1d) = e(tix1d) + 0.5 - d(tix1d,2);	% #1 down #2 hi
	e(tix2u) = e(tix2u) + 0.5 - d(tix2u,1);	% #2 up   #1 hi 
	e(tix2d) = e(tix2d) - 0.5 + d(tix2d,1);	% #2 down #1 low

	e        = cumsum(e) * (pi*diameter) / (2*nspokes);

	if (nargout>1) && ~isempty(SampleRate)
		% compute velocity based on window size twin
		twin = 0.025*5; % in seconds
		v = filter([1, zeros(1, round(twin*SampleRate)-1), -1], 1/twin, e);
		fprintf(1,'velocity estimation is delayed by %g milliseconds\n',twin*1000/2);
	end


	% State table for identifying a step forward (+1), 
	%    or backwards (-1), or no state change (oo). 

	% #1(t) #1(t+1) #2(t) #2(t+1)
	%
	%   0 0 0 0 	oo	
	%   0 0 0 1 	-1
	%   0 0 1 0 	+1	f
	%   0 0 1 1	oo	
	%
	%   0 1 0 0 	+1	f
	%   0 1 0 1	oo
	%   0 1 1 0 	oo
	%   0 1 1 1	-1
	%
	%   1 0 0 0 	-1
	%   1 0 0 1	oo
	%   1 0 1 0 	oo
	%   1 0 1 1	+1	f
	%
	%   1 1 0 0 	oo
	%   1 1 0 1	+1	f
	%   1 1 1 0 	-1
	%   1 1 1 1	oo

