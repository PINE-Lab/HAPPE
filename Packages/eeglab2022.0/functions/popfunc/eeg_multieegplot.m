% eeg_multieegplot() - Produce an eegplot() of a the average of an epoched dataset 
%                  (with optional pre-labelling of specific trials).
% Usage:
%   >> eeg_multieegplot( data,trialrej, elecrej, ...
%                                'key1', value, 'key2', value ... );
% Inputs:
%   data        - input data (channels x points or channels x points x trials).
%   trialrej    - array of 0s and 1s (depicting rejected trials) (size sweeps)
%   elecrej     - array of 0s and 1s (depicting electrodes rejected in 
%                 all trials) (size nbelectrodes x sweeps )
%   oldtrialrej - array of 0s and 1s (depicting rejected trials) (size sweeps)
%   oldelecrej  - array of 0s and 1s (depicting electrodes rejected in 
%                 all trials) (size nbelectrodes x sweeps )
%
% Note: 1) {'Key', value } Arguments are passed on to eegplot() 
%       2) To ignore previous rejections simply set 'oldtrialrej' and
%          'oldelecrej' arguments to empty ([]).
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eegplot(), eegplot2event(), eegplot2trial(), eeglab()

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

% 01-25-02 reformated help & license -ad 
% 03-07-02 corrected help -ad

function eeg_multieegplot( data, rej, rejE, oldrej, oldrejE, varargin);

if nargin < 1
   help eeg_multieegplot;
   return;
end
if nargin < 4
   oldrej = [];
end
if nargin < 5
   oldrejE = [];
end
if ~exist('command')
	command = '';
end
chans = size(data,1);
pnts  = size(data,2);

colnew = [0.8 0.8 1];
colold = [0.8 1 0.8];
 
if ndims(data) > 2 % --------------- considering epoched EEG data
   rejeegplot = [];
   % concatenate old rejection if not empty
   % --------------------------------------
   if ~isempty( oldrej )
   		oldrej 	 = find( oldrej > 0);
   		oldrejE  = oldrejE(:, oldrej)';
   		rejeegplot = trial2eegplot( oldrej, oldrejE, pnts, colold);
   end;	
   % convert for eegplot
   % -------------------
   if ~isempty(	rej )
   		rej    = find( rej > 0);
   		rejE   = rejE(:,rej)';
   		secondrejeegplot = trial2eegplot( rej, rejE, pnts, colnew); % see bottom of code for this function
   		rejeegplot = [ rejeegplot' secondrejeegplot' ]';

   		% remove duplicates
   		% -----------------
		%[tmp I] = unique_bc( rejeegplot(:,1) );
		%rejeegplot = rejeegplot(I,:);
   end
else % ---------------------------------------- considering continuous EEG
   % for continuous EEG, electrodes (rejE and oldrejE) are not considered yet
   % because there would be a format problem (these rejection are stored in
   % the event array).
   
   rejeegplot = [];
   if ~isempty(rej) %assuming nrejection x 2 (2 = begin end)
      s = size(rej, 1);
      rejeegplot = [ rej(3:4) colnew*ones(s,1)  zeros(s, chans) ];
   end;   

   % pooling previous rejections 
   if ~isempty(oldrej)
       tmp = [ oldrej(:, 3:4) colold*ones(s,1) zeros(size(oldrej,1), chans) ]; 
       rejeegplot = [rejeegplot; tmp];
   end
end

if isempty(varargin)   
	eegplot(data, 'winlength', 5, 'position', [100 300 800 500], 'winrej', rejeegplot, 'xgrid', 'off');
else
	eegplot(data, 'winlength', 5, 'position', [100 300 800 500], 'winrej', rejeegplot, 'xgrid', 'off', varargin{:} );
end;   	 	
return;

% convert eeglab format to eeplot format of rejection window
% ---------------------------------------------------------- 
function rejeegplot = trial2eegplot( rej, rejE, pnts, color)
   	rejeegplot = zeros(length(rej), size(rejE,2)+5);
   	rejeegplot(:, 6:end) = rejE;
   	rejeegplot(:, 1) = (rej(:)-1)*pnts;
   	rejeegplot(:, 2) = rej(:)*pnts-1;
   	rejeegplot(:, 3:5) = ones(size(rejeegplot,1),1)*color;
return
