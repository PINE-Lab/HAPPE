% pop_readedf() - Read a European data format .EDF data file. 
%
% Usage:
%   >> EEG = pop_readedf;             % an interactive window pops up
%   >> EEG = pop_readedf( filename ); % no pop-up window 
%
% Inputs:
%   filename       - European data format 16-bit EDF file 
% 
% Outputs:
%   EEG            - EEGLAB data structure
%
% Author: Arnaud Delorme, CNL / Salk Institute, 13 March 2002
%
% Notes: pop_readedf() uses the functions sdfopen() and sdfread(). Use
% the alternative function readedf() from the command line in 
% case of problem (note that this function does not recalibrate thte data).
%
% See also: sdfopen(), sdfread(), readedf()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 13 March 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% $Log: not supported by cvs2svn $
% Revision 1.20  2004/08/31 21:30:15  arno
% same
%
% Revision 1.19  2004/08/31 21:28:41  arno
% adapt to new BIOSIG
%
% Revision 1.18  2003/12/18 00:27:03  arno
% reading channel labels
%
% Revision 1.17  2003/12/16 17:00:37  arno
% remove warnings
%
% Revision 1.16  2003/06/05 15:38:18  arno
% only reading EDF format
%
% Revision 1.15  2003/04/10 18:06:43  arno
% default output argument
%
% Revision 1.14  2003/04/10 18:01:08  arno
% file filter
%
% Revision 1.13  2003/02/18 02:44:41  arno
% adding message to extract events
%
% Revision 1.12  2002/12/30 18:05:08  arno
% see also field updated
%
% Revision 1.11  2002/12/30 17:06:43  arno
% debugging for sdfread compatibility
%
% Revision 1.10  2002/12/27 18:06:20  scott
% header msg edit -sm
%
% Revision 1.9  2002/12/27 17:22:42  arno
% updating with new Alois version
%
% Revision 1.8  2002/11/14 23:42:02  arno
% reading both file formats EDF and BDF
%
% Revision 1.7  2002/11/12 01:30:53  arno
% update header
%
% Revision 1.6  2002/11/12 01:25:16  arno
% back to original JR
%
% Revision 1.5  2002/11/12 01:04:45  arno
% use new function of A. S.
%
% Revision 1.4  2002/10/15 17:04:06  arno
% drawnow
%
% Revision 1.3  2002/07/24 01:19:23  arno
% addind message
% ,
%
% Revision 1.2  2002/07/24 00:49:01  arno
% debugging
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

function [EEG, command] = pop_readedf(filename); 
EEG = [];
command = '';

if nargin < 1 
	% ask user
	[filename, filepath] = uigetfile('*.EDF;*.edf', 'Choose an EDF file -- pop_readedf()'); 
    drawnow;
	if filename == 0 return; end;
	filename = [filepath filename];
end;

% load datas
% ----------
EEG = eeg_emptyset;
fprintf('Reading EDF format using BIOSIG...\n');
EDF = sopen(filename, 'r');
[tmpdata EDF] = sread(EDF, Inf); tmpdata = tmpdata';
EEG.data = tmpdata;

%[EEG.data, header]  = readedf(filename);  
EEG.nbchan          = size(EEG.data,1);
EEG.pnts            = size(EEG.data,2);
EEG.trials          = 1;
EEG.srate           = EDF.SampleRate(1);
EEG.setname 		= 'EDF file';
disp('Event information might be encoded in the last channel');
disp('To extract these events, use menu File > Import event info > From data channel'); 
EEG.filename        = filename;
EEG.filepath        = '';
EEG.xmin            = 0; 
EEG.chanlocs        = struct('labels', cellstr(EDF.Label));

EEG = eeg_checkset(EEG);
command = sprintf('EEG = pop_readedf(''%s'');', filename); 

return;
