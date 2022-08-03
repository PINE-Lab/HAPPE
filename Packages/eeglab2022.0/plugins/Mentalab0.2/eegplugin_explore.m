% eegplugin_explore() - is the top-level Explore EEGLab plug-in function
%
% Usage:
%   >>  eegplugin_explore(fig,try_strings,catch_strings);
%
% Inputs:
%   fig - the handle of the main EEGLAB window
%    
% See also: 
%   pop_sample, eeglab 

% Copyright (C) 2021 Mentalab GmbH
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
function vers = eegplugin_explore(fig, try_strings, catch_strings)
vers = '0.2';

importmenu = findobj(fig, 'tag', 'import data');
uimenu(importmenu, 'label', 'Import Mentalab Explore data', 'callback', ...
    ['[EEG ORN LASTCOM] = pop_loadfile;' ...
    '  [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);' ......
    '  [ALLEEG ORN] = eeg_store(ALLEEG, ORN, (CURRENTSET + 1));' ...
    ' eeglab redraw']);
