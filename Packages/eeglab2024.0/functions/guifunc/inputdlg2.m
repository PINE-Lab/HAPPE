% INPUTDLG2 - inputdlg function clone with coloring and help for 
%               EEGLAB.
%
% Usage:
%   >> Answer = inputdlg2(Prompt,Title,LineNo,DefAns,funcname);
% 
% Inputs:
%   Same as inputdlg. Using the optional additional funcname parameter 
%   the function will create a help button. The help message will be
%   displayed using the POPHELP function.
%
% Output:
%   Same as inputdlg
%
% Note: The advantage of this function is that the color of the window
%       can be changed and that it displays an help button. Edit 
%       supergui to change window options. Also the parameter LineNo
%       can only be one.
%
% Author: Arnaud Delorme, CNL / Salk Institute, La Jolla, 11 August 2002
%
% See also: SUPERGUI, INPUTGUI

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, arno@salk.edu
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

function [result] = inputdlg2(Prompt,Title,LineNo,DefAns,funcname)

if nargin < 2
   help inputdlg2;
   return;
end
if nargin < 3
   LineNo = 1;
end
if nargin < 4
   DefAns = {};
end
if nargin < 2
   help inputdlg2;
   return;
end
if nargin < 5
	funcname = '';
end
	
if ~iscell(Prompt), Prompt = { Prompt }; end
if isempty(DefAns)
    DefAns = cell(1,length(Prompt));
    DefAns(:) = { '' }; 
end

if length(Prompt) ~= length(DefAns)
	error('inputdlg2: prompt and default answer cell array must have the same size');
end

geometry = {};
listgui = {};

% determine if vertical or horizontal
% -----------------------------------
geomvert = [];
for index = 1:length(Prompt)
	geomvert = [geomvert size(Prompt{index},1) 1];  % default is vertical geometry
end
if all(geomvert == 1) && length(Prompt) > 1
	geomvert = []; % horizontal
end

for index = 1:length(Prompt)
	if ~isempty(geomvert) % vertical
		geometry = { geometry{:} [ 1] [1 ]};
	else
		geometry = { geometry{:} [ 1 0.6 ]};
	end
	listgui = { listgui{:} { 'Style', 'text', 'string', Prompt{index}}  ...
				{ 'Style', 'edit', 'string', DefAns{index} } };
end

result = inputgui(geometry, listgui, ['pophelp(''' funcname ''');'], Title, [], 'normal', geomvert);
