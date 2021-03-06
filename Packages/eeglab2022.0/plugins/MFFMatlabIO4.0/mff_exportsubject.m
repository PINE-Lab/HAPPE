% mff_exportsubject - export subject information into MFF ressources
%
% Usage:
%   mff_exportsubject(EEG, mffFile);
%
% Inputs:
%  EEG     - EEGLAB structure
%  mffFile - filename/foldername for the MFF file

% This file is part of mffmatlabio.
%
% mffmatlabio is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% mffmatlabio is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with mffmatlabio.  If not, see <https://www.gnu.org/licenses/>.

function mff_exportsubject(EEG, mffFile)

if ~isfield(EEG.etc, 'subject') || isempty(EEG.etc.subject)
    return;
end

mff_path;

% Create an MFFFactory object.
mfffactorydelegate = javaObject('com.egi.services.mff.api.LocalMFFFactoryDelegate');
mfffactory = javaObject('com.egi.services.mff.api.MFFFactory', mfffactorydelegate);

% Create Signal object and read in event track file.
sURI = fullfile(mffFile, 'subject.xml');
objectType = javaObject('com.egi.services.mff.api.MFFResourceType', javaMethod('valueOf', 'com.egi.services.mff.api.MFFResourceType$MFFResourceTypes', 'kMFF_RT_Subject'));
mfffactory.createResourceAtURI(sURI, objectType);
sObject = mfffactory.openResourceAtURI(sURI, objectType);

variables = { 'Fields' 'array' { 'Name' 'char' {}; 'Data' 'char' {}; 'DataType' 'char' {} } };

sObject = mff_setobj(sObject, EEG.etc.subject, variables);
sObject.saveResource();
