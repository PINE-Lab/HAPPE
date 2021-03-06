% mff_exportinfo - export MFF 'info.xml' file from EEGLAB structure
%
% Usage:
%   mff_exportinfo(EEG, mffFile);
%
% Inputs:
%  EEG     - EEGLAB structure
%  mffFile - filename/foldername for the MFF file (MFF file/folder must
%            already exist)

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

function EEG = mff_exportinfo(EEG, mffFile)

mff_path;
mfffactorydelegate = javaObject('com.egi.services.mff.api.LocalMFFFactoryDelegate');
mfffactory         = javaObject('com.egi.services.mff.api.MFFFactory', mfffactorydelegate);

infoType = javaObject('com.egi.services.mff.api.MFFResourceType', javaMethod('valueOf', 'com.egi.services.mff.api.MFFResourceType$MFFResourceTypes', 'kMFF_RT_Info'));
if mfffactory.createResourceAtURI(fullfile(mffFile, 'info.xml'), infoType)
    fprintf('Info.xml file created successfully\n');
else
    fprintf('Info.xml ressource already exist, overwriting\n');
end
info = mfffactory.openResourceAtURI( fullfile(mffFile, 'info.xml'), infoType);

info.setMFFVersion(3);
if ~isfield(EEG.etc, 'recordingtime')
    error('Unknow record time');
end
begTime = mff_encodetime(EEG.etc.recordingtime, EEG.etc.timezone);
info.setRecordTime(begTime);
info.saveResource();
