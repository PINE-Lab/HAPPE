% mff_exportcategories - export MFF EEG event from EEGLAB structure to 
%                    'categories.xml' file
% Usage:
%   mff_exportcategories(EEG, mffFile);
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

function indtle = mff_exportcategories(EEG, mffFile)

indtle = []; % time locking event indices
if EEG.trials == 1
    return;
end
mff_path;

% Create a factory.
mfffactorydelegate = javaObject('com.egi.services.mff.api.LocalMFFFactoryDelegate');
mfffactory         = javaObject('com.egi.services.mff.api.MFFFactory', mfffactorydelegate);

categoriesRType = javaObject('com.egi.services.mff.api.MFFResourceType', javaMethod('valueOf', 'com.egi.services.mff.api.MFFResourceType$MFFResourceTypes', 'kMFF_RT_Categories'));
if mfffactory.createResourceAtURI(fullfile(mffFile, 'categories.xml'), categoriesRType)
    fprintf('categories.xml file created successfully\n');
else
    fprintf('categories.xml ressource already exist, overwriting\n');
end
catsResource = mfffactory.openResourceAtURI( fullfile(mffFile, 'categories.xml'), categoriesRType);

% find time-locking events
% ------------------------
events   = EEG.event;
eventlat = abs(eeg_point2lat( [ events.latency ], [ events.epoch ], EEG.srate, [EEG.xmin EEG.xmax]));
indtle    = find(eventlat == 0);
if length(indtle) < EEG.trials
    indtle    = find(eventlat < 0.02);
    if length(indtle) ~= EEG.trials
        if isfield(EEG.event, 'begintime')
            indtle = find(cellfun(@isempty, { EEG.event.begintime }));
            if length(indtle) ~= EEG.trials
                disp('WARNING: NOT THE RIGHT NUMBER OF TIME-LOCKING EVENTS');
            end
        end
    end
end

events = events(indtle);
uniqueType = unique( { events.type } );
jCatList = javaObject('java.util.ArrayList');
epochLen = [];
for iType = 1:length(uniqueType)
    trials = strmatch( uniqueType{iType}, { events.type }, 'exact'); % get time locking events
    
    catObj = javaObject('com.egi.services.mff.api.Category');
    catObj.setName( events(trials(1)).type );
    jList = javaObject('java.util.ArrayList');
    for iTrial = 1:length(trials)
        segmentObj = javaObject('com.egi.services.mff.api.Segment');
        epoch      = events(trials(iTrial)).epoch;
        
        segmentObj.setName( events(trials(iTrial)).type );
        segmentObj.setBeginTime(  round((epoch-1)*EEG.pnts/round(EEG.srate)*1000000) );
        segmentObj.setEndTime(    round((epoch  )*EEG.pnts/round(EEG.srate)*1000000) );
        segmentObj.setEventBegin( round((events(trials(iTrial)).latency-1)/round(EEG.srate)*1000000) );
        if isfield(events(trials(iTrial)), 'duration') && ~isempty(events(trials(iTrial)).duration)
             duration = events(trials(iTrial)).duration;
        else duration = 0;
        end
        segmentObj.setEventEnd(   round((events(trials(iTrial)).latency-1+duration)/round(EEG.srate)*1000000) );
        if isfield(events, 'status')
            segmentObj.setStatus(  events(trials(iTrial)).status );
        end
        if isfield(events, 'mffkeysbackup')
            % if additional keys are present -> average (not sure this heuristic is always correct)
            segmentObj.setName( 'Average' );
            
            if ~isempty(events(trials(iTrial)).mffkeysbackup)
                jListKeys = javaObject('java.util.ArrayList');
                tmpKeys = eval(events(trials(iTrial)).mffkeysbackup);
                for iKey = 1:length(tmpKeys);
                    keyObj = javaObject('com.egi.services.mff.api.Key');
                    keyObj.setCode(tmpKeys(iKey).code);
                    keyObj.setData(tmpKeys(iKey).data);
                    keyObj.setDataType(tmpKeys(iKey).datatype);
                    keyObj.setDescription(tmpKeys(iKey).description);
                    jListKeys.add(keyObj);
                end
                segmentObj.setKeys(jListKeys);
            end
        end
        
        % check epoch length
        epochLenTmp = segmentObj.getEndTime() - segmentObj.getBeginTime();
        if isempty(epochLen)
            epochLen = epochLenTmp;
        elseif epochLen ~= epochLenTmp
            error('Epoch length export error')
        end
        
        jList.add(segmentObj);
    end
    catObj.setSegments(jList);
    jCatList.add(catObj);
end
        
catsResource.setCategories(jCatList);
catsResource.saveResource();
