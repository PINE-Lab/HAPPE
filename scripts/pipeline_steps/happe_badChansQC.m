% happe_badChansQC() - A helper script for HAPPE that assesses the
%                      performance of bad channel rejection. Metrics
%                      include the number of channels included in
%                      processing (good channels), the percent of good
%                      channels remaining from the original user-selected
%                      set, and the list of channels (if any) marked as
%                      bad.
%
% Usage: 
%   >> dataQC = happe_badChansQC(EEG, dataQC, chanIDs, currFile)
%
% Inputs:
%   EEG      - The EEG, in EEGLAB (Delorme & Makeig, 2004) format, with bad
%              channels removed.
%   dataQC   - The matrix holding the data quality assessment metrics.
%   chanIDs  - The matrix holding the channels of interest for the current 
%              subject.
%   currFile - The index of the current file in the list of file names.
%              Used in updating the dataQC variable.
%
% Outputs:
%   dataQC   - The matrix containing the data quality assessment metrics
%              updated with the metrics for the current file.
%
% Author: A.D. Monachino, PINE Lab at Northeastern University, 2021
%
% This file is part of HAPPE.
% Copyright 2018, 2023 Alexa Monachino, Kelsie Lopez, Laurel Gabard-Durnam
%
% HAPPE is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
% 
% HAPPE is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
% details.
% 
% You should have received a copy of the GNU General Public License along
% with HAPPE. If not, see <https://www.gnu.org/licenses/>.

function dataQC = happe_badChansQC(EEG, dataQC, chanIDs, currFile)
fprintf('Evaluating bad channel detection...\n') ;
dataQC{currFile,3} = size(EEG.chanlocs,2) ;
dataQC{currFile, 4} = dataQC{currFile,3}/dataQC{currFile,2}*100 ;

if isempty(EEG.chanlocs); dataQC{currFile,5} = 'NA' ;
else
    badChans = setdiff(chanIDs, {EEG.chanlocs.labels}, 'stable') ;
    if isempty(badChans) && isempty(dataQC{currFile,5})
        dataQC{currFile,5} = 'None' ;
    elseif ~isempty(badChans)
        dataQC{currFile,5} = [sprintf('%s ', badChans{1:end-1}), ...
            badChans{end}] ;
    end
end
end