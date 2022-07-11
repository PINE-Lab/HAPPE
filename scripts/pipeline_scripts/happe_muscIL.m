% happe_muscIL() - A helper script for HAPPE that runs ICA on the data
%                  and uses the ICLabel (Pion-Tonachini et al., 2019)
%                  EEGLAB (Delorme & Makeig, 2004) plugin to reject
%                  artifacts with at least 25% probability of being muscle. 
%
% Usage: 
%   >> [outEEG, dataQC] = happe_muscIL(inEEG, dataQC, FileNames, currFile, inputExt)
%
% Inputs:
%   inEEG     - The EEG, in EEGLAB (Delorme & Makeig, 2004) format.
%   dataQC    - The data quality assessment matrix for updated values
%               regarding the number and percent of ICs rejected.
%   FileNames - A list of the raw file names. Used to save intermediate
%               outputs.
%   currFile  - The index of the current file in the list of files. Used to
%               properly input data QC values and to save intermediate
%               outputs.
%   inputExt  - A user-specified character array or empty character used to
%               properly save intermediate outputs.
%
% Outputs:
%   outEEG    - The EEG, in EEGLAB (Delorme & Makeig, 2004) format, after 
%               rejecting ICs.
%   dataQC    - The data quality assessment matrix updated with the number
%               and percent of ICs rejected.
%
% Authors: A.D. Monachino, PINE Lab at Northeastern University, 2022
%
% This file is part of HAPPE.
% Copyright 2018, 2021 Alexa Monachino, Kelsie Lopez, Laurel Gabard-Durnam
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

function [outEEG, dataQC] = happe_muscIL(inEEG, dataQC, FileNames, ...
    currFile, inputExt)
% RUN ICA: Use EEGLAB's (Delorme & Makeig, 2004) pop_runica function to 
% generate the ICA decomposition for the current file. The data will be 
% filtered at 1Hz for optimal ICA performance.
EEG_ICA = pop_runica(pop_eegfiltnew(inEEG, 1, [], [], 0, [], 0), ...
    'extended', 1, 'interupt', 'off') ;

% SAVE ICA DECOMPOSITION
pop_saveset(EEG_ICA, 'filename', strrep(FileNames{currFile}, inputExt, ...
    '_ICAdecomp.set')) ;

% LABEL AND FLAG COMPONENTS: Using ICLabel (Pion-Tonachini et al., 2019),
% first classify the independent components. Then, flag all components
% with at least 25% probability of being muscle artifact for rejection.
EEG_ICA = pop_iclabel(EEG_ICA, 'default') ;
EEG_ICA = pop_icflag(EEG_ICA, [NaN NaN; 0.25 1; NaN NaN; NaN NaN; ...
    NaN NaN; NaN NaN; NaN NaN]);
artICs = find(EEG_ICA.reject.gcompreject == 1) ;

% APPLY ICA WEIGHTS: Add the ICA weights from the decomposed data to the
% non-ICA-decomposed data. 
inEEG.icaweights = EEG_ICA.icaweights ;
inEEG.icasphere = EEG_ICA.icasphere ;
inEEG.icawinv = EEG_ICA.icawinv ;
inEEG.icaact = EEG_ICA.icaact ;
% Confirm that the weights, etc., transfered correctly.
inEEG = eeg_checkset(inEEG) ;

% REJECT ICs: If there are artifacts to reject, reject them from the data.
if ~isempty(artICs); outEEG = pop_subcomp(inEEG, artICs, 0) ;
elseif length(artICs) == length(EEG_ICA.reject.gcompreject)
    fprintf(2, ['All independent components rejected. No further processing' ...
        ' possible.']) ;
    error('HAPPE:allICsRej', ['All independent components rejected. No ' ...
        'further processing possible.']) ;
else; outEEG = inEEG ;
end

% ADD QC METRICS: Add the number and percent ICs rejected to the data
% quality assessment matrix.
dataQC{currFile, 7} = length(artICs) ;
dataQC{currFile, 8} = length(artICs)/length(EEG_ICA.reject.gcompreject)*100 ;

% SAVE DATA POST-IC REJECTION
pop_saveset(outEEG, 'filename', strrep(FileNames{currFile}, inputExt, ...
    '_muscILcleaned.set')) ;
end