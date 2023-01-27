% happe_interpChanBySeg() - An adapted version of FASTER code (Nolan et 
%                           al., 2010). For more information, refer to the
%                           relevant paper.
%
% Usage: 
%   >> [EEG, dataQC] = happe_interpChanBySeg(EEG, dataQC, currFile)
%
% Inputs:
%   EEG      - The EEG, in EEGLAB (Delorme & Makeig, 2004) format, prior to
%              filtering.
%   dataQC   - The matrix holding the data quality assessment metrics.
%   currFile - The index of the current file in the file names list.
%
% Outputs:
%   EEG      - The EEG, in EEGLAB (Delorme & Makeig, 2004) format, after 
%              filtering.
%   dataQC   - The matrix containing the data quality assessment metrics
%              updated with the metrics for the current file.
%
% Orignial Authors: Nolan, et al., 2010
% Copyright 2010 Hugh Nolan, Robert Whelan and Richard Reilly, Trinity 
% College Dublin, Ireland
%
% Adapted for HAPPE by: L.J. Gabard-Durnam, PINE Lab at Northeastern
%                       University, 2021
%                       A.D. Monachino, PINE Lab at Northeastern 
%                       University, 2021
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

function [EEG, dataQC] = happe_interpChanBySeg(EEG, dataQC, currFile, params)
fprintf('Interpolating bad data...\n') ;
if params.loadInfo.chanlocs.inc && params.reref.on && params.reref.flat
    EEGtemp = pop_select(EEG, 'channel', setdiff({EEG.chanlocs.labels}, ...
        params.reref.chan)) ;
else; EEGtemp = EEG ;
end

eegChans = 1:size(EEGtemp.chanlocs, 2) ;
chanNames = {EEGtemp.chanlocs.labels} ;
rejOps.measure = [1 1 1 1] ;
rejOps.z = [3 3 3 3] ;
if length(size(EEGtemp.data)) > 2
    status = '' ;
    lengthsEp = cell(1, size(EEGtemp.data, 3)) ;
    for v=1:size(EEGtemp.data, 3)
        listProps = single_epoch_channel_properties(EEGtemp, v, eegChans);
        lengthsEp{v} = eegChans(logical(min_z(listProps, rejOps)));
        badChanNames = chanNames(lengthsEp{v}) ;
        status = [status sprintf('[%d:',v) sprintf(' %s', ...
            badChanNames{:}) ']'] ; %#ok<AGROW> 
    end
    EEG = h_epoch_interp_spl(EEG, lengthsEp, EEG.nbchan) ;
    EEG.saved = 'no' ;

    % Add info about which channels were interpolated for each 
    % segment to the dataQM output.
    EEG.etc.epoch_interp_info = status ;
    dataQC{currFile,9} = cellstr(status) ;
end
end