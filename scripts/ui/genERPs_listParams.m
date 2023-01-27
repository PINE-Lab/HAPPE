% genERPs_listParams - A helper function for HAPPE's generateERPs script
%                      adapted from HAPPE's UI script 'listParams'.
%                      Lists out the parameters to the command window for
%                      user review.
%
% Usage: 
%   >> genERPs_listParams(params)
%
% Inputs:
%   params - A struct containing all the user-specified parameters to be
%            output in the command window.
%
% Outputs:
%
% Author: A.D. Monachino, PINE Lab at Northeastern University, 2022
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

function genERPs_listParams(params)
fprintf(['---------------------------------------------\nPARAMETER ' ...
    'SETTINGS:\n']) ;

%% AVERAGE OR INDIVIDUAL TRIALS
fprintf('Average or Individual Trials: ') ;
if params.indivTrials
    fprintf('Individual\n - Output format: ')
    if params.export.csv; fprintf('Multiple .csv files\n') ;
    else; fprintf('Single Excel file with multiple sheets\n') ;
    end
else
    fprintf('Average\n') ;
end

%% CHANNELS OF INTEREST
fprintf('Channels of Interest: ') ;
if strcmpi(params.chans.subset, 'all'); fprintf('All\n') ;
elseif strcmpi(params.chans.subset, 'coi_include')
    fprintf([sprintf('%s, ', params.chans.IDs{1:end-1}), ...
        params.chans.IDs{end} '\n']) ;
elseif strcmpi(params.chans.subset, 'coi_exclude')
    fprintf('All except ') ;
    fprintf([sprintf('%s, ', params.chans.IDs{1:end-1}), ...
        params.chans.IDs{end} '\n']) ;    
end

%% BAD CHANNEL INCLUSION/EXCLUSION
fprintf('Bad Channels: ') ;
if params.badChans.inc; fprintf('Included\n') ;
else; fprintf(['Excluded\n - File: ' params.badChans.file '\n']) ;
end

%% CALCULATING VALUES
fprintf('Calculating ERP Values: ') ;
if params.calcVals.on
    fprintf('On\n - Windows: ') ;
    for i=1:size(params.calcVals.windows,1)
        fprintf([params.calcVals.windows{i,3} ' in ' params.calcVals.windows{i,1} ...
            '-' params.calcVals.windows{i,2} '\n']);
        if i~=size(params.calcVals.windows,1)
            fprintf('            ') ;
        end
    end
    fprintf(' - Mean Amplitude: Calculate by ') ;
    if params.calcVals.meanAmpMethod(1) && ~params.calcVals.meanAmpMethod(2)
        fprintf('user-specified latency windows\n') ;
    elseif ~params.calcVals.meanAmpMethod(1) && params.calcVals.meanAmpMethod(2)
        fprintf('zero-bound script-generated latency windows\n') ;
    elseif params.calcVals.meanAmpMethod(1) && params.calcVals.meanAmpMethod(2)
        fprintf('both user-specified and zero-bound latency windows\n') ;
    end
    
    fprintf(' - Area Under the Curve: Calculate by ') ;
    if params.calcVals.aucMethod(1) && ~params.calcVals.aucMethod(2)
        fprintf('user-specified latency windows\n') ;
    elseif ~params.calcVals.aucMethod(1) && params.calcVals.aucMethod(2)
        fprintf('zero-bound script-generated latency windows\n') ;
    elseif params.calcVals.aucMethod(1) && params.calcVals.aucMethod(2)
        fprintf('both user-specified and zero-bound latency windows\n') ;
    end
else; fprintf('Off\n') ;
end

fprintf('---------------------------------------------\n') ;
end