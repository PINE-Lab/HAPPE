% isReprocessed() - A helper script for HAPPE. Determines whether the
%                   user is using a pre-set group of parameters.
%
% Usage: 
%   >> [preExist, params, changeParams] = isPreExist(reprocessing, ver)
%
% Inputs:
%   reprocessing - A binary number [0|1] indicating whether or not data is
%                  being reprocessed. 
%   ver          - A string representing the current version of HAPPE
%
% Outputs:
%   preExist     - A binary number [0|1] indicating whether or not the
%                  parameters are pre-existing.
%   params       - Either the previously created parameters, or an empty
%                  struct to be filled in later.
%   changeParams - A binary number [0|1] indicating whether or not the
%                  parameters should be changed.
%
% Author: A.D. Monachino, PINE Lab at Northeastern University, 2021
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

function [preExist, params, changeParams] = isPreExist(reprocessing, ver)
% DETERMINE IF LOADING PRE-EXISTING SET USING USER-INPUT VIA COMMAND WINDOW
fprintf('Load pre-existing set of input parameters [Y/N]\n') ;
preExist = choose2('n','y') ;

% IF LOADING PRE-EXISTING PARAMETERS...
if preExist
    while true
        % Use command window user input to collect the file name of the
        % pre-existing parameters. If the entered file is not an existing
        % file, repeats until a valid input is entered.
        fprintf(['Enter your parameter file, including the full path ' ...
            'and file extension:\n']) ;
        paramFile = input('> ', 's') ;
        if isfile(paramFile); break;
        else; fprintf('Invalid input: please enter the correct file\n') ;
        end
    end
    % LOAD THE PARAMETER FILE
    fprintf('Loading parameters...\n') ;
    load(paramFile) ;
    fprintf('Parameters loaded.\n') ;
    
    % VALIDATE THE PARAMETER SET
    % If the parameters were created in a different version of HAPPE, ask
    % the user if they would like to continue anyway. If the user indicates
    % not to continue, end this run via error.
    if isfield(params, 'HAPPEver') && ~strcmpi(params.HAPPEver, ver)
        fprintf(['These parameters were saved using a different version' ...
            ' of HAPPE\nand may be incompatible with the running pipeline.\n' ...
            'Proceed anyway? [Y/N]\n']) ;
        if ~choose2('n', 'y'); error(['Mismatch between HAPPE and input ' ...
                'parameter versions']) ; end
    % If there is no version associated with HAPPE, end the run via error
    % as it is unlikely that the script will be able to run without valid
    % parameters.
    elseif ~isfield(params, 'HAPPEver'); error(['Unable to identify the ' ...
            'HAPPE version of this parameter set.']) ;
    end
    
    % PROMPT THE USER ON WHETHER OR NOT TO CHANGE PARAMETERS:
    % First list out the loaded params in the command window for the user
    % to review.
    listParams(params) ;
    fprintf('Change an existing parameter? [Y/N]\n') ;
    changeParams = choose2('n', 'y') ;
    
% IF NOT LOADING PRE-EXISTING PARAMETERS...
else
    % IF REPROCESSING DATA, prompt the user to enter the correct set of
    % parameters for the data.
    if reprocessing; fprintf(['PLEASE ENSURE THE PARAMETERS (UP TO ' ...
            'SEGMENTATION) THAT YOU ENTER MATCH YOUR PREVIOUS RUN.\n']); end
    % Set parameters to an empty struct and indicate that the parameters
    % have not been changed.
    params = struct() ;
    changeParams = 0 ;
end
end