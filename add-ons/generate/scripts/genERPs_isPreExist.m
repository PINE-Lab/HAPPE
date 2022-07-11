% genERPs_isPreExist - A helper function for HAPPE's generateERPs script
%                      adapted from HAPPE's UI script 'isPreExist'.
%                      Determines if loading a pre-existing set of
%                      generateERP input parameters. If loading parameters,
%                      import the parameters and ask if changing.
%                      Otherwise, return an empty struct to be filled in by
%                      the user.
%
% Usage: 
%   >> [preExist, params, changeParams] = genERPs_isPreExist()
%
% Inputs:
%
% Outputs:
%   preExist     - A boolean value [0|1] indicating whether the parameters
%                  are pre-existing.
%   params       - Either an empty struct (if no existing params) or a
%                  loaded struct of existing parameter values.
%   changeParams - A boolean value [0|1] indicating whether the loaded
%                  parameter set should be changed.
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

function [preExist, params, changeParams] = genERPs_isPreExist()
% DETERMINE IF LOADING PRE-EXISTING SET USING USER-INPUT VIA COMMAND WINDOW
fprintf('Load pre-existing set of input parameters? [Y/N]\n') ;
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
    
    % PROMPT THE USER ON WHETHER OR NOT TO CHANGE PARAMETERS:
    % First list out the loaded params in the command window for the user
    % to review.
    genERPs_listParams(params) ;
    fprintf('Change an existing parameter? [Y/N]\n') ;
    changeParams = choose2('n', 'y') ;
    
% IF NOT LOADING PRE-EXISTING PARAMETERS...
% Set parameters to an empty struct and indicate that the parameters
% have not been changed.
else; params = struct() ; changeParams = 0 ;
end
end