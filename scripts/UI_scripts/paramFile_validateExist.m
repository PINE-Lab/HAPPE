% paramFile_validateExist() - A helper function for HAPPE that recursively
%                             checks to see if a given parameter file
%                             exists. Allows a user to create a new file or
%                             overwrite the existing file if it already 
%                             exists.
%
% Usage: 
%   >> param_file = paramFile_validateExist(param_file, indx)
%
% Inputs:
%   param_file - A string/character array forming the name of a parameter
%                set a user wants to load.
%   indx       - An integer value used to create a unique name for the
%                parameter set should the user choose the 'new' option.
%
% Outputs:
%   param_file - A string/character array containing a valid name for a
%                parameter file based on user inputs and existing files.
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

function paramFile = paramFile_validateExist(paramFile, defaultName, indx)
    % If the file already exists, prompt the user to overwrite the existing
    % file or save a new one.
    if isfile([pwd filesep paramFile]) 
        disp('A set of input parameters with this name already exists.') ;
        disp('  overwrite = Overwrite the existing file') ;
        disp('  new = Save a new file with a different name') ;
        if choose2('overwrite', 'new')
            % If the user indicates 'new' and has a custom name, try
            % validating the new custom name.
            if indx == 0
                disp('Enter a new name for the file:') ;
                paramFile = paramFile_validateExist(input('> ', 's'), ...
                    defaultName, 0) ;
            % If the user indicates 'new' and has a default name, try
            % validating the next consecutive default name.
            else; paramFile = paramFile_validateExist([defaultName ...
                    datestr(now, 'dd-mm-yyyy') '_' num2str(indx) '.mat'], ...
                    defaultName, indx+1) ;
            end
        end
    else; return ;
    end
end