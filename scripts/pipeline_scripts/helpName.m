% helpName() - A helper function for HAPPE that creates a save name for
%              HAPPE's output files (quality assessments and/or error log).
%              Searches for a file with an existing name and adds a number
%              incrementally until no file with that name is found.
%
% Usage: 
%   >> saveName = helpName(saveName)
%
% Inputs:
%   saveName - The original save name for the file.
%
% Outputs:
%   saveName - The updated save name for the file.
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

function saveName = helpName(saveName, ext)
saveNameTemp = saveName ;
i = 2 ;
while isfile(saveNameTemp)
    saveNameTemp = strrep(saveName, ext, ['_' num2str(i) ext]) ;
    i = i+1 ;
end
saveName = saveNameTemp ;
end