% userInput_cellArray() - A helper function for HAPPE that collects
%                         elements of a cell array through user input via
%                         the command line. The user enters each element
%                         and uses the newline key (Enter/Return) to
%                         seperate elements. Entering 'done', terminates
%                         the process and duplicates are removed from the
%                         list. Will continue until 'done' is entered,
%                         irrespective of user mistakes.
%
% Usage: 
%   >> cell_inputs = UI_cellArray(indx, cell_inputs)
%
% Inputs:
%  indx -        The starting index to begin entering elements into the
%                cell array.
%  cell_inputs - The cell array to which elements will be added.
%
% Outputs:
%   cell_inputs - A cell array containing each user input as an individual
%                 cell.
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

function cell_inputs = UI_cellArray(indx, cell_inputs)
    while true
        user_input = input('> ', 's') ;
        if strcmpi(user_input, 'done')
            cell_inputs = unique(cell_inputs, 'stable');
            break ;   
        else
            cell_inputs{indx} = user_input ;
            indx = indx + 1 ;
        end
    end
end
