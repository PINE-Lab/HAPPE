% choose2() - A helper function for HAPPE, allowing the user to make a
%             choice between two options via command-line input without
%             throwing an error or incorrectly accepting invalid input.
%             Is not case-sensitive.
%
% Usage: 
%   >> answ = choose2(c1, c2)
%
% Inputs:
%   c1   - The first user option in the form of a string/char array. This
%          option, when selected, is coded as a 0.
%   c2   - The second user option in the form of a string/char array. This
%          option, when selected, is coded as a 1.
%
% Outputs:
%   answ - A binary value of 0 or 1, reflecting the user's selection of c1
%          or c2, respectively.
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

function answ = choose2(c1, c2)
while true
    answ = input('> ', 's') ;
    if strcmpi(answ, c1); answ = 0 ; break ;
    elseif strcmpi(answ, c2); answ = 1 ; break ;
    else; fprintf("Invalid input: please enter %s or %s\n", c1, c2) ;
    end
end
end