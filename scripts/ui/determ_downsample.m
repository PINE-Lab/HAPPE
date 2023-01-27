% determ_downsample() - A helper script for setParams.m from HAPPE.
%                       Determines whether the user intends to downsample
%                       the data, and if so, to what frequency, through 
%                       user input via the command line.
%
% Usage: 
%   >> freq = determ_downsample()
%
% Inputs:
%
% Outputs:
%   freq - An integer representing the frequency to downsample the data to.
%          If the frequency is 0, that indicates that downsampling should 
%          not be performed.
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

function freq = determ_downsample()
fprintf('Resample data? [Y/N]\n') ;
if choose2("n", "y")
    fprintf(['HAPPE supports resampling to 250, 500, and 1000 Hz.\nResample' ...
        ' frequency:\n']) ;
    while true
        freq = input('> ') ;
        if ismember(freq, [250, 500, 1000]); break ;
        else; disp("Invalid input: please enter 250, 500, or 1000.") ;
        end
    end
else; freq = 0 ;
end
end