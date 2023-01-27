% qm_ERROR() - A helper function for HAPPE fills in rows of quality metrics
%              with "NaN" when HAPPE fails during the processing of any
%              particular file.
%
% Usage: 
%   >> out_metric = qm_ERROR(in_metric, i, current_file)
%
% Inputs:
%   in_metric    - The variable holding the quality metrics.
%   i            - An integer value of the length of the quality metrics.
%   current_file - The index of the current file.
%
% Outputs:
%   out_metric   - The variable holding the quality metrics including the
%                  newly entered NaN values.
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

function out_metric = qm_ERROR(in_metric, i, current_file)
if size(in_metric, 1) < current_file
    out_metric = [in_metric; NaN(1, i)] ;
else
    out_metric = in_metric ;
end
end