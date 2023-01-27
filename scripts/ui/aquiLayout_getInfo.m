% aquiLayout_getInfo() - Using the user-specified layout and HAPPE's
%                        netdata library, return the net info for that 
%                        layout from the netdata library.
%
% Usage: 
%   >> net_info = aquiLayout_getInfo(layout_type, netdata_lib)
%
% Inputs:
%   layout_type - Aquisition layout in the format [#,#], specifying the net
%                 type and the number of channels. Output from
%                 determ_aquiLayout.m.
%   netdata_lib - HAPPE's netdata library containing the relevant
%                 information about HAPPE's supported nets.
%
% Outputs:
%   net_info    - The struct containing the information about the selected 
%                 aquisition layout from HAPPE's netdata library
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

function net_info = aquiLayout_getInfo(layout_type, netdata_lib)
% GSN Nets:
if layout_type(1,1) == 1
    if layout_type(1,2) == 64; net_info = netdata_lib.GSN64 ;
    end
elseif layout_type(1,1) == 2
    if layout_type(1,2) == 32; net_info = netdata_lib.GSNHC32 ;
    elseif layout_type(1,2) == 64; net_info = netdata_lib.GSN64 ;
    elseif layout_type(1,2) == 128; net_info = netdata_lib.GSNHC128 ;
    elseif layout_type(1,2) == 256; net_info = netdata_lib.GSNHC256 ;
    else; disp("Invalid input: this layout type is not supported by this version of HAPPE.") ;
    end

% BioSemi Nets:
elseif layout_type(1,1) == 2
    if layout_type(1,2) == 16; net_info = netdata_lib.Biosemi16 ;
    elseif layout_type(1,2) == 32; net_info = netdata_lib.Biosemi32 ;
    elseif layout_type(1,2) == 64; net_info = netdata_lib.Biosemi64 ;
    elseif layout_type(1,2) == 128; net_info = netdata_lib.Biosemi128 ;
    else; disp("Invalid input: this layout type is not supported by this version of HAPPE.") ;
    end

% BrainProducts' Standard BrainCap:
elseif layout_type(1,1) == 3
    if layout_type(1,2) == 22; net_info = netdata_lib.BC22 ;
    elseif layout_type(1,2) == 32; net_info = netdata_lib.BC32 ;
    elseif layout_type(1,2) == 64; net_info = netdata_lib.BC64 ;
    elseif layout_type(1,2) == 96; net_info = netdata_lib.BC96 ;
    elseif layout_type(1,2) == 128; net_info = netdata_lib.BC128 ;
    else; disp("Invalid input: this layout type is not supported by this version of HAPPE.") ;
    end

% BrainProducts' Wet-Sponge R-Net for actiCHamp Plus:
elseif layout_type(1,1) == 4
    if layout_type(1,2) == 32; net_info = netdata_lib.RNPAC32 ;
    elseif layout_type(1,2) == 64; net_info = netdata_lib.RNPAC64 ;
    elseif layout_type(1,2) == 96; net_info = netdata_lib.RNPAC96 ;
    elseif layout_type(1,2) == 128; net_info = netdata_lib.RNPAC128 ;
    else; disp("Invalid input: this layout type is not supported by this version of HAPPE.") ;
    end

% Other Nets:
else
    error("Invalid input: this layout type is not supported by this version of HAPPE.") ;
end
end