% PURPOSE: subroutine 
%
% FORMAT:
%
% Seterplabfontsize(fontsizeGUI, fontunitsGUI)
%
% INPUT
%
% fontsizeGUI       - font size value
% fontunitsGUI      - font units ('pixels', 'points')
%
%
% *** This function is part of ERPLAB Toolbox ***
% Author: Javier Lopez-Calderon & Johanna Kreither
% Center for Mind and Brain
% University of California, Davis,
% Davis, CA
% 2012

%b8d3721ed219e65100184c6b95db209bb8d3721ed219e65100184c6b95db209b
%
% ERPLAB Toolbox
% Copyright � 2007 The Regents of the University of California
% Created by Javier Lopez-Calderon and Steven Luck
% Center for Mind and Brain, University of California, Davis,
% javlopez@ucdavis.edu, sjluck@ucdavis.edu
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function Seterplabfontsize(fontsizeGUI, fontunitsGUI)
if nargin<1
        fs = erpworkingmemory('fontsizeGUI');
        fu = erpworkingmemory('fontunitsGUI');
        
        %
        % Call GUI
        %
        answer = setvalueGUI({fs fu});
        
        if isempty(answer)
                disp('User selected Cancel')
                return
        end
        
        fontsizeGUI  = answer{1};
        fontunitsGUI = answer{2};
end
erpworkingmemory('fontsizeGUI', fontsizeGUI);
erpworkingmemory('fontunitsGUI', fontunitsGUI);