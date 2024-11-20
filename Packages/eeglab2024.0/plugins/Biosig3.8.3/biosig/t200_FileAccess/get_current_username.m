function [username] = get_current_username()
% GET_CURRENT_USERNAME returns the username of the user running
%    the current Matlab or Octave instance. This is used in biosig
%    when storing processed data in GDF format as the "Technician"
%    information.  

% Copyright (C) 2021 by Alois Schloegl <alois.schloegl@gmail.com>
%    This is part of the BIOSIG-toolbox https://biosig.sourceforge.io/
%
%    BioSig is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    BioSig is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with BioSig.  If not, see <http://www.gnu.org/licenses/>.


username = 'unknown';
if ispc()
	[status,username] = dos('echo %USERNAME%');
elseif ismac()
	[status,username] = system('id -F');
elseif isunix()
	[status,username] = system('id -un');
elseif exist('OCTAVE_VERSION','builtin') && ~ispc()
	t = getpwuid(getuid);
	if isfield(t,'gecos')
		username = strtok(t.gecos,',');
	end
end

