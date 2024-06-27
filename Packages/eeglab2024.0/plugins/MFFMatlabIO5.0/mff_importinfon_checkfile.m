% mff_importinfon_checkfile - support function for MFF_IMPORTINFON
%                             that checks a data file and fix it if
%                             necessary so it can be imported with the JAVA
%                             librairy
%
% Usage:
%   mff_exportsignal(mffFile);
%%
% Output:
%  File is read and overwritten if necessary. The old file is saved under
%  xxxxx_backup

% This file is part of mffmatlabio.
%
% mffmatlabio is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% mffmatlabio is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with mffmatlabio.  If not, see <https://www.gnu.org/licenses/>.

function mff_importinfon_checkfile(mffFile)

% Check if file exists
if ~exist(mffFile,'file')
    error('File %s does not exist',mffFile);
end

% make a copy of the file
try
    filebackup = [mffFile '_backup'];
    copyfile(mffFile, filebackup);
catch
    error('Could not make a backup of file %s; check that the folder is writable',mffFile);
end

% Open as text file and make a copy line by line
fid = fopen(filebackup,'r');
fid2 = fopen(mffFile,'w');
if fid == -1 || fid2 == -1
    error('Could not open file %s',mffFile);
end
while ~feof(fid)
    line = fgetl(fid);
    if contains(line,'<ch n="')
        line = strrep(line,'<ch n="','<ch n=');
        line = strrep(line,'">"','>');
    end
    fprintf(fid2,'%s\n',line);
end
fclose(fid);
fclose(fid2);

