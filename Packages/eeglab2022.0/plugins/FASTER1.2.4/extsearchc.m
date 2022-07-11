function [paths names] = extsearchc(startDir,extension,norecurse)

% Copyright (C) 2010 Hugh Nolan, Robert Whelan and Richard Reilly, Trinity College Dublin,
% Ireland
% nolanhu@tcd.ie, robert.whelan@tcd.ie
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    if (exist('norecurse','var')~=1)
        norecurse=0;
    end
    [paths names] = dirsearchc(startDir,extension,norecurse);
    
    function [paths names] = dirsearchc(currDir,searchstring,norecurse)
        if nargin < 2
            fprintf('Usage: [paths names] = dirsearch(currDir, searchstring)\n');
            paths = '';
            names = '';
            return;
        elseif (exist('norecurse','var')~=1)
            norecurse=0;
        end
        paths = {};
        names = {};
        list_currDir = dir(currDir);

        for u = 1:length(list_currDir)
            if (list_currDir(u).isdir==1 && strcmp(list_currDir(u).name,'.')~=1 && strcmp(list_currDir(u).name,'..')~=1 && norecurse==0)
                [temppaths tempnames] = dirsearchc(sprintf('%s%s%s',currDir,filesep,list_currDir(u).name),searchstring);
                paths = {paths{:} temppaths{:}};
                names = {names{:} tempnames{:}};
            elseif (length(list_currDir(u).name) > 4)
                extension = list_currDir(u).name(end-3 : end);
                    if strcmp(extension, searchstring) == 1
                        paths = {paths{:} currDir};
                        names = {names{:} list_currDir(u).name};
                    end
            end
        end

	end
end