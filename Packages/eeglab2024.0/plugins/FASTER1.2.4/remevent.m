function EEG = remevent(EEG,numbers)

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

    if ~iscell(numbers)
        str_numbers=cell(size(numbers));
        for u=1:length(numbers)
            str_numbers{u} = num2str(numbers(u));
        end
    else
        str_numbers = numbers;
    end
    
    if ischar(EEG.event(1).type)
        types={EEG.event.type};

        good_indices=true(size(EEG.event));
        for u=1:length(numbers)
            good_indices=good_indices & (~(strcmp(str_numbers{u},types)));
        end
        good_indices=find(good_indices);
    else
        types=[EEG.event.type];
        
        good_indices=true(size(EEG.event));
        for u=1:length(numbers)
            good_indices=good_indices & (types~=numbers(u));
        end
    end
    
    EEG.event = EEG.event(good_indices);