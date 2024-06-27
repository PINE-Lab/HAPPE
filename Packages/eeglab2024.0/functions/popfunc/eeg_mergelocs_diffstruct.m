% EEG_MERGELOCS - merge channel structure while preserving channel
%                   order
%
%      >> mergedlocs = eeg_mergelocs(loc1, loc2, loc3, ...);
%
% Inputs: 
%     loc1     - EEGLAB channel location structure
%     loc2     - second EEGLAB channel location structure
%
% Output: 
%     mergedlocs - merged channel location structure
%
% Author: Arnaud Delorme, August 2006

% Copyright (C) Arnaud Delorme, CERCO, 2006, arno@salk.edu
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function [alllocs warn] = eeg_mergelocs(varargin)

persistent warning_shown;
warn = 0;

% sort by length
% --------------
len = cellfun(@length, varargin);
[tmp so] = sort(len, 2, 'descend');
varargin = varargin(so);

alllocs = varargin{1};
for index = 2:length(varargin)
    
    % fuse while preserving order (assumes the same channel order)
    % ------------------------------------------------------------
    tmplocs = varargin{index};
    newlocs = myunion(alllocs, tmplocs);
    
    if length(newlocs) > length(union({ alllocs.labels }, { tmplocs.labels }))
        
        warn = 1;
        if isempty(warning_shown)
            disp('Warning: different channel montage order for the different datasets');
            warning_shown = 1;
        end
        
        % trying to preserve order of the longest array
        %----------------------------------------------
        if length(alllocs) < length(tmplocs)
            tmp     = alllocs;
            alllocs = tmplocs;
            tmplocs = tmp;
        end
        allchans = { alllocs.labels tmplocs.labels };
        [uniquechan ord1 ord2 ]  = unique_bc( allchans );
        
        [tmp rminds] = intersect_bc( uniquechan, { alllocs.labels });
        ord1(rminds) = [];
        tmplocsind = ord1-length(alllocs);
        
        newlocs = concatlocs(alllocs, tmplocs(tmplocsind));

    end
    alllocs = newlocs;
end

% union of two channel location structure
% without losing the order information
% ---------------------------------------
function alllocs = myunion(locs1, locs2)

   labs1 = { locs1.labels };
   labs2 = { locs2.labels };
   
   count1 = 1;
   count2 = 1;
   count3 = 1;
   alllocs = locs1; alllocs(:) = [];
   while count1 <= length(locs1) | count2 <= length(locs2)
       
       if count1 > length(locs1)
           alllocs = copyfields(alllocs, count3, locs2(count2));
           %alllocs(count3) = locs2(count2);
           count2 = count2 + 1;
           count3 = count3 + 1;
       elseif count2 > length(locs2)
           alllocs = copyfields(alllocs, count3, locs1(count1));
           %alllocs(count3) = locs1(count1);
           count1 = count1 + 1;
           count3 = count3 + 1;
       elseif strcmpi(labs1{count1}, labs2{count2})
           alllocs = copyfields(alllocs, count3, locs1(count1));
           %alllocs(count3) = locs1(count1);
           count1 = count1 + 1;
           count2 = count2 + 1;
           count3 = count3 + 1;
       elseif isempty(strmatch(labs1{count1}, labs2, 'exact'))
           alllocs = copyfields(alllocs, count3, locs1(count1));
           %alllocs(count3) = locs1(count1);
           count1 = count1 + 1;
           count3 = count3 + 1;
       else 
           alllocs = copyfields(alllocs, count3, locs2(count2));
           %alllocs(count3) = locs2(count2);
           count2 = count2 + 1;
           count3 = count3 + 1;
       end
       
   end
    
% concatenate channel structures with different fields   
function loc3 = concatlocs(loc1, loc2);

    fields1 = fieldnames(loc1);
    fields2 = fieldnames(loc2);
    if isequal(fields1, fields2)
        % the try, catch clause is necessary 
        % below seems to be a Matlab bug
        try, loc3 = [ loc1; loc2 ];
        catch
            try loc3 = [ loc1 loc2 ];
            catch, loc3 = [ loc1(:)' loc2(:)' ];
            end
        end
    else
        loc3 = loc1;
        for index = 1:length(loc2)
            loc3 = copyfields(loc3, length(loc1)+index, loc2(index));
        end
    end

% copy fields of a structure to another one
function [struct1] = copyfields(struct1, index1, struct2)

    fields1 = fieldnames(struct1);
    fields2 = fieldnames(struct2);

    if isequal(fields1, fields2)
        struct1(index1) = struct2;
    else
        for index = 1:length(fields2)
            struct1 = setfield(struct1, {index1}, fields2{index}, getfield(struct2, fields2{index}));
        end
    end
    
