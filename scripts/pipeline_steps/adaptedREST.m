% adapedREST - An adapted version of pop_REST_reref.m from the REST (Yao,
%              2001) EEGLAB (Delorme & Makieg, 2004) plugin. Strips the
%              necessity of using a GUI and reallocates some input
%              restrictions (e.g., overlapping channels in left and right
%              lists) to when the parameters are originally set in HAPPE's
%              setParams script. Otherwise uses code from the original
%              authors present in the plugin (provided in HAPPE's instance
%              of EEGLAB).
%
% Usage:
%   >> EEG = adaptedREST(EEG, cfg)
%
% Inputs:
%   EEG - The EEG struct, in EEGLAB format, prior to re-referencing using
%         REST.
%   cfg - The parameters needed to properly re-reference the data.
%
% Outputs:
%   EEG - The EEG struct, in EEGLAB format, after re-referencing using
%         REST.
%
% Original Authors: Li Dong and Shiang Hu, University of Electronic Science
%                   and Technology of China, 2019
% Copyright 2019 Key Laboratory for NeuroInformation of Ministry of 
% Education, School of Life Science and Technology, University of 
% Electronic Science and Technology of China
%
% Adapted for HAPPE by: A.D. Monachino, PINE Lab at Northeastern
%                       University, 2022
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

function EEG = adaptedREST(EEG, cfg)
%% SET LEADFIELD:
try
    % LOAD LEADFIELD
    if ~cfg.calc
        fprintf('----------------------------\nLoading Lead Field...\n') ;
        G = table2array(readtable(cfg.file)) ;
%         [~, ~, ext] = fileparts(cfg.file) ;
%         switch ext
%             case '.xlsx'; G = xlsread(cfg.file);
%             case '.xls'; G = xlsread(cfg.file);
%             case '.dat'; G = load(cfg.file);
%             case '.txt'; G = load(cfg.file,'-ascii');
%         end
        if sum(isnan(G(:)))>0
            error('NaN is contained in lead field matrix.');
        end
    
    % CALCULATE LEADFIELD
    else
        fprintf(['Calculating leadfield based on 3-concentric spheres' ...
            ' headmodel at once...\n']);
        % For Average/One Fixed Reference
        if strcmpi(cfg.ref, 'average') || strcmpi(cfg.ref, 'fixed')
            G = dong_getleadfield(EEG, cfg.chanlist);
            fprintf(['Lead Field Matrix: ' num2str(size(G,1)), ...
                ' sources X ' num2str(size(G,2)) ' channels\n']);
        
        % For IM/CM
        elseif strcmpi(cfg.ref, 'mastoid')
            % Check channel list for left and right channels
            cfg.chanlist = [cfg.chanL cfg.chanR] ;
            % Check XYZ coordinates of references for left and right channels
            if ~isempty(cfg.coordLeft) && all(isfinite(cfg.coordLeft))
                fprintf(['XYZ of Ref. for left channels: ' ...
                    num2str(cfg.coordLeft) '\n']);
            else
                error(['XYZ coordinates of reference for left channels' ...
                    ' are missing or contian an infinite value.']);
            end
            if ~isempty(cfg.coordRight) && all(isfinite(cfg.coordRight))
                fprintf(['XYZ of Ref. for right channels: ' ...
                    num2str(cfg.XYZ_rightRef) '\n']);
            else
                error(['XYZ coordinates of reference for right channels' ...
                    'are missing or contain an infinite value.']);
            end
            % Calculate leadfield
            G = dong_getleadfield(EEG, cfg.chanlist, [cfg.coordLeft; ...
                cfg.coordRight]) ;
            fprintf(['Lead Field Matrix: ' num2str(size(G,1)) ...
                ' sources X (' num2str(size(G,2)-2) ' channels + 2 ' ...
                'references)\n']);
        end
    end
catch ME
    rethrow(ME) ;
end

%% RE-REFERENCE:
try
    % SET ORIGINAL DATA, RESHAPING IF NECESSARY
    if ~isempty(cfg.chanlist)
        if ndims(EEG.data) == 3
            OrigData = EEG.data(cfg.chanlist, :);
            fprintf(['EEG data is in 3D format. Reshaping to channels x ' ...
                'timepoints...\n']) ;
        elseif ndims(EEG.data) == 2                                         %#ok<ISMAT> 
            OrigData = EEG.data(cfg.SelectChans.chanlist,:);
        end
    else; error('ERROR: No channels selected for re-referencing.') ;
    end
    fprintf(['EEG data: ' num2str(size(OrigData,1)) ' channels x ' ...
        num2str(size(OrigData,2)) ' time points.\n']) ;

    % REST ON AVERAGE OR FIXED REFERENCE DATA:
    if size(OrigData,1) == size(G,2) && (strcmpi(cfg.ref, 'average') || ...
            strcmpi(cfg.ref, 'fixed'))
        fprintf('Original reference is average or a fixed electrode...\n');
        OrigData = OrigData - repmat(mean(OrigData), size(OrigData,1),1);
%         switch cfg.OrigReferFlag
%             case 1 % average
%                 OrigData = OrigData - repmat(mean(OrigData),size(OrigData,1),1);
%             case 2 % re-refer to average
%                 OrigData = OrigData - repmat(mean(OrigData),size(OrigData,1),1);
%         end           NOTE: these appear to be the same, why split?

        % Refer to REST
        disp('Re-referencing to REST...');
        REST_EEG =  dong_rest_refer(OrigData, G) ;

    % REST ON MASTOID REFERENCE DATA:
    elseif size(OrigData,1) + 2 == size(G,2) && strcmpi(cfg.ref, 'mastoid')
        fprintf(['Original reference is ipsilateral mastoid or ' ...
            'contralateral mastoid...\n']);

        % Refer to REST
        disp('Re-referencing to REST...');
        REST_EEG = dong_rest_rerefer_cm(OrigData, G, ...
            1:length(cfg.chanL), (length(cfg.chanL)+1):(length(cfg.chanL)+ ...
            length(cfg.chanR)), [size(G,2)-1, size(G,2)]);
    else
        error(['Number of channels of the leadfield matrix and data are ' ...
            'NOT equal']) ;
    end

    % RESHAPE BACK TO 3D FORMAT IF NECESSARY:
    if ndims(EEG.data) == 3
        REST_EEG = reshape(REST_EEG, size(REST_EEG,1), size(EEG.data,2), ...
            size(EEG.data,3));
        fprintf('Reshaping to channels X timepoints X epochs...\n');
    end
    EEG.data = REST_EEG;       % EEG data refered to REST
    fprintf('Completed.\n') ;
catch ME
    rethrow(ME) ;
end
end