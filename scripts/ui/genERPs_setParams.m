% genERPs_isPreExist - A helper function for HAPPE's generateERPs script
%                      adapted from HAPPE's UI script 'isPreExist'.
%                      Determines if loading a pre-existing set of
%                      generateERP input parameters. If loading parameters,
%                      import the parameters and ask if changing.
%                      Otherwise, return an empty struct to be filled in by
%                      the user.
%
% Usage: 
%   >> params = genERPs_setParams(params, preExist, changedParams)
%
% Inputs:
%   params        - Either an empty struct or a loaded struct containing
%                   pre-existing parameters.
%   preExist      - A boolean value [0|1] indicating whether or not the
%                   parameters are pre-existing.
%   changedParams - A boolean value [0|1] indicating whether to change
%                   existing parameters.
%
% Outputs:
%   params        - A struct containing all the user-specified parameters.
%
% Author: A.D. Monachino, PINE Lab at Northeastern University, 2022
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

function params = genERPs_setParams(params, preExist, changedParams)

paramChoice = 'na' ;
while true
    %% BREAK IF NOT CHANGING PRE-EXISTING PARAMS
    if preExist && ~changedParams; break ; end
    
    %% IF CHANGING PARAMS...
    if changedParams
        fprintf(['Parameter to change: ave/indiv trials, channels of ' ...
            'interest, bad channel inclusion, calculating values, plotting.' ...
            '\nEnter "done" (without quotations) when finished changing ' ...
            'parameters.\n']) ;
        paramChoice = input('> ', 's') ;
    end
    
    %% DETERMINE IF RUNNING AVERAGE OR INDIVIDUAL TRIALS
    if ~preExist || strcmpi(paramChoice, 'ave/indiv trials')
        fprintf(['Trial type:\n  average = Average over trials\n  ' ...
            'individual = Individual trials\n']) ;
        params.indivTrials = choose2('average', 'individual') ;
        
        % SINGLE OR MULTIPLE FILES
        if params.indivTrials
            fprintf(['Choose your export format:\n  sheets = A single Excel file' ...
                ' with multiple sheets\n  files = Multiple .csv files\n']) ;
            params.export.csv = choose2('sheets', 'files') ;  
        else; params.export.csv = 1 ;
        end
    end 
    
    %% DETERMINE CHANNELS OF INTEREST
    % Ask the user if they want to include entered channels or exclude
    % entered channels.
    if ~preExist || strcmpi(paramChoice, 'channels of interest')
        [params.chans.subset, params.chans.IDs] = determ_chanIDs() ;
    end
    
    %% DETERMINE IF INCLUDING/EXCLUDING BAD CHANNELS
    % Ask the user whether or not to include bad/interpolated channels in
    % the analyses.
    if ~preExist || strcmpi(paramChoice, 'bad channel inclusion')
        fprintf(['Include bad channels in calculating ERP?\n  include = Keep bad ' ...
            'channels\n  exclude = Remove bad channels\n']) ;
        params.badChans.inc = choose2('exclude', 'include') ;
        
        % DETERMINE FILE CONTAINING THE BAD CHANNELS
        if ~params.badChans.inc
            while true
                fprintf(['Enter the file containing the bad channels, including ' ...
                    'the complete path.\nRefer to the HAPPE User Guide for ' ...
                    'instructions on creating this file and an example.\n']) ;
                params.badChans.file = input('> ', 's') ;
                if isfile(params.badChans.file); break ;
                else; fprintf('Invalid input: please enter an existing file.') ;
                end
            end
        end
    end
    
    %% DETERMINE IF CALCULATING VALUES
    % Determine, using user input from the command line, whether calculating
    % ERP values. If so, collect the latency windows, and the method(s) for
    % calculating area under the curve/50% area under the curve.
    if ~preExist || strcmpi(paramChoice, 'calculating values')
        fprintf('Calculate ERP values? [Y/N]\n') ;
        params.calcVals.on = choose2('n', 'y') ;
        if params.calcVals.on
            % COLLECT LATENCY WINDOWS
            params.calcVals.windows = [] ;
            fprintf(['Enter latency windows of interest with anticipated peak:\n' ...
                'Enter each window as two consecutive numbers followed by "max" or ' ...
                '"min" (without quotations).\nPress Enter/Return between entries.\n' ...
                'When you have entered all windows, input "done" (without quotations).' ...
                '\nExample: 10 100 max\n']) ;
            while true
                temp = split(input('> ', 's')) ;
                if length(temp) == 1 && strcmpi(temp, 'done'); break ;
                elseif size(temp, 1) == 3; temp = reshape(temp, 1, 3) ;
                else
                    fprintf(['Invalid input: please enter two numbers and "max"/"min"' ...
                        ' or "done" (without quotations).\n']) ;
                    continue ;
                end
                if str2num(temp{1}) > str2num(temp{2})
                    fprintf(['Invalid input: please make sure that the ' ...
                        'entries are consecutive.\n']) ;
                elseif str2num(temp{1}) == str2num(temp{2})
                    fprintf(['Invalid input: you cannot have a latency window of 0.\n']) ;
                elseif ~strcmpi(temp{3}, 'max') && ~strcmpi(temp{3}, 'min')
                    fprintf(['Invalid input: please specify "max" or "min" ' ...
                        '(without quotations).\n']) ;
                else
                    params.calcVals.windows = [params.calcVals.windows; temp] ;
                end
            end

            % DETERMINE METHOD FOR MEAN AMPLITUDE
            params.calcVals.meanAmpMethod = calcValMethods('mean amplitude') ;

            % DETERMINE METHOD FOR AUC/50% AUC
            params.calcVals.aucMethod = calcValMethods('area under the curve') ;
            
            if params.indivTrials
                % DETERMINE OUTPUT FORMAT FOR INDIVIDUAL TRIALS
                fprintf(['Choose an option:\n  subjects - Rows as trials, columns as ' ...
                    'values split by subject\n  values - Rows as subjects, columns as trials ' ...
                    'split by value.\nNOTE: See HAPPE User Guide for examples.\n']) ;
                while true
                    ui = input('> ', 's') ;
                    if strcmpi(ui, 'subjects'); params.export.catg = 'subject' ; break ;
                    elseif strcmpi(ui, 'values'); params.export.catg = 'value' ; break ;
                    else; fprintf('Invalid input: please enter "subjects" or "values".\n') ;
                    end
                end
            end
        end
    end
    
    %% DETERMINE IF PLOTTING
    if ~preExist || strcmpi(paramChoice,'plotting')
        fprintf('Plot the ERP waveforms? [Y/N]\n') ;
        params.plot = choose2('n','y') ;
    end

    %% DONE
   if ~preExist || strcmpi(paramChoice, 'done')
       fprintf('Please check your parameters before continuing.\n') ;
       genERPs_listParams(params) ;
       fprintf('Are the above parameters correct? [Y/N]\n') ;
       if choose2('n','y'); break ;
       elseif ~preExist; changedParams = 1 ; preExist = 1 ;
       end
   end 
end
end

function method = calcValMethods(value)
fprintf(['Choose a method for calculating ' value ':\n' ...
    '  windows = Restrict calculations to the user-specified latency window(s)\n' ...
    '  zeros = Calculate using points where the amplitude is 0\n' ...
    '  both = Calculate both by windows and by zeros\n']) ;
method = [0,0] ;
while true
    ui = input('> ', 's') ;
    if strcmpi(ui, 'windows'); method(1) = 1; break;
    elseif strcmpi(ui, 'zeros'); method(2) = 1; break ;
    elseif strcmpi(ui, 'both'); method = [1,1] ; break ;
    else; fprintf(['Invalid input: please enter "windows", "zeros", ' ...
            'or "both" (without quotations).\n']) ;
    end
end
end