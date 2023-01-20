%% SET FOLDERS AND PATH
clear ;
fprintf('Preparing HAPPE CheckCriteria_ParticipantInclusion add-on...\n') ;
happeDir = strrep(fileparts(which(mfilename('fullpath'))), '\add-ons', '') ;
addpath([happeDir filesep 'scripts'], ...
    [happeDir filesep 'scripts' filesep 'UI_scripts'], ...
    [happeDir filesep 'add-ons'], ...
    [happeDir filesep 'add-ons' filesep 'check'], ...
    [happeDir filesep 'add-ons' filesep 'scripts']) ;

%% DETERMINE AND SET PATH TO DATA
% Use input from the command line to set the path to the data. If an 
% invalid path is entered, repeat until a valid path is entered.
while true
    srcDir = input('Enter the path to the folder containing the dataset(s):\n> ','s') ;
    if exist(srcDir, 'dir') == 7; break ;
    else; fprintf(['Invalid input: please enter the complete path to the ' ...
            'folder containing the dataset(s).\n']) ;
    end
end
cd (srcDir) ;

%% COLLECT FILES TO RUN
fprintf('Gathering files...\n') ;
FileNames = [{dir('*.csv').name} {dir('*.xlsx').name}] ;
indxDataQC = [];
indxPipelineQC = [];
for i=1:size(FileNames,2)
    if contains(FileNames{i}, 'dataQC'); indxDataQC = [indxDataQC, i];                 %#ok<AGROW> 
    elseif contains(FileNames{i}, 'pipelineQC'); indxPipelineQC = [indxPipelineQC, i]; %#ok<AGROW> 
    end
end
if isempty(indxDataQC) && isempty(indxPipelineQC)
    error('ERROR: NO QC FILES DETECTED.') ;
end

%% DETERMINE IF LOADING PARAMETERS
fprintf('Load a pre-existing set of parameters? [Y/N]\n') ;
preExist = choose2('N','Y') ;
if preExist
    % Ask user for file
    while true
        % Use command window user input to collect the file name of the
        % pre-existing parameters. If the entered file is not an existing
        % file, repeats until a valid input is entered.
        fprintf(['Enter your parameter file, including the full path ' ...
            'and file extension:\n']) ;
        paramFile = input('> ', 's') ;
        if isfile(paramFile); break;
        else; fprintf('Invalid input: please enter the correct file\n') ;
        end
    end
    % Load file
    fprintf('Loading parameters...\n') ;
    load(paramFile) ;
    fprintf('Parameters loaded.\n') ;
    % List parameters
    genInc_listParams(params) ;
    % Ask if want to change parameters
    fprintf('Change an existing parameter? [Y/N]\n') ;
    changedParams = choose2('n', 'y') ;
else
    changedParams = 0 ;
    params.dataQC.single = {'Percent_Good_Chans_Selected', 0, NaN, ''; ...
        'Percent_Var_Retained_Post-Wav', 0, NaN, ''; ...
        'Percent_ICs_Rej', 0, NaN, ''; ...
        'Number_Segs_Post-Seg_Rej', 0, NaN, ''} ;
    params.dataQC.multiple = {'Number_TAG_Segs_Post-Seg_Rej', 0, {}; ...
        'Number_CONDITION_Segs_Post-Seg_Rej', 0, {}} ;
    params.dataQC.badChanROIs = {0, {NaN, NaN, NaN}} ;
    params.pipelineQC.linenoise = {'r FREQ hz pre/post linenoise reduction', 0, NaN, NaN} ;
    params.pipelineQC.ralldata = {'r alldata post/pre waveleting', 0, NaN} ;
    params.pipelineQC.rhz = {'r FREQ hz post/pre waveleting', 0, ...
        {'0.5', NaN; '1', NaN; '2', NaN; '5', NaN; '8', NaN; '12', NaN; ...
        '20', NaN; '30', NaN; '45', NaN; '70', NaN}} ;
end

%-------------------------------------------------------------------------%
%% SET PARAMETERS
while true
    userChoice = '' ;
    %% BREAK IF NOT CHANGING EXISTING PARAMETERS
    if preExist && ~changedParams; break; end

    %% IF CHANGING EXISTING PARAMETERS:
    if changedParams
        fprintf(['Parameter to change: Percent_Good_Chans_Selected,' ...
            ' Percent_Var_Retained_Post-Wav,\nPercent_ICs_Rej, ' ...
            'Number_Segs_Post-Seg_Rej, Number_TAG_Seg_Post-Seg_Rej,\n', ...
            'Number_CONDITION_Segs_Post-Seg_Rej, Bad Channel IDs, ' ...
            'line noise frequency,\nr alldata, r post/pre waveleting by ' ...
            'frequency\n']) ;
        userChoice = input('> ', 's') ;
    end

    %% DATA QC THRESHOLDS:
    % PERCENT GOOD CHANNELS, PERCENT VARIANCE POSTWAV, PERCENT ICS
    % REJECTED, NUMBER OF SEGMENTS RETAINED TOTAL%
    if ~preExist || any(ismember(params.dataQC.single(:,1), userChoice))
        for i=1:size(params.dataQC.single(:,1),1)
            if isempty(userChoice) || strcmpi(params.dataQC.single{i,1}, userChoice)
                fprintf(['Use ' params.dataQC.single{i,1} ' as inclusion' ...
                    ' criteria? [Y/N]\n']) ;
                params.dataQC.single{i,2} = choose2('N', 'Y') ;
                if params.dataQC.single{i,2}
                    fprintf(['Use a minimum or maximum value for exclusion?\n' ...
                        '  min = Any files with a value >= the threshold ' ...
                        'are included.\n  max = Any files with a value ' ...
                        ' < the threshold are included.\n']) ;
                    while true
                        params.dataQC.single{i,4} = input('> ','s') ;
                        if strcmpi(params.dataQC.single{i,4}, 'min') || ...
                            strcmpi(params.dataQC.single{i,4}, 'max')
                            break;
                        else; fprintf(['Invalid input: please enter "min" ' ...
                                'or "max" (without quotations)\n']) ;
                        end
                    end
                    params.dataQC.single{i,3} = str2double(input(['Enter ' ...
                        'minimum value for inclusion:\n> '], 's')) ;
                end
            end
        end
    end

    % NUMBER OF SEGMENTS RETAINED BY ONSET TAG, NUMBER OF SEGMENTS RETAINED
    % BY CONDITION
    if ~preExist || any(ismember(params.dataQC.multiple(:,1), userChoice))
        for i=1:size(params.dataQC.multiple(:,1),1)
            if isempty(userChoice) || strcmpi(params.dataQC.multiple{i,1}, userChoice)
                fprintf(['Use ' params.dataQC.multiple{i,1} ' as inclusion' ...
                    ' criteria? [Y/N]\n']) ;
                params.dataQC.multiple{i,2} = choose2('N','Y') ;
                if params.dataQC.multiple{i,2}
                    fprintf(['Enter your tag/condition name and the minimum ' ...
                        'number of segments\nfor inclusion, seperated by a ' ...
                        'space. Press enter/return between\neach entry. ' ...
                        'When you have finished entering tags/conditions ' ...
                        'and\nthresholds, enter "done" (without quotations)' ...
                        '.\nExample: VEP 15\n         faces 30\n']) ;
                    while true
                        temp = split(input('> ', 's'))' ;
                        if size(temp,2) == 1 && strcmpi(temp, 'done'); break ;
                        elseif size(temp,2) == 2
                            params.dataQC.multiple{i,3} = [params.dataQC.multiple{i,3}; temp] ;                                          
                        else
                            fprintf(['Invalid input: enter a tag/condition' ...
                                ' name and the inclusion threshold or ' ...
                                '"done" (without quotations).\n']) ;
                            continue ;
                        end
                    end
                    clear('temp') ;
                end
            end
        end
    end

    % REGIONS OF INTEREST Regions of interest and the max allowed to be missing per region. Can be
    % a single electrode
    if ~preExist || strcmpi('Bad Channel IDs', userChoice)
        fprintf('Use bad channel IDs as inclusion criteria? [Y/N]\n') ;
        params.dataQC.badChanROIs{1} = choose2('N', 'Y') ;
        if params.dataQC.badChanROIs{1}
            indx = 1 ;
            while true
                fprintf(['Enter the channel(s) to be used to assess ' ...
                    'inclusion, one at a time,\npressing enter/return' ...
                    ' after each entry. When you have finished\nentering ' ...
                    'all channels, input "done" (without quotations).\n']) ;
                params.dataQC.badChanROIs{2}{indx,1} = UI_cellArray(1, {}) ;
                params.dataQC.badChanROIs{2}{indx,2} = input(['Enter the ' ...
                    'maximum number of these channels that can be labeled ' ...
                    'as bad before the file is excluded:\n> ']) ;
                params.dataQC.badChanROIs{2}{indx,3} = input(['Enter a name' ...
                    ' for this ROI:\n> '], 's') ;
                fprintf('Enter another set of channels? [Y/N]\n') ;
                if choose2('Y','N'); break; end
                indx = indx + 1 ;
            end
        end
    else; params.dataQC.badChanROIs{2} = {} ;
    end
    
    %% PIPELINE QC THRESHOLDS:
    % PRIMARY LINE-NOISE FREQUENCY
    if ~preExist || strcmpi('line noise frequency', userChoice)
        fprintf(['Use r pre/post linenoise reduction as exclusion ' ...
            'criteria? [Y/N]\n']) ;
        params.pipelineQC.linenoise{1,2} = choose2('N', 'Y') ;
        if params.pipelineQC.linenoise{1,2}
            params.pipelineQC.linenoise{1,3} = input(['Primary line ' ...
                'noise frequency in Hz:\n> '],'s') ;
            params.pipelineQC.linenoise{1,4} = input(['Enter minimum value ' ...
                'for inclusion:\nExample: 0.5\n> ']) ;
        end
    end

    % R ALLDATA
    if ~preExist || strcmpi('r alldata', userChoice)
        fprintf(['Use r alldata post/pre waveleting as exclusion ' ...
            'criteria? [Y/N]\n'])
        params.pipelineQC.ralldata{1,2} = choose2('N','Y') ;
        if params.pipelineQC.ralldata{1,2}
            params.pipelineQC.ralldata{1,3} = input(['Enter minimum value ' ...
                'for inclusion:\nExample: 0.2\n> ']) ;
        end
    end

    % R FOR SPECIFIED FREQUENCIES
    if ~preExist || strcmpi('r post/pre waveleting by frequency', userChoice)
        % Determine if using this criteria
        fprintf(['Use r pre/post waveleting by frequency as exclusion ' ...
            'criteria? [Y/N]\n']) ;
        params.pipelineQC.rhz{1,2} = choose2('N','Y') ;
        if params.pipelineQC.rhz{1,2}
            % Choose custom set or all. If custom, enter the list of
            % frequencies.
            fprintf(['Use all included frequencies or a subset?\n  all = Use' ...
                ' all frequencies included in the QC sheet.\n  subset = Use ' ...
                'a subset of frequencies included in the QC sheet.\n']) ;
            if choose2('all', 'subset')
                fprintf(['Enter the frequencies of interest, one at a time, ' ...
                    'pressing enter/return between each entry.\nWhen you are ' ...
                    'finished, enter "done" (without quotations).\n']) ;
                freqs = {} ;
                while true
                    ui = input('> ','s') ;
                    if strcmpi(ui, 'done'); break;
                    elseif any(ismember(params.pipelineQC.rhz{1,3}(:,1),ui))
                        freqs = unique([freqs ui], 'stable') ;
                    else; fprintf([ui ' Hz is not a valid frequency.\n']) ;
                    end
                end
                params.pipelineQC.rhz{1,3} = [freqs' cell(size(freqs,2),1)] ;
            end
            % Use multiple values or a single cutoff value?
            fprintf(['Use a single value for exclusion for all frequencies or ' ...
                'specify by frequency?\n  single = Use a single ' ...
                ' value for every frequency\n  specify = Choose a value' ...
                ' for each frequency\n'])
            if choose2('single','specify')
                for i=1:size(params.pipelineQC.rhz{1,3},1)
                    fprintf(['Enter minimum value for inclusion at ' ...
                        params.pipelineQC.rhz{1,3}{i,1} ' Hz:\n']) ;
                    params.pipelineQC.rhz{1,3}{i,2} = str2double(input('> ','s')) ;
                end
            else
                val = input(['Enter minimum value for inclusion:\nExample: ' ...
                    '0.5\n> '], 's') ;
                for i=1:size(params.pipelineQC.rhz{1,3},1)
                    params.pipelineQC.rhz{1,3}{i,2} = str2double(val) ;
                end
            end
        end
    end

    %% CONFIRM PARAMETERS:
    if ~preExist || strcmpi('done', userChoice)
       fprintf('Please check your parameters before continuing.\n') ;
       genInc_listParams(params) ;
       fprintf('Are the above parameters correct? [Y/N]\n') ;
       if choose2('n','y'); break ;
       elseif ~preExist
           changedParams = 1 ;
           preExist = 1 ;
       end
   end
end
%-------------------------------------------------------------------------%
%% SAVE PARAMETERS
if ~preExist || changedParams
    % DETERMINE PARAMETER FILE NAME: Prompt to use a default or custom name
    % for parameter file. If file exists, ask to create new file with a 
    % different name or overwrite existing file.
    fprintf(['Parameter file save name:\n  default = Default name (input' ...
        'genInc_parameters_dd-mm-yyyy.mat).\n  custom = Create a custom file name' ...
        '\n']) ;
    if choose2('custom', 'default')
        paramFile = paramFile_validateExist(['genInc_parameters_' ...
            datestr(now, 'dd-mm-yyyy') '.mat'], 'genInc_parameters_', 2) ;
    else
        fprintf('File name (Do not include .mat):\n') ;
        paramFile = paramFile_validateExist([input('> ', 's') '.mat'], ...
            'inputParameters_', 0) ;
    end

    % SAVE PARAMETERS: Save the params variable to a .mat file using the
    % name created above.
    fprintf('Saving parameters...') ;
    save(paramFile, 'params') ;
    fprintf('Parameters saved.') ;
end


%% RUN ON DATAQC FILES
for currFile = 1:size(indxDataQC,2)
    try
        resultsTab = readtable(FileNames{indxDataQC(currFile)}, 'PreserveVariableNames', true) ;
        % DETERMINE INCLUSION BY PERCENT GOOD CHANNELS, PERCENT VARIANCE,
        % PERCENT ICS REJECTED, AND NUMBER OF SEGMENTS POST-REJECTION
        singles = cell(size(resultsTab,1), sum(cell2mat(params.dataQC.single(:,2)))) ;
        indx = 1;
        for currThresh=1:size(params.dataQC.single(:,1),1)
            if params.dataQC.single{currThresh,2}
                vals = resultsTab.(params.dataQC.single{currThresh,1}) ;
                for i=1:size(vals,1)
                    if (strcmpi(params.dataQC.single{currThresh,4}, 'min') && ...
                            vals(i) > params.dataQC.single{currThresh,3}) || ...
                            (strcmpi(params.dataQC.single{currThresh,4}, 'max') && ...
                            vals(i) < params.dataQC.single{currThresh,3})
                        singles{i, indx} = 'include' ;
                    else; singles{i, indx} = 'exclude' ;
                    end
                end
                indx = indx+1 ;
            end
        end
    
        % DETERMINE INCLUSION BY NUMBER OF SEGMENTS POST-REJECTION BY TAG
        tags = cell(size(resultsTab,1), size(params.dataQC.multiple{1,3}, 1)) ;
        if params.dataQC.multiple{1,2}
            for currThresh = 1:size(params.dataQC.multiple{1,3},1)
                vals = resultsTab.(strrep(params.dataQC.multiple{1,1}, ...
                    'TAG', params.dataQC.multiple{1,3}{currThresh,1})) ;
                for i=1:size(vals,1)
                    if vals(i) > params.dataQC.multiple{1,3}{currThresh,2}
                        tags{i, currThresh} = 'include' ;
                    else; tags{i, currThresh} = 'exclude' ;
                    end
                end
            end
        end
        % DETERMINE INCLUSION BY NUMBER OF SEGMENTS POST-REJECTION BY
        % CONDITION
        conds = cell(size(resultsTab,1), size(params.dataQC.multiple{2,3}, 1)) ;
        if params.dataQC.multiple{2,2}
            for currThresh = 1:size(params.dataQC.multiple{2,3},1)
                vals = resultsTab.(strrep(params.dataQC.multiple{2,1}, ...
                    'CONDITION', params.dataQC.multiple{2,3}{currThresh,1})) ;
                for i=1:size(vals,1)
                    if vals(i) > params.dataQC.multiple{2,3}{currThresh,2}
                        tags{i, currThresh} = 'include' ;
                    else; tags{i, currThresh} = 'exclude' ;
                    end
                end
            end
        end
        
        % DETERMINE INCLUSION BY THE BAD CHANNELS IN A ROI
        badChans = cell(size(resultsTab,1), params.dataQC.badChanROIs{1, ...
            1}*size(params.dataQC.badChanROIs{1,2},1)) ;
        if params.dataQC.badChanROIs{1}
            badChanIDs = table2cell(resultsTab(:,'Bad_Chan_IDs')) ;
            for currThresh=1:size(params.dataQC.badChanROIs{1,2},1)
                for i=1:size(badChanIDs,1)
                    currBadIDs = split(badChanIDs{i}) ;
                    if length(intersect(currBadIDs, params.dataQC.badChanROIs{2}{currThresh,1})) ...
                            >= params.dataQC.badChanROIs{2}{currThresh,2}
                        badChans{i,currThresh} = 'exclude' ;
                    else; badChans{i,currThresh} = 'include' ;
                    end
                end
            end
        end
    
        % DETERMINE INCLUSION ACROSS CRITERIA
        output = [singles, tags, conds, badChans, cell(size(resultsTab,1),1)] ;
        for i=1:size(resultsTab,1)
            if ismember('exclude', output(i,1:end-1)); output{i,end} = 'exclude' ;
            else; output{i,end} = 'include' ;
            end
        end
    
        % OUTPUT A .CSV WITH THE INCLUSION/EXCLUSION INFORMATION
        singlenames = cell(1,size(singles,2)) ;
        indx = 1;
        for i=1:size(params.dataQC.single,1)
            if params.dataQC.single{i,2} && strcmpi(params.dataQC.single{i,4}, 'min')
                singlenames{indx} = [params.dataQC.single{i,1} ' > ' ...
                    num2str(params.dataQC.single{i,3})] ;
                indx = indx+1 ;
            elseif params.dataQC.single{i,2} && strcmpi(params.dataQC.single{i,4}, 'max')
                singlenames{indx} = [params.dataQC.single{i,1} ' < ' ...
                    num2str(params.dataQC.single{i,3})] ;
                indx = indx + 1 ;
            end
        end
        tagnames = cell(1,size(params.dataQC.multiple{1,3},1)) ;
        for i=1:size(tagnames,2)
            tagnames{i} = [strrep(params.dataQC.multiple{1,1}, 'TAG', ...
                params.dataQC.multiple{1,3}{i,1}) ' > ' ...
                num2str(params.dataQC.multiple{1,3}{i,2})] ;
        end
        condnames = cell(1,size(params.dataQC.multiple{2,3},1)) ;
        for i=1:size(condnames,2)
            condnames{i} = [strrep(params.dataQC.multiple{2,1}, 'COND', ...
                params.dataQC.multiple{2,3}{i,1}) ' > ' ...
                num2str(params.dataQC.multiple{2,3}{i,2})];
        end
        roinames = cell(1,size(params.dataQC.badChanROIs{1,2},1)) ;
        if params.dataQC.badChanROIs{1}
                for i=1:size(roinames,2)
                    roinames{i} = [params.dataQC.badChanROIs{1,2}{i,3} ': ' ...
                        sprintf('%s, ', params.dataQC.badChanROIs{1,2}{i,1}{1:end-1}) ...
                        params.dataQC.badChanROIs{1,2}{i,1}{end} ' > ' ...
                        num2str(params.dataQC.badChanROIs{1,2}{i,2})] ;
                end
        end
        tabnames = [singlenames, tagnames, condnames, roinames, 'Across_Criteria'] ;
        savename = replace(FileNames{indxDataQC(currFile)}, {'.csv','.xlsx'}, '_ParticipantInclusion.csv') ;
        indx = 2 ;
        while isfile(savename)
            savename = strrep(savename, '.csv', ['_' num2str(indx) '.csv']) ;
            indx = indx+1 ;
        end
        writetable(cell2table(output, 'VariableNames', tabnames, 'RowNames', ...
            table2cell(resultsTab(:,'Row'))), savename, 'WriteRowNames', ...
            true, 'QuoteStrings', true) ;
    catch ME
        fprintf(['Unable to process ' FileNames{indxDataQC(currFile)} '\n']) ;
    end
end

%% RUN ON PIPELINEQC FILES
for currFile = 1:size(indxPipelineQC,2)
    try
        resultsTab = readtable(FileNames{indxPipelineQC(currFile)}, ...
            'PreserveVariableNames', true) ;
    
        linenoise = cell(size(resultsTab,1),params.pipelineQC.linenoise{1,2}) ;
        if params.pipelineQC.linenoise{1,2}
            vals = resultsTab.(strrep(params.pipelineQC.linenoise{1}, 'FREQ', ...
                params.pipelineQC.linenoise{3})) ;
            for i=1:size(vals,1)
                if vals(i) > params.pipelineQC.linenoise{4}; linenoise{i} = 'include' ;
                else; linenoise{i} = 'exclude' ;
                end
            end
        end
    
        ralldata = cell(size(resultsTab,1),params.pipelineQC.ralldata{1,2}) ;
        if params.pipelineQC.ralldata{1,2}
            vals = resultsTab.(params.pipelineQC.ralldata{1}) ;
            for i=1:size(vals,1)
                if vals(i) > params.pipelineQC.ralldata{3}; ralldata{i} = 'include' ;
                else; ralldata{i} = 'exclude' ;
                end
            end
        end
    
        rhz = cell(size(resultsTab,1), params.pipelineQC.rhz{2}*size(params.pipelineQC.rhz{3},1)) ;
        if params.pipelineQC.rhz{2}
            for currThresh=1:size(params.pipelineQC.rhz{3},1)
                vals = resultsTab.(strrep(params.pipelineQC.rhz{1}, 'FREQ', ...
                    params.pipelineQC.rhz{3}{currThresh,1})) ;
                for i=1:size(vals,1)
                    if vals(i) > params.pipelineQC.rhz{1,3}{currThresh,2}
                        rhz{i, currThresh} = 'include' ;
                    else; rhz{i, currThresh} = 'exclude' ;
                    end
                end
            end
        end
    
        % DETERMINE INCLUSION ACROSS CRITERIA
        output = [linenoise, ralldata, rhz, cell(size(resultsTab,1),1)] ;
        for i=1:size(resultsTab,1)
            if ismember('exclude', output(i,1:end-1)); output{i,end} = 'exclude' ;
            else; output{i,end} = 'include' ;
            end
        end
    
        % OUTPUT A .CSV WITH THE INCLUSION/EXCLUSION INFORMATION
        if params.pipelineQC.linenoise{2}
            lnnames = [strrep(params.pipelineQC.linenoise{1,1}, 'FREQ', ...
                params.pipelineQC.linenoise{1,3}) ' > ' ...
                num2str(params.pipelineQC.linenoise{4})] ;
        else; lnnames = {} ;
        end
        if params.pipelineQC.ralldata{2}
            rallnames = [params.pipelineQC.ralldata{1} ' > ' ...
                num2str(params.pipelineQC.ralldata{3})] ;
        else; rallnames = {} ;
        end
        if params.pipelineQC.rhz{2}
            rhznames = cell(1, size(params.pipelineQC.rhz{3},1)) ;
            for i=1:size(rhznames,2)
                rhznames{i} = [strrep(params.pipelineQC.rhz{1}, 'FREQ', ...
                    params.pipelineQC.rhz{1,3}{i,1}) ' > ' ...
                    num2str(params.pipelineQC.rhz{1,3}{i,2})] ;
            end
        else; rhznames = {} ;
        end
        tabnames = [lnnames, rallnames, rhznames, 'Across_Criteria'] ;
        
        savename = replace(FileNames{indxPipelineQC(currFile)}, {'.csv','.xlsx'}, '_ParticipantInclusion.csv') ;
        indx = 2 ;
        while isfile(savename)
            savename = strrep(savename, '.csv', ['_' num2str(indx) '.csv']) ;
            indx = indx+1 ;
        end
        writetable(cell2table(output, 'VariableNames', tabnames, 'RowNames', ...
            table2cell(resultsTab(:,'Row'))), savename, 'WriteRowNames', ...
            true, 'QuoteStrings', true) ;
    catch ME
        fprintf(['Unable to process ' FileNames{indxPipelineQC(currFile)} '\n']) ;
    end
end

%% SUPPORT FUNCTIONS
function genInc_listParams(params)
    for i=1:size(params.dataQC.single,1)
        fprintf([params.dataQC.single{i,1} ': ']) ;
        if params.dataQC.single{i,2}
            fprintf([params.dataQC.single{i,4} ' of ' ...
                num2str(params.dataQC.single{i,3}) '\n']) ;
        else; fprintf('NA\n') ;
        end
    end

    for i=1:size(params.dataQC.multiple,1)
        thresh = params.dataQC.multiple{i,3} ;
        if i==1; rep = 'TAG'; else; rep = 'CONDITION'; end
        for j=1:size(thresh,1)
            fprintf([strrep(params.dataQC.multiple{i,1}, rep, thresh{j,1}) ...
                ': min of ' num2str(thresh{j,2}) '\n']) ;
        end
    end
    
    fprintf('Bad Channels: ') ;
    if params.dataQC.badChanROIs{1}
        thresh = params.dataQC.badChanROIs{2} ;
        fprintf('\n') ;
        for i=1:size(thresh,1)
            fprintf([' - ' thresh{i,3} ': max of ' num2str(thresh{i,2}) ...
                ' channels from ' sprintf('%s, ', thresh{i,1}{1:end-1}) ...
                thresh{i,1}{end} '\n']) ;
        end
    else; fprintf('NA\n') ;
    end

    if params.pipelineQC.linenoise{2}
        fprintf([strrep(params.pipelineQC.linenoise{1}, 'FREQ', ...
            params.pipelineQC.linenoise{3}) ': min of ' ...
            num2str(params.pipelineQC.linenoise{4}) '\n']) ;
    end

    if params.pipelineQC.ralldata{1}
        fprintf([params.pipelineQC.ralldata{1} ': min of ' ...
            num2str(params.pipelineQC.ralldata{3}) '\n']) ;
    end

    if params.pipelineQC.rhz{2}
        for i=1:size(params.pipelineQC.rhz{3},1)
            fprintf([strrep(params.pipelineQC.rhz{1}, 'FREQ', ...
                params.pipelineQC.rhz{3}{i,1}) ': min of ' ...
                num2str(params.pipelineQC.rhz{3}{i,2}) '\n']) ;
        end
    end
end