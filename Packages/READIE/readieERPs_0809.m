% READIE ERPs - an adapted version of HAPPE+ER's (Monachino et al., 2022)
%               generateERPs script. Creates the ERP timeseries and
%               calculates common measures that are output in a format that
%               can be used with READIE (TBD).
%
% Developed at Northeastern University's PINE Lab and the BEAD Lab at the
% University of Southern California
%
% Author: A.D. Monachino, PINE Lab at Northeastern University, 2022
% Adapted for READIE by: A.D. Monachino, BEAD Lab at the University of
% Southern California, 2024

%% SET FOLDERS AND PATH
% Use input from the command line to set the path to the data. If an 
% invalid path is entered, repeat until a valid path is entered.
clear ;
fprintf('Preparing READIE_erps...\n') ;
scriptDir = fileparts(which(mfilename('fullpath'))) ;
addpath([scriptDir filesep 'support']) ;

%% DETERMINE AND SET PATH TO THE DATA
% Use input from the command line to set the path to the data. If an 
% invalid path is entered, repeat until a valid path is entered.
while true
    srcDir = input('Enter the path to the folder containing the processed dataset(s):\n> ','s') ;
    if exist(srcDir, 'dir') == 7; break ;
    else; disp("Invalid input: please enter the complete path to the folder containing the dataset(s).") ;
    end
end
cd(srcDir) ;

%% CREATE OUTPUT FOLDERS
% Create the folders in which to store outputs
fprintf('Creating output folders...\n') ;
if ~isfolder([srcDir filesep 'readieERPs']); mkdir([srcDir filesep 'readieERPs']) ;
end
addpath('readieERPs') ;
cd([srcDir filesep 'readieERPs']) ;
if ~isfolder([srcDir filesep 'readieERPs' filesep 'ERP_timeseries'])
    mkdir([srcDir filesep 'readieERPs' filesep 'ERP_timeseries']) ;
end
if ~isfolder([srcDir filesep 'readieERPs' filesep 'ERP_calculatedVals'])
    mkdir([srcDir filesep 'readieERPs' filesep 'ERP_calculatedVals']) ;
end
fprintf('Output folders created.\n') ;

%% SET PARAMETERS
% DETERMINE IF USING PRE-EXISTING SET OF PARAMETERS
[preExist, params, changedParams] = readieERPs_isPreExist() ;

% SET PARAMETERS THROUGH USER INPUT
params = readieERPs_setParams(params, preExist, changedParams) ;

% SAVE INPUT PARAMETERS
% If created a new or changed a parameter set, save as a new .mat file
if ~preExist || changedParams
    % Prompt to use a default or custom name for parameter the file. 
    % If the file exists, ask to create new file with a different name
    % or overwrite existing file.
    fprintf(['Parameter file save name:\n  default = Default name (readieERPs' ...
        '_parameters_dd-mm-yyyy.mat)\n  custom = Create a custom file name' ...
        '\n']) ;
    if choose2('custom', 'default')
        paramFile = readieERPs_validateExist(['readieERPs_parameters_' ...
            datestr(now, 'dd-mm-yyyy') '.mat'], 'readieERPs_parameters_', 2) ;
    else
        fprintf('File name (Do not include .mat):\n') ;
        paramFile = readieERPs_validateExist([input('> ', 's') '.mat'], ...
            'readieERPs_parameters_', 0) ;
    end

    % Save the params variable to a .mat file using the name created above.
    fprintf('Saving parameters...\n') ;
    save(paramFile, 'params') ;
    fprintf('Parameters saved.\n') ;
end
clear('preExist', 'changedParams', 'paramFile') ;


%% COLLECT FILES
cd(srcDir) ;
if str2double(params.format) == 1; ext = '.txt' ;
elseif str2double(params.format) == 2; ext = '.set' ;
end
FileNames = {dir(['*' ext]).name} ;
if size(FileNames,2) < 1; error('ERROR: No files detected') ; end

%% READ IN BAD CHANNELS
if ~params.badChans.inc
    % Load the File and convert table to cell
    badChans = table2cell(readtable(params.badChans.file)) ;

    % Split the string into a cell array
    for currtrial=1:size(badChans, 1)
        badChans{currtrial,2} = split(badChans{currtrial,2})' ;
    end
end

%% INITIALIZE OUTPUT TABLES
allSubsERP = cell(1,size(FileNames,2)+1) ;
MeanAmp_wind = cell(1,size(FileNames,2)+1) ;    % Mean Amplitude by User-Defined Windows
MeanAmp_zero = cell(1,size(FileNames,2)+1) ;    % Mean Amplitude by Zero-Bound Windows

%% PER FILE...
for currfile=1:size(FileNames, 2)+1
    % If processing a file...
    if currfile <= size(FileNames,2)
        try
            % LOAD FILE
            if strcmpi(params.format, '1'); currsub = importdata(FileNames{currfile}) ;
            elseif strcmpi(params.format, '2')
                currsub = pop_loadset('filename', FileNames{currfile}, ...
                    'filepath', srcDir) ;
            end
    
            % COMPILE LIST OF BAD CHANNELS
            if ~params.badChans.inc
                subBadChans = badChans{contains(badChans(:,1), ... 
                    strrep(FileNames{currfile}, ext, '')), 2} ;
            end
    
            % SET CHANNELS OF INTEREST
            if strcmpi(params.chans.subset, 'coi_exclude')
                if strcmpi(params.format, '1')
                    subChanIDs = setdiff(currsub.colheaders, [params.chans.IDs, 'Time']) ;
                elseif strcmpi(params.format, '2')
                    subChanIDs = setdiff({currsub.chanlocs.labels}, params.chans.IDs) ;
                end
                if ~params.badChans.inc; subChanIDs = setdiff(subChanIDs, subBadChans) ; 
                end
            elseif strcmpi(params.chans.subset, 'coi_include')
                subChanIDs = params.chans.IDs ;
                if ~params.badChans.inc; subChanIDs = setdiff(subChanIDs, ...
                        intersect(subBadChans, params.chans.IDs)) ; 
                end
            elseif strcmpi(params.chans.subset, 'all')
                if strcmpi(params.format, '1')
                    subChanIDs = setdiff(currsub.colheaders, 'Time') ;
                elseif strcmpi(params.format, '2')
                    subChanIDs = {currsub.chanlocs.labels} ;
                end
                if ~params.badChans.inc; subChanIDs = setdiff(subChanIDs, subBadChans) ; 
                end
            end
    
            % GET COLUMN NUMBERS
            chanCols = [] ;
            if strcmpi(params.format, '1')
                if ~isfield(currsub, 'colheaders')
                    currsub.colheaders = strsplit(currsub.textdata{1}) ;
                end
                for int=1:length(subChanIDs)
                    chanCols = [chanCols find(ismember(currsub.colheaders, ...
                        subChanIDs{int}))] ;
                end
            elseif strcmpi(params.format, '2')
                for int=1:length(subChanIDs)
                    chanCols = [chanCols find(ismember({currsub.chanlocs.labels}, ...
                        subChanIDs{int}))] ;
                end
            end
    
            % SPLIT BY TRIALS
            % Find the range of lats in the data
            if strcmpi(params.format, '1')
                lats = unique(currsub.data(:,1)) ;
                trialBounds = [find(currsub.data(:,1)==min(lats)), ...
                        find(currsub.data(:,1)==max(lats))] ;
            elseif strcmpi(params.format, '2')
                lats = (currsub.xmin:1/currsub.srate:currsub.xmax)'*1000 ;
            end
        catch e
            fprintf(['Failed to process ' FileNames{currfile} '...\n']) ;
            continue ;
        end
    % If processing the average of files...
    else
        temp = [] ;
        for i = 1:size(FileNames,2)
            if strcmpi(params.format, '1')
                if size(allSubsERP{i},2) > 3 
                    temp = [temp allSubsERP{i}(:,end-3)] ;
                else; temp = [temp allSubsERP{i}(:,end)] ;
                end
            elseif strcmpi(params.format, '2')
                temp = [temp allSubsERP{i}(:,end-3)] ;
            end
        end
        allSubsERP{end} = mean(temp, 2, 'omitnan') ;
        if strcmpi(params.format, '1')
            trialBounds = [1 size(allSubsERP{end},1)] ;
        elseif strcmpi(params.format, '2')
            currsub.trials = 1 ;
        end
    end

    % SET NUMBER OF TRIALS
    if strcmpi(params.format, '1')
        numTrials = size(trialBounds, 1) ;
    elseif strcmpi(params.format, '2')
        numTrials = currsub.trials ;
    end

    % CORRECT WINDOW LATENCIES
    currsub_winLats = params.calcVals.windows ;
    for i=1:size(params.calcVals.windows,1)
        for j=1:2
            if ~any(ismember(lats, ...
                    str2num(params.calcVals.windows{i,j}))) %#ok<*ST2NM> 
                [minVal, closestIndx] = min(abs(lats- ...
                    str2num(params.calcVals.windows{i,j}))) ;
                currsub_winLats{i,j} = ...
                    num2str(lats(closestIndx)) ;
            end
        end
    end
    
    % TABLES TO HOLD VALUES FOR ALL THIS SUBJECT'S TRIALS
    currSubERPs = NaN(size(lats,1), numTrials + 4) ;                % Current Subject's ERP Waveform(s)
    if params.calcVals.meanAmpMethod(1)                                     % Mean Amplitude by User-Defined Windows
        currSubMeanAmp_wind = cell(numTrials + 1, ...
            size(params.calcVals.windows, 1)) ;
    end
    if params.calcVals.meanAmpMethod(2)                                     % Mean Amplitude by Zero-Bound Windows
        currSubMeanAmp_zero = cell(numTrials + 1, 1) ;           
    end
    
    %% PER TRIAL
    for currtrial=1:numTrials+1
        try
            % FOR TRIALS WITHIN THE DATA...
            if currtrial <= numTrials && currfile <= size(FileNames,2)
                if strcmpi(params.format, '1')
                    currTrialData = currsub.data(trialBounds(currtrial, ...
                        1):trialBounds(currtrial, 2), :) ;
                    currTrialERP = mean(currTrialData(:, chanCols), 2, 'omitnan') ;
                elseif strcmpi(params.format, '2')
                    currTrialERP = mean(currsub.data(chanCols, :, currtrial), ...
                    1, 'omitnan')' ;
                end  
            % FOR THE AVERAGE OF TRIALS FOR A FILE...
            elseif currtrial > numTrials && currfile <= size(FileNames,2)
                % Calculate the average of all trials, omitting NaN
                currTrialERP = mean(currSubERPs(:,1:end-4),2, 'omitnan') ;
            % FOR THE AVERAGE OF TRIALS FOR THE AVERAGE OF FILES...
            else
                currTrialERP = allSubsERP{end} ;
            end
            
            % Add ERP to all ERPs for the current subject
            currSubERPs(:,currtrial) = currTrialERP ;
            
            % CALCULATE VALUES
            % REMOVE BASELINE FOR CALCULATIONS
            currTrialERP_noBL = currTrialERP(find(lats==0):end,:) ;
            
            % CREATE WINDOWS
            [currTrialWins, currTrialGlobal] = createWindows(lats, ...
                currTrialERP_noBL, currsub_winLats) ;
            if params.calcVals.meanAmpMethod(2)
               zeroCrossWins = createZeroWins(currTrialGlobal) ; 
            end
            
            % CALCULATE MEAN AMPLITUDE
            if params.calcVals.meanAmpMethod(1)
                for i=1:size(currTrialWins,2)
                    currSubMeanAmp_wind{currtrial,i} = mean(currTrialWins{i}(:,2), 'omitnan') ;
                end
            end
            if params.calcVals.meanAmpMethod(2)
                currSubMeanAmp_zero{currtrial} = calcAmpsZeros(zeroCrossWins) ;
            end
        catch e
            fprintf(['Failed to process trial ' num2str(currtrial) '...\n']) ;
            continue ;
        end
    end
    
    try
        % CALCULATE 95% CONFIDENCE INTERVAL AND STANDARD ERROR
        if numTrials > 1
            currSubERPs(:,end-2:end-1) = prctile(currSubERPs(:,1:end-4)', ...
                abs([0,100]-(100-95)/2))' ;
            stdError = [] ;
            for i=1:size(currSubERPs,1)
                data = currSubERPs(i,1:end-4)' ;
                stdError = [stdError; std(data)/sqrt(length(data))] ;
            end
            currSubERPs(:,end) = stdError ;
            clear('stdError') ;
        else
            if currfile <= size(FileNames,2)
                currSubERPs = currSubERPs(:,1) ;
            else
                if size(FileNames,2) > 1
                    currSubERPs(:,end-2:end-1) = prctile(temp(:,1:end)', ...
                        abs([0,100]-(100-95)/2))' ;
                    stdError = [] ;
                    for i = 1:size(temp,1)
                        data = temp(i,1:end)' ;
                        stdError = [stdError; std(data)/sqrt(length(data))] ;
                    end
                    currSubERPs(:,end) = stdError ;
                    currSubERPs = [currSubERPs(:,1) currSubERPs(:,end-2:end)] ;
                end
            end
        end
        
        % ADD THE VALUES TO THE SUBJECT-LEVEL ARRAY OF DATA:
        allSubsERP{currfile} = currSubERPs ;
        if params.calcVals.meanAmpMethod(1)
            MeanAmp_wind{currfile} = currSubMeanAmp_wind ;
        end
        if params.calcVals.meanAmpMethod(2)
            MeanAmp_zero{currfile} = currSubMeanAmp_zero ;
        end
    catch; continue ;
    end
end

%% PLOT AVERAGE ERP WAVEFORMS
aveToPlot = [] ;
for i=1:size(FileNames,2)
    if size(allSubsERP{i},2) > 1
        aveToPlot = [aveToPlot allSubsERP{i}(:,end-3)] ;
    else
        aveToPlot = [aveToPlot allSubsERP{i}(:,end)] ;%#ok<*AGROW> 
    end
end
aveToPlot = [aveToPlot allSubsERP{end}] ;
if params.plot
    if size(FileNames,2) > 1
        tiledlayout(1,3) ;
        tile1 = nexttile ;
        p1 = plot(tile1, lats, aveToPlot(:,1:end-4)) ;
        title(tile1, 'Average ERP Over Trials for Each File') ;
        
        tile2 = nexttile ;
        p2 = plot(tile2, lats, [aveToPlot(:,end-3), ...
            aveToPlot(:,end-3)-aveToPlot(:,end), ...
            aveToPlot(:,end-3)+aveToPlot(:,end)], '--k') ;
        p2(1).LineStyle = '-' ;
        p2(1).Color = 'r' ;
        title(tile2, 'Average ERP Across Files w/ Standard Error') ;
        
        tile3 = nexttile;
        p3 = plot(tile3, lats, aveToPlot(:,1:end-3), 'k') ;
        p3(end).Color = 'r' ;
        title(tile3, 'All Files Average + Average Across Files') ;
    else
        plot(lats, aveToPlot(:,1)) ;
        title('Average ERP for Individual File') ;
    end
    
    %% PLOT TRIAL WAVEFORMS BY FILE
    for i=1:size(FileNames,2)
        figure() ;
        currPlot = allSubsERP{i} ;
        if size(currPlot,2) > 1
            tiledlayout(1,3) ;
            tile1 = nexttile ;
            p1 = plot(tile1, lats, currPlot(:,1:end-4)) ;
            title(tile1, sprintf('ERP Per Trial for %s', FileNames{i})) ;

            tile2 = nexttile ;
            p2 = plot(tile2, lats, [currPlot(:,end-3), ...
                currPlot(:,end-3)-currPlot(:,end), ...
                currPlot(:,end-3)+currPlot(:,end)], '--k') ;
            p2(1).LineStyle = '-' ;
            p2(1).Color = 'r' ;
            title(tile2, 'Average ERP Across Trials w/ Standard Error') ;

            tile3 = nexttile;
            p3 = plot(tile3, lats, currPlot(:,1:end-3), 'k') ;
            p3(end).Color = 'r' ;
            title(tile3, 'All Trials + Average Over Trials') ;
        else
            plot(lats, currPlot(:,1)) ;
            title('ERP for Single Trial') ;
        end
    end
end

%% COMPILE STATS TABLES
cd([srcDir filesep 'readieERPs' filesep 'ERP_timeseries']) ;
% ABBREVIATE LONG FILE NAMES
outputFileNames = FileNames ;
for i=1:size(FileNames,2)
    if size(FileNames{i},2) > 63
        outputFileNames{i} = [FileNames{i}(1:60) '...'] ;
    end
end

% AVERAGE OVER TRIALS ALL FILES - ERP WAVEFORM
saveName = ['AllSubsAve_readieERPs_' datestr(now, 'dd-mm-yyyy') '.csv'] ;
indx = 2 ;
while isfile(saveName)
    saveName = strrep(saveName, '.csv', ['_' num2str(indx) '.csv']) ;
    indx = indx + 1 ;
end
if size(FileNames,2) > 1
    printTab = array2table(aveToPlot, 'VariableNames', [outputFileNames ...
        'Average ERP Waveform' '95% CI Lower Bound' '95% CI Upper Bound' ...
        'Standard Error']) ;
else
    printTab = array2table(aveToPlot(:,1), 'VariableNames', outputFileNames) ;
end
writetable(printTab, saveName) ;

% INDIVIDUAL TRIALS - ERP WAVEFORM
for i=1:size(FileNames,2)
    if size(allSubsERP{i},2) > 1
        varNames = cell(1, size(allSubsERP{i},2)-4) ;
        for j=1:size(allSubsERP{i},2)-4
            varNames{1,j} = ['Trial ' num2str(j)] ;
        end
        varNames = [varNames 'Average ERP Waveform' '95% CI Lower Bound' ...
            '95% CI Upper Bound' 'Standard Error'] ;
    else; varNames = {'Trial 1'} ;
    end
    printTab = array2table(allSubsERP{i}, 'VariableNames', varNames);
    saveName = strrep(FileNames{i}, ext, ['_' datestr(now, 'dd-mm-yyyy') '.csv']);
    indx = 2 ;
    while isfile(saveName)
        saveName = strrep(saveName, '.csv', ['_' num2str(indx) '.csv']) ;
        indx = indx + 1 ;
    end
    writetable(printTab, saveName) ;
end
 
% CALCULATED VALUES OUTPUT
cd([srcDir filesep 'readieERPs' filesep 'ERP_calculatedVals']) ;
% Mean Amplitude
meanAmpNames = cell(1, params.calcVals.meanAmpMethod(1)*size(params.calcVals.windows,1) + ...
    params.calcVals.meanAmpMethod(2)*1) ;
meanAmpTab = cell(size(FileNames,2)+1, params.calcVals.meanAmpMethod(1)*size(params.calcVals.windows,1) + ...
    params.calcVals.meanAmpMethod(2)*1) ;
if params.calcVals.meanAmpMethod(1)
    for i=1:size(params.calcVals.windows,1)
        meanAmpNames{i} = ['Mean Amplitude for Window ' ...
            params.calcVals.windows{i,1} '-' params.calcVals.windows{i,2}] ;
    end
    for i=1:size(FileNames,2)+1
        meanAmpTab(i,1:size(params.calcVals.windows)) = MeanAmp_wind{i}(end,:) ;
    end
end
if params.calcVals.meanAmpMethod(2)
    meanAmpNames{1,end} = 'Mean Amplitude for 0-Bound Windows' ;
    for i=1:size(FileNames,2)+1
        meanAmpTab{i,end} = MeanAmp_zero{i}(end,1) ;
    end
end
meanAmpTab = cell2table(meanAmpTab, 'VariableNames', meanAmpNames, ...
    'RowNames', [outputFileNames'; 'Average ERP Waveform']) ;


printValTab = meanAmpTab ;
saveValName = ['AllSubsAve_readieERPvals_' datestr(now, 'dd-mm-yyyy') '.csv'] ;
indx = 2 ;
while isfile(saveValName)
    saveValName = strrep(saveName, '.csv', ['_' num2str(indx) '.csv']) ;
    indx = indx + 1 ;
end
writetable(printValTab, saveValName, 'WriteRowNames', true, 'QuoteStrings', ...
    true);
    
measures = meanAmpNames ;
for i =1:size(FileNames,2)
    printValTab = [] ;
    if params.calcVals.meanAmpMethod(1)
        printValTab = [printValTab MeanAmp_wind{i}] ;
    end
    if params.calcVals.meanAmpMethod(2)
        printValTab = [printValTab MeanAmp_zero{i}] ;
    end
   
    varNames = cell(1,size(printValTab,1)) ;
    for j=1:size(printValTab,1)-1
       varNames{1,j} = ['Trial ' num2str(j)] ;
    end
    varNames{size(printValTab,1)} = 'Average of Trials' ;
    printValTab = cell2table(printValTab, 'VariableNames', ...
        measures, 'RowNames', varNames) ;
   
    saveValName = [strrep(FileNames{i}, ext, '') '_readieERPvals_' ...
        datestr(now, 'dd-mm-yyyy') '.csv'] ;
    indx = 2 ;
    while isfile(saveValName)
        saveValName = strrep(saveValName, '.csv', ['_' ...
            num2str(indx) '.csv']) ;
        indx = indx + 1 ;
    end
    
    writetable(printValTab, saveValName, 'WriteRowNames', ...
        true, 'QuoteStrings', true);
end
            