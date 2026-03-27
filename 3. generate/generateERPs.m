% Generate ERPs - a post-processing script in association with HAPPE+ER to
% create ERP waveforms and calculate common ERP measures.
%
% Developed at Northeastern University's PINE Lab
%
% For a detailed description of this script and user options, please see 
% the following manuscript(s):
%   Monachino, et al., (2022)
%
% Author: A.D. Monachino, PINE Lab at Northeastern University, 2022
%
% This file is part of HAPPE.
% Copyright 2018-2023 Alexa Monachino, Kelsie Lopez, Laurel Gabard-Durnam
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

%% SET FOLDERS AND PATH
% Use input from the command line to set the path to the data. If an 
% invalid path is entered, repeat until a valid path is entered.
clear ;
fprintf('Preparing HAPPE generateERPs...\n') ;
happeDir = strrep(fileparts(which(mfilename('fullpath'))), [filesep '3. ' ...
    'generate'], '') ;
addpath([happeDir filesep '3. generate'], ...
    [happeDir filesep 'files'], ...
    [happeDir filesep 'scripts'], ...
    [happeDir filesep 'scripts' filesep 'ui'], ...
    [happeDir filesep 'scripts' filesep 'support']) ;

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
if ~isfolder([srcDir filesep 'generateERPs']); mkdir([srcDir filesep 'generateERPs']) ;
end
addpath('generateERPs') ;
cd([srcDir filesep 'generateERPs']) ;
if ~isfolder([srcDir filesep 'generateERPs' filesep 'ERP_timeseries'])
    mkdir([srcDir filesep 'generateERPs' filesep 'ERP_timeseries']) ;
end
if ~isfolder([srcDir filesep 'generateERPs' filesep 'ERP_calculatedVals'])
    mkdir([srcDir filesep 'generateERPs' filesep 'ERP_calculatedVals']) ;
end
fprintf('Output folders created.\n') ;

%% SET PARAMETERS
% DETERMINE IF USING PRE-EXISTING SET OF PARAMETERS
[preExist, params, changedParams] = genERPs_isPreExist() ;

% SET PARAMETERS THROUGH USER INPUT
params = genERPs_setParams(params, preExist, changedParams) ;

% SAVE INPUT PARAMETERS
% If created a new or changed a parameter set, save as a new .mat file
if ~preExist || changedParams
    % Prompt to use a default or custom name for parameter the file. 
    % If the file exists, ask to create new file with a different name
    % or overwrite existing file.
    fprintf(['Parameter file save name:\n  default = Default name (genERP' ...
        '_parameters_dd-mm-yyyy.mat)\n  custom = Create a custom file name' ...
        '\n']) ;
    if choose2('custom', 'default')
        paramFile = paramFile_validateExist(['genERP_parameters_' ...
            datestr(now, 'dd-mm-yyyy') '.mat'], 'genERP_parameters_', 2) ;
    else
        fprintf('File name (Do not include .mat):\n') ;
        paramFile = paramFile_validateExist([input('> ', 's') '.mat'], ...
            'genERP_parameters_', 0) ;
    end

    % Save the params variable to a .mat file using the name created above.
    fprintf('Saving parameters...\n') ;
    save(paramFile, 'params') ;
    fprintf('Parameters saved.\n') ;
end
clear('preExist', 'changedParams', 'paramFile') ;


%% COLLECT FILES
cd(srcDir) ;
if params.indivTrials; ext = '_IndivTrial' ;
else; ext = '_AveOverTrials' ; 
end
fprintf(sprintf(['Enter the suffix used for this dataset, including stimulus tag' ...
    ' (if applicable).\nIf no extension beyond "%s", press ' ...
    'enter/return.\n'], ext)) ;
suffix = input('> ', 's') ;
FileNames = {dir(['*' ext suffix '.txt']).name} ;
if size(FileNames,2) < 1; error('ERROR: No files detected') ; end

%% READ IN BAD CHANNELS
if ~params.badChans.inc
    % Load the File and convert table to cell
    badChans = table2cell(readtable(params.badChans.file,Delimiter={','})) ;

    % Split the string into a cell array
    for currtrial=1:size(badChans, 1)
        badChans{currtrial,2} = split(badChans{currtrial,2})' ;
    end
end

%% INITIALIZE OUTPUT TABLES
allSubsERP = cell(1,size(FileNames,2)+1) ;
if params.calcVals.on
    allSubPeaks = cell(1,size(FileNames,2)+1) ;     % Peaks (Max and Mins)
    MeanAmp_wind = cell(1,size(FileNames,2)+1) ;    % Mean Amplitude by User-Defined Windows
    MeanAmp_zero = cell(1,size(FileNames,2)+1) ;    % Mean Amplitude by Zero-Bound Windows
    AUC_wind = cell(1,size(FileNames,2)+1) ;        % Area Under Curve by User-Defined Windows
    AUC50_wind = cell(1,size(FileNames,2)+1) ;      % 50% Area Under Curve by User-Defined Windows
    AUC_zero = cell(1,size(FileNames,2)+1) ;        % Area Under Curve by Zero-Bound Windows
    AUC50_zero = cell(1,size(FileNames,2)+1) ;      % 50% Area Under Curve by Zero-Bound Windows
end

%% PER FILE...
for currfile=1:size(FileNames, 2)+1
    % If processing a file...
    if currfile <= size(FileNames,2)
        try
            % LOAD FILE
            currsub = importdata(FileNames{currfile}) ;
    
            % COMPILE LIST OF BAD CHANNELS
            if ~params.badChans.inc
                subBadChans = badChans{contains(badChans(:,1), ... 
                    strrep(FileNames{currfile}, ['_processed' ext ...
                    suffix '.txt'], '')), 2} ;
            end
    
            % SET CHANNELS OF INTEREST
            if strcmpi(params.chans.subset, 'coi_exclude')
                subChanIDs = setdiff(currsub.colheaders, [params.chans.IDs, 'Time']) ;
                if ~params.badChans.inc
                    subChanIDs = setdiff(subChanIDs, subBadChans) ; 
                end
            elseif strcmpi(params.chans.subset, 'coi_include')
                subChanIDs = params.chans.IDs ;
                if ~params.badChans.inc
                    subChanIDs = setdiff(subChanIDs, ...
                        intersect(subBadChans, params.chans.IDs)) ; 
                end
            elseif strcmpi(params.chans.subset, 'all')
                subChanIDs = setdiff(currsub.colheaders, 'Time') ;
                if ~params.badChans.inc
                    subChanIDs = setdiff(subChanIDs, subBadChans) ; 
                end
            end
    
            % GET COLUMN NUMBERS
            chanCols = [] ;
            if ~isfield(currsub, 'colheaders')
                currsub.colheaders = strsplit(currsub.textdata{1}) ;
            end
            for int=1:length(subChanIDs)
                chanCols = [chanCols find(ismember(currsub.colheaders, ...
                    subChanIDs{int}))] ;
            end
    
            % SPLIT BY TRIALS
            % Find the range of lats in the data
            lats = unique(currsub.data(:,1)) ;
    
            trialBounds = [find(currsub.data(:,1)==min(lats)), ...
                    find(currsub.data(:,1)==max(lats))] ;
        catch e
            fprintf(['Failed to process ' FileNames{currfile} '...\n']) ;
            continue ;
        end
    % If processing the average of files...
    else
        temp = [] ;
        if params.indivTrials; endInd = 3;
        else; endInd = 0; 
        end
        for i = 1:size(FileNames,2)
            temp = [temp allSubsERP{i}(:,end-endInd)] ;
        end
        allSubsERP{end} = mean(temp,2,'omitnan') ;
        trialBounds = [1 size(allSubsERP{end},1)] ;
    end
    
    % TABLES TO HOLD VALUES FOR ALL THIS SUBJECT'S TRIALS
    currSubERPs = NaN(size(lats,1), size(trialBounds,1)+4) ;                % Current Subject's ERP Waveform(s)
    if params.calcVals.on
        currSubPeaks = cell(size(trialBounds,1)+1, ...                      % Peaks (Max and Min)
            2*size(params.calcVals.windows,1)+6) ;
        if params.calcVals.meanAmpMethod(1)                                 % Mean Amplitude by User-Defined Windows
            currSubMeanAmp_wind = cell(size(trialBounds,1)+1, ...
                size(params.calcVals.windows,1)) ;
        end
        if params.calcVals.meanAmpMethod(2)                                 % Mean Amplitude by Zero-Bound Windows
            currSubMeanAmp_zero = cell(size(trialBounds,1)+1,1) ;           
        end
        if params.calcVals.aucMethod(1)                                     % 50% & 100% Area Under Curve by User-Defined Windows
            currSubAUC_wind = cell(size(trialBounds,1)+1, ...
                size(params.calcVals.windows,1)+1) ;
            currSubAUC50_wind = cell(size(trialBounds,1)+1, ...
                2*size(params.calcVals.windows,1)+2) ;
        end
        if params.calcVals.aucMethod(2)                                     % 50% & 100% Area Under Curve by Zero-Bound Windows
            currSubAUC_zero = cell(size(trialBounds,1)+1,1) ;
            currSubAUC50_zero = cell(size(trialBounds,1)+1,1) ;
        end
    end
    
    %% PER TRIAL
    for currtrial=1:size(trialBounds,1)+1
        try
            % FOR TRIALS WITHIN THE DATA...
            if currtrial <= size(trialBounds,1) && currfile <= size(FileNames,2)
                % Pull out the current trial's data:
                currTrialData = currsub.data(trialBounds(currtrial, ...
                    1):trialBounds(currtrial, 2), :) ;
    
                % Calculate the average ERP of the specified channels for this
                % trial:
                trialChanData = [] ;
                for i=1:size(chanCols,2)
                    trialChanData = [trialChanData currTrialData(:, chanCols(i))] ;
                end
                currTrialERP = mean(trialChanData,2) ;
            % FOR THE AVERAGE OF TRIALS FOR A FILE...
            elseif currtrial > size(trialBounds,1) && currfile <= size(FileNames,2)
                % Calculate the average of all trials, omitting NaN
                currTrialERP = mean(currSubERPs(:,1:end-4),2, 'omitnan') ;
            % FOR THE AVERAGE OF TRIALS FOR THE AVERAGE OF FILES...
            else
                currTrialERP = allSubsERP{end} ;
            end
            
            % Add ERP to all ERPs for the current subject
            currSubERPs(:,currtrial) = currTrialERP ;
            
            % CALCULATE VALUES
            if params.calcVals.on
                % CORRECT WINDOW LATENCIES - only works if all files have same
                % latencies
                if currfile == 1
                    for i=1:size(params.calcVals.windows,1)
                        for j=1:2
                            if ~any(ismember(lats, ...
                                    str2num(params.calcVals.windows{i,j}))) %#ok<*ST2NM> 
                                [minVal, closestIndx] = min(abs(lats- ...
                                    str2num(params.calcVals.windows{i,j}))) ;
                                params.calcVals.windows{i,j} = ...
                                    num2str(lats(closestIndx)) ;
                            end
                        end
                    end
                end
                
                % REMOVE BASELINE FOR CALCULATIONS
                currTrialERP_noBL = currTrialERP(find(lats==0):end,:) ;
                
                % CREATE WINDOWS
                [currTrialWins, currTrialGlobal] = createWindows(lats, ...
                    currTrialERP_noBL, params.calcVals.windows) ;
                if params.calcVals.meanAmpMethod(2) || params.calcVals.aucMethod(2)
                   zeroCrossWins = createZeroWins(currTrialGlobal) ; 
                end
                
                % CALCULATE PEAKS
                currSubPeaks(currtrial,1:end-2) = calcPeaks(params.calcVals.windows, ...
                    currTrialWins, currTrialGlobal) ;
                [currSubPeaks{currtrial,end-1}, currSubPeaks{currtrial,end}] = ...
                    getAllPeaks(currTrialGlobal) ;
                
                % CALCULATE MEAN AMPLITUDE
                if params.calcVals.meanAmpMethod(1)
                    for i=1:size(currTrialWins,2)
                        currSubMeanAmp_wind{currtrial,i} = mean(currTrialWins{i}(:,2), 'omitnan') ;
                    end
                end
                if params.calcVals.meanAmpMethod(2)
                    currSubMeanAmp_zero{currtrial} = calcAmpsZeros(zeroCrossWins) ;
                end
                
                % CALCULATE AUC AND 50% AUC
                if params.calcVals.aucMethod(1)
                    [currSubAUC_wind(currtrial,1:end), currSubAUC50_wind(currtrial,1:end)] = ...
                        calcAUC(currTrialWins, currTrialGlobal) ;
                end
                if params.calcVals.aucMethod(2)
                    [currSubAUC_zero{currtrial}, currSubAUC50_zero{currtrial}] = ...
                        calcAUCZeros(zeroCrossWins) ;
                end
            end
        catch
            fprintf(['Failed to process trial ' num2str(currtrial) '...\n']) ;
            continue ;
        end
    end
    
    try
        % CALCULATE 95% CONFIDENCE INTERVAL AND STANDARD ERROR
        if size(trialBounds,1) > 1
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
        if params.calcVals.on
            allSubPeaks{currfile} = currSubPeaks ;
            if params.calcVals.meanAmpMethod(1)
                MeanAmp_wind{currfile} = currSubMeanAmp_wind ;
            end
            if params.calcVals.meanAmpMethod(2)
                MeanAmp_zero{currfile} = currSubMeanAmp_zero ;
            end
            if params.calcVals.aucMethod(1)
                AUC_wind{currfile} = currSubAUC_wind ;
                AUC50_wind{currfile} = currSubAUC50_wind ;
            end
            if params.calcVals.aucMethod(2)
                AUC_zero{currfile} = currSubAUC_zero ;
                AUC50_zero{currfile} = currSubAUC50_zero ;
            end
        end
    catch; continue ;
    end
end

%% PLOT AVERAGE ERP WAVEFORMS
aveToPlot = [] ;
if params.indivTrials
    for i=1:size(FileNames,2)
        aveToPlot = [aveToPlot allSubsERP{i}(:,end-3)] ;                    %#ok<*AGROW> 
    end
    aveToPlot = [aveToPlot allSubsERP{end}] ;
else
    for i=1:size(allSubsERP,2)
        aveToPlot = [aveToPlot allSubsERP{i}] ;
    end
end
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
    if params.indivTrials
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
end

%% COMPILE STATS TABLES
cd([srcDir filesep 'generateERPs' filesep 'ERP_timeseries']) ;
% ABBREVIATE LONG FILE NAMES
outputFileNames = FileNames ;
for i=1:size(FileNames,2)
    if size(FileNames{i},2) > 30
        outputFileNames{i} = [FileNames{i}(1:30-3) '...'] ;
    end
end

% AVERAGE OVER TRIALS ALL FILES - ERP WAVEFORM
saveName = ['AllSubsAve_generatedERPs' suffix '_' datestr(now, 'dd-mm-yyyy') '.csv'] ;
indx = 2 ;
while isfile(saveName)
    saveName = strrep(saveName, '.csv', ['_' num2str(indx) '.csv']) ;
    indx = indx + 1 ;
end
if length(FileNames) > 1
    printTab = array2table(aveToPlot, 'VariableNames', [outputFileNames ...
        'Average ERP Waveform' '95% CI Lower Bound' '95% CI Upper Bound' ...
        'Standard Error']) ;
else
    printTab = array2table(aveToPlot(:,1), 'VariableNames', outputFileNames) ;
end
writetable(printTab, saveName) ;

% INDIVIDUAL TRIALS - ERP WAVEFORM
if params.indivTrials
    if ~params.export.csv
        saveName = ['IndivSubs_generatedERPs' suffix '_' datestr(now, ...
            'dd-mm-yyyy') '.xlsx'] ;
        indx = 2 ;
        while isfile(saveName)
            saveName = strrep(saveName, '.xlsx', ['_' num2str(indx) ...
                '.xlsx']) ;
            indx = indx + 1 ;
        end
    end
    for i=1:size(FileNames,2)
        varNames = cell(1, size(allSubsERP{i},2)-4) ;
        for j=1:size(allSubsERP{i},2)-4
            varNames{1,j} = ['Trial ' num2str(j)] ;
        end
        printTab = array2table(allSubsERP{i}, 'VariableNames', [varNames ...
            'Average ERP Waveform' '95% CI Lower Bound' '95% CI Upper Bound' ...
            'Standard Error']);
        if params.export.csv
            saveName = strrep(FileNames{i}, ['_processed_IndivTrial' suffix '.txt'], ...
                ['_generatedERPs' suffix '_' datestr(now, 'dd-mm-yyyy') '.csv']);
            indx = 2 ;
            while isfile(saveName)
                saveName = strrep(saveName, '.csv', ['_' num2str(indx) '.csv']) ;
                indx = indx + 1 ;
            end
            writetable(printTab, saveName) ;
        else
            writetable(printTab, saveName, 'Sheet', outputFileNames{i}, ...
            	'WriteRowNames', true) ;
        end
    end
end

if params.calcVals.on
    cd([srcDir filesep 'generateERPs' filesep 'ERP_calculatedVals']) ;
    % AVERAGE OVER TRIALS ALL FILES - ERP MEASURE VALUES
    % Peaks:
    peakNames = cell(1, 2*size(params.calcVals.windows,1)+6) ;
    for i=1:size(params.calcVals.windows,1)
        peakNames{i*2-1} = [params.calcVals.windows{i,3} ...
            ' Value for Window ' params.calcVals.windows{i,1} '-' ...
            params.calcVals.windows{i,2}] ;
        peakNames{i*2} = ['Latency at ' params.calcVals.windows{i,3} ...
            ' for Window ' params.calcVals.windows{i,1} '-' ...
            params.calcVals.windows{i,2}] ;
    end
    peakNames(end-5:end) = {'Global Max Value', ...
        'Latency at Global Max', 'Global Min Value', 'Latency at Global Min', ...
        'All Maxes (with Values)' 'All Mins (with Values)'} ;
    peakTab = cell(size(FileNames,2)+1,size(peakNames,2)) ;
    for i=1:size(allSubPeaks,2)
        peakTab(i,:) = allSubPeaks{i}(end,:) ;
    end
    peakTab = cell2table(peakTab, 'VariableNames', peakNames, 'RowNames', ...
        [outputFileNames'; 'Average ERP Waveform']) ;
    
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
    
    % Area Under the Curve
    AUCNames = cell(1, params.calcVals.aucMethod(1)*(size(params.calcVals.windows,1)+1) + ...
        params.calcVals.aucMethod(2)*1) ;
    AUC50Names = cell(1, params.calcVals.aucMethod(1)*(2*size(params.calcVals.windows,1)+2) + ...
        params.calcVals.aucMethod(2)*1) ;
    AUCTab = cell(size(FileNames,2)+1, params.calcVals.aucMethod(1)*(size(params.calcVals.windows,1)+1) + ...
        params.calcVals.aucMethod(2)*1) ;
    AUC50Tab = cell(size(FileNames,2)+1, params.calcVals.aucMethod(1)*(2*size(params.calcVals.windows,1)+2) + ...
        params.calcVals.aucMethod(2)*1) ;
    if params.calcVals.aucMethod(1)
        for i=1:size(params.calcVals.windows,1)
            AUCNames{i} = ['Area Under Curve for Window ' ...
                params.calcVals.windows{i,1} '-' params.calcVals.windows{i,2}] ;
            AUC50Names{i*2-1} = ['50% Area Under Curve for Window' ...
                params.calcVals.windows{i,1} '-' params.calcVals.windows{i,2}] ;
            AUC50Names{i*2} = ['Latency at 50% Area Under Curve for Window' ...
                params.calcVals.windows{i,1} '-' params.calcVals.windows{i,2}] ;
        end
        AUCNames{size(params.calcVals.windows,1)+1} = 'Global Area Under the Curve' ;
        AUC50Names(2*size(params.calcVals.windows,1)+1:2*size(params.calcVals.windows,1)+2) = ...
            {'Global 50% Area Under Curve', 'Latency at 50% Area Under the Curve'} ;
        for i=1:size(FileNames,2)+1
            AUCTab(i,1:size(params.calcVals.windows,1)+1) = AUC_wind{i}(end,:) ;
            AUC50Tab(i,1:(2*size(params.calcVals.windows,1)+2)) = AUC50_wind{i}(end,:) ;
        end
    end
    if params.calcVals.aucMethod(2)
        AUCNames{1,end} = 'Area Under Curve for 0-Bound Windows' ;
        for i=1:size(FileNames,2)+1
            AUCTab{i,end} = AUC_zero{i}(end,1) ;
        end
        AUC50Names{1,end} = '50% Area Under Curve for 0-Bound Windows' ;
        for i=1:size(FileNames,2)+1
            AUC50Tab{i,end} = AUC50_zero{i}(end,1) ;
        end
    end
    AUCbothTab = cell2table([AUCTab AUC50Tab], 'VariableNames', [AUCNames AUC50Names], ...
        'RowNames', [outputFileNames'; 'Average ERP Waveform']) ;
    
    printValTab = [peakTab meanAmpTab AUCbothTab] ;
    saveValName = ['AllSubsAve_generatedERPvals' suffix '_' ...
        datestr(now, 'dd-mm-yyyy') '.csv'] ;
    indx = 2 ;
    while isfile(saveValName)
        saveValName = strrep(saveName, '.csv', ['_' num2str(indx) '.csv']) ;
        indx = indx + 1 ;
    end
    writetable(printValTab, saveValName, 'WriteRowNames', true, 'QuoteStrings', ...
        true);
    
    if params.indivTrials
        measures = [peakNames meanAmpNames AUCNames AUC50Names] ;
        % EXPORT BY SUBJECT
        if strcmpi(params.export.catg, 'subject')
            if ~params.export.csv
                saveValName = ['IndivSubs_generatedERPvals' suffix '_' ...
                    datestr(now, 'dd-mm-yyyy') '.xlsx'] ;
                indx = 2 ;
                while isfile(saveValName)
                    saveValName = strrep(saveValName, '.xlsx', ['_' ...
                        num2str(indx) '.xlsx']) ;
                    indx = indx + 1 ;
                end
            end
            for i =1:size(FileNames,2)
                printValTab = allSubPeaks{i} ;
                if params.calcVals.meanAmpMethod(1)
                    printValTab = [printValTab MeanAmp_wind{i}] ;
                end
                if params.calcVals.meanAmpMethod(2)
                    printValTab = [printValTab MeanAmp_zero{i}] ;
                end
                if params.calcVals.aucMethod(1)
                    printValTab = [printValTab AUC_wind{i} AUC50_wind{i}] ;
                end
                if params.calcVals.aucMethod(2)
                    printValTab = [printValTab AUC_zero{i} AUC50_zero{i}] ;
                end
               
                varNames = cell(1,size(printValTab,1)) ;
                for j=1:size(printValTab,1)-1
                   varNames{1,j} = ['Trial ' num2str(j)] ;
                end
                varNames{size(printValTab,1)} = 'Average of Trials' ;
                printValTab = cell2table(printValTab, 'VariableNames', ...
                    measures, 'RowNames', varNames) ;
               
                if params.export.csv
                    saveValName = [strrep(FileNames{i}, ['_processed_IndivTrial' ...
                        suffix '.txt'], '') '_generatedERPvals' suffix '_' ...
                        datestr(now, 'dd-mm-yyyy') '.csv'] ;
                    indx = 2 ;
                    while isfile(saveValName)
                        saveValName = strrep(saveValName, '.csv', ['_' ...
                            num2str(indx) '.csv']) ;
                        indx = indx + 1 ;
                    end
                    
                    writetable(printValTab, saveValName, 'WriteRowNames', ...
                        true, 'QuoteStrings', true);
                else
                   writetable(printValTab, saveValName, 'Sheet', outputFileNames{i}, ...
                       'WriteRowNames', true) ; % 'QuoteStrings', true) ;
                end
            end
            
        % EXPORT BY METRIC
        else
            byMeasure = cell(1,size(measures,2)) ;
            % Determine Max Number of Trials
            maxTrials = zeros(1, size(FileNames,2)) ;
            for i=1:size(FileNames,2)
                maxTrials(i) = size(allSubsERP{i},2)-3 ;
            end
            maxTrials = max(maxTrials) ;
            indx = 0 ;
            
            % Peaks:
            for i=1:size(peakNames,2)
                byMeasure{i} = cell(size(FileNames,2)+1, maxTrials) ;
                for j=1:size(FileNames,2)+1
                    curr = allSubPeaks{1,j}(:,i)' ;
                    byMeasure{i}(j,1:size(curr,2)-1) = curr(1,1:end-1) ;
                    byMeasure{i}(j,end) = curr(1,end) ;
                end
                indx = indx+1 ;
            end
            
            % Mean Amplitude:
            if params.calcVals.meanAmpMethod(1)
                for i=1:size(MeanAmp_wind{1},2)
                    byMeasure{i+indx} = cell(size(FileNames,2)+1, maxTrials) ;
                    for j=1:size(FileNames,2)+1
                        curr = MeanAmp_wind{1,j}(:,i)' ;
                        byMeasure{i+indx}(j,1:size(curr,2)-1) = curr(1,1:end-1) ;
                        byMeasure{i+indx}(j,end) = curr(1,end) ;
                    end
                end
                indx = indx+i ;
            end
            if params.calcVals.meanAmpMethod(2)
                indx = indx+1 ;
                byMeasure{indx} = cell(size(FileNames,2)+1, maxTrials) ;
                for j=1:size(FileNames,2)+1
                    curr = MeanAmp_zero{1,j}(:,1)' ;
                    byMeasure{indx}(j,1:size(curr,2)-1) = curr(1,1:end-1) ;
                    byMeasure{indx}(j,end) = curr(1,end) ;
                end
            end
            
            % Area Under Curve:
            if params.calcVals.aucMethod(1)
                for i=1:size(AUC_wind{1},2)
                    byMeasure{i+indx} = cell(size(FileNames,2)+1, maxTrials) ;
                    for j=1:size(FileNames,2)+1
                        curr = AUC_wind{1,j}(:,i)' ;
                        byMeasure{i+indx}(j,1:size(curr,2)-1) = curr(1,1:end-1) ;
                        byMeasure{i+indx}(j,end) = curr(1,end) ;
                    end
                end
                indx = indx+i ;
            end
            if params.calcVals.aucMethod(2)
                indx = indx+1 ;
                byMeasure{indx} = cell(size(FileNames,2)+1, maxTrials) ;
                for j=1:size(FileNames,2)+1
                    curr = AUC_zero{1,j}(:,1)' ;
                    byMeasure{indx}(j,1:size(curr,2)-1) = curr(1,1:end-1) ;
                    byMeasure{indx}(j,end) = curr(1,end) ;
                end
            end
            
            % 50% Area Under Curve
            if params.calcVals.aucMethod(1)
                for i=1:size(AUC50_wind{1},2)
                    byMeasure{i+indx} = cell(size(FileNames,2)+1, maxTrials) ;
                    for j=1:size(FileNames,2)+1
                        curr = AUC50_wind{1,j}(:,i)' ;
                        byMeasure{i+indx}(j,1:size(curr,2)-1) = curr(1,1:end-1) ;
                        byMeasure{i+indx}(j,end) = curr(1,end) ;
                    end
                end
                indx = indx+i ;
            end
            if params.calcVals.aucMethod(2)
                indx = indx+1 ;
                byMeasure{indx} = cell(size(FileNames,2)+1, maxTrials) ;
                for j=1:size(FileNames,2)+1
                    curr = AUC50_zero{1,j}(:,1)' ;
                    byMeasure{indx}(j,1:size(curr,2)-1) = curr(1,1:end-1) ;
                    byMeasure{indx}(j,end) = curr(1,end) ;
                end
            end
            
            if ~params.export.csv
                saveValName = ['ByMeasure_generatedERPvals' suffix '_' ...
                    datestr(now, 'dd-mm-yyyy') '.xlsx'] ;
                indx = 2 ;
                while isfile(saveValName)
                    saveValName = strrep(saveValName, '.xlsx', ['_' ...
                        num2str(indx) '.xlsx']) ;
                    indx = indx + 1 ;
                end
            end
            
            varNames = cell(1, maxTrials) ;
            for i=1:maxTrials-1
                varNames{i} = ['Trial ' num2str(i)] ;
            end
            varNames{end} = 'Average ERP Waveform' ;
            
            for i=1:size(measures,2)
                sheetName = strrep(measures{i}, ' ', '_') ;
                if params.export.csv
                    saveValName = [sheetName '_generatedERPvals' suffix ...
                        '_' datestr(now, 'dd-mm-yyyy') '.csv'] ;
                    indx = 2 ;
                    while isfile(saveValName)
                        saveValName = strrep(saveValName, '.csv', ['_' ...
                            num2str(indx) '.csv']) ;
                        indx = indx + 1 ;
                    end
                    printValTab = cell2table(byMeasure{i}, 'VariableNames', ...
                        varNames, 'RowNames', [FileNames'; ...
                        'Average ERP Waveform']) ;
                    writetable(printValTab, saveValName, 'WriteRowNames', ...
                        true, 'QuoteStrings', true);
                else
                    if size(sheetName,2) >= 30; sheetName = sheetName(1:30) ; end
                    printValTab = cell2table(byMeasure{i}, 'VariableNames', ...
                        varNames, 'RowNames', [outputFileNames'; ...
                        'Average ERP Waveform']) ;
                    writetable(printValTab, saveValName, 'Sheet', sheetName, ...
                        'WriteRowNames', true) ;
                end
            end
        end
    end
end