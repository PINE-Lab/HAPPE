% Generate Functional Connectivity - dwPLI - A post-processing script to
% calculate dwPLI on HAPPE-processed data. Relies on functions from
% FieldTrip (Oostenveld et al., 2011), which can be found here:
% https://www.fieldtriptoolbox.org/
%
% If using this script, please ensure you have FieldTrip downloaded and
% stored seperately from your HAPPE folder, and that it is not added to the
% path at the same time as CleanLine (included in your HAPPE download).
%
% If using this script for analyses in publications, cite both HAPPE and
% FieldTrip.
%
% Developed at Northeastern University's PINE Lab
%
% For a detailed description of this script and user options, please see 
% the following manuscript(s):
%   Monachino, et al., (----)
%
% Authors: A.D. Monachino, PINE Lab at Northeastern University, 2022
%          Laurel J. Gabard-Durnam, PINE Lab at Northestern University, 2022
%          Adela Desowska, ..., 2022
%          Lizzie Shepherd, ..., 2022
%
% This file is part of HAPPE.
% Copyright 2018-2022 Alexa Monachino, Kelsie Lopez, Laurel Gabard-Durnam
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
clear ;
fprintf('Preparing HAPPE generatedwPLI add-on...\n') ;
happeDir = strrep(fileparts(which(mfilename('fullpath'))), '\add-ons', '') ;
addpath([happeDir filesep 'scripts'], ...
    [happeDir filesep 'scripts' filesep 'UI_scripts'], ...
    [happeDir filesep 'add-ons'], ...
    [happeDir filesep 'add-ons' filesep 'generate'], ...
    [happeDir filesep 'add-ons' filesep 'generate' filesep 'scripts']) ;

%% DETERMINE AND SET FIELDTRIP PATH
while true
    fieldTripDir = input('Enter your FieldTrip folder location, including full file path:\n> ','s') ;
    if exist(fieldTripDir, 'dir') == 7; break ;
    else; disp("Invalid input: please enter the complete path to the folder.") ;
    end
end
addpath(fieldTripDir) ;
addpath([fieldTripDir filesep 'external' filesep 'eeglab']) ;

%% DETERMINE AND SET PATH TO THE DATA
while true
    srcDir = input('Enter the path to the folder containing the processed dataset(s):\n> ','s') ;
    if exist(srcDir, 'dir') == 7; break ;
    else; disp("Invalid input: please enter the complete path to the folder containing the dataset(s).") ;
    end
end
cd(srcDir) ;

%% DETERMINE OUTPUT FORMAT
fprintf(['Choose your export format:\n  sheets = A single Excel file' ...
    ' with multiple sheets\n  files = Multiple .csv files\n']) ;
exportCSV = choose2('sheets', 'files') ;

%% DEFINE BANDS
fprintf(['Choose how to define bands:\n  default - Delta (2-4), Theta (4-6),' ...
    ' LowAlpha (6-9), HighAlpha (9-12), Beta (13-30), Gamma (30-45)\n  ' ...
    'custom - Create custom bands using user input\n']) ;
customBands = choose2('default', 'custom') ;
if customBands
    bands = [] ;
    fprintf(['Enter the frequency bands of interest:\n' ...
    'Enter each band as a list, each element seperated by a blank' ...
    ' space,\nwith the band name as the first item in the list. ' ...
    'Press enter/return between entries.\nWhen you ' ...
    'have entered all bands, input "done" (without quotations).\n' ...
    'For optimization purposes, we do not recommend crossing 30Hz within' ...
    'a single band.\nExample: LowAlpha 6 9\n']) ;
    while true
        temp = split(input('> ', 's'))' ;
        if size(temp,2) == 1 && strcmpi(temp, 'done'); break ;
        elseif size(temp,2) == 3
            bands = [bands; temp] ;                                         %#ok<AGROW> 
        else
            fprintf(['Invalid input: enter a band name and the two frequency' ...
                ' limits or "done" (without quotations).\n']) ;
            continue ;
        end
    end
else
    bands = {'Delta', '2', '4'; 'Theta', '4', '6'; 'LowAlpha', '6', '9'; 'HighAlpha', ...
        '9', '12'; 'Beta', '13', '20'; 'Gamma', '30', '45'} ;
end

fprintf('Calculate dwPLI across all frequencies? [Y/N]\n') ;
calcAll = choose2('N', 'Y') ;
if calcAll; bands = [bands; {'All', NaN, NaN}] ; end

%% SUBSET
fprintf('Perform analyses between subsets of electrodes? [Y/N]\n') ;
compSubset = choose2('N','Y') ;
if compSubset
    subsetList = [];
    indx = 1 ;
    while true
        fprintf(['Enter the electrode(s) in the first subset, one at a time, ' ...
        'pressing enter/return between entries.\nExample: E71\nWhen you ' ...
        'have entered all channels in this subset, input "done" (without' ...
        ' quotations)\n']) ;
        temp = UI_cellArray(1, []) ;
        subsetList{indx,1} = unique(temp(~cellfun('isempty',temp))) ;
        fprintf(['Enter the electrode(s) in the second subset, the same' ...
            ' way as the first subset.\n']) ;
        temp = UI_cellArray(1, []) ;
        subsetList{indx,2} = unique(temp(~cellfun('isempty',temp))) ;
        indx = indx+1 ;
        fprintf('Compare another pair of subsets? [Y/N]\n') ;
        if ~choose2('N','Y'); break; end
    end  
end

%% CREATE OUTPUT FOLDERS
% Create the folders in which to store outputs
fprintf('Creating output folder...\n') ;
if ~isfolder([srcDir filesep 'generate_dwPLI'])
    mkdir([srcDir filesep 'generate_dwPLI']) ;
end
addpath('generate_dwPLI') ;
fprintf('Output folder created.\n') ;

%% COLLECT FILES
fprintf(['Enter any text, including spaces, between "processed" and ".set"' ...
    ' (if applicable).\nIf no text exists press enter/return.\nExample: ' ...
    'For "##_processed_rerun-2022.set", enter "_rerun-2022" (without quotations).\n']) ;
suffix = input('> ', 's') ;
FileNames = {dir(['*_processed' suffix '.set']).name} ;
if size(FileNames,2) < 1; error('ERROR: No files detected') ; end

%% INITIALIZE VARIABLES
FC_mat = cell(size(FileNames,2), size(bands,1)) ; % Functional Connectivity Output Matrix
FC_ave = cell(size(FileNames,2), size(bands,1)) ;
FC_abs = cell(size(FileNames,2), size(bands,1)) ;
varNames = cell(size(FileNames,2),1) ;
if compSubset; FC_subsets = cell(size(FileNames,2), size(subsetList,1), size(bands,1)) ; end

%% SET CFG SETTINGS
% Laplacian Filter:
cfg_laplace = struct('channel', 'all', 'reref', 'yes', 'refmethod', 'laplace', ...
    'refchannel', 'all') ;
% Gamma FFT Processing:
cfg_TF_high = struct('method', 'mtmfft', 'keeptrials', 'yes', 'output', ...
    'fourier', 'taper', 'dpss', 'foi', 30:0.5:100, 'tapsmofrq', 3) ;
% Non-Gamma FFT Processing:
cfg_TF_low = struct('method', 'mtmfft', 'keeptrials', 'yes', 'output', ...
    'fourier', 'taper', 'hanning', 'foi', 0.1:1:30) ;
% All Freqs FFT Processing:
cfg_TF_all = struct('method', 'mtmfft', 'keeptrials', 'yes', 'output', ...
    'fourier', 'taper', 'hanning', 'foi', 0.1:1:100) ;
% Band Average dwPLI:
cfg_band = struct('frequency', [], 'avgoverfreq', 'yes') ;

%% RUN OVER FILES
for currfile=1:size(FileNames,2)
    % LOAD THE DATA
    EEGraw = load('-mat', FileNames{currfile}) ;

    % CONVERT TO FIELDTRIP FORMAT
    EEG = eeglab2fieldtrip(EEGraw, 'preprocessing') ;

    % PRE-PROCESS TO ADD LAPLACIAN FILTER
    EEG = ft_preprocessing(cfg_laplace, EEG) ;
    varNames{currfile} = EEG.label ;

    % IF SUBSETS, COLLECT THE CHANNEL INDEXES
    if compSubset
        setIndxs = cell(size(subsetList,1), size(subsetList,2)) ;
        for currSet=1:size(subsetList,1)
            for j=1:2
                set = cell(1,size(subsetList{currSet,j},2)) ;
                for i=1:size(subsetList{currSet,j},2)
                    set{i} = find(ismember(EEG.label, subsetList{currSet,j}{i})) ;
                end
                setIndxs{currSet,j} = set ;
            end
        end
    end

    % CALCULATE dwPLI FC
    dwPLI_high = ft_connectivityanalysis(struct('method', 'wpli_debiased'), ...
        ft_freqanalysis(cfg_TF_high, EEG)) ;
    dwPLI_low = ft_connectivityanalysis(struct('method', 'wpli_debiased'), ...
        ft_freqanalysis(cfg_TF_low, EEG)) ;
    dwPLI_all = ft_connectivityanalysis(struct('method', 'wpli_debiased'), ...
        ft_freqanalysis(cfg_TF_all, EEG)) ;

    % AVERAGE dwPLI OVER FREQUENCY BANDS
    for currband=1:size(bands,1)
        cfg_band.frequency = [str2double(bands{currband,2}) str2double(bands{currband,3})] ;
        if all(cfg_band.frequency >= 30)
            FC_mat{currfile, currband} = ft_selectdata(cfg_band, dwPLI_high).wpli_debiasedspctrm ;
            FC_ave{currfile, currband} = mean(FC_mat{currfile,currband}, ...
                'all', 'omitnan') ;
            FC_abs{currfile, currband} = mean(abs(FC_mat{currfile,currband}), ...
                'all', 'omitnan') ;
        elseif ~any(isnan(cfg_band.frequency))
            FC_mat{currfile, currband} = ft_selectdata(cfg_band, dwPLI_low).wpli_debiasedspctrm ;
            FC_ave{currfile, currband} = mean(FC_mat{currfile,currband}, ...
                'all', 'omitnan') ;
            FC_abs{currfile, currband} = mean(abs(FC_mat{currfile,currband}), ...
                'all', 'omitnan') ;
        else
            cfg_band.frequency = 'all' ;
            cfg_band.nanmean = 'yes' ;
            FC_mat{currfile, currband} = ft_selectdata(cfg_band, dwPLI_all).wpli_debiasedspctrm ;
        end
        if compSubset
            for i=1:size(subsetList,1)
                temp = FC_mat{currfile, currband} ;
                FC_subsets{currfile,i,currband} = mean(temp(cell2mat(setIndxs{i,1}), ...
                    cell2mat(setIndxs{i,2})), 'all', 'omitnan') ;
            end
        end
    end
end

%% PRINT OUT TABLE
cd([srcDir filesep 'generate_dwPLI']) ;
if customBands
    bandNames = cell(1, size(bands,1)) ;
    for i=1:stop
        bandNames{i} = [bands{i,1} ' (' bands{i,2} '-' bands{i,3} ')'] ;
    end
    if calcAll; bandNames{end} = 'All' ; end
else; bandNames = bands(:,1)' ;
end

if compSubset
    subNames = cell(1, size(subsetList,1)) ;
    for i=1:size(subsetList,1)
        set1 = subsetList{i,1} ;
        set2 = subsetList{i,2} ;
        str = [sprintf('%s, ', set1{1:end-1}) set1{end} ' x ' ...
            sprintf('%s, ', set2{1:end-1}) set2{end}] ;
        if length(str) >= 30; str = [str(1:27) '...'] ; end
        subNames{i} = str ;
    end
end

% Print out Average dwPLI values
fc_ave_savename = helpName(['generatedwPLI_AVE' suffix '_' datestr(now, 'dd-mm-yyyy') ...
    '.csv']) ;
writetable(cell2table(FC_ave, 'VariableNames', bandNames, 'RowNames', ...
    FileNames), fc_ave_savename, 'WriteRowNames', true, 'QuoteStrings', true) ;

% Print out Absolute Value Average dwPLI values
fc_abs_savename = helpName(['generatedwPLI_ABS' suffix '_' datestr(now, 'dd-mm-yyyy') ...
    '.csv']) ;
writetable(cell2table(FC_abs, 'VariableNames', bandNames, 'RowNames', ...
    FileNames), fc_abs_savename, 'WriteRowNames', true, 'QuoteStrings', true) ;

% Export the matrices and subsets (if applicable)
if exportCSV
    for currfile=1:size(FileNames,2)
        for currband=1:size(bands,1)
            saveName = helpName(strrep(FileNames{currfile}, ...
                ['_processed' suffix '.set'], ['_dwPLI_' bandNames{currband} ...
                '_' datestr(now, 'dd-mm-yyyy') '.csv']), '.csv') ;
            writetable(array2table(FC_mat{currfile, currband}, 'VariableNames', ...
                varNames{currfile}, 'RowNames', varNames{currfile}), saveName, ...
                'WriteRowNames', true, 'QuoteStrings', true) ;
         end
    end

    if compSubset
        for currband=1:size(bands,1)
            subsetSaveName = helpName(['dwPLI_' bandNames{currband} ...
                '_subsets_' datestr(now, 'dd-mm-yyyy') '.csv'], '.csv') ;
            writetable(array2table(squeeze(FC_subsets(:, :, currband)), 'VariableNames', ...
                subNames, 'RowNames', FileNames), subsetSaveName, ...
                'WriteRowNames', true, 'QuoteStrings', true) ;
        end
    end
else
    if compSubset
        subsetSaveName = helpName(['dwPLI_subsets_' datestr(now, ...
            'dd-mm-yyyy') '.xlsx'], '.xlsx') ;
    end
    for currfile=1:size(FileNames,2)
        saveName = helpName(strrep(FileNames{currfile}, ['_processed' ...
            suffix '.set'], ['_dwPLI_' datestr(now, 'dd-mm-yyyy') ...
            '.xlsx']), '.xlsx') ;
        for currband=1:size(bands,1)
            writetable(array2table(FC_mat{currfile, currband}, 'VariableNames', ...
                varNames{currfile}, 'RowNames', varNames{currfile}), ...
                saveName, 'WriteRowNames', true, 'Sheet', bandNames{currband}) ;
            if compSubset
                writetable(array2table(squeeze(FC_subsets(:, :, currband)), ...
                    'VariableNames', subNames, 'RowNames', FileNames), ...
                    subsetSaveName, 'WriteRowNames', true, 'Sheet', bandNames{currband}) ;
            end
        end
    end
end