% inEEG = EEG data matrix chans x samps
% srate = EEG sampling rate
% EEGchan_inc = whether there is an ECG chan in the data
% ECGchan_ID = the ECG chan ID OR the channels in which ECG is present
% peakWinSize = window size for peak detection in fraction of a second (.3 for 300 ms)

function outEEG = ECGone(inEEG, params)
%[outEEG, template, peakIndxs] = ECGone(inEEG, params, varargin)
%% VALIDATE VARARGIN
% vis = 0 ;
tolerance = 10 ;
zsqCutOff = 1 ;

%% DETERMINE ECG CHANNEL/CHANNELS WITH ARTIFACT
% Get THE index/indicies for ECG based on the EEG's channel locations and
% user-specified channel IDs
chanLabels = {inEEG.chanlocs.labels} ;
chanIndxs = [] ;
if ~params.ECGchan.inc; stop = size(params.ECGchan.ID, 2);
else; stop = 1;
end
for i=1:stop
    chanIndxs = [chanIndxs find(strcmpi(chanLabels, params.ECGchan.ID{i}))] ;  %#ok<AGROW> 
end

%% PREPARE EEG FOR TEMPLATE CREATION
% If the ECG channel is not included, filter the EEG to enable proxy
% creation
if ~params.ECGchan.inc
    EEGfilt = eeg_regepochs(pop_eegfiltnew(inEEG, [], 100, [], 0, [], 0), ...
        'recurrence', params.procEpoch, 'limits', [0 params.procEpoch], ...
        'rmbase', NaN) ;
end

% Segment the EEG into epochs with a length defined by the user (30 or 60
% seconds recommended)
EEG = eeg_regepochs(inEEG, 'recurrence', params.procEpoch, 'limits', ...
    [0 params.procEpoch], 'rmbase', NaN) ;

%% PROCCESS EACH EPOCH...
for currEpoch=1:size(EEG.data,3)
    %% EXTRACT CURRENT EPOCH FROM EEG
    epochData = squeeze(EEG.data(:,:,currEpoch)) ;
    cleanedEpoch = epochData ;
    
    %% GATHER ECG CHANNEL
    % If the ECG channel is included in the recording, simply pull it out
    % from the data. If the channel is not included, create a proxy using
    % channels designated by the user to contain ECG artifact.
    if params.ECGchan.inc; ECGchan = epochData(chanIndxs,:) ;
    else
        % Use the filtered data to emphasize the peaks through
        % multiplication
        ECGchan = squeeze(EEGfilt.data(:,:,currEpoch)).*squeeze(EEGfilt.data(:,:,currEpoch)) ;
        % Set ECGchan to be the average of all all the non-outlier
        % artifact-laden channels
        ECGchan = mean(ECGchan(chanIndxs(~isnan(nanOutliers(mean(ECGchan(chanIndxs, ...
            :), 2, 'omitnan')))), :), 1, 'omitnan') ;
    end
    
    %% FIND PEAKS IN THE ECG CHANNEL
    [peaks, peakIndxs] = findpeaks(ECGchan, 'MINPEAKDISTANCE', ...
        round(inEEG.srate*params.peakWinSize)) ;

    % Ensure the peak is not on the edges
    peakIndxs = peakIndxs(2:end-1) ;

    % If approximating the channel from artifact, confirm the data is
    % sufficiently "peaky". If not, end the run via error.
    if ~params.ECGchan.inc
%         outlierIndxs = find(peaks > (1.5*iqr(peaks) + median(peaks, 'omitnan')) ...
%             | peaks < (-1.5*iqr(peaks) + median(peaks, 'omitnan'))) ;

        if median((peaks - median(ECGchan))./median(ECGchan)) < params.peaky
            fprintf('Insufficient ECG artifact to create template.\n') ;
            continue ;
        end
    end

    %% PLOT THE ECG CHANNEL AND ITS PEAKS, if visualizations are enabled
%     if vis % ***
%         figure ;
%         hold off ;
%         plot(ECGchan) ;
%         hold on ;
%         plot(peakIndxs, ECGchan(peakIndxs), 'r*') ;
%     end

    %% CREATE TEMPLATE
    % Calculate the number of samples in 120 milliseconds (the estimated
    % window for an ECG signal)
    numSamps = floor(0.12*EEG.srate) ;
    % Initialize matricies to hold the template
    tempTemplate = zeros(numSamps, size(epochData, 1)) ;
    template = zeros(length(peakIndxs), numSamps, size(epochData, 1)) ;
    
    % For each peak...
    for currPeak=1:size(peakIndxs,2)
        % Create a window of indicies centered around the peak using the 
        % number of samples per 120 milliseconds
        peakWindIndxs = (peakIndxs(currPeak)-numSamps/2):(peakIndxs(currPeak)-1+numSamps/2) ;

        % For each channel, collect the values in the generated window
        for currChan=1:size(epochData, 1)
            template(currPeak,:,currChan) = epochData(currChan, peakWindIndxs) ;
        end
    end

    % Find the max absolute value for each peak (rows) for each channel
    % (cols)
    chanMax = squeeze(max(abs(template), [], 2, 'omitnan')) ;

    % For each channel...
    for currChan=1:size(epochData, 1)
        % Set each window with an outlier max value to NaN in the template
        template(isnan(nanOutliers(chanMax(:,currChan))),:,currChan) = NaN ;
        % Calculate the mean values across peaks for the current channel
        tempTemplate(:,currChan) = squeeze(mean(template(:,:,currChan),1, 'omitnan')) ;
    end
    template = tempTemplate ; % rows = each point in window, cols = channels

    %% SET BAD CHANNELS TO 0
%     template(:,badChansInds) = 0 ;

    %% ENSURE THE PEAKS ARE CENTERED WITH A USER-DEFINED TOLERANCE
    % Get the index of the max value for each channel, accounting for
    % multiple points with the same value (e.g., a flat peak across two
    % points)
    indxMax = nan(size(epochData, 1),1) ;
    for currChan = 1:size(epochData, 1)
        indxMax(currChan) = median(find(abs(template(:,currChan)) == ...
            max(abs(template(:,currChan)))), 'omitnan') ;                  
    end
    indxMax = round(indxMax) ;

    % Ensure that the peak is within some tolerance of the midpoint. If
    % not, set that column to zeros.
    template(:, (indxMax < (numSamps/2)-tolerance) | ...
        (indxMax > (numSamps/2)+tolerance)) = 0 ;

    %% CONVERT THE TEMPLATE TO A HANNING WINDOW
    % Create a hanning window using the number of samples in 120
    % milliseconds
    wind = hann(numSamps) ;
    % For each channel subtract the mean of the channel from each value 
    % and multiply it by the hanning window - What is the purpose of this?
    for currChan = 1:size(epochData, 1)
        template(:, currChan) = (template(:, currChan) - ...
            mean(template(:,currChan), 'omitnan')).*wind ;
    end
    
    %% COMPLETE A SECOND PASS, if signal was approximated from artifact
    if ~params.ECGchan.inc
        % CALCULATE RESIDUALS
        % Need to build in safety for if +/- tolerance is out of bounds
        resids = zeros(length(peakIndxs), size(epochData,1)) ;
        for currIndx = 1:length(peakIndxs)
            resids(currIndx,:) = var(zscore(epochData(:, peakIndxs(currIndx) - ...
                tolerance:peakIndxs(currIndx)+tolerance)') - ...
                zscore(template(round(0.5*numSamps)-tolerance: ...
                round(0.5*numSamps)+tolerance, :))) ;
        end
        resids = resids' ;

        % PLOT RESIDUALS, if visualizations are enabled
%         if vis
%             figure ;
%             hold off;
%             imagesc(resids, [0, 100]) ;
%         end

        % REMOVE OUTLIER INDICIES using residuals and z-square cutoffs
        peakIndxs(median(resids) > zsqCutOff) = [] ;                        %#ok<UDIM> 

        % END THE RUN IF ALL INDICIES ARE REMOVED
        if isempty(peakIndxs)
            fprintf('All incidies removed.\n') ;
            continue ;
        end

        % PLOT THE ECG CHANNEL AND ITS PEAKS, if visualizations are enabled
%         if vis
%             figure ;
%             hold off ;
%             plot(ECGchan) ;
%             hold on ;
%             plot(peakIndxs, ECGchan(peakIndxs), 'r*') ;
%         end

        % REMAKE TEMPLATE EXCUDING REMOVED INDICIES
        template = zeros(length(peakIndxs), numSamps, size(epochData, 1)) ;

        for currPeak=1:length(peakIndxs)
            peakWindIndxs = (peakIndxs(currPeak)-numSamps/2):(peakIndxs(currPeak)-1+numSamps/2) ;
            for i=1:size(epochData,1)
                template(currPeak,:,i) = epochData(i, peakWindIndxs) ;
            end
        end

        if size(template,1) > 1; template = squeeze(mean(template,'omitnan')) ;
        else; template = squeeze(template) ;
        end

        for i=1:size(epochData,1)
            template(:,i) = (template(:,i) - mean(template(:,i), ...
                'omitnan')).*wind ;
        end
%         Remove the bad channels from the template
%         template(:, global_badchs) = 0 ;
    end

    % SUBTRACT THE TEMPLATE FROM THE DATA MATRIX
    for currPeak = 1:length(peakIndxs)
        peakWindIndxs = (peakIndxs(currPeak)-numSamps/2):(peakIndxs(currPeak)-1+numSamps/2) ;
        cleanedEpoch(:, peakWindIndxs) = epochData(:, peakWindIndxs) - template' ;
    end

%     % Plot the EEG before and after
%     if vis
%     end
%     % Plot the difference between before and after
%     if vis
%     end

    EEG.data(:,:,currEpoch) = cleanedEpoch ;
end
outEEG = inEEG ;
outEEG.data = reshape(EEG.data, size(inEEG.data,1), []) ;
end

function x = nanOutliers(x)
    x(x > (1.5*iqr(x) + median(x, 'omitnan')) | x < (-1.5*iqr(x) + ...
        median(x, 'omitnan'))) = NaN ;
end