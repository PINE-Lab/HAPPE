function [EEG, command] = loadcurryTwoOne(fullfilename, varargin)
%   Import a Neuroscan Curry file into EEGLAB. Currently supports Curry6,
%   Curry7, and Curry8 data files (both continuous and epoched). Epoched
%   datasets are loaded in as continuous files with boundary events. Data
%   can be re-epoched using EEGLAB/ERPLAB functions.
%
%   Input Parameters:
%        1    Specify the filename of the Curry file (extension should be either .cdt, .dap, .dat, or .rs3). 
%
%   Example Code:
%
%       >> EEG = pop_loadcurry;   % an interactive uigetfile window
%       >> EEG = loadcurry;   % an interactive uigetfile window
%       >> EEG = loadcurry('C:\Studies\File1.cdt');    % no pop-up window 
%
%   Optional Parameters:
%       1     'CurryLocations' - Boolean parameter to determine if the sensor
%               locations are carried forward from Curry [1, 'True'] or if the channel
%               locations from EEGLAB should be used [0, 'False' Default].
%
%   Author for reading into Matlab: Neuroscan 
%   Author for translating to EEGLAB: Matthew B. Pontifex, Health Behaviors and Cognition Laboratory, Michigan State University, August 26, 2015
%
%
%   revision 2.1 - 
%     Updated to make sure event latencies are in double format.
%
%   revision 2.0 - 
%     Updated for Curry8 compatibility and compatibility with epoched
%     datasets. Note that Curry only carries forward trigger events used in
%     the epoching process.
%
%   revision 1.3 -
%     Revised to be backward compatible through r2010a - older versions may work but have not been tested.
%
%   revision 1.2 -
%     Fixed an issue related to validating the trigger markers.
%
%   revision 1.1 - 
%     Fixed a problem with reading files in older versions of matlab.
%     Added import of impedance check information for the most recent check
%        as well as the median of the last 10 checks. Data is available in
%        EEG.chanlocs
%     Moved function to loadcurry() and setup pop_loadcurry() as the pop up shell. 
%     Created catch for user cancelling file selection dialog. Now throws
%        an error to not overwrite what is in the EEG variable
%
%   If there is an error with this code, please email pontifex@msu.edu with the issue and I'll see what I can do.


    command = '';
    if nargin < 1 % No file was identified in the call
        try
            % flip to pop_loadcurry()
            [EEG, command] = pop_loadcurry();
        catch
            % only error that should occur is user cancelling prompt
            error('loadcurry(): File selection cancelled. Error thrown to avoid overwriting data in EEG.')
        end
    else
        
        if ~isempty(varargin)
             r=struct(varargin{:});
        end
        try, r.CurryLocations; catch, r.CurryLocations = 'False'; end
        try, r.Force; catch, r.Force = 'False'; end
        try, r.AltFile; catch, r.AltFile = 'False'; end
        if strcmpi(r.CurryLocations, 'True') | (r.CurryLocations == 0)
            r.CurryLocations = 'True';
        end

        EEG = [];
        EEG = eeg_emptyset;
        [pathstr,name,ext] = fileparts(fullfilename);
        filename = [name,ext];
        filepath = [pathstr, filesep];
        file = [pathstr, filesep, name];

        % Ensure the appropriate file types exist under that name
        boolfiles = 1;
        curryvers = 0;
        if strcmpi(ext, '.cdt')
            curryvers = 8;
            if (exist([file '.cdt'], 'file') == 0) || (exist([file '.cdt.dpa'], 'file') == 0)
                boolfiles = 0;
                error('Error in pop_loadcurry(): The requested filename "%s" in "%s" does not have both file components (.cdt, .cdt.dpa) created by Curry 8.', name, filepath)
            end
        else
            curryvers = 7;
            if (exist([file '.dap'], 'file') == 0) || (exist([file '.dat'], 'file') == 0) || (exist([file '.rs3'], 'file') == 0)
                boolfiles = 0;
                error('Error in pop_loadcurry(): The requested filename "%s" in "%s" does not have all three file components (.dap, .dat, .rs3) created by Curry 6 and 7.', name, filepath)
            end
        end

        if (boolfiles == 1)

            %% Provided by Neuroscan enclosed within Program Files Folder for Curry7 (likely to be included in Curry8)
            % Received updated version on 6-19-2016 from Michael Wagner, Ph.D., Senior Scientist, Compumedics Germany GmbH, Heußweg 25, 20255 Hamburg, Germany
            % Modified to retain compatibility with earlier versions of Matlab and Older Computers by Pontifex
            
            if (curryvers == 7)
                datafileextension = '.dap';
            elseif (curryvers == 8)
                datafileextension = '.cdt.dpa';
            end
            
            % Open parameter file
            fid = fopen([file, datafileextension],'rt');
            if (fid == -1)
               error('Error in loadcurry(): Unable to open file.') 
            end
            try
                cell = textscan(fid,'%s','whitespace','','endofline','§');
            catch
                % In case of earlier versions of Matlab or Older Computers
                fclose(fid); 
                fid = fopen([file, datafileextension],'rt');
                f = dir([file, datafileextension]);
                try
                    cell = textscan(fid,'%s','whitespace','','endofline','§','BufSize',round(f.bytes+(f.bytes*0.2)));
                catch
                    fclose(fid); 
                    fid = fopen([file, datafileextension],'rt');
                    cell = textscan(fid,'%s','whitespace','','BufSize',round(f.bytes+(f.bytes*0.2)));
                end
            end
            fclose(fid);            
            cont = cell2mat(cell{1});

            % read parameters from file
            % tokens (second line is for Curry 6 notation)
            tok = { 'NumSamples'; 'NumChannels'; 'NumTrials'; 'SampleFreqHz';  'TriggerOffsetUsec';  'DataFormat'; 'DataSampOrder';   'SampleTimeUsec'; 
                    'NUM_SAMPLES';'NUM_CHANNELS';'NUM_TRIALS';'SAMPLE_FREQ_HZ';'TRIGGER_OFFSET_USEC';'DATA_FORMAT';'DATA_SAMP_ORDER'; 'SAMPLE_TIME_USEC' };

            % scan in cell 1 for keywords - all keywords must exist!
            nt = size(tok,1);
            a = zeros(nt,1);
            for i = 1:nt
                 ctok = tok{i,1};
                 ix = strfind(cont,ctok);
                 if ~isempty ( ix )
                     text = sscanf(cont(ix+numel(ctok):end),' = %s');     % skip =
                     if strcmp ( text,'ASCII' ) || strcmp ( text,'CHAN' ) % test for alphanumeric values
                         a(i) = 1;
                     else 
                         c = sscanf(text,'%f');         % try to read a number
                         if ~isempty ( c )
                             a(i) = c;                  % assign if it was a number
                         end
                     end
                 end 
            end

            % derived variables. numbers (1) (2) etc are the token numbers
            nSamples    = a(1)+a(1+nt/2);
            nChannels   = a(2)+a(2+nt/2);
            nTrials     = a(3)+a(3+nt/2);
            fFrequency  = a(4)+a(4+nt/2);
            fOffsetUsec = a(5)+a(5+nt/2);
            nASCII      = a(6)+a(6+nt/2);
            nMultiplex  = a(7)+a(7+nt/2);
            fSampleTime = a(8)+a(8+nt/2);

            if ( fFrequency == 0 && fSampleTime ~= 0 )
                fFrequency = 1000000 / fSampleTime;
            end   
            
            %Search for Impedance Values
            tixstar = strfind(cont,'IMPEDANCE_VALUES START_LIST');
            tixstop = strfind(cont,'IMPEDANCE_VALUES END_LIST');

            impedancelist = []; 
            impedancematrix = [];

            if (~isempty(tixstar)) && (~isempty(tixstop))
                text = cont(tixstar:tixstop-1);
                tcell = textscan(text,'%s');
                tcell = tcell{1,1};
                for tcC = 1:size(tcell,1)
                   tcell{tcC} = str2num(tcell{tcC}); % data was read in as strings - force to numbers
                   if ~isempty(tcell{tcC}) % skip if it is not a number
                       impedancelist(end+1) = tcell{tcC};
                   end
                end

                % Curry records last 10 impedances
                impedancematrix = reshape(impedancelist,[(size(impedancelist,2)/10),10])';
                impedancematrix(impedancematrix == -1) = NaN; % screen for missing
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % open file containing labels
            if (curryvers == 7)
                datafileextension = '.rs3';
            elseif (curryvers == 8)
                datafileextension = '.cdt.dpa';
            end
            
            fid = fopen([file, datafileextension],'rt');
            if (fid == -1)
               error('Error in loadcurry(): Unable to open file.') 
            end
            try
                cell = textscan(fid,'%s','whitespace','','endofline','§');
            catch
                fclose(fid);
                fid = fopen([file, datafileextension],'rt');
                f = dir([file, datafileextension]);
                try
                    cell = textscan(fid,'%s','whitespace','','endofline','§','BufSize',round(f.bytes+(f.bytes*0.2)));
                catch
                    fclose(fid);
                    fid = fopen([file, datafileextension],'rt');
                    cell = textscan(fid,'%s','whitespace','','BufSize',round(f.bytes+(f.bytes*0.2)));
                end
            end
            fclose(fid);
            cont = cell2mat(cell{1});

            % read labels from rs3 file
            % initialize labels
            labels = num2cell(1:nChannels);

            for i = 1:nChannels
                text = sprintf('EEG%d',i);
                labels(i) = cellstr(text);
            end
                        
            % scan in cell 1 for LABELS (occurs four times per channel group)
            ix = strfind(cont,[char(10),'LABELS']);
            nt = size(ix,2);
            nc = 0;
            
            for i = 4:4:nt                                                      % loop over channel groups
                newlines = ix(i-1) + strfind(cont(ix(i-1)+1:ix(i)),char(10));   % newline
                last = nChannels - nc;
                for j = 1:min(last,size(newlines,2)-1)                          % loop over labels
                    text = cont(newlines(j)+1:newlines(j+1)-1);
                    if isempty(strfind(text,'END_LIST'))
                        nc = nc + 1;
                        labels(nc) = cellstr(text);
                    else 
                        break
                    end
                end 
            end

            % read sensor locations from rs3 file
            % initialize sensor locations
            sensorpos = zeros(3,0);

            % scan in cell 1 for SENSORS (occurs four times per channel group)
            ix = strfind(cont,[char(10),'SENSORS']);
            nt = size(ix,2);
            nc = 0;

            for i = 4:4:nt                                                      % loop over channel groups
                newlines = ix(i-1) + strfind(cont(ix(i-1)+1:ix(i)),char(10));   % newline
                last = nChannels - nc;
                for j = 1:min(last,size(newlines,2)-1)                          % loop over labels
                    text = cont(newlines(j)+1:newlines(j+1)-1);
                    if isempty(strfind(text,'END_LIST'))
                        nc = nc + 1;
                        tcell = textscan(text,'%f');                           
                        posx = tcell{1}(1);
                        posy = tcell{1}(2);
                        posz = tcell{1}(3);
                        sensorpos = cat ( 2, sensorpos, [ posx; posy; posz ] );
                    else 
                        break
                    end
                end 
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % read events from cef/ceo file
            % initialize events
            ne = 0;                                                             % number of events
            events = zeros(4,0);
            annotations = cellstr('empty');

             % open file containing labels
            if (curryvers == 7)
                datafileextension = '.cef';
                datafileextensionalt = '.ceo';
            elseif (curryvers == 8)
                datafileextension = '.cdt.cef';
                datafileextensionalt = '.cdt.ceo';
            end
            
            % find appropriate file
            fid = fopen([file, datafileextension],'rt');
            if fid < 0
                fid = fopen([file, datafileextensionalt],'rt');
                f = dir([file, datafileextensionalt]);
            else
                f = dir([file, datafileextension]);
            end

            if fid >= 0                
                try
                    cell = textscan(fid,'%s','whitespace','','endofline','§');
                catch
                    fclose(fid);
                    fid = fopen([file, datafileextension],'rt');
                    if fid < 0
                        fid = fopen([file, datafileextensionalt],'rt');
                        f = dir([file, datafileextensionalt]);
                    else
                        f = dir([file, datafileextension]);
                    end
                    try
                        cell = textscan(fid,'%s','whitespace','','endofline','§','BufSize',round(f.bytes+(f.bytes*0.2)));
                    catch
                        fclose(fid);
                        fid = fopen([file, datafileextension],'rt');
                        if fid < 0
                            fid = fopen([file, datafileextensionalt],'rt');
                            f = dir([file, datafileextensionalt]);
                        else
                            f = dir([file, datafileextension]);
                        end
                        cell = textscan(fid,'%s','whitespace','','BufSize',round(f.bytes+(f.bytes*0.2)));
                    end
                end
                fclose(fid);
                cont = cell2mat(cell{1});

                % scan in cell 1 for NUMBER_LIST (occurs five times)
                ix = strfind(cont,'NUMBER_LIST');

                newlines = ix(4) - 1 + strfind(cont(ix(4):ix(5)),char(10));     % newline
                last = size(newlines,2)-1;
                for j = 1:last                                                  % loop over labels
                    text = cont(newlines(j)+1:newlines(j+1)-1);
                    tcell = textscan(text,'%d');                           
                    sample = tcell{1}(1);                                       % access more content using different columns
                    type = tcell{1}(3);
                    startsample = tcell{1}(5);
                    endsample = tcell{1}(6);
                    ne = ne + 1;
                    events = cat ( 2, events, [ sample; type; startsample; endsample ] );
                end

                % scan in cell 1 for REMARK_LIST (occurs five times)
                ix = strfind(cont,'REMARK_LIST');
                na = 0;

                newlines = ix(4) - 1 + strfind(cont(ix(4):ix(5)),char(10));     % newline
                last = size(newlines,2)-1;
                for j = 1:last                                                  % loop over labels
                    text = cont(newlines(j)+1:newlines(j+1)-1);
                    na = na + 1;
                    annotations(na) = cellstr(text);
                end    
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % read dat file
            if (curryvers == 7)
                datafileextension = '.dat';
            elseif (curryvers == 8)
                datafileextension = '.cdt';
            end
            
            if nASCII == 1
                fid = fopen([file, datafileextension],'rt');
                if (fid == -1)
                   error('Error in loadcurry(): Unable to open file.') 
                end
                f = dir([file, datafileextension]);
                try
                    fclose(fid);
                    fid = fopen([file, datafileextension],'rt');
                    cell = textscan(fid,'%f',nChannels*nSamples*nTrials);
                catch
                    fclose(fid);
                    fid = fopen([file, datafileextension],'rt');
                    cell = textscan(fid,'%f',nChannels*nSamples*nTrials, 'BufSize',round(f.bytes+(f.bytes*0.2)));
                end
                fclose(fid);
                data = reshape([cell{1}],nChannels,nSamples*nTrials);
            else
                fid = fopen([file, datafileextension],'rb');
                if (fid == -1)
                   error('Error in loadcurry(): Unable to open file.') 
                end
                data = fread(fid,[nChannels,nSamples*nTrials],'float32');
                fclose(fid);
            end

            % transpose?
            if nMultiplex == 1
                data = data';
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % time axis
            time = linspace(fOffsetUsec/1000,fOffsetUsec/1000+(nSamples*nTrials-1)*1000/fFrequency,nSamples*nTrials);

            %% Created to take this data and place it into EEGLAB format (v13.4.4b)
            
            EEG.setname = 'Neuroscan Curry file';
            if (curryvers == 7)
                datafileextension = '.dap';
            elseif (curryvers == 8)
                datafileextension = '.cdt';
            end
            EEG.filename = [name, datafileextension];
            EEG.filepath = filepath;
            EEG.comments = sprintf('Original file: %s%s', filepath, [name, datafileextension]);
            EEG.ref = 'Common';
            EEG.trials = nTrials;
            EEG.pnts = nSamples;
            EEG.srate = fFrequency;
            EEG.times = time;
            EEG.data = double(data);
            EEG.xmin = 0;
            EEG.xmax = (EEG.pnts-1)/EEG.srate+EEG.xmin;
            EEG.nbchan = size(EEG.data,1);
            EEG.urchanlocs = [];
            EEG.chaninfo.plotrad = [];
            EEG.chaninfo.shrink = [];
            EEG.chaninfo.nosedir = '+X';
            EEG.chaninfo.nodatchans = [];
            EEG.chaninfo.icachansind = [];
            
            % Populate channel labels
            EEG.chanlocs = struct('labels', [], 'ref', [], 'theta', [], 'radius', [], 'X', [], 'Y', [], 'Z', [],'sph_theta', [], 'sph_phi', [], 'sph_radius', [], 'type', [], 'urchan', []);
            for cC = 1:(numel(labels))
                EEG.chanlocs(cC).labels = char(upper(labels(cC))); % Convert labels to uppercase and store as character array string
                EEG.chanlocs(cC).urchan = cC;
            end

            if strcmpi(r.CurryLocations, 'True')
                % Populate channel locations
                % LPS sensor system:
                % from right towards left, 
                % from anterior towards posterior, 
                % from inferior towards superior

                % MATLAB/EEGLAB system:
                % x is towards the nose, 
                % y is towards the left ear, 
                % z towards the vertex.

                try, sensorpos; booler = 0; catch; booler = 1; end
                if (booler == 0)
                    for cC = 1:size(sensorpos,2)
                       EEG.chanlocs(cC).Y = sensorpos(1,cC); 
                       EEG.chanlocs(cC).X = sensorpos(2,cC)*-1; 
                       EEG.chanlocs(cC).Z = sensorpos(3,cC); 
                    end
                    % Populate other systems based upon these values
                    EEG.chanlocs = convertlocs( EEG.chanlocs, 'auto');
                end
                EEG.history = sprintf('%s\nEEG = loadcurry(''%s%s'', ''CurryLocations'', ''True'');', EEG.history, filepath, [name, datafileextension]); 
            else
                % Use default channel locations
                try
                    tempEEG = EEG; % for dipfitdefs
                    dipfitdefs;
                    tmpp = which('eeglab.m');
                    tmpp = fullfile(fileparts(tmpp), 'functions', 'resources', 'Standard-10-5-Cap385_witheog.elp');
                    userdatatmp = { template_models(1).chanfile template_models(2).chanfile  tmpp };
                    try
                        [T, tempEEG] = evalc('pop_chanedit(tempEEG, ''lookup'', userdatatmp{1})');
                    catch
                        try
                            [T, tempEEG] = evalc('pop_chanedit(tempEEG, ''lookup'', userdatatmp{3})');
                        catch
                            booler = 1;
                        end
                    end
                    EEG.chanlocs = tempEEG.chanlocs;
                catch
                    booler = 1;
                end
                EEG.history = sprintf('%s\nEEG = loadcurry(''%s%s'');', EEG.history, filepath, [name, datafileextension]);
            end
            
            % Place impedance values within the chanlocs structure
            try
                if ~isempty(impedancematrix)
                    if (size(impedancematrix,2) == size({EEG.chanlocs.labels},2)) % number of channels matches number of impedances
                        impedancematrix(impedancematrix == -1) = NaN; % screen for missing values
                        impedancelist = nanmedian(impedancematrix);
                        for cC = 1:size({EEG.chanlocs.labels},2)
                           EEG.chanlocs(cC).impedance = impedancematrix(1,cC)/1000; 
                           EEG.chanlocs(cC).median_impedance = impedancelist(1,cC)/1000;
                        end
                    end
                end
            catch
                booler = 1;
            end
            
            % Populate Event List
            if ~isempty(events)
                EEG.event = struct('type', [], 'latency', [], 'urevent', []);
                EEG.urevent = struct('type', [], 'latency', []);
                for cC = 1:size(events,2)
                    EEG.event(cC).urevent = cC;
                    EEG.event(cC).type = events(2,cC);
                    EEG.event(cC).latency = double(events(1,cC));
                    EEG.urevent(cC).type = EEG.event(cC).type;
                    EEG.urevent(cC).latency = double(EEG.event(cC).latency);
                end
            else
                % Event list is empty
                % Determine if Triggers are present
                if ~isempty(find(strcmpi(labels,'Trigger')))

                    % Remove baseline from trigger channel
                    EEG.data(find(strcmpi(labels,'Trigger')),:) = EEG.data(find(strcmpi(labels,'Trigger')),:)-EEG.data(find(strcmpi(labels,'Trigger')),1); 

                    % Populate list based on values above 0, triggers may last more than one sample
                    templat = find(EEG.data(find(strcmpi(labels,'Trigger')),:)>0);
                    templatrem = [];
                    for cC = 2:numel(templat)
                        % If the sampling point is one off
                        if ((templat(cC)-1) == templat(cC-1))
                           templatrem(end+1) = templat(cC);
                        end
                    end
                    templat = setdiff(templat,templatrem);
                    if ~isempty(templat)
                        EEG.event = struct('type', [], 'latency', [], 'urevent', []);
                        EEG.urevent = struct('type', [], 'latency', []);
                        % Populate event list
                        for cC = 1:numel(templat)
                            try
                                EEG.event(cC).urevent = cC;
                                EEG.event(cC).type = EEG.data(find(strcmpi(labels,'Trigger')),templat(cC));
                                EEG.urevent(cC).type = EEG.event(cC).type;
                                EEG.event(cC).latency = double((templat(cC)-1));
                                EEG.urevent(cC).latency = double(EEG.event(cC).latency);
                            catch
                                continue
                            end
                        end
                    end
                end
            end
            
            % Remove Trigger Channel
%             if ~isempty(find(strcmpi(labels,'Trigger')))
%                 EEG.data(find(strcmpi(labels,'TRIGGER')),:) = [];
%                 EEG.chanlocs(find(strcmpi(labels,'TRIGGER'))) = [];
%             end
            EEG.nbchan = size(EEG.data,1);

            % Handle Epoched Datasets
            if (nTrials > 1)
                % Data has epochs
                [T, EEG] = evalc('eeg_epoch2continuous(EEG)'); % Force to continous and add boundary events
                
            else
                % Data is continuous
            
                % Check to See if Behavioral Data is Available
                % Neurscan STIM2 Data is stored with the same file type as the  EEG data in Curry 6 & 7 which is problematic
                % Check for PsychoPy .PSYDAT file
                % 'Trial','Event','Duration','ISI','ITI','Type','Resp','Correct','Latency','ClockLatency','Trigger','MinRespWin','MaxRespWin','Stimulus'
                AltFile = [file '.psydat'];
%                 if ~strcmpi(r.AltFile, 'False')
%                     AltFile = r.AltFile;
%                 end
%                 if ~(exist(AltFile, 'file') == 0) 
% 
%                     fid = fopen(AltFile,'rt');
%                     if (fid ~= -1)
%                         cell = textscan(fid,'%s');
%                         fclose(fid);
%                         cont = cell{1};
% 
%                         % Check file version
%                         if strcmpi(cont(1,1),'gentask.....=') % Could be Neuroscan Stim2 or modified Psychopy formats
%                             if strcmpi(cont(2,1),'PsychoPy_Engine_3')
%                                 EEG = importengine3psychopy(EEG, AltFile, 'Force', r.Force);
%                             end
%                         end
%                     end
%                 end
            end
            
            EEG = eeg_checkset(EEG);
            EEG.history = sprintf('%s\nEEG = eeg_checkset(EEG);', EEG.history);

        end
    end
   
end

