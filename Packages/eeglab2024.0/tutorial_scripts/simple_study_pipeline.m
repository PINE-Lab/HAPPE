% check folder
eeglab;
if ~exist('task-P300_events.json', 'file')
    error('Download the data from https://openneuro.org/datasets/ds003061/ and go to the downloaded folder');
else
    filepath = fileparts(which('task-P300_events.json'));
end

% import data
pop_editoptions( 'option_storedisk', 1); % only one dataset in memory at a time
[STUDY, ALLEEG] = pop_importbids(filepath, 'studyName','Oddball', 'subjects', [1:2]); % when using all subjects, one subject is truncated and cause the pipeline to return an error

% remove non-ALLEEG channels (it is also possible to process ALLEEG data with non-ALLEEG data
ALLEEG = pop_select( ALLEEG,'nochannel',{'EXG1','EXG2','EXG3','EXG4','EXG5','EXG6','EXG7','EXG8', 'GSR1', 'GSR2', 'Erg1', 'Erg2', 'Resp', 'Plet', 'Temp'});

% compute average reference
ALLEEG = pop_reref( ALLEEG,[]);

% clean data using the clean_rawdata plugin
ALLEEG = pop_clean_rawdata( ALLEEG,'FlatlineCriterion',5,'ChannelCriterion',0.87, ...
    'LineNoiseCriterion',4,'Highpass',[0.25 0.75] ,'BurstCriterion',20, ...
    'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian', ...
    'WindowCriterionTolerances',[-Inf 7] ,'fusechanrej',1);

% recompute average reference interpolating missing channels (and removing
% them again after average reference - STUDY functions handle them automatically)
ALLEEG = pop_reref( ALLEEG,[],'interpchan',[]);

% run ICA reducing the dimention by 1 to account for average reference 
plugin_askinstall('picard', 'picard', 1); % install Picard plugin
ALLEEG = pop_runica(ALLEEG, 'icatype','picard','concatcond','on','options',{'pca',-1});

% run ICLabel and flag artifactual components
ALLEEG = pop_iclabel(ALLEEG, 'default');
ALLEEG = pop_icflag( ALLEEG,[NaN NaN;0.9 1;0.9 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN]);

% extract data epochs
ALLEEG = pop_epoch( ALLEEG,{'oddball_with_reponse','standard'},[-1 2] ,'epochinfo','yes');
ALLEEG = eeg_checkset( ALLEEG );
ALLEEG = pop_rmbase( ALLEEG,[-1000 0] ,[]);

% create STUDY design
STUDY = std_maketrialinfo(STUDY, ALLEEG);
STUDY = std_makedesign(STUDY, ALLEEG, 1, 'name','STUDY.design 1','delfiles','off', ...
    'defaultdesign','off','variable1','type','values1',{'oddball_with_reponse','standard'},...
    'vartype1','categorical','subjselect', STUDY.subject);

% precompute ERPs at the STUDY level
[STUDY, ALLEEG] = std_precomp(STUDY, ALLEEG, {},'savetrials','on','rmicacomps','on','interp','on','recompute','on','erp','on');

% plot ERPS
STUDY = pop_erpparams(STUDY, 'topotime',350);
chanlocs = eeg_mergelocs(ALLEEG.chanlocs); % get all channels from all datasets
STUDY = std_erpplot(STUDY,ALLEEG,'channels', {chanlocs.labels}, 'design', 1);

% revert default option
pop_editoptions( 'option_storedisk', 0);
