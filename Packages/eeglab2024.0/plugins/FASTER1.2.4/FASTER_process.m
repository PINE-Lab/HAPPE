function EEG=FASTER_process(option_wrapper,log_file)

% Copyright (C) 2010 Hugh Nolan, Robert Whelan and Richard Reilly, Trinity College Dublin,
% Ireland
% nolanhu@tcd.ie, robert.whelan@tcd.ie
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

EEG=[];
try
    tic;
    o=option_wrapper.options;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % File options %
    %%%%%%%%%%%%%%%%
    % 1 File name including full path (string)
    % 2 Reference channel (integer > 0)
    % 3 Number of data channels (integer > 0)
    % 4 Number of extra channels (integer > 0)
    % 5 Channel locations file including full path (string)
    % 6 Save options (cell)
    %%%%%%%%%%%%%%%%

    using_ALLEEG=o.file_options.using_ALLEEG;
    prefix=o.file_options.file_prefix;
    %prefix_ALLEEG=o.file_options.prefix_ALLEEG;

    fullfilename = o.file_options.current_file;
    ref_chan = o.channel_options.ref_chan;
    eeg_chans = o.channel_options.eeg_chans;
    if (eeg_chans==0)
        eeg_chans=[];
    end
    ext_chans = o.channel_options.ext_chans;
    if (ext_chans==0)
        ext_chans=[];
    end
    channel_locations_file = o.file_options.channel_locations;
    save_options = o.save_options;
    cutoff_markers = o.file_options.cutoff_markers;
    do_reref = o.channel_options.do_reref;

    if (~do_reref)
        ref_chan=[];
    end
    [filepath,filename,extension] = fileparts(fullfilename);

    %log_file = fopen([filepath filesep filename '.log'],'a');

    c=clock;
    months={'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'};
    fprintf(log_file,'\n%d/%s/%d %d:%d:%d\n',c(3),months{c(2)},c(1),c(4),c(5),round(c(6)));
    fprintf(log_file,'%.2f - Opened log file.\n',toc);

    %%%%%%%%%%%%%%%%%%%%%%
    % File setup section %
    %%%%%%%%%%%%%%%%%%%%%%

    % Import .bdf file or load .set file
    % Note: import all channels and then remove the unnecessary ones, as
    % otherwise the event channel gets removed and we have no event data.
    if strcmpi(extension,'.bdf') && ~using_ALLEEG
        fprintf('Importing %s.\n',fullfilename);
        EEG = pop_biosig(fullfilename);
        EEG.setname = filename;
        EEG = pop_select(EEG, 'nochannel',length(eeg_chans)+length(ext_chans)+1:size(EEG.data,1));
        if (do_reref)
            if (max(EEG.data(ref_chan,:))==0 && min(EEG.data(ref_chan,:))==0)
                fprintf(log_file,'%.2f - Reference channel %d is already zeroed. Data was not re-referenced.\n',toc,ref_chan);
            elseif (o.ica_options.keep_ICA && ~isempty(EEG.icaweights))
                fprintf(log_file,'%.2f - Data was not re-referenced to maintain existing ICA weights. Bad channel detection may be ineffective.\n',toc,ref_chan);
            else
                EEG = h_pop_reref( EEG, ref_chan, 'exclude', ext_chans, 'keepref', 'on');
            end
        end

        filename = [o.file_options.file_prefix filename];
        filepath=o.file_options.oplist{o.file_options.current_file_num};
        mkdir([filepath filesep 'Intermediate']);
        EEG = pop_saveset(EEG,'filename',[filename '.set'],'filepath',filepath,'savemode','onefile');
        fprintf(log_file,'%.2f - Imported and converted file %s.\n',toc,fullfilename);
    elseif strcmpi(extension,'.set') && ~using_ALLEEG
        fprintf('Loading %s.\n',fullfilename);
        EEG = pop_loadset('filename',[filename '.set'],'filepath',filepath);
        fprintf(log_file,'%.2f - Loaded file %s.\n',toc,fullfilename);
        if ~isempty(o.file_options.output_folder_name)
            filepath=o.file_options.oplist{o.file_options.current_file_num};
            mkdir([filepath filesep 'Intermediate']);
        else
            filepath=o.file_options.oplist{o.file_options.current_file_num};            
            mkdir([filepath filesep 'Intermediate']);
            pop_saveset(EEG,'filename',['Original_' filename '.set'],'filepath',[filepath filesep 'Intermediate']);
            delete(fullfilename);
            if exist([fullfilename(1:end-4) '.fdt'],'file')
                delete([fullfilename(1:end-4) '.fdt']);
            end
            if exist([fullfilename(1:end-4) '.dat'],'file')
                delete([fullfilename(1:end-4) '.dat']);
            end
        end
        filename = [o.file_options.file_prefix filename];
        EEG.filename = [filename '.set'];
    elseif using_ALLEEG
        EEG=evalin('base',sprintf('ALLEEG(%d);',o.file_options.plist{o.file_options.current_file_num}));
        filepath=o.file_options.oplist{o.file_options.current_file_num};

        if ~isempty(EEG.filename)
            filename=sprintf('%s%s.set',prefix,EEG.filename);
        elseif ~isempty(EEG.setname)
            filename=sprintf('%sALLEEG(%d)_%s.set',prefix,o.file_options.current_file_num,EEG.setname);
        else
            filename=sprintf('%sALLEEG(%d).set',prefix,o.file_options.current_file_num);
        end
        EEG.filepath=filepath;
        EEG.filename=filename;
        mkdir([filepath filesep 'Intermediate']);
        EEG = pop_select(EEG, 'nochannel',length(eeg_chans)+length(ext_chans)+1:size(EEG.data,1));
        if (do_reref)
            if (max(EEG.data(ref_chan,:))==0 && min(EEG.data(ref_chan,:))==0)
                fprintf(log_file,'%.2f - Reference channel %d is already zeroed. Data was not re-referenced.\n',toc,ref_chan);
            elseif (o.ica_options.keep_ICA && ~isempty(EEG.icaweights))
                fprintf(log_file,'%.2f - Data was not re-referenced to maintain existing ICA weights. Bad channel detection may be ineffective.\n',toc,ref_chan);
            else
                EEG = h_pop_reref( EEG, ref_chan, 'exclude', ext_chans, 'keepref', 'on');
            end
        end
    else
        EEG=[];
        fprintf('Unknown file format.\n');
        fprintf(log_file,'%.2f - Unknown file format. Cannot process.\n',toc);
        return;
    end
    EEG = eeg_checkset(EEG);

    % Check if channel locations exist, and if not load them from disk.
    if (~isfield(EEG.chanlocs,'X') || ~isfield(EEG.chanlocs,'Y') || ~isfield(EEG.chanlocs,'Z') || isempty(EEG.chanlocs)) || isempty([EEG.chanlocs(:).X]) || isempty([EEG.chanlocs(:).Y]) || isempty([EEG.chanlocs(:).Z])
        EEG = pop_chanedit(EEG, 'load', {channel_locations_file});
        EEG.saved='no';
        fprintf(log_file,'%.2f - Loaded channel locations file from %s.\n',toc,channel_locations_file);
    end
    %EEG = pop_saveset(EEG,'savemode','resave');

    %%%%%%%%%%%%%%%%
    % Save options %
    %%%%%%%%%%%%%%%%
    do_saves=(~using_ALLEEG || (o.file_options.save_ALLEEG && ~isempty(EEG.filename)) || ~isempty(o.file_options.output_folder_name));
    if (~do_saves)
        save_options = zeros(size(save_options));
    else
        EEG = pop_saveset(EEG,'filename',[filename '.set'],'filepath',filepath,'savemode','onefile');
    end
    save_before_filter = save_options(1);
    save_before_interp = save_options(2);
    save_before_epoch = save_options(3);
    save_before_ica_rej = save_options(4);
    save_before_epoch_interp = save_options(5);

    if save_before_filter
        EEGBAK=EEG;
        EEGBAK.setname = ['pre_filt_' EEG.setname];
        pop_saveset(EEGBAK,'filename',['1_pre_filt_' EEG.filename],'filepath',[filepath filesep 'Intermediate'],'savemode','onefile');
        clear EEGBAK;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Filtering %
    %%%%%%%%%%%%%

    resample_frequency=o.filter_options.resample_freq;
    do_resample=o.filter_options.resample_on;
    % Downsampling is done later (shouldn't really be done at all).

    do_hipass=o.filter_options.hpf_on;
    do_lopass=o.filter_options.lpf_on;
    do_notch=o.filter_options.notch_on;

    if any(any(isnan(EEG.data)))
        fprintf('NaN in EEG data before filtering.\n');
    end

    if do_hipass
        w_h=o.filter_options.hpf_freq;
        t_h=o.filter_options.hpf_bandwidth;
        r_h=o.filter_options.hpf_ripple;
        a_h=o.filter_options.hpf_attenuation;

        [m, wtpass, wtstop] = pop_firpmord([w_h-(t_h) w_h+(t_h)], [0 1], [10^(-1*abs(a_h)/20) (10^(r_h/20)-1)/(10^(r_h/20)+1)], EEG.srate);
        if mod(m,2);m=m+1;end;
        EEG = pop_firpm(EEG, 'fcutoff', w_h, 'ftrans', t_h, 'ftype', 'highpass', 'wtpass', wtpass, 'wtstop', wtstop, 'forder', m);
        EEG.saved='no';
        fprintf(log_file,'%.2f - Highpass filter: %.3fHz, transition band: %.2f, order: %d.\n',toc,w_h,t_h,m);
    end

    if do_lopass
        w_l=o.filter_options.lpf_freq;
        t_l=o.filter_options.lpf_bandwidth;
        r_l=o.filter_options.lpf_ripple;
        a_l=o.filter_options.lpf_attenuation;

        [m, wtpass, wtstop] = pop_firpmord([w_l-(t_l) w_l+(t_l)], [1 0], [(10^(r_l/20)-1)/(10^(r_l/20)+1) 10^(-1*abs(a_l)/20)], EEG.srate);
        if mod(m,2);m=m+1;end;
        EEG = pop_firpm(EEG, 'fcutoff', w_l, 'ftrans', t_l, 'ftype', 'lowpass', 'wtpass', wtpass, 'wtstop', wtstop, 'forder', m);
        EEG.saved='no';
        fprintf(log_file,'%.2f - Lowpass filter: %.3fHz, transition band: %.2f, order: %d.\n',toc,w_l,t_l,m);
    end

    if do_notch
        for n=1:length(o.filter_options.notch_freq)
            w_n=[o.filter_options.notch_freq(n)-o.filter_options.notch_bandwidth1/2 o.filter_options.notch_freq(n)+o.filter_options.notch_bandwidth1/2];
            t_n=o.filter_options.notch_bandwidth2;
            r_n=o.filter_options.notch_ripple;
            a_n=o.filter_options.notch_attenuation;

            [m, wtpass, wtstop] = pop_firpmord([w_n(1)-(t_n) w_n(1)+(t_n) w_n(2)-(t_n) w_n(2)+(t_n)], [0 1 0], [10^(-1*abs(a_n)/20) (10^(r_n/20)-1)/(10^(r_n/20)+1) 10^(-1*abs(a_n)/20)], EEG.srate);
            if mod(m,2);m=m+1;end;
            EEG = pop_firpm(EEG, 'fcutoff', w_n, 'ftrans', t_n, 'ftype', 'bandstop', 'wtpass', wtpass, 'wtstop', wtstop, 'forder', m);
            EEG.saved='no';
            fprintf(log_file,'%.2f - Notch filter: %.3f to %.3fHz, transition band: %.2f, order: %d.\n',toc,w_n(1),w_n(2),t_n,m);
        end
    end

    if (do_saves), EEG = pop_saveset(EEG,'savemode','resave'); end

    if save_before_interp
        EEGBAK=EEG;
        EEGBAK.setname = ['pre_interp_' EEG.setname];
        pop_saveset(EEGBAK,'filename',['2_pre_interp_' EEG.filename],'filepath',[filepath filesep 'Intermediate'],'savemode','onefile');
        clear EEGBAK;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data cutoff point      %
    % Will be re-implemented %
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % if ~isempty(cutoff_markers) && any(cutoff_markers)
    % 	cutoff_point=[0 size(EEG.data,2)+1];
    % 	for u=1:length(EEG.event)
    % 		if EEG.event(u).type == cutoff_markers(1) || strcmp(EEG.event(u).type,cutoff_markers(1))
    % 			cutoff_point(1)=EEG.event(u).latency; % Finds the last 255 (check this one)
    % 		end
    % 		if EEG.event(u).type == cutoff_markers(2) || strcmp(EEG.event(u).type,cutoff_markers(2))
    % 			cutoff_point(2)=EEG.event(u).latency; % Finds the last 255 (check this one)
    % 		end
    % 	end
    % 	if cutoff_point(1) > 1
    % 		EEG = pop_select( EEG, 'nopoint',[1 cutoff_point(1)] );
    % 	end
    % 	if cutoff_point(2) < size(EEG.data,2)
    % 		EEG = pop_select( EEG, 'nopoint',[cutoff_point(2) size(EEG.data(:,:),2)] );
    % 	end
    % end

    % %New cutoff points for VESPA
    %
    % EEG = remevent(EEG,768);EEG=remevent(EEG,33536);
    %
    %
    % first_real_event = -1;
    % last_real_event = -1;
    %
    % for u=1:length(EEG.event)-2
    %
    % 	if ((EEG.event(u).latency - EEG.event(u+1).latency) * (1000/EEG.srate) < 100 && (EEG.event(u+1).latency - EEG.event(u+2).latency) * (1000/EEG.srate) < 100 && first_real_event == -1)
    % 		first_real_event = u;
    % 	end
    %
    % 	if (first_real_event ~= -1 && (EEG.event(u).latency - EEG.event(u+1).latency) * (1000/EEG.srate) > 100 && (EEG.event(u+1).latency - EEG.event(u+2).latency) * (1000/EEG.srate) > 100 && last_real_event == -1)
    % 		last_real_event = u;
    % 	end
    %
    % end
    %
    % first_real_time=max(EEG.event(first_real_event).latency - EEG.srate,1);
    %
    % if (last_real_event==-1)
    % 	last_real_time=min(EEG.event(end).latency + EEG.srate,size(EEG.data(:,:),2));
    % else
    % 	last_real_time=min(EEG.event(last_real_event).latency + EEG.srate,1);
    % end
    %
    % EEG = pop_select( EEG, 'point',[first_real_time:last_real_time] );
    % EEG.saved='no';
    % EEG = pop_saveset(EEG,'savemode','resave');
    %
    % fprintf(log_file,'Cropped between %.2f and %.2f seconds.\n',first_real_time/EEG.srate,last_real_time/EEG.srate);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Channel interpolation options %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1 Automatic interpolation of bad channels on or off (1 / 0)
    % 2 Radius for channel interpolation hypersphere (integer > 0)
    % 3 Automatic interpolation of channels per single epoch at end of process (1 / 0)
    % 4 Radius for epoch interpolation hypersphere (integer > 0)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    chans_to_interp=[];
    do_auto_interp = o.channel_options.channel_rejection_on;

    if do_auto_interp
        list_properties = channel_properties(EEG,eeg_chans,ref_chan);
        lengths = min_z(list_properties,o.channel_options.rejection_options); % Need to edit to make rejection_options.measure a vector, instead of multiple fields
        chans_to_interp = union(eeg_chans(logical(lengths)),o.channel_options.bad_channels);
        chans_to_interp = setdiff(chans_to_interp,ref_chan); % Ref chan may appear bad, but we shouldn't interpolate it!
        if (o.channel_options.exclude_EOG_chans)
            chans_to_interp = setdiff(chans_to_interp,o.ica_options.EOG_channels);
        end
        if ~o.channel_options.interp_after_ica
            if ~isempty(chans_to_interp)
                fprintf('Interpolating channel(s)');
                fprintf(' %d',chans_to_interp);
                fprintf('.\n');
                EEG = h_eeg_interp_spl(EEG,chans_to_interp,ext_chans);
                EEG.saved='no';
                fprintf(log_file,'%.2f - Interpolated channels',toc); fprintf(log_file,' %d',chans_to_interp); fprintf(log_file,'.\n');
            end
        end
    end

    if (do_saves), EEG = pop_saveset(EEG,'savemode','resave'); end

    if save_before_epoch
        EEGBAK=EEG;
        EEGBAK.setname = ['pre_epoch_' EEG.setname];
        pop_saveset(EEGBAK,'filename',['3_pre_epoch_' EEG.filename],'filepath',[filepath filesep 'Intermediate'],'savemode','onefile');
        clear EEGBAK;
    end

    %%% Do resampling here (if done pre-filtering, it creates problems). %%%
    %%% It does anyway, it seems. %%%
    if do_resample
        old_name = EEG.setname;
        old_srate = EEG.srate;
        EEG = pop_resample( EEG, resample_frequency);
        EEG.setname = old_name;
        fprintf(log_file,'%.2f - Resampled from %dHz to %dHz.\n',toc,old_srate,resample_frequency);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Epoch options %
    %%%%%%%%%%%%%%%%%
    % 1 Epoching on or off (1 / 0)
    % 2 Markers to epoch from (array of integers or cell of strings)
    % 3 Epoch length (vector of 2 floats, 1 negative, 1 positive) - seconds
    % 4 Baseline length for mean subtraction (vector of 2 integers) (0 => baseline subtraction off) - milliseconds
    % 5 Auto epoch rejection on or off (1 / 0)
    % 6 Radius for epoch rejection hypersphere (integer > 0)
    %%%%%%%%%%%%%%%%%
    markers = o.epoch_options.epoch_markers;
    epoch_length = o.epoch_options.epoch_limits;
    baseline_time = o.epoch_options.baseline_sub * 1000;
    do_epoch_rejection = o.epoch_options.epoch_rejection_on;
    do_epoching = ((~isempty(markers) && o.epoch_options.markered_epoch) || o.epoch_options.unmarkered_epoch) && any(o.epoch_options.epoch_limits) && length(o.epoch_options.epoch_limits)==2;

    %%%%%%%%%%%%%%
    % Epoch data %
    %%%%%%%%%%%%%%
    if do_epoching
        oldname = EEG.setname;
        if ~o.epoch_options.unmarkered_epoch
            EEGt = h_epoch(EEG,markers,epoch_length);
            EEG.setname = oldname;
            EEG.saved='no';
            if isnumeric(markers)
                fprintf(log_file,'%.2f - Epoched data on markers',toc);
                fprintf(log_file,' %d',markers);
                fprintf(log_file,'.\n');
            else
                fprintf(log_file,'%.2f - Epoched data on markers',toc);
                fprintf(log_file,' %s',markers{:});
                fprintf(log_file,'.\n');
            end
            if size(EEG.data,3)==0
                fprintf(log_file,'Epoch length too short, no epochs were generated.\n');
            else
                EEG=EEGt;
                clear EEGt;
            end
        else
            EEG = eeg_regepochs(EEG,o.epoch_options.unmarkered_epoch_interval,epoch_length,NaN);
            EEG.setname = oldname;
            EEG.saved='no';
            fprintf(log_file,'%.2f - Epoched data every %.2f seconds.\n',toc,o.epoch_options.unmarkered_epoch_interval);
        end

        % Remove epoch baselines after epoching:
        if any(baseline_time)
            EEG = pop_rmbase( EEG, baseline_time);
        end
    end
    if (size(EEG.data,3)>1)
        % Rereference just to print baseline variance, as otherwise the initial
        % BL variance is with a single reference, and the final in average
        % reference
        EEGtemp = h_pop_reref(EEG, [], 'exclude',ext_chans, 'refstate', ref_chan);
        fprintf(log_file,'Initial baseline variance: %.2f.\n',median(var(mean(EEGtemp.data(:,1:round(EEGtemp.srate*-1*EEGtemp.xmin),:),3),[],2)));
        clear EEGtemp;
    end

    if (do_saves), EEG = pop_saveset(EEG,'savemode','resave'); end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Epoch rejection section %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if do_epoch_rejection && size(EEG.data,3)>1
        if (o.channel_options.interp_after_ica)
            list_properties = epoch_properties(EEG,setdiff(eeg_chans,chans_to_interp));
        else
            list_properties = epoch_properties(EEG,eeg_chans);
        end
        [lengths] = min_z(list_properties,o.epoch_options.rejection_options);
        EEG=pop_rejepoch(EEG, find(lengths),0);
        fprintf(log_file,'%.2f - Rejected %d epochs',toc,length(find(lengths)));
        fprintf(log_file,' %d',find(lengths));
        fprintf(log_file,'.\n');
        EEG.saved='no';
    end

    if (do_saves), EEG = pop_saveset(EEG,'savemode','resave'); end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Average reference %
    %%%%%%%%%%%%%%%%%%%%%
    if (do_reref && ~o.ica_options.keep_ICA)
        if ~o.channel_options.interp_after_ica
            EEG = h_pop_reref(EEG, [], 'exclude',ext_chans, 'refstate', ref_chan);
        else
            EEG = h_pop_reref(EEG, [], 'exclude',[ext_chans chans_to_interp], 'refstate', ref_chan);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ICA options %
    %%%%%%%%%%%%%%%
    % 1 ICA on or off (1 / 0)
    % 2 Auto component rejection on or off (1 / 0)
    % 3 Radius for component rejection hypersphere (integer > 0)
    % 4 EOG channels (vector of integers)
    %%%%%%%%%%%%%%%
    do_ica = o.ica_options.run_ica;
    k_value = o.ica_options.k_value;
    do_component_rejection = o.ica_options.component_rejection_on;
    EOG_chans = o.ica_options.EOG_channels;
    ica_chans = o.ica_options.ica_channels;

    %%%%%%%%%%
    % Do ICA %
    %%%%%%%%%%
    if do_ica && (~o.ica_options.keep_ICA || isempty(EEG.icaweights))
        num_pca = min(floor(sqrt(size(EEG.data(:,:),2) / k_value)),(size(EEG.data,1) - length(chans_to_interp) - 1));
        num_pca = min(num_pca,length(setdiff(ica_chans,chans_to_interp)));
        if (o.channel_options.interp_after_ica)
            %EEG = pop_runica(EEG,  'icatype', 'runica', 'dataset',1, 'chanind',setdiff(ica_chans,chans_to_interp),'options',{'extended',1,'pca',num_pca});
            ica_chans=intersect(setdiff(ica_chans,chans_to_interp),union(eeg_chans,ext_chans));
            EEG = pop_runica(EEG,  'dataset',1, 'chanind',setdiff(ica_chans,chans_to_interp),'options',{'extended',1,'pca',num_pca});
        else
            %EEG = pop_runica(EEG,  'icatype', 'runica', 'dataset',1, 'chanind',ica_chans,'options',{'extended',1,'pca',num_pca});
            ica_chans=intersect(ica_chans,union(eeg_chans,ext_chans));
            EEG = pop_runica(EEG,  'dataset',1, 'chanind',ica_chans,'options',{'extended',1,'pca',num_pca});
        end
        EEG.saved='no';
        fprintf(log_file,'%.2f - Ran ICA.\n',toc);
    end

    if (do_saves), EEG = pop_saveset(EEG,'savemode','resave'); end

    if save_before_ica_rej
        EEGBAK=EEG;
        EEGBAK.setname = ['pre_comp_rej_' EEG.setname];
        pop_saveset(EEGBAK,'filename',['4_pre_comp_rej_' EEG.filename],'filepath',[filepath filesep 'Intermediate'],'savemode','onefile');
        clear EEGBAK;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Component rejection section %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Also includes topoplots   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if do_component_rejection && ~isempty(EEG.icaweights)
        EEG = eeg_checkset(EEG);
        original_name=EEG.setname;
        if do_lopass
            list_properties = component_properties(EEG,EOG_chans,[w_l-(t_l/2) w_l+(t_l/2)]);
        elseif ~isempty(o.ica_options.lopass_freq) && o.ica_options.lopass_freq~=0
            list_properties = component_properties(EEG,EOG_chans,[o.ica_options.lopass_freq-5 o.ica_options.lopass_freq+5]);
        else
            list_properties = component_properties(EEG,EOG_chans);
            o.ica_options.rejection_options.measure(2)=0;
        end
        [lengths] = min_z(list_properties,o.ica_options.rejection_options);
        bad_comps=find(lengths);

        % Plot stuff
        if (o.ica_options.IC_images)
            p=1;
            activations=eeg_getica(EEG);
            perc_vars = var(activations(:,:),[],2);
            perc_vars = 100*perc_vars./sum(perc_vars);
            for u=1:size(EEG.icawinv,2)
                if ~mod(u-1,16)
                    if (u~=1)
                        saveas(h,sprintf('%s%sIntermediate%sComponents_%d.png',filepath,filesep,filesep,p));
                        p=p+1;
                        close(h);
                    end
                    h=figure;
                end
                subplot(4,4,1+mod(u-1,16));
                %                 if (size(EEG.icawinv,1)~=length(EEG.chanlocs))
                %                    topoplot(EEG.icawinv(:,u),EEG.chanlocs(setdiff(1:length(EEG.chanlocs),chans_to_interp)));
                topoplot(EEG.icawinv(:,u),EEG.chanlocs(EEG.icachansind));
                %                 else
                %                     topoplot(EEG.icawinv(:,u),EEG.chanlocs);
                %end
                title(sprintf('Component %d\n%.1f%% variance',u,perc_vars(u)));
                if ~isempty(find(bad_comps==u, 1))
                    c=get(h,'Children');
                    c2=get(c(1),'Children');
                    set(c2(5),'FaceColor',[0.6 0 0]);
                    x=get(c2(5),'XData');
                    x(1:end/2)=1.5*(x(1:end/2));
                    set(c2(5),'XData',x);
                    y=get(c2(5),'YData');
                    y(1:end/2)=1.5*(y(1:end/2));
                    set(c2(5),'YData',y);
                end
            end
            %p=p+1;
            saveas(h,sprintf('%s%sIntermediate%sComponents_%d.png',filepath,filesep,filesep,p));
            if ~isempty(h)
                close(h);
            end
        end

        % Reject
        if ~isempty(find(lengths,1))
            fprintf('Rejecting components');
            fprintf(' %d',find(lengths));
            fprintf('.\n');
            EEG = pop_subcomp(EEG, find(lengths), 0);
            fprintf(log_file,'%.2f - Rejected %d components',toc,length(find(lengths)));
            fprintf(log_file,' %d',find(lengths));
            fprintf(log_file,'.\n');
        else
            fprintf('Rejected no components.\n');
            fprintf(log_file,'%.2f - Rejected no components.\n',toc);
        end
        EEG.setname=original_name;
        EEG.saved='no';
    elseif ~isempty(EEG.icawinv) && o.ica_options.IC_images
        activations=eeg_getica(EEG);
        perc_vars = var(activations(:,:),[],2);
        perc_vars = 100*perc_vars./sum(perc_vars);
        p=1;
        for u=1:size(EEG.icawinv,2)
            if ~mod(u-1,16)
                if (u~=1)
                    saveas(h,sprintf('%s%sIntermediate%sComponents_%d.png',filepath,filesep,filesep,p));
                    p=p+1;
                    close(h);
                end
                h=figure;
            end
            subplot(4,4,1+mod(u-1,16));
            topoplot(EEG.icawinv(:,u),EEG.chanlocs);
            title(sprintf('Component %d\n%.1f%% variance',u,perc_vars(u)));
        end
        %p=p+1;
        saveas(h,sprintf('%s%sIntermediate%sComponents_%d.png',filepath,filesep,filesep,p));
        if ~isempty(h)
            close(h);
        end
    end

    if (do_saves), EEG = pop_saveset(EEG,'savemode','resave'); end

    if save_before_epoch_interp
        EEGBAK=EEG;
        EEGBAK.setname = ['pre_epoch_interp_' EEG.setname];
        pop_saveset(EEGBAK,'filename',['5_pre_epoch_interp_' EEG.filename],'filepath',[filepath filesep 'Intermediate'],'savemode','onefile');
        clear EEGBAK;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Interpolation section part 2 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if o.channel_options.interp_after_ica
        if ~isempty(chans_to_interp)
            fprintf('Interpolating channel(s)');
            fprintf(' %d',chans_to_interp);
            fprintf('.\n');
            EEG = h_eeg_interp_spl(EEG,chans_to_interp,ext_chans);
            EEG.saved='no';
            fprintf(log_file,'%.2f - Interpolated channels',toc); fprintf(log_file,' %d',chans_to_interp); fprintf(log_file,'.\n');
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Epoch interpolation section %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    do_epoch_interp=o.epoch_interp_options.epoch_interpolation_on;
    if do_epoch_interp && length(size(EEG.data)) > 2
        status = '';
        lengths_ep=cell(1,size(EEG.data,3));
        for v=1:size(EEG.data,3)
            list_properties = single_epoch_channel_properties(EEG,v,eeg_chans);
            lengths_ep{v}=eeg_chans(logical(min_z(list_properties,o.epoch_interp_options.rejection_options)));
            status = [status sprintf('%d: ',v) sprintf('%d ',lengths_ep{v}) sprintf('\n')];
        end
        EEG=h_epoch_interp_spl(EEG,lengths_ep,ext_chans);
        EEG.saved='no';
        epoch_interps_log_file=fopen([filepath filesep filename '_epoch_interpolations.txt'],'a');
        fprintf(epoch_interps_log_file,'%s',status);
        fclose(epoch_interps_log_file);
        fprintf(log_file,'%.2f - Did per-epoch interpolation cleanup.\n',toc);
        fprintf(log_file,['See ' filename(1:end-4) '_epoch_interpolations.txt for details.\n']);
    end

    if ~isempty(o.channel_options.op_ref_chan)
        EEG = h_pop_reref(EEG, o.channel_options.op_ref_chan, 'exclude',ext_chans, 'refstate', [], 'keepref', 'on');
    end
    if (do_saves), EEG = pop_saveset(EEG,'savemode','resave'); end

    if using_ALLEEG
        fprintf('Done with ALLEEG(%d) - %s.\nTook %d seconds.\n',o.file_options.current_file_num,EEG.setname,toc);
    else
        fprintf('Done with file %s.\nTook %d seconds.\n',[filepath filesep filename extension],toc);
    end

    fprintf(log_file,'%.2f - Finished.\n',toc);
    if (size(EEG.data,3>1))
        fprintf(log_file,'Final baseline variance: %.2f.\n',median(var(mean(EEG.data(:,1:round(EEG.srate*-1*EEG.xmin),:),3),[],2)));
        % More stats here!
    end
    fclose(log_file);

    if (using_ALLEEG)
        assignin('base','FASTER_TMP_EEG',EEG);
        if o.file_options.overwrite_ALLEEG
            evalin('base',sprintf('ALLEEG(%d)=FASTER_TMP_EEG; clear FASTER_TMP_EEG',o.file_options.current_file_num));
        else
            evalin('base','[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, FASTER_TMP_EEG);clear FASTER_TMP_EEG;');
        end
    end
catch
    m=lasterror;
    EEG_state{1}=evalc('disp(EEG)');
    try
        if ~isempty(fopen(log_file))
            frewind(log_file);
            EEG_state{2}=fscanf(log_file,'%c',inf);
            
            try fclose(log_file); catch; end;
        end
    catch
    end
    EEG_state{3}=option_wrapper;
    EEG_state{4}=builtin('version');
    if exist('eeg_getversion','file')
        EEG_state{5}=eeg_getversion;
    else
        EEG_state{5}=which('eeglab');
    end
    
    assignin('caller','EEG_state',EEG_state);
    rethrow(m);
end
end
