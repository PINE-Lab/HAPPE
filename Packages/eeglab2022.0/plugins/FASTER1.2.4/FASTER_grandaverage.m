function FASTER_grandaverage(startDir,option_wrapper,all_errors,top_log,plist,nlist)

if ~option_wrapper.options.averaging_options.make_GA
    return;
end

using_ALLEEG=option_wrapper.options.file_options.using_ALLEEG;
save_ALLEEG=option_wrapper.options.file_options.save_ALLEEG;
outDir=option_wrapper.options.file_options.output_folder_name;

if using_ALLEEG
    length_ALLEEG=evalin('base','length(ALLEEG);');
    if ~(length(plist)==length_ALLEEG)
        plist=cell(1,length_ALLEEG);
    end
end

try
    if (option_wrapper.options.averaging_options.make_GA && (~isempty(plist) || using_ALLEEG))
        GA_markers=option_wrapper.options.averaging_options.GA_markers;

        if ~isempty(GA_markers)
            GAs=cell(size(GA_markers));
            p=ones(size(GA_markers));
        else
            p=1;
        end
        GA_epoch_length=option_wrapper.options.averaging_options.GA_epoch_length;

        %% Do per-file or per-ALLEEG-dataset averaging

        for i=1:length(plist)
            % If the output folder option is set, get the proper subdirectory
            % under that, otherwise use the original directory
            if ~using_ALLEEG && ~isempty(outDir)
                dir_structure=get_dir_structure(plist{i},option_wrapper.options.file_options.folder_name);
                filepath=outDir;
                for v=1:length(dir_structure)
                    filepath=[filepath filesep dir_structure{v}];
                end
            else
                filepath=plist{i};
            end

            if (~isempty(option_wrapper.options.file_options.searchstring))
                searchstring2=option_wrapper.options.file_options.searchstring;
            else
                searchstring2=nlist{i};
            end
            if using_ALLEEG || (~isempty(strfind(nlist{i},searchstring2)) || ~isempty(strfind(filepath,searchstring2))) && exist([filepath filesep option_wrapper.options.file_options.file_prefix nlist{i}(1:end-4) '.set'],'file')
                if ~using_ALLEEG
                    EEGt=pop_loadset('filepath',filepath,'filename',[option_wrapper.options.file_options.file_prefix nlist{i}(1:end-4) '.set']);
                elseif ~evalin('base',sprintf('isempty(ALLEEG(%d).data);',i))
                    EEGt=evalin('base',sprintf('ALLEEG(%d);',i));
                end
                if ~isempty(option_wrapper.options.averaging_options.GA_markers) && ~isempty(EEGt) && ~isempty(EEGt.data)
                    if iscell(GA_markers)
                        for v=1:length(GA_markers)
                            [EEGt1 did_epoch]=h_epoch(EEGt,GA_markers{v},GA_epoch_length + [2/EEGt.srate -2/EEGt.srate]);
                            if ~isempty(EEGt1.data) && size(EEGt1.data,3)>1 && did_epoch
                                if exist('trimmean','file')==2
                                    trim_mean_on=option_wrapper.options.averaging_options.GA_trimmed_mean;
                                    trim_mean_perc=option_wrapper.options.averaging_options.GA_trimmed_mean_perc;
                                    GAs{v}(:,:,p(v))=trimmean(EEGt1.data,trim_mean_on*trim_mean_perc,3);
                                else
                                    GAs{v}(:,:,p(v))=mean(EEGt1.data,3);
                                end
                                p(v)=p(v)+1;
                            end
                        end
                        %p=p+1;
                    else
                        for v=1:length(GA_markers)
                            [EEGt1 did_epoch]=h_epoch(EEGt,GA_markers(v),GA_epoch_length + [2/EEGt.srate -2/EEGt.srate]);
                            if ~isempty(EEGt1.data) && size(EEGt1.data,3)>1 && did_epoch
                                if exist('trimmean','file')==2
                                    trim_mean_on=option_wrapper.options.averaging_options.GA_trimmed_mean;
                                    trim_mean_perc=option_wrapper.options.averaging_options.GA_trimmed_mean_perc;
                                    GAs{v}(:,:,p(v))=trimmean(EEGt1.data,trim_mean_on*trim_mean_perc,3);
                                else
                                    GAs{v}(:,:,p(v))=mean(EEGt1.data,3);
                                end
                                p(v)=p(v)+1;
                            else

                            end
                        end
                        %p=p+1;
                    end
                elseif ~isempty(EEGt.data) && size(EEGt.data,3)>1
                    if exist('trimmean','file')==2
                        trim_mean_on=option_wrapper.options.averaging_options.GA_trimmed_mean;
                        trim_mean_perc=option_wrapper.options.averaging_options.GA_trimmed_mean_perc;
                        GAs(:,:,p)=trimmean(EEGt.data,trim_mean_on*trim_mean_perc,3);
                    else
                        GAs(:,:,p)=mean(EEGt.data,3);
                    end
                    p=p+1;
                elseif size(EEGt.data,3)==1
                    warning('Continuous dataset not included in grand average. A blank epoch may be present.');
                else
                    warning('Empty dataset not included in grand average. A blank epoch may be present.');
                end
            end
        end

        %% Do rejection

        if ~isempty(option_wrapper.options.averaging_options.GA_markers) && exist('GAs','var') && ~isempty(GAs)
            if iscell(GA_markers)
                for v=1:length(GA_markers)
                    if ~isempty(GAs{v})
                        cl=EEGt.chanlocs; ci=EEGt.chaninfo;
                        EEGt=make_EEG(GAs{v},GA_markers{v},EEGt.srate,GA_epoch_length + [2/EEGt.srate -2/EEGt.srate]);
                        EEGt.chanlocs=cl; EEGt.chaninfo=ci;
                        if (option_wrapper.options.averaging_options.subject_removal_on) && size(EEGt.data,3)>1
                            list_properties=GA_properties(EEGt,option_wrapper.options.channel_options.eeg_chans,option_wrapper.options.ica_options.EOG_channels);
                            lengths=min_z(list_properties,option_wrapper.options.averaging_options.rejection_options);
                            if ~isempty(find(lengths, 1))
                                bad_subjs=find(lengths);
                                EEGt=pop_rejepoch(EEGt,bad_subjs,0);
                                fprintf(top_log,'In the grand average for marker %s, the following files were removed:',GA_markers{v});
                                fprintf(top_log,'%s%s%s\n',plist{bad_subjs},filesep,nlist{bad_subjs});
                                fprintf(top_log,'\n');
                            end
                        end
                        EEGt=eeg_checkset(EEGt);
                        if ~isempty(outDir)
                            pop_saveset(EEGt,'filename',[option_wrapper.options.file_options.file_prefix EEGt.setname '.set'],'filepath',outDir);
                        elseif ~using_ALLEEG || save_ALLEEG
                            pop_saveset(EEGt,'filename',[option_wrapper.options.file_options.file_prefix EEGt.setname '.set'],'filepath',startDir);
                        end
                        if using_ALLEEG
                            assignin('base','FASTER_Temp_EEG',EEGt);
                            evalin('base','[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, FASTER_Temp_EEG);clear FASTER_Temp_EEG;eeglab redraw;');
                        end
                    else
                        fprintf(top_log,'No markers %s detected in any files.\n',GA_markers{v})
                    end
                end
            else
                for v=1:length(GA_markers)
                    if ~isempty(GAs{v})
                        cl=EEGt.chanlocs; ci=EEGt.chaninfo;
                        EEGt=make_EEG(GAs{v},GA_markers(v),EEGt.srate,GA_epoch_length + [2/EEGt.srate -2/EEGt.srate]);
                        EEGt.chanlocs=cl; EEGt.chaninfo=ci;
                        if (option_wrapper.options.averaging_options.subject_removal_on) && size(EEGt.data,3)>1
                            list_properties=GA_properties(EEGt,option_wrapper.options.channel_options.eeg_chans,option_wrapper.options.ica_options.EOG_channels);
                            lengths=min_z(list_properties,option_wrapper.options.averaging_options.rejection_options);
                            if ~isempty(find(lengths, 1))
                                bad_subjs=find(lengths);
                                EEGt=pop_rejepoch(EEGt,bad_subjs,0);
                                fprintf(top_log,'In the grand average for marker %d, the following files were removed:',GA_markers(v));
                                fprintf(top_log,'%s%s%s\n',plist{bad_subjs},filesep,nlist{bad_subjs});
                                fprintf(top_log,'\n');
                            end
                        end
                        EEGt=eeg_checkset(EEGt);
                        if ~isempty(outDir)
                            pop_saveset(EEGt,'filename',[option_wrapper.options.file_options.file_prefix EEGt.setname '.set'],'filepath',outDir);
                        elseif ~using_ALLEEG || save_ALLEEG
                            pop_saveset(EEGt,'filename',[option_wrapper.options.file_options.file_prefix EEGt.setname '.set'],'filepath',startDir);
                        end
                        if using_ALLEEG
                            assignin('base','FASTER_Temp_EEG',EEGt);
                            evalin('base','[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, FASTER_Temp_EEG);clear FASTER_Temp_EEG;eeglab redraw;');
                        end
                    else
                        fprintf(top_log,'No markers %d detected in any files.\n',GA_markers(v))
                    end
                end
            end
        elseif exist('GAs','var') && ~isempty(GAs)
            cl=EEGt.chanlocs; ci=EEGt.chaninfo;
            EEGt=make_EEG(GAs,1,EEGt.srate,GA_epoch_length + [2/EEGt.srate -2/EEGt.srate]);
            EEGt.chanlocs=cl; EEGt.chaninfo=ci;
            if (option_wrapper.options.averaging_options.subject_removal_on) && size(EEGt.data,3)>1
                list_properties=GA_properties(EEGt,option_wrapper.options.channel_options.eeg_chans,option_wrapper.options.ica_options.EOG_channels);
                lengths=min_z(list_properties,option_wrapper.options.averaging_options.rejection_options);
                if ~isempty(find(lengths, 1))
                    bad_subjs=find(lengths);
                    EEGt=pop_rejepoch(EEGt,bad_subjs,0);
                    fprintf(top_log,'In the grand average, the following files were removed:');
                    fprintf(top_log,'%s%s%s\n',plist{bad_subjs},filesep,nlist{bad_subjs});
                    fprintf(top_log,'\n');
                end
            end
            EEGt=eeg_checkset(EEGt);

            if ~isempty(outDir)
                pop_saveset(EEGt,'filename',[option_wrapper.options.file_options.file_prefix 'GA.set'],'filepath',outDir);
            elseif ~using_ALLEEG || save_ALLEEG
                pop_saveset(EEGt,'filename',[option_wrapper.options.file_options.file_prefix 'GA.set'],'filepath',startDir);
            end
            if using_ALLEEG
                assignin('base','FASTER_TMP_EEG',EEGt);
                evalin('base','[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, FASTER_TMP_EEG);clear FASTER_TMP_EEG;eeglab redraw;');
            end
        else
            fprintf(top_log,'No epoched files were found.\n')
        end
    end
catch
    m=lasterror;
    fprintf('Error - %s.\n',m.message);
    fprintf(top_log,'Error in grand averaging - %s.\n',m.message);
    if option_wrapper.debug
        if exist('EEGt','var')
            EEG_state{1}=evalc('disp(EEGt)');
        end
        EEG_state{2}=option_wrapper;
        EEG_state{3}=builtin('version');
        if exist('eeg_getversion','file')
            EEG_state{4}=eeg_getversion;
        else
            EEG_state{4}=which('eeglab');
        end
        all_errors{end+1,1}=m;
        all_errors{end,2}=EEG_state;
        save([startDir filesep option_wrapper.options.file_options.file_prefix  'FASTER_errors.mat'],'all_errors','-mat');
    end
end
