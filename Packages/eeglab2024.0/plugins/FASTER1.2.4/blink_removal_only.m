function blink_removal_only(startDir,blink_chans)
[paths names]=extsearchc(startDir,'.bdf');
prefix = 'EOG_removal_';

for u=1:length(paths)
    try
        EEG = pop_loadset('filename',['6_pre_comp_rej_' names{u}(1:end-4) '.set'],'filepath',[paths{u} filesep 'Intermediate']);
        
        for w=1:size(EEG.icaact,1)
            for v=1:length(blink_chans)
                x=zeros(1,length(blink_chans));
                if ~(max(EEG.data(blink_chans(v),:))==0 && min(EEG.data(blink_chans(v),:))==0);
                    f = corrcoef(EEG.icaact(w,:),EEG.data(blink_chans(v),:));
                    x(v) = abs(f(1,2));
                else
                    x(v) = 0;
                end
            end
            list_properties(w,1)=max(x);
        end

        [lengths] = min_z(list_properties);
        EEG = pop_subcomp(EEG, find(lengths >= 1), 0);
        EEG = eeg_checkset(EEG);

        lengths_ep=cell(1,size(EEG.data,3));
        for v=1:size(EEG.data,3)
            list_properties = single_epoch_channel_properties(EEG,v,64);
            lengths_ep{v}=find(min_z(list_properties));
        end

        EEG=h_epoch_interp_spl(EEG,lengths_ep,65:72);

        pop_saveset(EEG,'filename',[prefix names{u}(1:end-4) '.set'],'filepath',[paths{u}]);
    catch
        fprintf('\nSkipped file %s%s%s.\n\n',paths{u},filesep,names{u});
    end
end