function EEG = make_EEG(data,event_num,srate,epoch_lims)
EEG.data=data;
EEG.srate=srate;
EEG.xmin = epoch_lims(1);
EEG.xmax = epoch_lims(2);
EEG.nbchan=size(EEG.data,1);
EEG.pnts = size(EEG.data,2);
EEG.trials=size(EEG.data,3);

EEG.icawinv = [];
EEG.icasphere = [];
EEG.icaweights = [];
EEG.icaact = [];
EEG.ref = [];
%EEG.setname = [];
if ~ischar(event_num)
    EEG.setname = ['GA' sprintf('_%d',event_num)];
else
    EEG.setname = ['GA_' event_num];
end
EEG.filename = [];
EEG.filepath = [];
EEG.saved = 'no';
EEG.chanlocs = '';
EEG.chaninfo = '';
EEG.comments = '';
EEG=eeg_checkset(EEG);

if size(EEG.data,3 > 1)
    for u=1:size(EEG.data,3)
        if ~ischar(event_num)&&isscalar(event_num)
            EEG.event(u).type=sprintf('%d',event_num);
        else
            EEG.event(u).type=1;
        end
        EEG.event(u).latency=EEG.xmin*-1*EEG.srate + (EEG.xmax - EEG.xmin)*EEG.srate*(u-1);
        EEG.urevent(u).type=EEG.event(u).type;
        EEG.urevent(u).latency=EEG.event(u).latency;
        EEG.event(u).epoch=u;
        EEG.event(u).urevent=u;
        EEG.epoch(u).event=u;
        EEG.epoch(u).eventlatency=0;
        EEG.epoch(u).eventtype={[1]};
        EEG.epoch(u).urevent={[u]};
    end
else
    EEG.event=[];
    EEG.urevent=[];
end

EEG=eeg_checkset(EEG);