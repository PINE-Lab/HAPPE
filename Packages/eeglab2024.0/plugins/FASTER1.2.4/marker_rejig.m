function EEG=marker_rejig(EEG,marker1,N,operation,marker2,out_marker)
types = [EEG.event(:).type];

for u=1:length(types)
    switch(operation)
        case 'after'
            if (u>N && types(u)==marker1 && types(u-N)==marker2)
                EEG.event(u).type=out_marker;
            end
        case 'before'
            if (u<=(length(types)-N) && types(u)==marker1 && types(u+N)==marker2)
                EEG.event(u).type=out_marker;
            end
        case 'not_after'
            if (u>N && types(u)==marker1 && types(u-N)~=marker2)
                EEG.event(u).type=out_marker;
            end
        case 'not_before'
            if (u<=(length(types)-N) && types(u)==marker1 && types(u+N)~=marker2)
                EEG.event(u).type=out_marker;
            end
    end
end
