function [EEG did_epoch] = h_epoch(EEG,markers,epoch_length)
tempevents={EEG.event(:).type};
if isnumeric(tempevents{1});
	tempevents2=textscan(sprintf('%d ',tempevents{:}),'%s ');
	tempevents2=tempevents2{1}';
	events = unique(tempevents2);
else
	events = unique(tempevents);
end
if ~iscell(markers)
	epoch_markers=textscan(num2str(markers),'%s ');
	epoch_markers=epoch_markers{1}';
else
	epoch_markers=markers;
end
[cell_markers ia]=intersect(events,epoch_markers);
if ~isempty(ia)
    EEG=pop_epoch(EEG,cell_markers,epoch_length);
    did_epoch=1;
else
    warning('No epochs that were specified are present in the file. Data was not epoched.')
    did_epoch=0;
end