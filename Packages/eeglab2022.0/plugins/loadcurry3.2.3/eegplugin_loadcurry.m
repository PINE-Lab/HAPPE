function eegplugin_loadcurry(fig,try_strings,catch_strings)

% create menu
toolsmenu1 = findobj(fig, 'tag', 'import data');

uimenu( toolsmenu1, 'label', 'From Neuroscan Curry files', 'separator','on','callback', '[EEG LASTCOM]=pop_loadcurry;  [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); eeglab redraw');
