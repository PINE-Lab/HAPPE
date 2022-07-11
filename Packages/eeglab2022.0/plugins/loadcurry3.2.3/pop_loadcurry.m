function [EEG, command] = pop_loadcurry(fullfilename, varargin)
%   Import a Neuroscan Curry file into EEGLAB. Currently supports Curry version 6, 7, 8,
%   and 9 data files (both continuous and epoched). Epoched
%   datasets are loaded in as continuous files with boundary events. Data
%   can be re-epoched using EEGLAB/ERPLAB functions.
%
%   Input Parameters:
%        1    Specify the filename of the Curry file (extension should be either .cdt, .dap, .dat, or .rs3). 
%
%   Example Code:
%
%       >> EEG = pop_loadcurry;   % an interactive uigetfile window
%       >> EEG = pop_loadcurry('C:\Studies\File1.cdt');    % no pop-up window 
%       >> EEG = loadcurry('C:\Studies\File1.cdt');    % no pop-up window 
%
%   Author for reading into Matlab: Neuroscan 
%   Author for translating to EEGLAB: Matthew B. Pontifex, Health Behaviors and Cognition Laboratory, Michigan State University, February 17, 2021
%   Github: https://github.com/mattpontifex/loadcurry
%
%   revision 3.0 - Curry9 compatibility. Increased the verbosity of the
%     equivalent call so that it provides the default parameters.
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
%   If there is an error with this code, please email pontifex@msu.edu with
%   the issue and I'll see what I can do.

    command = '';
    if nargin < 1 % No file was identified in the call

        [filename, filepath] = uigetfile({'*.dap;*.DAP;*.rs3;*.RS3;*.cdt;*.CDT' 'All Curry data files'; '*.cdt;*.CDT' 'Curry 8 and 9 data files'; '*.dap;*.DAP;*.rs3;*.RS3' 'Curry 6 and 7 data files'}, 'Choose a Neuroscan Curry file -- pop_loadcurry()'); 
        
        drawnow;
        if filename == 0 % the user aborted
            error('pop_loadcurry(): File selection cancelled. Error thrown to avoid overwriting data in EEG.')
            
        else % The user selected something
            
            % Get filename components
            [pathstr,name,ext] = fileparts([filepath, filename]);
            fullfilename = [filepath, filename];
            file = [pathstr, filesep, name];
            
            % Ensure the appropriate file types exist under that name
            boolfiles = 1;
            curryvers = 0;
            if strcmpi(ext, '.cdt')
                curryvers = 9;
                if (exist([file '.cdt'], 'file') == 0) || ((exist([file '.cdt.dpa'], 'file') == 0) && (exist([file '.cdt.dpo'], 'file') == 0))
                    boolfiles = 0;
                    
                    if (exist([file '.cdt'], 'file') == 0)
                        error('Error in pop_loadcurry(): The requested filename "%s" in "%s" does not have a .cdt file created by Curry 8 and 9.', name, filepath)
                    end
                    if ((exist([file '.cdt.dpa'], 'file') == 0) && (exist([file '.cdt.dpo'], 'file') == 0))
                        error('Error in pop_loadcurry(): The requested filename "%s" in "%s" does not have a .cdt.dpa/o file created by Curry 8 and 9.', name, filepath)
                    end
                end
            else
                curryvers = 7;
                if (exist([file '.dap'], 'file') == 0) || (exist([file '.dat'], 'file') == 0) || (exist([file '.rs3'], 'file') == 0)
                    boolfiles = 0;
                    error('Error in pop_loadcurry(): The requested filename "%s" in "%s" does not have all three file components (.dap, .dat, .rs3) created by Curry 6 and 7.', name, filepath)
                end
            end
            
            if (boolfiles == 1) & (curryvers ~= 0) % All files exist
                EEG = [];
                EEG = eeg_emptyset;
                [EEG, command] = loadcurry(fullfilename);
            else
                return % should be impossible to occur
            end
            
            % Show user the direct call equivalent
            if (curryvers > 7)
                command = sprintf('\nEquivalent Command:\n\tEEG = loadcurry(''%s%s'', ''KeepTriggerChannel'', ''True'', ''CurryLocations'', ''False'');',filepath, [name '.cdt']);
                disp(command)
            elseif (curryvers == 7)
                command = sprintf('\nEquivalent Command:\n\tEEG = loadcurry(''%s%s'', ''KeepTriggerChannel'', ''True'', ''CurryLocations'', ''False'');',filepath, [name '.dap']);
                disp(command)
            end
        end

    else % file was specified in the call
        
        EEG = [];
        EEG = eeg_emptyset;
        [EEG, command] = loadcurry(fullfilename);
    end
        
end

