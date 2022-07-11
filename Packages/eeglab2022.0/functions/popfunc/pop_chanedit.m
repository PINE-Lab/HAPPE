% pop_chanedit() - Edit the channel locations structure of an EEGLAB dataset,
%                  EEG.chanlocs. For structure location and file formats, 
%                  see >> help readlocs  
% 
%                  EEG.chanlocs. For structure location and file formats,
%                  see >> help readlocs
%
% Usage:    >> EEG = pop_chanedit( EEG, 'key1', value1, 'key2', value2, ... ); 
%           >> [ chanlocs options ] = pop_chanedit( chanlocs, 'key1', value1);
%           >> [ chanlocs chaninfo options ] = pop_chanedit( chanlocs, chaninfo, ...
%                        'key1', value1, 'key2', value2, ... );
%
% Graphic interface:
%   "Channel information ('field name')" - [edit boxes] display channel field
%                   contents for the current channel. Command line equivalent
%                   to modify these fields: 'transform'
%   "Opt. 3D center" - [button] optimally re-center 3-D channel coordinates. Uses
%                 chancenter(). Command line equivalent: 'convert', { 'chancenter'
%                 [xc yc zc] }, [xc yc zc] being the center of the sphere. Use []
%                 to find the center of the best fitting sphere.
%   "Rotate axis" - [button] force one electrode to one position and rotate the other
%                 electrodes accordingly. Command line equivalent: 'forcelocs'.
%   "Transform axis" - [button] perform any operation on channel fields. Command
%                 line equivalent: 'transform'.
%   "Xyz->polar & sph." - [button] convert 3-D cartesian coordinates to polar and
%                 3-D spherical coordinates. This is useful when you edit the
%                 coordinates manually. Command line equivalent: 'convert', 'cart2all'.
%   "Sph.->polar & xyz" - [button] convert 3-D spherical coordinates to polar and
%                 3-D cartesian coordinates. Command line equivalent: 'convert', 'sph2all'.
%   "Polar->sph & xyz" - [button] convert 2-D polar coordinates to 3-D spherical and
%                 3-D cartesian coordinates. Command line equivalent: 'convert', 'topo2all'.
%                 Note that if spherical radii are absent, they are forced to 1.
%   "Set head radius" - [button] change head size radius. This is useful
%                 to make channels location compatible with a specified spherical model.
%                 Command line equivalent: 'headrad'.
%   "Set channel types" - [button] set channel type names for a range of data channels.
%   "Delete chan" - [button] delete channel. Command line equivalent: 'delete'.
%   "Insert chan" - [button] insert channel before current channel.
%                 Command line equivalent: 'insert'.
%   "<<" - [button] scroll channel backward by 10.
%   "<" - [button] scroll channel backward by 1.
%   ">" - [button] scroll channel forward by 1.
%   ">>" - [button] scroll channel forward by 10.
%   "Append chan" - [button] append channel after the current channel.
%                 Command line equivalent: 'append'.
%   "Plot 2D"     - [button] plot channel locations in 2-D using topoplot()
%   "Plot radius [value (0.2-1.0), []=auto)" - [edit box] default plotting radius
%                 in 2-D polar views. This does NOT affect channel locations; it
%                 is only used for visualization. This parameter is attached to the
%                 chanlocs structure and is then used in all 2-D scalp topoplots.
%                 Default -> to data limits. Command line equivalent: 'plotrad'.
%   "Nose along +X" - [list] Indicate the direction of the nose. This information
%                 is used in functions like topoplot(), headplot() and dipplot().
%                 Command line equivalent: 'nosedir'.
%   "Plot 3D"     - [button] plot channel positions in 3-D using plotchans3d()
%   "Read locations" - [button] read location file using readlocs()
%                 Command line equivalent: 'load'.
%   "Read help"   - [button] display readlocs() function help.
%   "Save .ced"   - [button] save channel locations in native EEGLAB ".ced" format.
%                 Command line equivalent: 'save'.
%   "Save others" - [button] save channel locations in other formats using
%                 pop_writelocs() (see readlocs() for available channel formats).
%   "Cancel"      - [button] cancel all editing.
%   "Help"        - [button] display this help message.
%   "OK"          - [button] save edits and propagate to parent.
%
% Inputs:
%   EEG      - EEG dataset
%   chanlocs - EEG.chanlocs structure
%
% Optional inputs:
%   'convert'     - {conversion_type [args]} Conversion type may be: 'cart2topo'
%                   'sph2topo', 'topo2sph', 'sph2cart', 'cart2sph', or 'chancenter'.
%                   See help messages for these functions. Args are only relevant
%                   for 'chancenter'. More info is given in the graphic interface
%                   description above.
%   'transform'   - String command for manipulating arrays. 'chan' is full channel
%                   info. Fields that can be manipulated are 'labels', 'theta'
%                   'radius' (polar angle and radius), 'X', 'Y', 'Z' (cartesian
%                   3-D) or 'sph_theta', 'sph_phi', 'sph_radius' for spherical
%                   horizontal angle, azimuth and radius.
%                   Ex: 'chans(3) = chans(14)', 'X = -X' or a multi-step transform
%                   with steps separated by ';': Ex. 'TMP = X; X = Y; Y = TMP'
%   'changechan'  - {number value1 value2 value3 ...} Change the values of all fields
%                   for the given channel number, mimimally {num label theta radius}.
%                   Ex: 'changechan' {12 'PXz' -90 0.30}
%   'changefield' - {number field value} Change field value for channel number number.
%                   Ex: {34 'theta' 320.4}.
%   'insert'      - {number label theta radius X Y Z sph_theta sph_phi sph_radius }
%                   Insert new channel and specified values before the current channel
%                    number. If the number of values is less than 10, remaining
%                   fields will be 0. (Previously, this parameter was termed 'add').
%   'append'      - {num label theta radius X Y Z sph_theta sph_phi sph_radius }
%                   same as 'insert' (above) but insert the the new channel after
%                   the current channel number.
%   'delete'      - [int_vector] Vector of channel numbers to delete.
%   'forcelocs'   - [cell] call forcelocs() to force a particular channel to be at a
%                   particular location on the head sphere; rotate other channels
%                   accordingly.
%   'skirt'       - Topographical polar skirt factor (see >> help topoplot)
%   'shrink'      - Topographical polar shrink factor (see >> help topoplot)
%   'load'        - [filename|{filename, 'key', 'val'}] Load channel location file
%                   optional arguments (such as file format) to the function
%                   readlocs() can be specified if the input is a cell array.
%   'save'        - 'filename' Save text file with channel info.
%   'eval'        - [string] evaluate string ('chantmp' is the name of the channel
%                   location structure).
%   'headrad'     - [float] change head radius.
%   'lookup'      - [string] look-up channel numbers for standard locations in the
%                   channel location file given as input.
%  'rplurchanloc' - [1,0] [1] Replace EEG.urchanlocs.If EEG.urchanlocs is empty this
%                   option will be ignored and will be set to [1]
%
% Outputs:
%   EEG        - new EEGLAB dataset with updated channel location structures 
%                EEG.chanlocs, EEG.urchanlocs, EEG.chaninfo
%   chanlocs   - updated channel location structure
%   chaninfo   - updated chaninfo structure
%   options    - structure containing plotting options (equivalent to EEG.chaninfo)
%
% Ex:    EEG = pop_chanedit(EEG,'load', { 'dummy.elp' 'elp' }, 'delete', [3 4], ...
%                       'convert', { 'xyz->polar' [] -1 1 }, 'save', 'mychans.loc' )
%        % Load polhemus file, delete two channels, convert to polar (see
%        % cart2topo() for arguments) and save into 'mychans.loc'.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 20 April 2002
%
% See also: readlocs()

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 15 March 2002, arno@salk.edu
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

% hidden parameter
%   'gui' - [figure value], allow to process the same dialog box several times

function [chansout, chaninfo, urchans, com] = pop_chanedit(chans, orichaninfo, varargin);

% process multiple datasets
% -------------------------
if length(chans) > 1 && isfield(chans, 'setname') && nargin > 2
    [ chansout, com ] = eeg_eval( 'pop_chanedit', chans, 'params', { orichaninfo, varargin{:} } );
    return;
end

urchans  = [];
com ='';
if nargin < 1
    help pop_chanedit;
    return;
end
chansout = chans;
chaninfo   = [];
fig      = [];

if nargin < 2
    orichaninfo = [];
end

if isempty(chans) || all(~ishandle(chans))
    % in case an EEG structure was given as input
    % -------------------------------------------
    if isfield(chans, 'chanlocs')
        
        % process multiple datasets
        if length(chans) > 1
            sameAsFirst = arrayfun(@(x)isequaln(chans(1).chanlocs, x.chanlocs), chans(2:end));
            if ~all(sameAsFirst)
                error( [ 'All datasets need to have the exact same channel structure.' 10 'If you want to look up channel location for all datasets,' 10 'do so for the first one and write a loop' ] );
            end
            
            % pop up GUI for first dataset
            EEG = chans(1);
            [chansout, chaninfo, urchans, com] = pop_chanedit(EEG, orichaninfo, varargin{:});
            chans(1) = chansout;
            
            % pop up GUI for first dataset
            for iDat = 2:length(chans)
                EEG = chans(iDat);
                eval(com);
                chans(iDat) = EEG;
            end
            chansout = chans;
            return;
        end
            
        dataset_input = 1;
        EEG           = chans;
        chans         = EEG(1).chanlocs;
        nchansori     = EEG.nbchan;
        if isfield(EEG, 'chaninfo')
            chaninfo = EEG(1).chaninfo;
        else chaninfo = [];
        end
        if isfield(EEG, 'urchanlocs')
            urchans = EEG(1).urchanlocs;
        end
    else
        nchansori     = 0;
        dataset_input = 0;
        chaninfo      = orichaninfo;
    end

    % dealing with additional parameters
    % ----------------------------------
    if nargin > 1 && ~ischar(orichaninfo), % nothing
        if nargin > 2
            if ~ischar(varargin{1})
                urchans  = varargin{1};
                varargin = varargin(2:end);
            end
        end
    elseif nargin > 1 && ~isempty(orichaninfo) && ischar(orichaninfo)
        varargin = { orichaninfo varargin{:} };
        if isequal(orichaninfo, chaninfo)
            chaninfo    = [];
        end
        orichaninfo = [];
    end
    
    % insert "no data channels" in channel structure
    % ----------------------------------------------
    nbchan = length(chans);
    [tmp chaninfo chans] = eeg_checkchanlocs(chans, chaninfo);

    if isfield(chaninfo, 'shrink') && ~isempty(chaninfo.shrink)
        icadefs;
        if SHRINKWARNING
            warndlg2( [ 'You are currently shrinking channel locations for display.' 10 ...
                'A new option (more anatomically correct) is to plot channels' 10 ...
                'outside head limits so the shrink option has been disabled.' 10 ...
                '(Edit the icadefs file to disable this message)' ], 'Shrink factor warning');
        end
    end

    oldchaninfo = chaninfo;
end

% Checking for flag to replace 'urchanloc' field
indx_tmp = find(strcmp(varargin,'rplurchanloc'));
flag_replurchan = 0;
if ~isempty(indx_tmp)
    flag_replurchan = varargin{indx_tmp+1};
end

if nargin < 3 && isstruct(chans)

    totaluserdat = {};
    % lookup channel locations if necessary
    % -------------------------------------
    if ~all(cellfun('isempty', {chans.labels})) && all(cellfun('isempty', {chans.theta}))
        [chans chaninfo urchans com] = pop_chanedit(chans, chaninfo, 'lookupgui', []);
        for index = 1:length(chans)
            chans(index).ref       = '';
            chans(index).datachan  = 1;
        end
        if ~isempty(com)
            totaluserdat = com;
            %[chans chaninfo urchans com] = pop_chanedit(chans, chaninfo, com{:});
        end
    end

    commentfields = { 'Channel label ("label")', ...
                      'Polar angle ("theta")', 'Polar radius ("radius")', ...
                      'Cartesian X ("X")', ...
                      'Cartesian Y ("Y")', ...
                      'Cartesian Z ("Z")', ...
                      'Spherical horiz. angle ("sph_theta")', ...
                      'Spherical azimuth angle ("sph_phi")', ...
                      'Spherical radius ("sph_radius")' ...
                      'Channel type' 'Reference' ...
                      'Index in backup ''urchanlocs'' structure' ...
                      'Channel in data array (set=yes)' };

    % add field values
    % ----------------
    geometry = { 1 };
    tmpstr = sprintf('Channel information ("field_name"):');
    uilist = { { 'Style', 'text', 'string', tmpstr, 'fontweight', 'bold'  } };

    uiconvert = { ...
        { 'Style', 'pushbutton', 'string', 'Opt. head center', 'callback', 'pop_chanedit(gcbf, [], ''chancenter'', []);' } ...
        { 'Style', 'pushbutton', 'string', 'Rotate axis'     , 'callback', 'pop_chanedit(gcbf, [], ''forcelocs'', []);' } ...
        { 'Style', 'pushbutton', 'string', 'Transform axes'  , 'callback', 'pop_chanedit(gcbf, [], ''transform'', []);' } ...
        {  }, ...
        { 'Style', 'pushbutton', 'string', 'xyz -> polar & sph.', 'callback', 'pop_chanedit(gcbf, [], ''convert'', {''cart2all''});' }, ...
        { 'Style', 'pushbutton', 'string', 'sph. -> polar & xyz', 'callback', 'pop_chanedit(gcbf, [], ''convert'', {''sph2all'' });' }, ...
        { 'Style', 'pushbutton', 'string', 'polar -> sph. & xyz', 'callback', 'pop_chanedit(gcbf, [], ''convert'', {''topo2all''});' }, ...
        {  }, ...
        { 'Style', 'pushbutton', 'string', 'Set head radius', 'callback',   'pop_chanedit(gcbf, [], ''headrad'', []);' } ...
        { 'Style', 'pushbutton', 'string', 'Set channel types', 'callback', 'pop_chanedit(gcbf, [], ''settype'', []);' } ...
        { 'Style', 'pushbutton', 'string', 'Set reference', 'callback',     'pop_chanedit(gcbf, [], ''setref'' , []);' } ... 
        { } { } };

    % create text and edit for each field
    % -----------------------------------
    allfields = { 'labels' 'theta' 'radius' 'X' 'Y' 'Z' 'sph_theta' 'sph_phi' 'sph_radius' 'type' 'ref' 'urchan' 'datachan' };
    for index = 1:length(allfields)-1
        cbfield = [ 'valnumtmp   = str2num(get(findobj(gcbf, ''tag'', ''chaneditnumval''), ''string''));' ...
                    'pop_chanedit(gcbf, [], ''changefield'', { valnumtmp ''' allfields{index} ''' get(gcbo, ''string'') });' ...
                    'clear valnumtmp;' ];
        geometry = { geometry{:} [1.5 1 0.2 1] };
        uilist   = { uilist{:}, ...
            { 'Style', 'text', 'string', commentfields{index} }, ...
            { 'Style', 'edit', 'tag', [ 'chanedit' allfields{index} ], 'string', ...
            num2str(getfield(chans,{1}, allfields{index})), 'horizontalalignment', 'center', 'callback', cbfield } ...
              { } uiconvert{index} };
    end

    % special checkbox for chandata field
    % -----------------------------------
    geometry = { geometry{:} [2 0.35 0.5 1] };
    cbfield = [ 'valnumtmp   = str2num(get(findobj(gcbf, ''tag'', ''chaneditnumval''), ''string''));' ...
                'pop_chanedit(gcbf, [], ''changefield'', { valnumtmp ''' allfields{end} ''' get(gcbo, ''value'') });' ...
                'clear valnumtmp;' ];
    uilist   = { uilist{:}, ...
        { 'Style', 'text', 'string', commentfields{end} }, ...
        { 'Style', 'checkbox', 'tag', [ 'chanedit' allfields{end}], 'string', '' 'value', 1 'callback', cbfield } { } uiconvert{end} };

    % add buttons
    % -----------
    geometry =  { geometry{:} [1] [1.15 0.5 0.6 1.9 0.4 0.4 1.15] [1.15 0.7 0.7 1 0.7 0.7 1.15] };
    cb_del    = [ 'valnum   = str2num(char(get(findobj(gcbf,''tag'', ''chaneditnumval''), ''string'')));' ...
                  'pop_chanedit(gcbf, [], ''deletegui'', valnum);' ];
    cb_insert = [ 'valnum   = str2num(char(get(findobj(gcbf,''tag'', ''chaneditnumval''), ''string'')));' ...
                  'pop_chanedit(gcbf, [], ''insert'', valnum);' ];
    cb_append = [ 'valnum   = str2num(char(get(findobj(gcbf,''tag'', ''chaneditnumval''), ''string'')));' ...
                  'pop_chanedit(gcbf, [], ''append'', valnum);' ];

    uilist   = { uilist{:}, ...
        { }, ...
        { 'Style', 'pushbutton', 'string', 'Delete chan',  'callback', cb_del }, ...
        { },{ }, ...
        { 'Style', 'text'      , 'string', ['Channel number (of ' int2str(length(chans)) ')'], ...
        'fontweight', 'bold', 'tag', 'chaneditscantitle' }, { },{ },{ }, ...
        { 'Style', 'pushbutton', 'string', 'Insert chan',  'callback', cb_insert } ...
        { 'Style', 'pushbutton', 'string', '<<', 'callback', [ 'pop_chanedit(gcbf, [], ''movecursor'', -10);' ] } ...
        { 'Style', 'pushbutton', 'string', '<',  'callback', [ 'pop_chanedit(gcbf, [], ''movecursor'', -1);' ] } ...
        { 'Style', 'edit'      , 'string', '1', 'tag', 'chaneditnumval', 'callback', [ 'pop_chanedit(gcbf, []);' ] } ...
        { 'Style', 'pushbutton', 'string', '>',  'callback', [ 'pop_chanedit(gcbf, [], ''movecursor'', 1);' ] } ...
        { 'Style', 'pushbutton', 'string', '>>', 'callback', [ 'pop_chanedit(gcbf, [], ''movecursor'', 10);' ] } ...
        { 'Style', 'pushbutton', 'string', 'Append chan',  'callback', cb_append }, ...
        };

    % add sorting options
    % -------------------
    cb_rplurchan = [ 'valnumtmp   = get(findobj(gcbf, ''tag'', ''rplurchan''), ''value'');' ...
                     'pop_chanedit(gcbf, [], ''rplurchanloc'', valnumtmp);' ...
                     'clear valnumtmp;' ];
    noseparam = strmatch(upper(chaninfo.nosedir), { '+X' '-X' '+Y' '-Y' });
    if isempty(noseparam), error('Wrong value for nose direction'); end
    geometry = { geometry{:} [1] [0.9 1.3 0.6 1.1 0.9] [1] [1 1 1 1 1] [1]};
    uilist   = { uilist{:},...
        { } ...
        { 'Style', 'pushbutton', 'string', 'Plot 2-D', 'callback', 'pop_chanedit(gcbf, [], ''plot2d'', []);' },...
        { 'Style', 'text', 'string', 'Plot radius (0.2-1, []=auto)'} ...
        { 'Style', 'edit', 'string', char(chaninfo.plotrad), 'tag', 'plotrad' 'callback' 'pop_chanedit(gcbf, [], ''plotrad'', []);' } ...
        { 'Style', 'popupmenu',  'string', 'Nose along +X|Nose along -X|Nose along +Y|Nose along -Y', ...
        'tag' 'nosedir' 'value',noseparam, 'callback' 'pop_chanedit(gcbf,[],''nosedir'',[]);' 'listboxtop' noseparam } ...
        { 'Style', 'pushbutton', 'string', 'Plot 3-D (xyz)',     'callback', 'pop_chanedit(gcbf, [], ''plot3d'', []);' } ...
        {}, ...
        { 'Style', 'pushbutton', 'string', 'Read locations',              'callback', 'pop_chanedit(gcbf,[],''load'',[]);' }, ...
        { 'Style', 'pushbutton', 'string', 'Read locs help',              'callback', 'pophelp(''readlocs.m'');' }, ...
        { 'Style', 'pushbutton', 'string', 'Look up locs',                'callback', 'pop_chanedit(gcbf,[], ''lookupgui'', []);' }, ...
        { 'Style', 'pushbutton', 'string', 'Save (as .ced)',              'callback', 'pop_chanedit(gcbf,[], ''save'',[]);' } ...
        { 'Style', 'pushbutton', 'string', 'Save (other types)',          'callback', 'pop_chanedit(gcbf,[], ''saveothers'',[]);' } ...
        { 'Style', 'checkbox'  , 'string', 'Overwrite Original Channels', 'callback', cb_rplurchan, 'tag'   , 'rplurchan',  'value', flag_replurchan }...
        };

    % evaluation of command below is required to center text (if
    % declared a text instead of edit, the uicontrol is not centered)
    comeval = [ 'set(findobj( ''tag'', ''chanediturchan''), ''style'', ''text'', ''backgroundcolor'', [.66 .76 1] );' ...
                'set(findobj( ''tag'', ''chaneditref''),    ''style'', ''text'', ''backgroundcolor'', [.66 .76 1] );' ...
                'set(findobj( ''tag'', ''ok''), ''callback'', ''valnumtmp   = get(findobj(gcbf, ''''tag'''', ''''rplurchan''''), ''''value''''); pop_chanedit(gcbf, [],''''rplurchanloc'''',valnumtmp, ''''return'''', []);'')' ];

    userdata.chans     = chans;
    userdata.nchansori = nchansori;
    userdata.chaninfo  = chaninfo;
    userdata.urchans   = urchans ;
    userdata.commands  = totaluserdat;

    [results userdata returnmode] = inputgui( 'geometry', geometry, 'uilist', uilist, 'helpcom', ...
        'pophelp(''pop_chanedit'');', 'title', 'Edit channel info -- pop_chanedit()', ...
        'userdata', userdata, 'eval' , comeval );

    if length(results) == 0, 
        com = ''; 
        if dataset_input, chansout = EEG; end; 
        return; 
    end

    % transfer events back from global workspace
    chans      = userdata.chans;
    chaninfo   = userdata.chaninfo;
    urchans    = userdata.urchans; 
    if ~isempty(userdata.commands)
        com = sprintf('%s=pop_chanedit(%s, %s);', inputname(1), inputname(1), vararg2str(userdata.commands));
    end
    % Updating flag
    flag_replurchan = results{15};
else
    
    % call from command line or from a figure
    % ---------------------------------------
    currentpos = 0;
    if ishandle(chans)
        fig         = chans;
        userdata    = get(fig, 'userdata');
        chans       = userdata.chans;
        nchansori   = userdata.nchansori;
        chaninfo    = userdata.chaninfo;
        urchans     = userdata.urchans;
        currentpos  = str2num(get(findobj(fig, 'tag', 'chaneditnumval'), 'string'));
    end
    
    args = varargin;
    % no interactive inputs
    % scan all the fields of g
    % ------------------------
    for curfield = 1:2:length(args)
        switch lower(args{curfield})
            case 'rplurchanloc'
                if flag_replurchan, urchans = eeg_checkchanlocs(chans, chaninfo); end
                args{curfield}     = 'rplurchanloc';
                args{ curfield+1 } =  flag_replurchan;
            case 'return'
                [tmpchans] = eeg_checkchanlocs(chans);
                if nchansori ~= 0 && nchansori ~= length(tmpchans)
                    if ~popask(strvcat(['The number of data channels (' int2str(length(tmpchans)) ') not including fiducials does not'], ...
                            ['correspond to the initial number of channels (' int2str(nchansori) '), so for consistency purposes'], ...
                            'new channel information will be ignored if this function was called from EEGLAB', ...
                            'If you have added a reference channel manually, check the "Data channel" checkbox is off'))
                    else
                        set(findobj(fig, 'tag', 'ok'), 'userdata', 'stop');
                    end
                else
                    set(findobj(fig, 'tag', 'ok'), 'userdata', 'stop');
                end
                args = {};
            case 'plot3d', % GUI only
                tmpind = find(~cellfun('isempty', { chans.X }));
                if ~isempty(tmpind),
                    plotchans3d([ [ chans(tmpind).X ]' [ chans(tmpind).Y ]' [ chans(tmpind).Z ]'], { chans(tmpind).labels }); 
                else disp('cannot plot: no XYZ coordinates');
                end
                args = {};
            case 'plot2d', % GUI only
                plotrad = str2num(get(findobj(fig, 'tag', 'plotrad'), 'string'));
                figure; topoplot([],chans, 'style', 'blank', 'drawaxis', 'on', 'electrodes', ...
                                 'labelpoint', 'plotrad', plotrad, 'chaninfo', chaninfo);
                args = {};
            case 'movecursor', % GUI only
                currentpos = max(1,min(currentpos+args{curfield+1},length(chans)));
                args = {};
            case 'plotrad',
                if isempty( args{curfield+1} )
                     args{curfield+1} = str2num(get(findobj(fig, 'tag', 'plotrad'), 'string'));
                end
                chaninfo.plotrad = args{curfield+1};
            case 'forcelocs',
                if ~isempty(fig) % GUI BASED
                    [ comtmp tmpforce ] = forcelocs(chans); 
                    if ~isempty(tmpforce), 
                        args{curfield+1} = tmpforce{1};
                    end
                end
                if ~isempty(args{curfield+1})
                    chans = forcelocs(chans,args{curfield+1});
                    disp('Convert XYZ coordinates to spherical and polar');
                end
            case 'chancenter',
                if ~isempty(fig)
                    [chans newcenter tmpcom] = pop_chancenter(chans);
                    args{curfield  } = 'eval';
                    args{curfield+1} = tmpcom;
                end
            case 'convert',
                if iscell(args{curfield+1})
                    method=args{curfield+1}{1};
                    extraargs = args{curfield+1}(2:end);
                else
                    method=args{curfield+1};
                    extraargs = {''};
                end
                if ~isempty(fig) && ~strcmp(method, 'chancenter')
                    tmpButtonName=questdlg2( strvcat('This will modify fields in the channel structure', ...
                        'Are you sure you want to apply this function ?'), 'Confirmation', 'Cancel', 'Yes','Yes');
                    if ~strcmpi(tmpButtonName, 'Yes'), return; end
                end
                switch method
                    case 'chancenter',
                        if isempty(extraargs)
                            [X Y Z]=chancenter( [chans.X ]', [ chans.Y ]', [ chans.Z ]',[]);
                        else
                            [X Y Z]=chancenter( [chans.X ]', [ chans.Y ]', [ chans.Z ]', extraargs{:});
                        end
                        if isempty(X), return; end
                        for index = 1:length(chans)
                            chans(index).X  = X(index);
                            chans(index).Y  = Y(index);
                            chans(index).Z  = Z(index);
                        end
                        disp('Note: automatically convert XYZ coordinates to spherical and polar');
                        chans = convertlocs(chans, 'cart2all');
                    otherwise
                        chans = convertlocs(chans, method, 'verbose', 'on');
                end
            case 'settype'                
                if ~isempty(fig)
                    args{curfield+1} = inputdlg2({'Channel indices' 'Type (e.g. EEG)' }, ...
                                             'Set channel type', 1, { '' '' }, 'pop_chanedit');
                end
                try, tmpchans = args{curfield+1}{1}; tmptype = args{curfield+1}{2};catch, return; end
                if isempty(tmpchans) && isempty(tmptype), return; end
                if ischar(tmpchans)
                    tmpchans = eval( [ '[' tmpchans ']' ], 'settype: error in channel indices');
                end
                if ~ischar(tmptype), tmptype = num2str(tmptype); end
                for index = 1:length(tmpchans)
                    if tmpchans(index) > 0 && tmpchans(index) <= length(chans)
                        chans( tmpchans(index) ).type = tmptype;
                    end
                end
            case 'setref'
                if ~isempty(fig)
                    disp('Note that setting the reference only changes the reference labels');
                    disp('Use the re-referencing menu to change the reference');
                    args{curfield+1} = inputdlg2({'Channel indices' 'Reference (e.g. Cz)' }, ...
                                                 'Set channel reference', 1, { '' '' }, 'pop_chanedit');
                end
                try, tmpchans = args{curfield+1}{1}; tmpref = args{curfield+1}{2};catch, return; end
                if isempty(tmpchans) && isempty(tmpref), return; end
                if ischar(tmpchans)
                    tmpchans = eval( [ '[' tmpchans ']' ], 'settype: error in channel indices');
                end
                if ~ischar(tmpref), tmpref = num2str(tmpref); end
                for index = 1:length(tmpchans)
                    if tmpchans(index) > 0 && tmpchans(index) <= length(chans)
                        chans( tmpchans(index) ).ref = tmpref;
                    end
                end
            case 'transform'
                if ~isempty(fig)
                    args{curfield+1} = inputdlg2({'Enter transform: (Ex: TMP=X; X=-Y; Y=TMP or Y(3) = X(2), etc.' }, ...
                                                  'Transform', 1, { '' }, 'pop_chanedit');
                end
                try, tmpoper = args{curfield+1}; catch, return; end
                if isempty(deblank(tmpoper)), return; end
                if iscell(tmpoper), tmpoper = tmpoper{1}; end
                tmpoper = [ tmpoper ';' ];
                [eloc, labels, theta, radius, indices] = readlocs(chans);
                if isempty(findstr(tmpoper, 'chans'))
                    try,
                        X          = [ chans(indices).X ];
                        Y          = [ chans(indices).Y ];
                        Z          = [ chans(indices).Z ];
                        sph_theta  = [ chans(indices).sph_theta ];
                        sph_phi    = [ chans(indices).sph_phi ];
                        sph_radius = [ chans(indices).sph_radius ];
                        eval(tmpoper);
                        
                        for ind = 1:length(indices)
                            chans(indices(ind)).X = X(min(length(X),ind));
                            chans(indices(ind)).Y = Y(min(length(Y),ind));
                            chans(indices(ind)).Z = Z(min(length(Z),ind));
                            chans(indices(ind)).theta  = theta(min(length(theta),ind));
                            chans(indices(ind)).radius = radius(min(length(radius),ind));
                            chans(indices(ind)).sph_theta  = sph_theta(min(length(sph_theta),ind));
                            chans(indices(ind)).sph_phi    = sph_phi(min(length(sph_phi),ind));
                            chans(indices(ind)).sph_radius = sph_radius(min(length(sph_radius),ind));
                        end
                        
                        if     ~isempty(findstr(tmpoper, 'X')),          chans = convertlocs(chans, 'cart2all'); end
                        if     ~isempty(findstr(tmpoper, 'Y')),          chans = convertlocs(chans, 'cart2all'); end
                        if     ~isempty(findstr(tmpoper, 'Z')),          chans = convertlocs(chans, 'cart2all'); end
                        if     ~isempty(findstr(tmpoper, 'sph_theta')),  chans = convertlocs(chans, 'sph2all');
                        elseif ~isempty(findstr(tmpoper, 'theta')),      chans = convertlocs(chans, 'topo2all'); end
                        if     ~isempty(findstr(tmpoper, 'sph_phi')),    chans = convertlocs(chans, 'sph2all');  end
                        if     ~isempty(findstr(tmpoper, 'sph_radius')), chans = convertlocs(chans, 'sph2all');
                        elseif ~isempty(findstr(tmpoper, 'radius')),     chans = convertlocs(chans, 'topo2all'); end
                    catch, disp('Unknown error when applying transform'); end
                else
                    eval(tmpoper);
                end
                
            case 'headrad'
                if ~isempty(fig) % GUI
                     tmpres = inputdlg2({'Enter new head radius (same unit as DIPFIT head model):' }, ...
                                        'Head radius', 1, { '' }, 'pop_chanedit'); 
                     if ~isempty(tmpres),
                          args{ curfield+1 } = str2num(tmpres{1});
                     else return;
                     end
                end
                if ~isempty( args{ curfield+1 } )
                    allrad = [ chans.sph_radius ];
                    if length(unique(allrad)) == 1 % already spherical
                        chans = pop_chanedit(chans, 'transform', [ 'sph_radius = ' num2str( args{ curfield+1 } ) ';' ]);
                    else % non-spherical, finding best match
                        factor = args{ curfield+1 } / mean(allrad);
                        chans = pop_chanedit(chans, 'transform', [ 'sph_radius = sph_radius*' num2str( factor ) ';' ]);
                        disp('Warning: electrodes do not lie on a sphere. Sphere model fitting for');
                        disp('         dipole localization will work but generate many warnings');
                    end
                    chans = convertlocs(chans, 'sph2all');
                end
                
            case 'shrink'
                chans(1).shrink = args{ curfield+1 };
                
            case 'plotrad'
                chans(1).plotrad = args{ curfield+1 };
                
            case 'deletegui'
                chans(args{ curfield+1 })=[];
                currentpos = min(length(chans), currentpos);
                args{ curfield   } = 'delete';
                
            case 'delete'
                chans(args{ curfield+1 })=[];
                
            case 'changefield'
                tmpargs = args{ curfield+1 };
                if length( tmpargs ) < 3
                    error('pop_chanedit: not enough arguments to change field value');
                end
                if ~isempty(strmatch( tmpargs{2}, { 'X' 'Y' 'Z' 'theta' 'radius' 'sph_theta' 'sph_phi' 'sph_radius'}))
                     if ~isnumeric(tmpargs{3}), tmpargs{3} = str2num(tmpargs{3}); end
                end
                for paramChan = 2:2:length(tmpargs)
                    eval([ 'chans(' int2str(tmpargs{1}) ').'  tmpargs{paramChan} '=' reformat(tmpargs{paramChan+1} ) ';' ]);
                end
            case { 'insert' 'add' 'append' }
                tmpargs = args{ curfield+1 };
                allfields = fieldnames(chans);
                if isnumeric(tmpargs)
                    tmpargs2    = cell(1, length(allfields)+1);
                    tmpargs2{1} = tmpargs;
                    tmpargs     = tmpargs2;
                    if strcmpi(allfields{end}, 'datachan'), tmpargs{end} = 0; end
                end
                if length( tmpargs ) < length(allfields)+1
                    error('pop_chanedit: not enough arguments to change all field values');
                end
                num = tmpargs{1};
                if strcmpi(lower(args{curfield}), 'append'), num=num+1; currentpos = currentpos+1; end
                chans(end+1) = chans(end);
                chans(num+1:end) = chans(num:end-1);
                for index = 1:length( allfields )
                    chans = setfield(chans, {num}, allfields{index}, tmpargs{index+1});
                end
                if isfield(chans, 'datachan')
                    if isempty(chans(num).datachan)
                        chans(num).datachan = 0;
                    end
                end
            case 'changechan'
                tmpargs = args{ curfield+1 };
                num = tmpargs{1};
                allfields = fieldnames(chans);
                if length( tmpargs ) < length(allfields)+1
                    error('pop_chanedit: not enough arguments to change all field values');
                end
                for index = 1:length( allfields )
                    eval([ 'chans(' int2str(num) ').' allfields{index} '=' reformat(tmpargs{index+1}) ';' ]);
                end
                
            case 'load'
                if ~isempty(fig) % GUI
                    [tmpf tmpp] = uigetfile('*.*', 'Load a channel location file'); 
                    drawnow;
                    if ~isequal(tmpf, 0),
                        tmpformats = readlocs('getinfos');
                        tmpformattype = { 'autodetect' tmpformats(1:end-1).type };
                        tmpformatstr  = { 'autodetect' tmpformats(1:end-1).typestring };
                        tmpformatdesc = { 'Autodetect file format from file extension' tmpformats(1:end-1).description };
                        %cb_listbox    = 'tmpdesc=get(gcbf, ''userdata''); set(findobj(gcbf, ''tag'', ''strdesc''), ''string'',  strmultiline([ ''File format: '' tmpdesc{get(gcbo, ''value'')} ], 30, 10)); clear tmpdesc;'' } }, ''pophelp(''''readlocs'''')'',' ...
                	    %        'Read electrode file'', tmpformatdesc, ''normal'', 4);
                        %txtgui = [ strmultiline([ 'File format: Autodetect file format from file extension'], 20, 10) 10 10 ];
                        %tmpfmt = inputgui( 'geometry', {[1 1]}, ...
                        %                   'uilist'  , { { 'style', 'text', 'string', txtgui 'tag' 'strdesc' }, ...
                        %                                 { 'style', 'listbox', 'string', strvcat(tmpformatstr) 'callback' '' } }, ...
                        %                   'geomvert', [10], ...
                        %                   'helpcom' , 'pophelp(''readlocs'');');
                        tmpfmt = inputgui( 'geometry', {[1 1 1] [1]}, ...
                                           'uilist'  , { { 'style', 'text', 'string', 'File format:' 'tag' 'strdesc' } {} {}, ...
                                                         { 'style', 'listbox', 'string', strvcat(tmpformatstr) 'callback' '' } }, ...
                                           'geomvert', [1 8], ...
                                           'helpcom' , 'pophelp(''readlocs'');');
                        if isempty(tmpfmt), 
                             args{ curfield+1 } = [];
                        else args{ curfield+1 } = { fullfile(tmpp, tmpf) 'filetype' tmpformattype{tmpfmt{1}} };
                        end
                    else args{ curfield+1 } = [];
                	end
                end
                
                tmpargs = args{ curfield+1 };
                if ~isempty(tmpargs),
                    if ischar(tmpargs)
                        [chans] = readlocs(tmpargs);
                        [tmp tmp2 chans]  = eeg_checkchanlocs(chans);
                        chaninfo          = [];
                        chaninfo.filename = tmpargs;
                    else
                        [chans] = readlocs(tmpargs{:});
                        [tmp tmp2 chans]  = eeg_checkchanlocs(chans);
                        chaninfo          = [];
                        chaninfo.filename = tmpargs{1};
                    end
                    
                    % backup file content etc...
                    % --------------------------
                    tmptext         = loadtxt( chaninfo.filename, 'delim', [], 'verbose', 'off', 'convert', 'off');
                    chaninfo.filecontent = strvcat(tmptext{:});
                    
                    % set urchan structure
                    % --------------------
                    urchans = chans;
                    for index = 1:length(chans)
                        chans(index).urchan = index;
                    end
                end
                if ~isfield(chans, 'datachan')
                    chans(1).datachan = [];
                end
                for index = 1:length(chans)
                    if isempty(chans(index).datachan)
                        chans(index).datachan = 1;
                    end
                end
                
            case 'eval'
                tmpargs = args{ curfield+1 };
                eval(tmpargs);
                
            case 'saveothers'
                com = pop_writelocs(chans);
                args{ curfield }   = 'eval';
                args{ curfield+1 } = com;

            case 'save'
                if ~isempty(fig)
                     [tmpf tmpp] = uiputfile('*.ced', 'Save channel locs in EEGLAB .ced format');
                     drawnow;
                     args{ curfield+1 } = fullfile(tmpp, tmpf);
                end
                tmpargs = args{ curfield+1 };
                [tmp, tmp, filext] = fileparts(tmpargs);
                if isempty(tmpargs), return; end
                fid = fopen(tmpargs, 'w');
                if fid ==-1, error('Cannot open file'); end
                
                allfields = fieldnames(chans);
                fields = { 'labels' 'theta' 'radius' 'X' 'Y' 'Z' 'sph_theta' 'sph_phi' 'sph_radius' 'type' };
                tmpdiff = setdiff(fields, allfields);
                if ~isempty(tmpdiff), error(sprintf('Field "%s" missing in channel location structure', tmpdiff{1})); end
                if strcmpi(filext,'.ced')
                    origfield = fields;  % Setting the order as expected in readlocs.m
                else
                    origfield = intersect(allfields,fields,'stable');  % Getting the original order from file
                end
                fprintf(fid, 'Number\t');
                for field = 1:length(fields)
                    fprintf(fid, '%s\t', origfield{field});
                end
                fprintf(fid, '\n');
                for index=1:length(chans)
                    fprintf(fid, '%d\t',  index);
                    for field = 1:length(fields)
                        tmpval = getfield(chans, {index}, origfield{field});
                        if ischar(tmpval)
                            fprintf(fid, '%s\t',  tmpval);
                        else
                            fprintf(fid, '%3.3g\t',  tmpval);
                        end
                    end
                    fprintf(fid, '\n');
                end
                if isempty(tmpargs), chantmp = readlocs(tmpargs); end
                
            case 'nosedir'
                nosevals = { '+X' '-X' '+Y' '-Y' };
                if ~isempty(fig)
                    tmpval = get(findobj(gcbf, 'tag', 'nosedir'), 'value');
                    args{ curfield+1 } = nosevals{tmpval};
                    warndlg2( [ 'Changing the nose direction will force EEGLAB to physically rotate ' 10 ...
                                'electrodes, so next time you call this interface, nose direction will' 10 ...
                                'be +X. If your electrodes are currently aligned with a specific' 10 ...
                                'head model, you will have to rotate them in the model coregistration' 10 ... 
                                'interface to realign them with the model.'], 'My Warn Dialog');
                end
                chaninfo.nosedir = args{ curfield+1 };
                if isempty(strmatch(chaninfo.nosedir, nosevals))
                    error('Wrong value for nose direction');
                end
                
            case { 'lookup' 'lookupgui' }
                if strcmpi(lower(args{curfield}), 'lookupgui')
                    standardchans = { 'Fp1' 'Fpz' 'Fp2' 'Nz' 'AF9' 'AF7' 'AF3' 'AFz' 'AF4' 'AF8' 'AF10' 'F9' 'F7' 'F5' ...
                        'F3' 'F1' 'Fz' 'F2' 'F4' 'F6' 'F8' 'F10' 'FT9' 'FT7' 'FC5' 'FC3' 'FC1' 'FCz' 'FC2' ...
                        'FC4' 'FC6' 'FT8' 'FT10' 'T9' 'T7' 'C5' 'C3' 'C1' 'Cz' 'C2' 'C4' 'C6' 'T8' 'T10' ...
                        'TP9' 'TP7' 'CP5' 'CP3' 'CP1' 'CPz' 'CP2' 'CP4' 'CP6' 'TP8' 'TP10' 'P9' 'P7' 'P5' ...
                        'P3' 'P1' 'Pz' 'P2' 'P4' 'P6' 'P8' 'P10' 'PO9' 'PO7' 'PO3' 'POz' 'PO4' 'PO8' 'PO10' ...
                        'O1' 'Oz' 'O2' 'O9' 'O10' 'CB1' 'CB2' 'Iz' };
                    for indexchan = 1:length(chans)
                        if isempty(chans(indexchan).labels), chans(indexchan).labels = ''; end
                    end
                    [tmp1 ind1 ind2] = intersect_bc( lower(standardchans), {chans.labels});
                    if ~isempty(tmp1) || isfield(chans, 'theta')

                        % finding template location files
                        % -------------------------------
                        setmodel = [ 'tmpdat = get(gcbf, ''userdata'');' ...
                            'tmpval = get(gcbo, ''value'');' ...
                            'set(findobj(gcbf, ''tag'', ''elec''), ''string'', tmpdat{tmpval});' ...
                            'clear tmpval tmpdat;' ];
                        try
                            EEG = eeg_emptyset; % for dipfitdefs
                            dipfitdefs;
                            userdatatmp = { template_models(1).chanfile template_models(2).chanfile 'Standard-10-5-Cap385_witheog.elp' }; % last file in the path (see eeglab.m)
                            clear EEG;
                        catch, userdatatmp = { 'Standard-10-5-Cap385.sfp' 'Standard-10-5-Cap385.sfp' 'Standard-10-5-Cap385_witheog.elp' }; % files are in the path (see eeglab.m)
                        end

                        % other commands for help/load
                        % ----------------------------
                        comhelp = [ 'warndlg2(strvcat(''The template file depends on the model'',' ...
                            '''you intend to use for dipole fitting. The default file is fine for'',' ...
                            '''spherical model.'');' ];
                        commandload = [ '[filename, filepath] = uigetfile(''*'', ''Select a text file'');' ...
                            'if filename ~=0,' ...
                            '   set(findobj(''parent'', gcbf, ''tag'', ''elec''), ''string'', [ filepath filename ]);' ...
                            'end;' ...
                            'clear filename filepath tagtest;' ];
                        if ~isfield(chans, 'theta'),                    message =1;
                        elseif all(cellfun('isempty', {chans.theta })), message =1;
                        else                                            message =2;
                        end
                        if message == 1
                            textcomment = strvcat('Only channel labels are present currently, but some of these labels have known', ...
                                'positions. Do you want to look up coordinates for these channels using the electrode', ...
                                'file below? If you have a channel location file for this dataset, press cancel, then', ...
                                'use button "Read location" in the following gui. If you do not know, just press OK.');
                        else
                            textcomment = strvcat('Some channel labels may have known locations.', ...
                                'Do you want to look up coordinates for these channels using the electrode', ...
                                'file below? If you do not know, press OK.');
                        end
                        uilist = { { 'style' 'text' 'string' textcomment } ...
                            { 'style' 'popupmenu'  'string' [ 'use BESA file for 4-shell dipfit spherical model' ...
                            '|use MNI coordinate file for BEM dipfit model|Use spherical file with eye channels' ] ...
                            'callback' setmodel 'value' 2 } ...
                            { } ...
                            { 'style' 'edit'       'string' userdatatmp{2} 'tag' 'elec' } ...
                            { 'style' 'pushbutton' 'string' '...' 'callback' commandload } ...
                            { } };
%                             { 'Style', 'checkbox', 'value', 0, 'string','Overwrite Original Channels' } };

                        res = inputgui( { 1 [1 0.3] [1 0.3] 1 }, uilist, 'pophelp(''pop_chanedit'')', 'Look up channel locations?', userdatatmp, 'normal', [4 1 1 1] );
                        if ~isempty(res)
                            chaninfo.filename = res{2};
                            args{ curfield   } = 'lookup';
                            args{ curfield+1 } = res{2};
                            com = args;
                        else
                            return;
                        end
                    end
                else    
                    chaninfo.filename = args{ curfield+1 };
                end
                if strcmpi(chaninfo.filename, 'standard-10-5-cap385.elp')
                    dipfitdefs;
                    chaninfo.filename = template_models(1).chanfile;
                elseif strcmpi(chaninfo.filename, 'standard_1005.elc')
                    dipfitdefs;
                    chaninfo.filename = template_models(2).chanfile;
                elseif strcmpi(chaninfo.filename, 'standard_1005.ced')
                    dipfitdefs;
                    chaninfo.filename = template_models(2).chanfile;
                end
                tmplocs = readlocs( chaninfo.filename, 'defaultelp', 'BESA' );                
                for indexchan = 1:length(chans)
                    if isempty(chans(indexchan).labels), chans(indexchan).labels = ''; end
                end
                [tmp ind1 ind2] = intersect_bc(lower({ tmplocs.labels }), lower({ chans.labels }));
                if ~isempty(tmp)
                    chans = struct('labels', { chans.labels }, 'datachan', { chans.datachan }, 'type', { chans.type });
                    [ind2 ind3] = sort(ind2);
                    ind1 = ind1(ind3);
                    
                    for index = 1:length(ind2)
                        chans(ind2(index)).theta      = tmplocs(ind1(index)).theta;
                        chans(ind2(index)).radius     = tmplocs(ind1(index)).radius;
                        chans(ind2(index)).X          = tmplocs(ind1(index)).X;
                        chans(ind2(index)).Y          = tmplocs(ind1(index)).Y;
                        chans(ind2(index)).Z          = tmplocs(ind1(index)).Z;
                        chans(ind2(index)).sph_theta  = tmplocs(ind1(index)).sph_theta;
                        chans(ind2(index)).sph_phi    = tmplocs(ind1(index)).sph_phi;
                        chans(ind2(index)).sph_radius = tmplocs(ind1(index)).sph_radius;
                    end
                    if isfield(tmplocs, 'type')
                        for index = 1:length(ind2)
                            chans(ind2(index)).type   = tmplocs(ind1(index)).type;
                        end
                    end
                    
                    tmpdiff = setdiff_bc([1:length(chans)], ind2);
                    if ~isempty(tmpdiff)
                        fprintf('Channel lookup: no location for ');
                        for index = 1:(length(tmpdiff)-1)
                            fprintf('%s,', chans(tmpdiff(index)).labels);
                        end
                        fprintf('%s\nSend us standard location for your channels at eeglab@sccn.ucsd.edu\n', ...
                            chans(tmpdiff(end)).labels);
                    end
                    if ~isfield(chans, 'type'), chans(1).type = []; end
                end
                if ~isempty(findstr(args{ curfield+1 }, 'standard_10')) && ...
                        ~isempty(findstr(args{ curfield+1 }, '.elc'))
                    chaninfo.nosedir = '+Y';
                else
                    chaninfo.nosedir = '+X';
                end
                if flag_replurchan, urchans = eeg_checkchanlocs(chans, chaninfo); end
                for index = 1:length(chans)
                    chans(index).urchan    = index;
                    chans(index).ref       = '';
                end
        end
    end
    
end

% call from a figure
% ------------------
if ~isempty(fig)
    userdata.chans    = chans;
    userdata.chaninfo = chaninfo;
    userdata.commands = { userdata.commands{:} args{:} };
    userdata.urchans  = urchans;
    set(fig, 'userdata', userdata);

    set(findobj(fig, 'tag', 'chaneditnumval'), 'string', num2str(currentpos));
    set(findobj(fig, 'tag', 'chaneditscantitle'), 'string', ['Channel number (of ' int2str(length(chans)) ')']);

    % update GUI with current channel info
    allfields = fieldnames(chans);
    if ~isempty(chans)
        for index = 1:length(allfields)
            obj = findobj(fig, 'tag', [ 'chanedit' allfields{index}]);
            if strcmpi(allfields{index}, 'datachan') 
                set(obj, 'value', getfield(chans(currentpos), allfields{index}));
            else
                tmpval = getfield(chans(currentpos), allfields{index});
                if ischar(tmpval) && strcmpi(tmpval, '[]'), tmpval = ''; end
                set(obj, 'string', num2str(tmpval));
            end
        end
    else
        for index = 1:length(allfields)
            obj = findobj(fig, 'tag', [ 'chanedit' allfields{index}]);
            if strcmpi(allfields{index}, 'datachan') 
                set(obj, 'value', 0);
            else
                set(obj, 'string', '');
            end
        end
    end
else
    [chans chaninfo] = eeg_checkchanlocs(chans, chaninfo);
    if dataset_input,
         if nchansori == length(chans)
             for index = 1:length(EEG)
                 EEG(index).chanlocs = chans;
                 EEG(index).chaninfo = chaninfo;
             end
             % Updating urchanlocs            
             if flag_replurchan && ~isempty(urchans), EEG.urchanlocs = urchans; end
             EEG = eeg_checkset(EEG); % for channel orientation
         else
             disp('Channel structure size not consistent with the data so changes will be ignored');
             disp('Use the function pop_select(EEG, ''nochannel'', [x]); if you wish the remove data channels');
         end
         try chansout = EEG; catch, end
    else chansout = chans;
    end
end

return;

% format the output field
% -----------------------
function strval = reformat( val )
if isnumeric(val) && isempty(val), val = '[]'; end
if ischar(val), strval = [ '''' val '''' ];
else           strval = num2str(val);
end

% extract text using tokens (not used)
% ------------------------------------
function txt = inserttxt( txt, tokins, tokfind);
locfind = findstr(txt, tokfind);
for index = length(locfind):-1:1
    txt = [txt(1:locfind(index)-1) tokins txt(locfind(index):end)];
end

% ask for confirmation
% --------------------
function num = popask( text )
ButtonName=questdlg2( text, ...
    'Confirmation', 'Cancel', 'Yes','Yes');
switch lower(ButtonName),
    case 'cancel', num = 0;
    case 'yes',    num = 1;
end
