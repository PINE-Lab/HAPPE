% pop_dipfit_headmodel() - use Fieldtrip functions to generate headmodel 
%                          from anatomical T1 MRI.
%
% Usage: 
%  >> OUTEEG = pop_dipfit_headmodel( EEG ); % pop up window 
%  >> OUTEEG = pop_dipfit_headmodel( EEG, mriFile, 'key', 'val' ); 
%
% Inputs:
%   INEEG     - input dataset
%   mriFile   - name of a T1 MRI file
%
% Optional inputs:
%  'datatype' - ['eeg'|'meg'] data type. Default is automatically determined
%  'nasion'   - [X Y Z] nasion location in voxel
%  'lpa'      - [X Y Z] left ear fiducial location (lpa) in voxel
%  'rpa'      - [X Y Z] right ear fiducial location (rpa) in voxel
%  'nasion'   - [X Y Z] nasion location in voxel
%  'plotfiducial' - ['nasion'|'lpa'|'rpa'] plot fiducial position on MRI.
%                   Default is empty (nothing plotted). Note that the function
%                   abord after plotting. May also be a cell array.
%  'plotmesh'     - ['brain'|'skull'|'scalp'] plot mesh on top of MRI.
%                   Default is empty (nothiong plotted). 
%  'nasion'   - [X Y Z] nasion location in voxel
%  'ft_volumerealign'  - [struct] option for ft_volumerealign function
%  'ft_volumesegment'  - [struct] option for ft_volumesegment function
%  'ft_prepare_headmodel'  - [struct] option for ft_prepare_headmodel function
%                       Use struct('method', 'dipoli') to use the "dipoli" 
%                       method to extract mesh (Windows and Linux only)
%                       Use struct('numvertices', [800, 1600, 2400]) to set
%                       the number of vertices for the "brain", "skull" and
%                       "skin" surfaces.
% Outputs:
%   OUTEEG - output dataset with an updated dipfit structure.
%
% Authors: Arnaud Delorme, SCCN, 2022

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [EEG, com] = pop_dipfit_headmodel( EEG, mriFile, varargin )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1
  help pop_dipfit_headmodel;
  return
end

meeg = {'EEG' 'MEG'};
if isfield(EEG(1).chaninfo, 'type') && ~isempty(strfind(lower(EEG(1).chaninfo.type, 'meg')))
    meegFlag = 1;
else
    meegFlag = 0;
end

com = '';
if nargin < 2
    str = strvcat('Fisrt, select the subject''s anatomical T1 MRI.', ...
                  'Call the DIPFIT settings after this menu item to align', ...
                  'electrodes with the newly created head model.');

    questdlg2(str, 'Subject''s MRI', 'Continue', 'Continue');
    [filename, filepath] = uigetfile('*', 'Select the MRI file');
    mriFile = fullfile(filepath, filename);

    commandload = [ '[filename, filepath] = uigetfile(''*'', ''Select a text file'');' ...
                    'if filename ~=0,' ...
                    '   set(findobj(''parent'', gcbf, ''tag'', ''mrifile''), ''string'', [ filepath filename ]);' ...
                    'end' ...
                    'clear filename filepath tagtest;' ];

    % sidecar file
    [mriPath,mriFile2] = fileparts(mriFile);
    [~,mriFile2] = fileparts(mriFile2);
    fileSideCar = fullfile(mriPath, [ mriFile2 '.json' ]);
    if exist(fileSideCar, 'file')
        fidStr = 'Use BIDS coordsystem file';
        fidEnable = 'off';
    else
        fidStr = '';
        fidEnable = 'off';
    end
    geometry = { [3 0.5 0.5] [2 1.5 0.5] [2 1.5 0.5] [2 1.5 0.5] [2 1.5 0.5] [2 1.5 0.5] };
    uilist = { { 'style' 'text' 'string' 'Parameters to calculate head model from MRI' 'fontweight' 'bold' } {} ...
               { 'style' 'text' 'string' 'plot' } ...
               { 'style' 'text' 'string' 'MRI file' } ...
               { 'style' 'edit' 'string' mriFile 'tag' 'mrifile' } {} ...
               { 'style' 'text' 'string' 'Nasion X, Y, Z (MRI voxels space)' } ...
               { 'style' 'edit' 'string' fidStr 'enable' fidEnable 'tag' 'nasion'} ...
               { 'style' 'checkbox' 'string' '' 'tag' 'plotnasion' } ...
               { 'style' 'text' 'string' 'Left ear (LPA) X, Y, Z (MRI voxels space)' } ...
               { 'style' 'edit' 'string' fidStr 'enable' fidEnable 'tag' 'lpa'} ...
               { 'style' 'checkbox' 'string' '' 'tag' 'plotlpa' } ...
               { 'style' 'text' 'string' 'Right ear (RPA) X, Y, Z (MRI voxels space)' } ...
               { 'style' 'edit' 'string' fidStr 'enable' fidEnable 'tag' 'rpa'} ...
               { 'style' 'checkbox' 'string' '' 'tag' 'plotrpa' } ...
               { 'style' 'text' 'string' 'Select EEG (3 surfaces) or MEG (1 surface)' } ...
               { 'style' 'popupmenu' 'string' meeg 'value' meegFlag+1 'tag' 'meeg' } ...
               { 'style' 'checkbox' 'string' '' 'tag' 'plotsurf' } };
     
	[ result,~,~,restag] = inputgui( geometry, uilist, 'pophelp(''pop_dipfit_headmodel'')', 'Create headmodel from MRI - pop_dipfit_headmodel');
	if length(result) == 0 return; end
    
    mriFile = restag.mrifile;
    
    options = { 'datatype' meeg{restag.meeg} };
    if ~isequal(restag.nasion, fidStr)
        options = [ options { 'nasion' str2num(restag.nasion) }];
    end
    if ~isequal(restag.lpa, fidStr)
        options = [ options { 'lpa' str2num(restag.lpa) }];
    end
    if ~isequal(restag.rpa, fidStr)
        options = [ options { 'rpa' str2num(restag.rpa) }];
    end

    % plotting
    fiducials = {};
    if restag.plotnasion
        fiducials = [ fiducials { 'nasion' }];
    end
    if restag.plotlpa
        fiducials = [ fiducials { 'lpa' }];
    end
    if restag.plotrpa
        fiducials = [ fiducials { 'rpa' }];
    end
    if ~isempty(fiducials)
        options = [ options { 'plotfiducial' fiducials }];
    end
    if restag.plotsurf
        disp('Only the scalp mesh will be plotted; you may plot others (skull & brain) from the command line')
        options = [ options { 'plotmesh' 'scalp' }];
    end
else
    options = varargin;
end

g = finputcheck(options, { 'datatype' 'string'  {'eeg' 'meg'} meeg{meegFlag+1};
                           'nasion'   'integer' []            [];
                           'lpa'      'integer' []            [];
                           'rpa'      'integer' []            [];
                           'ft_volumerealign'      ''  []            [];
                           'ft_volumesegment'      ''  []            [];
                           'ft_prepare_headmodel'  ''  []            [];
                           'plotfiducial' '' [] {};
                           'plotmesh'     'string'  {'brain' 'skull' 'scalp' ''} '' }, 'pop_dipfit_headmodel');
if ischar(g)
    error(g);
end

mri = ft_read_mri(mriFile);
[mriPath,mriFile2] = fileparts(mriFile);
[~,mriFile2] = fileparts(mriFile2);
fileSideCar = fullfile(mriPath, [ mriFile2 '.json' ]);

% read JSON file
if ~isempty(g.nasion) && ~isempty(g.lpa) && ~isempty(g.rpa)
    if length(g.nasion) ~= 3 || length(g.lpa) ~= 3 || length(g.rpa) ~= 3
        error('Fiducial coordinates must have 3 values (x, y, z)')
    end
    fiducials = [...
        g.nasion;
        g.lpa;
        g.rpa];
elseif exist(fileSideCar, 'file')
    disp('JSON sidecar file found');
    fid = fopen(fileSideCar, 'r');
    if fid == -1
        error('Cannot open file %s', fileSideCar)
    end
    raw = fread(fid,inf);
    fclose(fid);
    coordinates = jsondecode(char(raw'));
    
    if isfield(coordinates, 'AnatomicalLandmarkCoordinates')
        if ~isfield(coordinates.AnatomicalLandmarkCoordinates, 'Nasion') || ...
                ~isfield(coordinates.AnatomicalLandmarkCoordinates, 'LPA') || ...
                ~isfield(coordinates.AnatomicalLandmarkCoordinates, 'RPA')
            error('Some anatomical coordinate missing')
        end

        cfg              = g.ft_volumerealign;
        cfg.method       = 'fiducial';
        % this information has been obtained from the .json associated with the anatomical image
        cfg.fiducial.nas = coordinates.AnatomicalLandmarkCoordinates.Nasion';
        cfg.fiducial.lpa = coordinates.AnatomicalLandmarkCoordinates.LPA';
        cfg.fiducial.rpa = coordinates.AnatomicalLandmarkCoordinates.RPA';
        fiducials = [...
            cfg.fiducial.nas;
            cfg.fiducial.lpa;
            cfg.fiducial.rpa];

        % plot fiducial
        if ~isempty(g.plotfiducial)
            if ~iscell(g.plotfiducial)
                g.plotfiducial = { g.plotfiducial };
            end

            for iFid = 1:length(g.plotfiducial)
                ind = strmatch(lower(g.plotfiducial{iFid}), {'nasion' 'lpa' 'rpa'} );
                if isempty(ind)
                    error('Fiducial not found')
                end
    
                cfg2 = [];
                cfg2.locationcoordinates = 'voxel'; % treat the location as voxel coordinates
                cfg2.location = fiducials(ind,:);
                mri2 = mri;
                mri2.transform = eye(4);
                mri2.transform(:,4) = 1;
                ft_sourceplot(cfg2, mri2);
            end
        end

        mri = ft_volumerealign(cfg, mri);        

    else
        disp('JSON sidecar file found, but it does not contain the AnatomicalLandmarkCoordinates field');
    end
else
    if ~isempty(g.nasion) || ~isempty(g.lpa) || ~isempty(g.rpa)
        error('You must provide 3 fiducials')
    else
        error('Fiducials must be provided')
    end
end

%mri  = ft_convert_coordsys(mri, 'acpc');
% extract volume
meshes        = {'brain','skull','scalp'};
cfg           = g.ft_volumesegment;
cfg.output    = meshes;
segmentedmri  = ft_volumesegment(cfg, mri);

% convert fiducials
fiducials = segmentedmri.transform * [ fiducials ones(size(fiducials,1),1)]';
fiducials = fiducials(1:3,:)';
chanlocs = [];
chanlocs.labels = 'Nasion';
chanlocs.X      = fiducials(1,1);
chanlocs.Y      = fiducials(1,2);
chanlocs.Z      = fiducials(1,3);
chanlocs(end+1).labels = 'LPA';
chanlocs(end).X = fiducials(2,1);
chanlocs(end).Y = fiducials(2,2);
chanlocs(end).Z = fiducials(2,3);
chanlocs(end+1).labels = 'RPA';
chanlocs(end).X = fiducials(3,1);
chanlocs(end).Y = fiducials(3,2);
chanlocs(end).Z = fiducials(3,3);
chanlocs = convertlocs(chanlocs, 'cart2all');

if strcmpi(g.datatype, 'eeg')
    % create head model
    % use ft_prepare_mesh to create boundary element
    % model with specific number of vertices
    cfg        = g.ft_prepare_headmodel;
    if ~isfield(cfg, 'method')
        warning backtrace off
        warning('Using ''bemcp'' method to extract mesh (not the default, but the most portable method)');
        warning backtrace on
        cfg.method ='bemcp'; % You can also specify 'openmeeg', 'bemcp', or another method.
    end
    cfg.tissue={'brain','skull','scalp'};
    headmodel  = ft_prepare_headmodel(cfg, segmentedmri);
else
    segmentedmri2 = segmentedmri;
    segmentedmri2 = rmfield(segmentedmri2, 'skull');
    segmentedmri2 = rmfield(segmentedmri2, 'scalp');
    cfg        = [];
    cfg.method = 'singleshell';
    headmodel  = ft_prepare_headmodel(cfg, segmentedmri2);
end

% plot mesh on MRI
if ~isempty(g.plotmesh)
    ind = strmatch(lower(g.plotmesh), meshes, 'exact' );
    if isempty(ind)
        error('Mesh not found')
    end
    cfg = [];
    cfg.intersectmesh = headmodel.bnd(ind);
    ft_sourceplot(cfg, mri);
end

% save data in DIPFIT structure
EEG.dipfit.mrifile  = mri;
EEG.dipfit.hdmfile  = headmodel;
EEG.dipfit.chanfile = chanlocs;
EEG.dipfit.coordformat = mri.coordsys;
EEG.dipfit.coord_transform = []; %'meg';

com = sprintf('EEG = pop_dipfit_headmodel(EEG, ''%s'', %s);', mriFile, vararg2str(options));

