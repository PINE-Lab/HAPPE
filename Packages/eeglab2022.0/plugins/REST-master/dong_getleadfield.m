function G = dong_getleadfield(EEG,chanlist,XYZmastoids)
%  calculating leadfield matrix at once
%    Input:
%       EEG: EEG structure loaded by EEGLAB. XYZ coordinated should be saved
%            in EEG.chanlocs.
%       chanlist: indices of selected channels re-referencing to REST.
%            channels X 1, e.g. 62 channels X 1
%       XYZmastoids: coordinates of left and right mastoids, e.g. 2 X 3.
%            The first row is xyz coordinates of left mastoids, the second
%            row is xyz coordinates of right mastoids. It is optional.
%            Default is empty.
%    Output:
%       G:   leadfield matrix calculated for dipoles in concentric spheres
%            based on associated Legendre fucntion.sources X channels, 
%            e.g. 3000 sources X 62 channels.
% -------------
% Copyright (C) 2020.8, Li Dong (Lidong@uestc.edu.cn)
% -------------
if nargin < 2
  error('Reqiured 2 inputs at least!!!');
elseif nargin == 2
  XYZmastoids = [];
end

% use xyz coordinates in the EEG.chanlocs
if ~isempty(EEG) || ~isempty(EEG.chanlocs)
  if isfield(EEG.chanlocs,'X') && isfield(EEG.chanlocs,'Y') && isfield(EEG.chanlocs,'Z')
    if ~isempty(EEG.chanlocs(1).X) && ~isempty(EEG.chanlocs(1).Y) &&~isempty(EEG.chanlocs(1).Z)
      channs = chanlist; % selected channs
      xyz_elec = zeros(length(channs),3);
      for nc = 1:length(channs)
        xyz_elec(nc,1) = EEG.chanlocs(channs(nc)).X;
        xyz_elec(nc,2) = EEG.chanlocs(channs(nc)).Y;
        xyz_elec(nc,3) = EEG.chanlocs(channs(nc)).Z;
      end
      
      if ~isempty(XYZmastoids) && all(isfinite(XYZmastoids(:)))
        xyz_elec = [xyz_elec;XYZmastoids];
      end
      
    else
      errordlg('EEG coordinates (EEG.chanlocs.X/Y/Z) are empty, please select lead field file OR load channel locations in EEGLAB first!!!!','Data Error');
      return
    end
  else
    errordlg('EEG coordinates (EEG.chanlocs.X/Y/Z) are empty, please select lead field file OR load channel locations in EEGLAB first!!!!','Data Error');
    return
  end
else
  errordlg('EEG or EEG.chanlocs are empty, please select lead field file OR load channel locations in EEGLAB first!!!!','Data Error');
  return
end
% -------------------
% load fixed dipoles and define their oritations
% it can be defined by a file with dipole coordinates
[ProgramPath, ~, ~] = fileparts(which('pop_REST_reref.m'));
xyz_dipoles = load([ProgramPath,filesep,'corti869-3000dipoles.dat']);

% Calculate the dipole orientations.
xyz_dipOri           = bsxfun ( @rdivide, xyz_dipoles, sqrt ( sum ( xyz_dipoles .^ 2, 2 ) ) );
xyz_dipOri ( 2601: 3000, 1 ) = 0;
xyz_dipOri ( 2601: 3000, 2 ) = 0;
xyz_dipOri ( 2601: 3000, 3 ) = 1;
% ------------------
% define headmodel
headmodel        = [];
headmodel.type   = 'concentricspheres';
headmodel.o      = [ 0.0000 0.0000 0.0000 ];
headmodel.r      = [ 0.8700,0.9200,1];
headmodel.cond   = [ 1.0000,0.0125,1];
headmodel.tissue = { 'brain' 'skull' 'scalp' };
% -------------------
% calculate leadfield
[G,~] = dong_calc_leadfield3(xyz_elec,xyz_dipoles,xyz_dipOri,headmodel);
G = G';