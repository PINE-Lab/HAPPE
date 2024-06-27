function [leadfield,headmodel] = dong_calc_leadfield3(xyz_elec,xyz_dipoles,xyz_dipOri,headmodel)
% Forward program:  calculate electric leadfield matrix for dipoles in concentric spheres
%  based on associated Legendre fucntion.
% -------------------------------------------------------------------------
% Input:
%    xyz_elec: xyz coordinates of electrodes. channels X 3, e.g. 129 channels X 3 xyz coordinates.
%    xyz_dipoles: xyz coordinates of dipoles. dipoles X 3, e.g. 3000 sources X 3 xyz coordinates.
%    xyz_dipOri: oritations of dipoles. if it is emtpy, the leadfield with x, y and z oritations 
%                will be generated. [] Or dipoles X 3, e.g. 3000 sources X 3 xyz coordinates
%    headmodel: head model of 1-3 concentric spheres.
%               headmodel.r: the radiis of spheres which are normailzed by
%                            the largest sphere.E.g. [0.87,0.92,1].
%               headmodel.type: head model type. it is concentric spheres
%               headmodel.cond (optional): the resistivities (Brain, CSF,
%                            Skull and Scalp) which are normalized by the largest
%                            resistivity. E.g. For 3 concentri spheres, default is [1,0.0125,1];
%               headmodel.tissue(optional): tissues. E.g. { 'brain' 'skull' 'scalp' };
%               headmodel.o: Center of the spheres (optional if
%                            origin).e.g. [0,0,0].
%               headmodel.terms: the constants for the Legendre expansion for EEG leadfields.
%                            The returning value is the constant for the third layer on the outer
%                            surface, e.g. K(3,r=R3).
% Output: leadfield: leadfield matrix of concentric sphere model.
%               No. of channels X No. of dipoles,e.g. 129 channels X 3000 sources
%               OR 
%               No. of channels X No. of dipoles*3 (correponding to x,y,
%               and z oritations),  leadfield.X->129 channels X 3000 sources
% -------------------------------------------------------------------------
% Yao, D. (2000). "High-resolution EEG mappings: a spherical harmonic spectra
%    theory and simulation results." Clin Neurophysiol 111(1): 81-92.

% Dong, L., et al. (2017). "MATLAB Toolboxes for Reference Electrode Standardization 
%                  Technique (REST) of Scalp EEG." Front Neurosci 11.
% Cuffin & Cohen 1979 Electroencephalogr. Clin. Neurophysiol., 47:131-146.

% Based on FieldTrip 20160222 functions: * eeg_leadfield4 by Robert Oostenveld
% Helped by Ricardo Bruna

% Copyright (C) 2019.8, Li Dong (Lidong@uestc.edu.cn)
% -------------------------------------------------------------------------
if ~isfield(headmodel,'r')
    error('radii values should be specified for each tissue type!');
end

% Defines the center of the spheres, if needed.
if ~isfield ( headmodel, 'o' )
    headmodel.o     = zeros ( 1, 3 );
end
% Centers the sphere in the origin.
xyz_dipoles = bsxfun ( @minus, xyz_dipoles,  headmodel.o );
xyz_elec    = bsxfun ( @minus, xyz_elec,  headmodel.o );
headmodel.o = zeros ( 1, 3 );

% check radii of spheres
headmodel.r = (headmodel.r(:))';
Nr = length(headmodel.r); % No. of spheres
if Nr == 3
    disp('3 concentric spheres are used');
elseif Nr == 1
    headmodel.r = headmodel.r ( [ 1 1 1 ] ); % If only one sphere creates the other two.
    disp('1 concentric spheres are used');
elseif Nr == 2
    headmodel.r(1,3) = headmodel.r(1,2);
    disp('2 concentric spheres are used');
elseif Nr > 3
    error('only 3 concentric spheres are supported');
end

[r_temp, indx] = sort(headmodel.r(:)); % sort the spheres from the smallest to the largest
headmodel.r = (r_temp./max(r_temp))'; % normailzed by the largest sphere.

% check xyz coordinates of electrodes
if isempty(xyz_elec)
    error('xyz coordinates of electrodes is empty');
end
Nelec = size(xyz_elec,1); % No. of electrodes
disp(['No. of electrodes: ',num2str(Nelec)]);

% check whether the electrode ly on the sphere, allowing 0.5% tolerance
dist = sqrt(sum(xyz_elec.^2,2));
if any(abs(dist-max(headmodel.r))>max(headmodel.r)*0.005)
    disp('electrodes do not ly on sphere surface -> using projection')
end

xyz_elec = max(headmodel.r) * xyz_elec ./ repmat(dist,1,3); % normalize to unit sphere 把坐标归一化到单位球表面，避免原始坐标不是在单位表面的

% check dipole coordinates
Ndips = size(xyz_dipoles, 1); % No. of dipoles
disp(['No. of diples: ',num2str(Ndips)]);

radii_dipoles = sqrt ( sum ( xyz_dipoles .^ 2, 2 ) );
if any(radii_dipoles >= min(headmodel.r))
    disp('some dipoles are not within sphere 1 -> rescale xyz coordinates of all dipoles');
    xyz_dipoles = ((min(headmodel.r)-0.001)./max(radii_dipoles))*xyz_dipoles;
end

% check conductivities 
if isfield(headmodel,'cond')
    if isempty(headmodel.cond)
        % it being empty indicates that the user did not specify a conductivity,
        % use a default instead
        if Nr == 1
            disp('default conductivities are used: 1-->[sphere1]')
            headmodel.cond  = [ 1,1,1 ];
            headmodel.tissue = {'sphere1'};
        elseif Nr == 2
            disp('default normailzed conductivities are used: [0.0125,1]-->[sphere1,sphere2]')
            headmodel.cond = [0.0125,1,1];
            headmodel.tissue = { 'sphere1','sphere2'};
        elseif Nr == 3
            disp('default normailzed conductivities are used: [1,0.0125,1]-->[brain,skull,scalp]')
            headmodel.cond = [1,0.0125,1];   % Brain,Skull and Scalp
            headmodel.tissue = { 'brain','skull','scalp' };
        elseif Nr >3
            error('more than 3 concentric spheres are not supported');
        end
    else
        headmodel.cond = (headmodel.cond(:))';
        if length(headmodel.cond) == 1
           headmodel.cond = headmodel.cond ( [ 1,1,1,] );
           headmodel.tissue = { 'sphere1'};
        elseif length(headmodel.cond) == 2
           headmodel.cond(1,3) = headmodel.cond (1,2);
           headmodel.tissue = { 'sphere1','sphere2'};
        elseif length(headmodel.cond) == 3 && ~isfield(headmodel,'tissue')
            headmodel.tissue = { 'brain','skull','scalp' };
        elseif length(headmodel.cond) > 3
            error('more than 3 conductivity values of concentric spheres are not supported');
        end
    end
else
    % it is not in the field indicates that the user did not specify a conductivity,
    % use a default instead
    if Nr == 1
        disp('default conductivities are used: 1-->[sphere1]')
        headmodel.cond  = [1,1,1];
        headmodel.tissue = {'sphere1'};
    elseif Nr == 2
        disp('default normailzed conductivities are used: [0.0125,1]-->[sphere1,sphere2]')
        headmodel.cond = [0.0125,1,1];
        headmodel.tissue = { 'sphere1','sphere2'};
    elseif Nr == 3
        disp('default normailzed conductivities are used: [1,0.0125,1]-->[brain,skull,scalp]')
        headmodel.cond = [1,0.0125,1];   % Brain,Skull and Scalp
        headmodel.tissue = { 'brain','skull','scalp' };
    elseif Nr > 3
        error('more than 3 concentric spheres are not supported');
    end
end

cond = headmodel.cond(indx);
headmodel.cond = cond./max(cond); % normalized by the largest resistivity.

% -------------------------------------------------------------------------
% use more convenient names for the radii and conductivities
r1 = headmodel.r(1); c1 = headmodel.cond(1);
try r2 = headmodel.r(2); c2 = headmodel.cond(2);catch;end;
try r3 = headmodel.r(3); c3 = headmodel.cond(3);catch;end;
% try r4 = headmodel.r(4); c4 = headmodel.cond(4);catch;end;

% Computes the constant factors for the sphere, if needed.
if isfield ( headmodel, 'terms' )
    if isempty(headmodel.terms)
        nterms = 60; % cut off tail at Lmax times loop (the dgree of associated Legendre function);
        degrees = 1:nterms;
        % calculating coefficents of Alm, Blm, Elm in equtions 7-10 (Yao 2000)
        % equations (6)-(11)
        s1 = c2/c1;
        s2 = c2/c3;
        f1 = r2/r3;
        chi1 = (degrees./(degrees+1)).*r3.^(2*degrees+1);
        
        alpha1 = f1.^(2*degrees+1)*(1-s2) - (1+degrees*s2./(1+degrees));
        beta1 = f1.^(-2*degrees-1)*(1-s2) - (1 + s2*(degrees+1)./degrees);
        lamda = alpha1./beta1;
        
        Alm = ((2*degrees+1).*((1+lamda.*r1.^(-2*degrees-1))))./(degrees*(1-s1).*r1.^(2*degrees+1)+lamda.*(degrees+s1*(degrees+1))) - r1.^(-2*degrees-1);
        Blm = (2*degrees+1)./(degrees*(1-s1).*r1.^(2*degrees+1) + lamda.*(degrees+s1*(degrees+1)));
        Elm = Blm.*(r2.^(2*degrees+1)+lamda).*(degrees+1)./((degrees+1).*r2.^(2*degrees+1)+degrees.*r3.^(2*degrees+1));
        
        K1 = Alm.*r1.^(2*degrees+1) + 1;
        K2 = Blm.*(r2.^(2*degrees+1) + lamda);
        K3 = Elm.*(r3.^(2*degrees+1) + chi1);
        if Nr == 1
            headmodel.terms = K1;
            headmodel.r = headmodel.r(1);
            headmodel.cond = headmodel.cond(1);
        elseif Nr == 2
            headmodel.terms = K2;
            headmodel.r = headmodel.r(1:end-1);
            headmodel.cond = headmodel.cond(1:end-1);
        elseif Nr == 3
            headmodel.terms = K3;
        end
    end
else
    nterms = 60; % cut off tail at Lmax times loop (the dgree of associated Legendre function);
    degrees = 1:nterms;
    % calculating coefficents of Alm, Blm, Elm in equtions 7-10 (Yao 2000)
    % equations (6)-(11)
    s1 = c2/c1;
    s2 = c2/c3;
    f1 = r2/r3;
    chi1 = (degrees./(degrees+1)).*r3.^(2*degrees+1);
    
    alpha1 = f1.^(2*degrees+1)*(1-s2) - (1+degrees*s2./(1+degrees));
    beta1 = f1.^(-2*degrees-1)*(1-s2) - (1 + s2*(degrees+1)./degrees);
    lamda = alpha1./beta1;
    
    Alm = ((2*degrees+1).*((1+lamda.*r1.^(-2*degrees-1))))./(degrees*(1-s1).*r1.^(2*degrees+1)+lamda.*(degrees+s1*(degrees+1))) - r1.^(-2*degrees-1);
    Blm = (2*degrees+1)./(degrees*(1-s1).*r1.^(2*degrees+1) + lamda.*(degrees+s1*(degrees+1)));
    Elm = Blm.*(r2.^(2*degrees+1)+lamda).*(degrees+1)./((degrees+1).*r2.^(2*degrees+1)+degrees.*r3.^(2*degrees+1));
    
    K1 = Alm.*r1.^(2*degrees+1) + 1;
    K2 = Blm.*(r2.^(2*degrees+1) + lamda);
    K3 = Elm.*(r3.^(2*degrees+1) + chi1);
    if Nr == 1
        headmodel.terms = K1;
        headmodel.r = headmodel.r(1);
        headmodel.cond = headmodel.cond(1);
    elseif Nr == 2
        headmodel.terms = K2;
        headmodel.r = headmodel.r(1:end-1);
        headmodel.cond = headmodel.cond(1:end-1);
    elseif Nr == 3
        headmodel.terms = K3;
    end
end
% -------------------------------------------------------------------------
% This implementation is adapted from:
% Cuffin & Cohen 1979 Electroencephalogr. Clin. Neurophysiol., 47:131-146.
% Based on FieldTrip 20160222 functions: * eeg_leadfield4 by Robert Oostenveld

% Gets the radius and the conductivity of the outter sphere.
oradius     = headmodel.r    (end);
ocond       = headmodel.cond (end);

% Gets the terms of the expansion.
terms       = headmodel.terms;

% Extracts the number of terms of the expansion.
nterms      = numel ( terms );
degrees     = 1: nterms;

% Calculates Gamma (Eqs. A2, A3, Cuffin & Cohen 1979).
gamma       = 1 ./ ( degrees .* terms ./ ( 2 .* degrees + 1 ) .^ 4 );

% % Extracts the number of dipoles and sensors.
% Ndips       = size ( xyz_dipoles, 1 );
% Nelec       = size ( xyz_elec, 1 );

% Calculates the radius of each dipole relative to the outter sphere.
radii = sqrt ( sum ( xyz_dipoles .^ 2, 2 ) );

% Generates the costant factors for each dipole.
% Cuffin & Cohen 1979. Eq A2. Part 1.
const1      = ( 2 * degrees + 1 ) .^ 4;
const2      = bsxfun ( @power, radii / oradius, degrees - 1 );
const3      = gamma .* 4 * pi * ocond * oradius ^ 2;
const       = bsxfun ( @rdivide, bsxfun ( @times, const1, const2 ), const3 );

% Initializes the leadfield matrix.
leadfield = zeros ( Nelec, 3, Ndips);

% Goes through each dipole.
for nd = 1: Ndips
    % Gets the current dipole.
    dippos      = xyz_dipoles ( nd, : );

    % Rotates the system in order to put the dipole in the positive z-axis.
    % This separates the radial (z) and tangential (xy) components.
    
    % Initializes the rotation matrix.
    rot         = eye ( 3, 'single' );
    
    % If dipole not in the z-axis rotates it.
    if dippos (1) ~= 0 || dippos (2) ~= 0
        val1         = norm ( dippos );
        val2         = norm ( dippos ( 1: 2 ) );
        rot ( 1, 1 ) = dippos (1) * dippos (3) / ( val1 * val2 );
        rot ( 1, 2 ) = dippos (2) * dippos (3) / ( val1 * val2 );
        rot ( 1, 3 ) = -1.0 * val2 / val1;
        rot ( 2, 1 ) = -1.0 * dippos (2) / val2;
        rot ( 2, 2 ) =        dippos (1) / val2;
        rot ( 2, 3 ) =                  0;
        rot ( 3, : ) = dippos / val1;
        
    % If dipole in negative z rotates it around the x-axis.
    elseif dippos (3) < 0
        rot ( 2, 2 ) = -1;
        rot ( 3, 3 ) = -1;
    end
    
    % Creates a rotated version of the electrodes definition.
    dipsens     = xyz_elec * rot';
    
    
    % Gets the electrode's position in spherical coordinates.
    [ phi, theta ] = cart2sph ( dipsens ( :, 1 ), dipsens ( :, 2 ), dipsens ( :, 3 ) );
    cos_theta   = double ( cos ( pi / 2 - theta ) );
    
    % Calculates the 0th and 1st order Legendre polynomials.
    P0          = my_plgndr ( nterms, 0, cos_theta );
    P1          = my_plgndr ( nterms, 1, cos_theta );
    
    % Discards PX_0.
    P0 ( :, 1 ) = [];
    P1 ( :, 1 ) = [];
    
    % Corrects P1.
    % Cuffin & Cohen 1979. Eq A2. Part 2.
    P1          = bsxfun ( @rdivide, P1, 1: nterms );
    
    % Creates the radial and tangential components.
    % Cuffin & Cohen 1979. Eq A2. Part 3.
    s_r         = bsxfun ( @times, const ( nd, : ), P0 );
    s_t         = bsxfun ( @times, const ( nd, : ), P1 );
    
    % Creates the leadfield for each direction of the dipole.
    % Cuffin & Cohen 1979. Eq A2. Part 4.
    lf          = zeros ( Nelec, 3);
    lf ( :, 1 ) = sum ( s_t, 2 ) .* -cos ( phi );
    lf ( :, 2 ) = sum ( s_t, 2 ) .* -sin ( phi );
    lf ( :, 3 ) = sum ( s_r, 2 );
    
    % Rotates back the leadfield and stores it.
    leadfield ( :, :, nd ) = lf * rot;
end

% % plot check
% figure;
% k3 = 2800;
% subplot(2,3,1)
% lf_temp1 = leadfield(:,1,k3);
% quiver3(xyz_dipoles(k3,1),xyz_dipoles(k3,2),xyz_dipoles(k3,3), ...
%     xyz_dipOri(k3,1),xyz_dipOri(k3,2),xyz_dipOri(k3,3),1, 'color','r');
% ft_plot_topo3d(xyz_elec,lf_temp1);
% title('leadfield-X');
% 
% subplot(2,3,2)
% lf_temp2 = leadfield(:,2,k3);
% quiver3(xyz_dipoles(k3,1),xyz_dipoles(k3,2),xyz_dipoles(k3,3), ...
%     xyz_dipOri(k3,1),xyz_dipOri(k3,2),xyz_dipOri(k3,3),1, 'color','r');
% ft_plot_topo3d(xyz_elec,lf_temp2);
% title('leadfield-Y');
% 
% subplot(2,3,3)
% lf_temp3 = leadfield(:,3,k3);
% quiver3(xyz_dipoles(k3,1),xyz_dipoles(k3,2),xyz_dipoles(k3,3), ...
%     xyz_dipOri(k3,1),xyz_dipOri(k3,2),xyz_dipOri(k3,3),1, 'color','r');
% ft_plot_topo3d(xyz_elec,lf_temp3);
% title('leadfield-Z');
% 
% subplot(2,3,4)
% lf_temp1 = lf_cuffin(:,1,k3);
% quiver3(xyz_dipoles(k3,1),xyz_dipoles(k3,2),xyz_dipoles(k3,3), ...
%     xyz_dipOri(k3,1),xyz_dipOri(k3,2),xyz_dipOri(k3,3),1, 'color','r');
% ft_plot_topo3d(xyz_elec,lf_temp1);
% title('leadfield-cuffin-X');
% 
% subplot(2,3,5)
% lf_temp2 = lf_cuffin(:,2,k3);
% quiver3(xyz_dipoles(k3,1),xyz_dipoles(k3,2),xyz_dipoles(k3,3), ...
%     xyz_dipOri(k3,1),xyz_dipOri(k3,2),xyz_dipOri(k3,3),1, 'color','r');
% ft_plot_topo3d(xyz_elec,lf_temp2);
% title('leadfield-cuffin-Y');
% 
% subplot(2,3,6)
% lf_temp3 = lf_cuffin(:,3,k3);
% quiver3(xyz_dipoles(k3,1),xyz_dipoles(k3,2),xyz_dipoles(k3,3), ...
%     xyz_dipOri(k3,1),xyz_dipOri(k3,2),xyz_dipOri(k3,3),1, 'color','r');
% ft_plot_topo3d(xyz_elec,lf_temp3);
% title('leadfield-cuffin-Z');

% Reshapes the leadfield in matrix form.
% leadfield = leadfield ( :, : );
% -------------------------------------------------------------------------
if isempty(xyz_dipOri)
    % Reshapes the leadfield in matrix form with x, y, z oritations.
    leadfield.X = reshape(leadfield ( :, 1, : ),Nelec,Ndips);
    leadfield.Y = reshape(leadfield ( :, 2, : ),Nelec,Ndips);
    leadfield.Z = reshape(leadfield ( :, 3, : ),Nelec,Ndips);
else
    if Ndips == size(xyz_dipOri,1)
        % Projects the source momentum over the radial(oritations defined by user) direction.
        leadfield   = permute ( leadfield, [ 1 3 2 ] );
        leadori     = permute ( xyz_dipOri, [ 3 1 2 ] );
        leadfield   = bsxfun ( @times, leadfield, leadori );
        leadfield   = sum ( leadfield, 3 );
    else
        warning('The No. of normals(oritations) are NOT equal to No. of dipoles, output the X, Y, Z leadfield, respectively.');
        % Reshapes the leadfield in matrix form with x, y, z oritations.
        leadfield.X = reshape(leadfield ( :, 1, : ),Nelec,Ndips);
        leadfield.Y = reshape(leadfield ( :, 2, : ),Nelec,Ndips);
        leadfield.Z = reshape(leadfield ( :, 3, : ),Nelec,Ndips);
    end
end
% -------------------------------------------------------------------------
% % plot check
% lf_exe = load('Lead_Field.dat'); % Leadfield calcualed by leadfield.exe
% 
% k3 = 2950; % dipole i
% 
% figure;
% subplot(2,2,1)
% plot3(xyz_dipoles(2601:3000,1),xyz_dipoles(2601:3000,2),xyz_dipoles(2601:3000,3),'.k');
% hold on;
% quiver3(xyz_dipoles(k3,1),xyz_dipoles(k3,2),xyz_dipoles(k3,3), ...
%     xyz_dipOri(k3,1),xyz_dipOri(k3,2),xyz_dipOri(k3,3),1, 'color','r');
% grid on;
% ft_plot_topo3d(xyz_elec,leadfield(:,k3));
% title('My')
% hold off;
% 
% subplot(2,2,2)
% quiver3(xyz_dipoles(k3,1),xyz_dipoles(k3,2),xyz_dipoles(k3,3), ...
%     xyz_dipOri(k3,1),xyz_dipOri(k3,2),xyz_dipOri(k3,3),1, 'color','r');
% ft_plot_topo3d(xyz_elec,lf_exe(k3,:));
% title('leadfield-exe');
% 
% subplot(2,2,3)
% plot(leadfield(:,k3),'-*b');hold on;
% plot(lf_exe(k3,:),'-*r');
% hold off;
% legend('lf-Dong','lf-exe');
% grid on;
% 
% subplot(2,2,4)
% plot(leadfield(:,k3)-(lf_exe(k3,:))','-g');
% grid on;