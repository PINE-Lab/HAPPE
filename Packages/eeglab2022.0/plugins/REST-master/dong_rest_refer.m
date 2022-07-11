function [data_z] = dong_rest_refer(data,G)
%   Main function of Reference Electrode Standardization Technique
%   Corrdinate system: Three-concentric-sphere head volume conductor 
%       model. The triangle shows the nose. The centre of the spheres is 
%       defined as the coordinate origin. The axis directed away from the
%       origin toward the left ear is defined as the -x axis, and that 
%       from the origin to the nasion is the +y axis. The +z axis is defined
%       as the axis that is perpendicular to both these axes and directed 
%       from the origin to the vertex.
%   Input: 
%         data:  The EEG potentials with average reference,channels X time points, 
%                e.g. 62 channels X 10000 time points.The original reference must be
%                average reference.
%         G: Lead Field matrix, sources X channels, e.g. 3000 sources X 62 channels.
%   Output:
%         data_z: The EEG potentials with zero reference, 
%                channels X time points.
%
%  Author: Shiang Hu 
%      Date: Sep. 22, 2016
%  Edit by Li Dong (Sep.24,2016)
%  Edit by Li Dong (Aug. 28, 2016)
%  Copyright (C) 2019.8, Li Dong (Lidong@uestc.edu.cn)
%  For more see http://www.neuro.uestc.edu.cn/rest/
%  Reference: Yao D (2001) A method to standardize a reference of scalp EEG recordings to a point at infinity.
%                       Physiol Meas 22:693?11. doi: 10.1088/0967-3334/22/4/305

if nargin < 2
    error('Please input the Lead Field matrix!');
end
G = G';
if size(data,1)~=size(G,1)
    error('No. of Channels of lead field matrix and data are NOT equal!');
end

if size(data,1)>size(data,2)
    warning('No. of channels > No. of time points in data???');
end

Gar = G - repmat(mean(G),size(G,1),1);
data_z = G * pinv(Gar,0.05) * data;  % the value 0.05 is for real data; 
                                     % for simulated data, it may be set as zero.
data_z = data + repmat(mean(data_z),size(G,1),1); % V = V_avg + AVG(V_0)

