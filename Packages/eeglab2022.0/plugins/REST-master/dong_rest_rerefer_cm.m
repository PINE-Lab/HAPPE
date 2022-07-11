function [data_z] = dong_rest_rerefer_cm(data,G,leftelec,rightelec,Aelec)
%   Main function of Reference Electrode Standardization Technique
%   re-refering EEG data from contralateral/ipsilateral mastoids to REST
%   Input: 
%         data:  The EEG potentials with contralateral mastoids reference 
%                (e.g. F1-A2, F3-A2,...F2-A1,F4-A1,...),channels X time points, 
%                e.g. 16 channels X 10000 time points.
%         G:     Lead Field matrix, sources X (EEG channels + 2 mastoids), 
%                e.g. 3000 sources X 18 channels.
%         leftelec: index of left electrodes excluding left mastoid in G. e.g. [1:2:16]
%         rightelec: index of right electrodes excluding right mastoid in G. e.g. [2:2:16]
%         Aelec: index of mastoids [A1,A2] or [A2,A1] in G. e.g. [17,18] 
%                For contralateral mastoids, Aelec should be [right mastoid,left mastoid],i.e. [A2,A1].
%                For ipsilateral mastoids, Aelec should be [left mastoid,right mastoid], i.e. [A1,A2].
%   Output:
%         data_z: The EEG potentials with zero reference, 
%                channels X time points.
%
%   Written by Li Dong (July 20, 2018)
%   For more see http://www.neuro.uestc.edu.cn/rest/
%   Reference: Yao D (2001) A method to standardize a reference of scalp EEG recordings to a point at infinity.
%                       Physiol Meas 22:693?11. doi: 10.1088/0967-3334/22/4/305
%  Li Dong*, Fali Li, Qiang Liu, Xin Wen, Yongxiu Lai, Peng Xu and Dezhong Yao*. 
%              MATLAB Toolboxes for Reference Electrode Standardization Technique (REST) 
%              of Scalp EEG. Frontiers in Neuroscience,  2017:11(601).

if nargin < 5
    error('5 inputs at least!');
end
G = G';

if size(data,1) ~= size(G,1)-2
    error('No. of Channels of lead field matrix and data are NOT matched!');
end

if size(data,1) > size(data,2)
    warning('No. of channels > No. of time points in data???');
end

Gar = G;
Gar(leftelec,:) = G(leftelec,:) - repmat(G(Aelec(1),:),length(leftelec),1); % left electrodes - Aelec(1)
Gar(rightelec,:) = G(rightelec,:) - repmat(G(Aelec(2),:),length(rightelec),1); % right electrodes - Aelec(2)
Gar(Aelec,:) = [];

data_z = G * pinv(Gar,0.05) * data;  % the value 0.05 is for real data; 
                                     % for simulated data, it may be set as zero.