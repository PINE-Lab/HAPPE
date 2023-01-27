% calcSNR_PSNR() - Calculate the Signal to Noise Ratio and Peak Signal to 
%                  Noise Ratio in a dataset pre/post- any particular 
%                  processing step.
%
% Usage: 
%   >> [SNR_alldata, PeakSNR_alldata] = calcSNR_PSNR(preEEG, postEEG, order)
%
% Inputs:
%   preEEG  - EEG signal pre-processing step in channels x samples format
%   postEEG - EEG signal post-processing step in channels x samples format
%   order   - The order in which to subtract the pre and post EEGs.
%                 {For 1, preEEG - postEEG; For 2, postEEG - preEEG} 
%
% Outputs:
%   SNR     - Signal to noise ratio over all channels and timepoints
%   PeakSNR - Peak signal to noise ratio over all channels and
%                     timepoints
%
% Author: A.D. Monachino, PINE Lab at Northeastern University, 2021
%         L.J. Gabard-Durnam, PINE Lab at Northeastern University, 2021
%
% This file is part of HAPPE.
% Copyright 2018, 2021 Alexa Monachino, Kelsie Lopez, Laurel Gabard-Durnam
%
% HAPPE is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
% 
% HAPPE is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
% details.
% 
% You should have received a copy of the GNU General Public License along
% with HAPPE. If not, see <https://www.gnu.org/licenses/>.

function [SNR, PeakSNR] = calcSNR_PSNR(preEEG, postEEG, order)
orderMod = (order == 1)*preEEG + (order == 2)*postEEG ;

% Calculate mean squared error per channel as intermediate variable:
squared_differences = (preEEG - postEEG).^2 ;
MSE = mean(squared_differences, 2) ;

% Complete SNR and PSNR calculations:
DEN = sum(squared_differences, 2) ;
NUM = sum(orderMod.^2, 2) ;
SNR = mean(20*log10(realsqrt(NUM)./realsqrt(DEN))) ;
PeakSNR = mean(20 * log10(max(orderMod, [], 2)./realsqrt(MSE))) ;   
end