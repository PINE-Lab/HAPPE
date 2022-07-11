function [cleanOxysignals, cleanDeoxysignals] = remNoiseCAR(oxy_signals,deoxy_signals,ExCh)
%
% remNoiseCAR removes respiration an blood pressure related noise from
% [(de)oxy-Hb] signals by using CAR. 
%
%   [cleanSignals] = remNoiseCAR(Signals,fs,ExCh)
%
% Input:
%   oxy_signals      ... Oxy-Signals with noise
%   deoxy_signals    ... Deoxy-Signals with noise
%   ExCh             ... Channels not used
%
%
% Output:
%   cleanOxysignals    ... corrected Oxy-signal without influence of noise
%   cleanDeoxysignals  ... corrected Deoxy-signal without influence of noise



% Copyright (C) 2012 by the Institute for Knowledge Discovery, Graz University of Technology 
% Guenther Bauernfeind <g.bauernfeind@tugraz.at>
% WWW: http://bci.tugraz.at/
% $Id: remNoiseCAR.m v0.1 2012-02-13 10:00:00 ISD$
%
% This software was implemented in the framework of the Styrian government project "Einfluss von
% Herz-Kreislauf-Parametern auf das Nah-Infrarot-Spektroskopie (NIRS) Signal" (A3-22.N-13/2009-8) 
% and is part of the BIOSIG-toolbox http://biosig.sf.net/
%
% LICENSE:
% 
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>. 


    
        CAR_ch_f=ExCh{1,1};
        
        CAR_ch=[1:1:size(oxy_signals,2)];
        CAR_ch(CAR_ch_f)=[];

        mean_oxy_signal=mean(oxy_signals(:,CAR_ch),2);
        mean_deoxy_signal=mean(deoxy_signals(:,CAR_ch),2);
        
        for tt=1:size(oxy_signals,2)
        cleanOxysignals(:,tt)=oxy_signals(:,tt)-mean_oxy_signal(:,1);
        cleanDeoxysignals(:,tt)=deoxy_signals(:,tt)-mean_deoxy_signal(:,1);
        end
        
   

end