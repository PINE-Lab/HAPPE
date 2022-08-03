function [influence] = calcInfluence(dirtySignal, Noise) 
% calcInfluence calculates the unknown influence of the Puls-Noise using
% the LMS algorithm
%
% [influence] = calcInfluence(dirtySignal, Noise)
%
% Input:
%   dirtySignal     ... dirty signal (either [oxy-Hb] or [deoxy-Hb]) with puls influence
%   Noise           ... preprocessed artificial Noise signal (either BP signal or signal from a
%                       fingerpuls sensor
%     
%
% Output:
%   influence       ... calculated unknown influence


% Copyright (C) 2012 by the Institute for Knowledge Discovery, Graz University of Technology 
% Guenther Bauernfeind <g.bauernfeind@tugraz.at>
% WWW: http://bci.tugraz.at/
% $Id: calcInfluence.m v0.1 2012-02-13 10:00:00 ISD$
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
 
M=length(dirtySignal); 
 
% Parameters 
N=50;
y=zeros(M,1); 
w = zeros(N,1); 
u1 = zeros(N,1); 



for i=1:1:M 
    if i<N 
        u1(N:-1:N-i+1)=Noise(i:-1:1); 
    else  
        u1(N:-1:1) = Noise(i:-1:i-N+1); 
    end 
    y(i) = w'*u1; 
 
    %LMS 
    k = 0.05*u1; 
    E = dirtySignal(i) - w'*u1; 
    w = w + k*E; 
   
    
end
 
influence=y;