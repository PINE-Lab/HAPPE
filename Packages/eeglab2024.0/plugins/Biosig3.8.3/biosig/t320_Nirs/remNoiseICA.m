function [cleanSignals] = remNoiseICA( Signals,bpNoise,respNoise,fs,ExCh)
%
% remNoiseICA removes respiration an blood pressure related noise from
% [(de)oxy-Hb] signals by using ICA. The [(de)oxy-Hb] signal is decomposed 
% into independent components (ICs) via SOBI ICA [1]. The coherence between 
% each IC and the noise signals is then calculated. ICs for which the coherence 
% with one of the artefact signals is higher than the mean of all the coherence 
% scores with that artefact signal plus 1 standard deviation are flagged for 
% removal. 
%
% The code uses the runica function from:
% Makeig, Scott et al. "EEGLAB: ICA Toolbox for Psychophysiological Research". 
% WWW Site, Swartz Center for Computational Neuroscience, Institute of Neural
% Computation, University of San Diego California
% <www.sccn.ucsd.edu/eeglab/>, 2000.
%
% Please be sure to have included
%
%   [cleanSignals]=remNoiseICA( Signals,bpNoise,respNoise,fs,ExCh)
%
% Input:
%   Signals    ... Signals with noise (either [oxy-Hb] or [deoxy-Hb])
%   bpNoise    ... Noise signal BP
%   respNoise  ... Noise signal respiration
%   fs         ... Sampling frequency
%   ExCh       ... Channels not used
%
%
% Output:
%   cleanSignals  ... corrected signal without influence of noise
%
%[1] Belouchrani A, Abed-Meraim K, Cardoso J, Moulines E. A blind source
%    seperation technique using second-order statistics. IEEE Transactions 
%    on signal processing, 45(2): 434-444, 1997


% Copyright (C) 2012 by the Institute for Knowledge Discovery, Graz University of Technology 
% Ian Daly <ian.daly@tugraz.at> and Guenther Bauernfeind <g.bauernfeind@tugraz.at>
% WWW: http://bci.tugraz.at/
% $Id: remNoiseICA.m v0.1 2012-02-13 10:00:00 ISD$
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

    artSig(1,:) = bpNoise;
    artSig(2,:) = respNoise;
    
    chsUse=[1:1:size(Signals,2)];
    chsUse(ExCh{1,1})=[];
        
    dirtySignals = Signals(:,chsUse);

    
    [ICweights ICsphere,compvars,bias,signs,lrates,S] = runica( dirtySignals','maxsteps',40 );
    
    % 2. Get activations.
    S = ICweights * ICsphere * dirtySignals(1:size(artSig,2),:)';

    allICsArt = [S; artSig];
    % Normalise signals.
    for i = 1:size( allICsArt,1 ),
        allICsArt(i,:) = allICsArt(i,:) - mean( allICsArt(i,:) );
        allICsArt(i,:) = allICsArt(i,:) ./ std( allICsArt(i,:) );
    end

    
    %artIndex  = [zeros(1,size(S,1)) ones(1,size(artSig,1))];
    
    % Get coherence.
    for k = 1:size(allICsArt,1)-3,
        [a b] =  mscohere(allICsArt(end,:),allICsArt(k,:),[],[],fs);%hanning(fs),fs/2,fs);
            
        cohVal(1,k)=mean( a );
    end
        
    for k = 1:size(allICsArt,1)-3,
        [a b] =  mscohere(allICsArt(end-1,:),allICsArt(k,:),[],[],fs);%hanning(fs),fs/2,fs);
            
        cohVal(2,k)=mean( a );
    end
    for k = 1:size(allICsArt,1)-3,
        [a b] =  mscohere(allICsArt(end-2,:),allICsArt(k,:),[],[],fs);%hanning(fs),fs/2,fs);
            
        cohVal(3,k)=mean( a );
    end
    
    indI = mean( cohVal(1,:) ) + std( cohVal(1,:) );
    [removeICs1] = find( cohVal(1,:) > indI );
    indI = mean( cohVal(2,:) ) + (std(cohVal(2,:)).*0.5);
    [removeICs2] = find( cohVal(2,:) > indI );
    removeICs = unique( [removeICs1 removeICs2] );
    removeICs = removeICs( find(removeICs <= size(ICweights,1) ) );

    %noRemovedICs = length( removeICs );
    
    % RunICA.
    ICweights2 = ICweights;
    ICsphere2 = ICsphere;
    S(removeICs,:) = 0;
    
    % 5. Translate from ICs back to original data.
    data = (ICweights2*ICsphere2)^-1 * S;
    
    cleanSignals = zeros(size(Signals));

    cleanSignals(1:size(data,2),chsUse) = data';
    

end