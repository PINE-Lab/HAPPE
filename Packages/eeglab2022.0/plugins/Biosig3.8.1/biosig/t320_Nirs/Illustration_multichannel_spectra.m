function [rOxy rDeoxy]=Illustration_multichannel_spectra(ChNr, oxy_Data, deoxy_Data, fs, dispFreq, ExCh)
% Illustration_multichannel_spectra 
% uses the calcNIRSspectra function to calculate and illustrates the [(de)oxy-Hb]spectra 
% of each used NIRS channel from a 3*11 measurement grid
%
% [rOxy rDeoxy]=Illustration_multichannel_spectra(ChNr, oxy_Data, deoxy_Data, fs, dispFreq, ExCh)
%
% Input:
%   ChNr          ... Number of the current channel
%   oxy_Data      ... Oxy data of the current channel
%   deoxy_Data    ... Deoxy data of the current channel
%   fs            ... Sampling frequency
%   dispFreq      ... Frequencies up to this value are plotted
%   ExCh          ... Channels not used
%
% Output:
%   rOxy          ... Structure containing the [oxy-Hb] spectrum
%   rDeoxy        ... Calculated [deoxy-Hb] spectrum

% Copyright (C) 2012 by the Institute for Knowledge Discovery, Graz University of Technology 
% Guenther Bauernfeind <g.bauernfeind@tugraz.at>
% WWW: http://bci.tugraz.at/
% $Id: Illustration_multichannel_spectra.m v0.1 2012-02-13 10:00:00 ISD$
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


   % check if current channel is used
        
        for i=1:1:size(ExCh{1,1},2)
            if ChNr==ExCh{1,1}(i)
                rOxy =[];
                rDeoxy=[];
                return
            end
        end
        
       

        
  %% Illustration definition for 3*11 grid
                    %first row
                    if ChNr==1
                        %set(gca,'Visible','off')
                        g1=axes('position',[0.05300 0.76667 0.07576 0.125]);axes(g1);
                    end
                    if ChNr==2
                        g2=axes('position',[0.14394 0.76667 0.07576 0.125]);axes(g2);end
                    if ChNr==3
                        g3=axes('position',[0.23485 0.76667 0.07576 0.125]);axes(g3);end
                    if ChNr==4
                        g4=axes('position',[0.32576 0.76667 0.07576 0.125]);axes(g4);end
                    if ChNr==5
                        g5=axes('position',[0.41667 0.76667 0.07576 0.125]);axes(g5);end
                    if ChNr==6
                        g6=axes('position',[0.50758 0.76667 0.07576 0.125]);axes(g6);end
                    if ChNr==7
                        g7=axes('position',[0.59848 0.76667 0.07576 0.125]);axes(g7);end
                    if ChNr==8
                        g8=axes('position',[0.68939 0.76667 0.07576 0.125]);axes(g8);end
                    if ChNr==9
                        g9=axes('position',[0.78030 0.76667 0.07576 0.125]);axes(g9);end
                    if ChNr==10
                        g10=axes('position',[0.87121 0.76667 0.07576 0.125]);axes(g10);end
                    
                    %second row
                    if ChNr==11
                        g11=axes('position',[0.00758 0.60000 0.07576 0.125]);axes(g11);end
                    if ChNr==12
                        g12=axes('position',[0.09848 0.60000 0.07576 0.125]);axes(g12);end
                    if ChNr==13
                        g13=axes('position',[0.18939 0.60000 0.07576 0.125]);axes(g13);end
                    if ChNr==14
                        g14=axes('position',[0.28030 0.60000 0.07576 0.125]);axes(g14);end
                    if ChNr==15
                        g15=axes('position',[0.37121 0.60000 0.07576 0.125]);axes(g15);end
                    if ChNr==16
                        g16=axes('position',[0.46212 0.60000 0.07576 0.125]);axes(g16);end
                    if ChNr==17
                        g17=axes('position',[0.55303 0.60000 0.07576 0.125]);axes(g17);end
                    if ChNr==18
                        g18=axes('position',[0.64394 0.60000 0.07576 0.125]);axes(g18);end
                    if ChNr==19
                        g19=axes('position',[0.73485 0.60000 0.07576 0.125]);axes(g19);end
                    if ChNr==20
                        g20=axes('position',[0.82575 0.60000 0.07576 0.125]);axes(g20);end
                    if ChNr==21
                        g21=axes('position',[0.91667 0.60000 0.07576 0.125]);axes(g21);end
                    
                    %third row
                    if ChNr==22
                        g22=axes('position',[0.05300 0.43334 0.07576 0.125]);axes(g22);end
                    if ChNr==23
                        g23=axes('position',[0.14394 0.43334 0.07576 0.125]);axes(g23);end
                    if ChNr==24
                        g24=axes('position',[0.23485 0.43334 0.07576 0.125]);axes(g24);end
                    if ChNr==25
                        g25=axes('position',[0.32576 0.43334 0.07576 0.125]);axes(g25);end
                    if ChNr==26
                        g26=axes('position',[0.41667 0.43334 0.07576 0.125]);axes(g26);end
                    if ChNr==27
                        g27=axes('position',[0.50758 0.43334 0.07576 0.125]);axes(g27);end
                    if ChNr==28
                        g28=axes('position',[0.59848 0.43334 0.07576 0.125]);axes(g28);end
                    if ChNr==29
                        g29=axes('position',[0.68939 0.43334 0.07576 0.125]);axes(g29);end
                    if ChNr==30
                        g30=axes('position',[0.78030 0.43334 0.07576 0.125]);axes(g30);end
                    if ChNr==31
                        g31=axes('position',[0.87121 0.43334 0.07576 0.125]);axes(g31);end
                    
                    %fourth row
                    if ChNr==32
                        g32=axes('position',[0.00758 0.26667 0.07576 0.125]);axes(g32);end
                    if ChNr==33
                        g33=axes('position',[0.09848 0.26667 0.07576 0.125]);axes(g33);end
                    if ChNr==34
                        g34=axes('position',[0.18939 0.26667 0.07576 0.125]);axes(g34);end
                    if ChNr==35
                        g35=axes('position',[0.28030 0.26667 0.07576 0.125]);axes(g35);end
                    if ChNr==36
                        g36=axes('position',[0.37121 0.26667 0.07576 0.125]);axes(g36);end
                    if ChNr==37
                        g37=axes('position',[0.46212 0.26667 0.07576 0.125]);axes(g37);end
                    if ChNr==38
                        g38=axes('position',[0.55303 0.26667 0.07576 0.125]);axes(g38);end
                    if ChNr==39
                        g39=axes('position',[0.64394 0.26667 0.07576 0.125]);axes(g39);end
                    if ChNr==40
                        g40=axes('position',[0.73485 0.26667 0.07576 0.125]);axes(g40);end
                    if ChNr==41
                        g41=axes('position',[0.82575 0.26667 0.07576 0.125]);axes(g41);end
                    if ChNr==42
                        g42=axes('position',[0.91667 0.26667 0.07576 0.125]);axes(g42);end
                    
                    %fifth row
                    if ChNr==43
                        g43=axes('position',[0.05300 0.10000 0.07576 0.125]);axes(g43);end
                    if ChNr==44
                        g44=axes('position',[0.14394 0.10000 0.07576 0.125]);axes(g44);end
                    if ChNr==45
                        g45=axes('position',[0.23485 0.10000 0.07576 0.125]);axes(g45);end
                    if ChNr==46
                        g46=axes('position',[0.32576 0.10000 0.07576 0.125]);axes(g46);end
                    if ChNr==47
                        g47=axes('position',[0.41667 0.10000 0.07576 0.125]);axes(g47);end
                    if ChNr==48
                        g48=axes('position',[0.50758 0.10000 0.07576 0.125]);axes(g48);end
                    if ChNr==49
                        g49=axes('position',[0.59848 0.10000 0.07576 0.125]);axes(g49);end
                    if ChNr==50
                        g50=axes('position',[0.68939 0.10000 0.07576 0.125]);axes(g50);end
                    if ChNr==51
                        g51=axes('position',[0.78030 0.10000 0.07576 0.125]);axes(g51);end
                    if ChNr==52
                        g52=axes('position',[0.87121 0.10000 0.07576 0.125]);axes(g52);end
        
        
%%   Calculate and display the [(de)oxy-Hb]spectra
                    

        %calculate [(de)oxy-Hb] spectra 
        
        rOxy=calcNIRSspectra(oxy_Data,fs);
        rDeoxy=calcNIRSspectra(deoxy_Data,fs);
         
        %display [oxy-Hb] spectra
        p1 = rOxy{1}.p;
        f1 = rOxy{1}.f;
        Color='r';
        idx = find(f1<=dispFreq);
        plot(f1(idx), 10*log10(p1(idx)),'Color', Color,'LineWidth',1.5);
        hold on
        
        %display [deoxy-Hb] spectra
        p1 = rDeoxy{1}.p;
        f1 = rDeoxy{1}.f;
        Color='b';
        idx = find(f1<=dispFreq);
        plot(f1(idx), 10*log10(p1(idx)),'Color', Color,'LineWidth',1.5);
        xlabel('f [Hz]','FontSize',6);
        



        

        
        
        
        
        
                 
                
 end
        