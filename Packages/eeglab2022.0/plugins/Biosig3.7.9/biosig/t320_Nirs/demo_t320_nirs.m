% DEMO_t320_nirs 
% Demonstrates Different Approaches for the 
% Removal of Physiological Artefacts  from NIRS Signals

% Copyright (C) 2012 by the Institute for Knowledge Discovery, Graz University of Technology 
% Guenther Bauernfeind <g.bauernfeind@tugraz.at>
% WWW: http://bci.tugraz.at/
% $Id: demo_t320_nirs.m v0.1 2012-02-13 10:00:00 ISD$
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




clear all
close all

% load file

system('wget http://pub.ist.ac.at/~schloegl/download/Demo_Daten_t320.mat');
load('Demo_Daten_t320.mat')


%% Settings
dispFreq=2;                                     %Display frequency for spectra calculation 
fs=Demo.Setting.fs_NIRS;                        %Sampling frequency NIRS data
fs_Physio=Demo.Setting.fs_Physio;               %Sampling frequency Physio data
ExCh=Demo.Setting.ExcludedSurroundingChannels;  %Channels not used for analysis
%% Signals
Oxy_signals=Demo.Data.oxy_signals;              %Oxy-Signals
Deoxy_signals=Demo.Data.deoxy_signals;          %Deoxy-Signals
BP_high=Demo.Data.bp;                           %BP signal
Resp_high=Demo.Data.resp;                       %Respiration Signal
ECG_high=Demo.Data.ecg;                         %ECG Signal
BPsys_high = calcBPlin(BP_high,fs_Physio,1);         %BPsys calculation
BPdia_high = calcBPlin(BP_high,fs_Physio,2);         %BPdia calculation

% Downsampling to 10Hz (Necessary for elimination of the physiolgical
% artifacts)
BP= resample(BP_high,fs,fs_Physio);
BPsys= resample(BPsys_high,fs,fs_Physio);
BPdia= resample(BPdia_high,fs,fs_Physio);
Resp= resample(Resp_high,fs,fs_Physio);
ECG= resample(ECG_high,fs,fs_Physio);

%% Raw [(de)oxy-Hb] data for comparison

% Illustration of [(de)oxy-Hb] raw spectra and
% calculation of averaged raw spectra

spect=[];
count=1;

for ChNr=1:1:size(Oxy_signals,2) 
    figure(1);
    orient landscape
    [rOxy rDeoxy]=Illustration_multichannel_spectra(ChNr, Oxy_signals(:,ChNr), ...
        Deoxy_signals(:,ChNr),fs, dispFreq, ExCh);
    
    if ~isempty(rOxy)
        spect(:,1,count)=rOxy{1}.p;
        spect(:,2,count)=rDeoxy{1}.p;
        count=count+1;
        Base=rOxy{1}.f;
    end
    
end

close

%Average [(de)oxy-Hb] spectra over all used channels

dat_spec_raw=[];

    for k = 1 : size(spect,2)
        dat_spec_raw(:,k)=mean(spect(:,k,:),3); 
    end
    
    
    
%% adaptPulsremove

% In this implementation the BP signal is used to define the puls inluence.
% Beside the BP signal also the signal from a fingerpuls sensor is applicable

for  ChNr=1:1:size(Oxy_signals,2) %beginn der Kanalschleife
        Data.oxy_signals_withoutpuls(:,ChNr)=adaptPulseremove(Oxy_signals(:,ChNr),BP,fs);
        Data.deoxy_signals_withoutpuls(:,ChNr)=adaptPulseremove(Deoxy_signals(:,ChNr),BP,fs);
end

% Illustration of [(de)oxy-Hb] puls cleaned spectra and
% calculation of averaged puls cleaned spectra

spect=[];
count=1;

for ChNr=1:1:size(Oxy_signals,2) 
    figure(2);
    orient landscape
    [rOxy rDeoxy]=Illustration_multichannel_spectra(ChNr, Data.oxy_signals_withoutpuls(:,ChNr), ...
        Data.deoxy_signals_withoutpuls(:,ChNr),fs, dispFreq, ExCh);
    
    if ~isempty(rOxy)
        spect(:,1,count)=rOxy{1}.p;
        spect(:,2,count)=rDeoxy{1}.p;
        count=count+1;
        Base=rOxy{1}.f;
    end
    
end

close

%Average [(de)oxy-Hb] spectra over all used channels

dat_spec_pulscleaned=[];

    for k = 1 : size(spect,2)
        dat_spec_pulscleaned(:,k)=mean(spect(:,k,:),3); 
    end

    
% Comparison raw and cleaned
                    figure(3)
                    % [oxy-Hb] raw
                    p1 = dat_spec_raw(:,1);
                    f1 = Base;
                    Color=[1 0.6 0.6];
                    idx = find(f1<=dispFreq);
                    plot(f1(idx), 10*log10(p1(idx)),'Color', Color,'LineWidth',1.5);
                    hold on
                    % [deoxy-Hb] raw
                    p1 = dat_spec_raw(:,2);
                    Color=[0.8 0.8 1];
                    plot(f1(idx), 10*log10(p1(idx)),'Color', Color,'LineWidth',1.5);
                    % [oxy-Hb] puls clean
                    p1 = dat_spec_pulscleaned(:,1);
                    Color='r';
                    idx = find(f1<=dispFreq);
                    plot(f1(idx), 10*log10(p1(idx)),'Color', Color,'LineWidth',1.5);
                    % [deoxy-Hb] puls clean
                    p1 = dat_spec_pulscleaned(:,2);
                    Color='b';
                    plot(f1(idx), 10*log10(p1(idx)),'Color', Color,'LineWidth',1.5);        
                    ylabel('Power spectrum ','FontSize',16,'FontWeight','bold')
                    xlabel('f (Hz)','FontSize',16,'FontWeight','bold')
                    title('Avg. Spectrum [(de)oxy-Hb] before and after puls removal' ,'FontSize',12,'FontWeight','demi','Interpreter','none')
                    legend('[oxy-Hb]_r_a_w','[deoxy-Hb]_r_a_w','[oxy-Hb]_c_l_e_a_n','[deoxy-Hb]_c_l_e_a_n')
                    
%% Respiration and Mayerwave (BP influence) removed by remNoiseTF

dispFreq=0.8;

%preprocessing Mayer wave influence 
Wnmayer=[0.07 0.13];                    
N=200;
b = fir1(N,Wnmayer/(fs/2),'bandpass');
Mayer=filtfilt(b,1,BPdia);  %Beside the BPdia also the BPsys or the HR signal is applicable

Mayer=Mayer-mean(Mayer);



for  ChNr=1:1:size(Oxy_signals,2) 
        Data.oxy_signals_withoutmayer(:,ChNr)=remNoiseTF(Oxy_signals(:,ChNr),Mayer,fs);
        Data.deoxy_signals_withoutmayer(:,ChNr)=remNoiseTF(Deoxy_signals(:,ChNr),Mayer,fs);
end


%preprocessing respiratory influence
 Wnresp=[0.20 0.40];                     
 N=200;
 b = fir1(N,Wnresp/(fs/2),'bandpass');
 Resp=filtfilt(b,1,Resp); 
 
 Resp=Resp-mean(Resp);
  
for  ChNr=1:1:size(Oxy_signals,2) 
        Data.oxy_signals_withoutmayer_resp(:,ChNr)=remNoiseTF(Data.oxy_signals_withoutmayer(:,ChNr),Resp,fs);
        Data.deoxy_signals_withoutmayer_resp(:,ChNr)=remNoiseTF(Data.deoxy_signals_withoutmayer(:,ChNr),Resp,fs);
end


% Illustration of [(de)oxy-Hb] Mayer and respiration cleaned spectra and
% calculation of averaged cleaned spectra

spect=[];
count=1;

for ChNr=1:1:size(Oxy_signals,2) 
    figure(4);
    orient landscape
    [rOxy rDeoxy]=Illustration_multichannel_spectra(ChNr, Data.oxy_signals_withoutmayer_resp(:,ChNr), ...
        Data.deoxy_signals_withoutmayer_resp(:,ChNr),fs, dispFreq, ExCh);
    
    if ~isempty(rOxy)
        spect(:,1,count)=rOxy{1}.p;
        spect(:,2,count)=rDeoxy{1}.p;
        count=count+1;
        Base=rOxy{1}.f;
    end
    
end

close

%Average [(de)oxy-Hb] spectra over all used channels

dat_spec_mayer_respcleaned=[];

    for k = 1 : size(spect,2)
        dat_spec_mayer_respcleaned(:,k)=mean(spect(:,k,:),3); 
    end

    
% Comparison raw and cleaned
                    figure(5)
                    % [oxy-Hb] raw
                    p1 = dat_spec_raw(:,1);
                    f1 = Base;
                    Color=[1 0.6 0.6];
                    idx = find(f1<=dispFreq);
                    plot(f1(idx), 10*log10(p1(idx)),'Color', Color,'LineWidth',1.5);
                    hold on
                    % [deoxy-Hb] raw
                    p1 = dat_spec_raw(:,2);
                    Color=[0.8 0.8 1];
                    plot(f1(idx), 10*log10(p1(idx)),'Color', Color,'LineWidth',1.5);
                    % [oxy-Hb] puls clean
                    p1 = dat_spec_mayer_respcleaned(:,1);
                    Color='r';
                    idx = find(f1<=dispFreq);
                    plot(f1(idx), 10*log10(p1(idx)),'Color', Color,'LineWidth',1.5);
                    % [deoxy-Hb] puls clean
                    p1 = dat_spec_mayer_respcleaned(:,2);
                    Color='b';
                    plot(f1(idx), 10*log10(p1(idx)),'Color', Color,'LineWidth',1.5);        
                    ylabel('Power spectrum ','FontSize',16,'FontWeight','bold')
                    xlabel('f (Hz)','FontSize',16,'FontWeight','bold')
                    title('Avg. Spectrum [(de)oxy-Hb] before and after BP (Mayer wave) and Respiration removement with TF' ,'FontSize',12,'FontWeight','demi','Interpreter','none')
                    legend('[oxy-Hb]_r_a_w','[deoxy-Hb]_r_a_w','[oxy-Hb]_c_l_e_a_n','[deoxy-Hb]_c_l_e_a_n')
                    
%% Respiration and Mayerwave (BP influence) removed by remNoiseICA

%preprocessing Mayer wave influence 
Wnmayer=[0.07 0.13];                    
N=200;
b = fir1(N,Wnmayer/(fs/2),'bandpass');
Mayer=filtfilt(b,1,BPdia);  %Beside the BPdia also the BPsys or the HR signal is applicable

% Mayer=Mayer-mean(Mayer);

%preprocessing respiratory influence
Wnresp=[0.20 0.40];                     
N=200;
b = fir1(N,Wnresp/(fs/2),'bandpass');
Resp=filtfilt(b,1,Resp); 

%  Resp=Resp-mean(Resp);

Data.oxy_signals_ICA=remNoiseICA(Oxy_signals,Mayer,Resp,fs,ExCh);
Data.deoxy_signals_ICA=remNoiseICA(Deoxy_signals,Mayer,Resp,fs,ExCh);

% Illustration of [(de)oxy-Hb] Mayer and respiration cleaned spectra and
% calculation of averaged cleaned spectra

spect=[];
count=1;

for ChNr=1:1:size(Oxy_signals,2) 
    figure(6);
    orient landscape
    [rOxy rDeoxy]=Illustration_multichannel_spectra(ChNr, Data.oxy_signals_ICA(:,ChNr), ...
        Data.deoxy_signals_ICA(:,ChNr),fs, dispFreq, ExCh);
    
    if ~isempty(rOxy)
        spect(:,1,count)=rOxy{1}.p;
        spect(:,2,count)=rDeoxy{1}.p;
        count=count+1;
        Base=rOxy{1}.f;
    end
    
end

close

%Average [(de)oxy-Hb] spectra over all used channels

dat_spec_ICA=[];

    for k = 1 : size(spect,2)
        dat_spec_ICA(:,k)=mean(spect(:,k,:),3); 
    end

    
% Comparison raw and cleaned
                    figure(7)
                    % [oxy-Hb] raw
                    p1 = dat_spec_raw(:,1);
                    f1 = Base;
                    Color=[1 0.6 0.6];
                    idx = find(f1<=dispFreq);
                    plot(f1(idx), 10*log10(p1(idx)),'Color', Color,'LineWidth',1.5);
                    hold on
                    % [deoxy-Hb] raw
                    p1 = dat_spec_raw(:,2);
                    Color=[0.8 0.8 1];
                    plot(f1(idx), 10*log10(p1(idx)),'Color', Color,'LineWidth',1.5);
                    % [oxy-Hb] puls clean
                    p1 = dat_spec_ICA(:,1);
                    Color='r';
                    idx = find(f1<=dispFreq);
                    plot(f1(idx), 10*log10(p1(idx)),'Color', Color,'LineWidth',1.5);
                    % [deoxy-Hb] puls clean
                    p1 = dat_spec_ICA(:,2);
                    Color='b';
                    plot(f1(idx), 10*log10(p1(idx)),'Color', Color,'LineWidth',1.5);        
                    ylabel('Power spectrum ','FontSize',16,'FontWeight','bold')
                    xlabel('f (Hz)','FontSize',16,'FontWeight','bold')
                    title('Avg. Spectrum [(de)oxy-Hb] before and after BP (Mayer wave) and Respiration removement with ICA' ,'FontSize',12,'FontWeight','demi','Interpreter','none')
                    legend('[oxy-Hb]_r_a_w','[deoxy-Hb]_r_a_w','[oxy-Hb]_c_l_e_a_n','[deoxy-Hb]_c_l_e_a_n')
%% Respiration and Mayerwave (BP influence) removed by remNoiseCAR
                    
[Data.oxy_signals_CAR, Data.deoxy_signals_CAR] = remNoiseCAR(Oxy_signals,Deoxy_signals,ExCh);

% Illustration of [(de)oxy-Hb] Mayer and respiration cleaned spectra and
% calculation of averaged cleaned spectra

spect=[];
count=1;

for ChNr=1:1:size(Oxy_signals,2) 
    figure(8);
    orient landscape
    [rOxy rDeoxy]=Illustration_multichannel_spectra(ChNr, Data.oxy_signals_CAR(:,ChNr), ...
        Data.deoxy_signals_CAR(:,ChNr),fs, dispFreq, ExCh);
    
    if ~isempty(rOxy)
        spect(:,1,count)=rOxy{1}.p;
        spect(:,2,count)=rDeoxy{1}.p;
        count=count+1;
        Base=rOxy{1}.f;
    end
    
end

close

%Average [(de)oxy-Hb] spectra over all used channels

dat_spec_ICA=[];

    for k = 1 : size(spect,2)
        dat_spec_CAR(:,k)=mean(spect(:,k,:),3); 
    end

    
% Comparison raw and cleaned
                    figure(9)
                    % [oxy-Hb] raw
                    p1 = dat_spec_raw(:,1);
                    f1 = Base;
                    Color=[1 0.6 0.6];
                    idx = find(f1<=dispFreq);
                    plot(f1(idx), 10*log10(p1(idx)),'Color', Color,'LineWidth',1.5);
                    hold on
                    % [deoxy-Hb] raw
                    p1 = dat_spec_raw(:,2);
                    Color=[0.8 0.8 1];
                    plot(f1(idx), 10*log10(p1(idx)),'Color', Color,'LineWidth',1.5);
                    % [oxy-Hb] puls clean
                    p1 = dat_spec_CAR(:,1);
                    Color='r';
                    idx = find(f1<=dispFreq);
                    plot(f1(idx), 10*log10(p1(idx)),'Color', Color,'LineWidth',1.5);
                    % [deoxy-Hb] puls clean
                    p1 = dat_spec_CAR(:,2);
                    Color='b';
                    plot(f1(idx), 10*log10(p1(idx)),'Color', Color,'LineWidth',1.5);        
                    ylabel('Power spectrum ','FontSize',16,'FontWeight','bold')
                    xlabel('f (Hz)','FontSize',16,'FontWeight','bold')
                    title('Avg. Spectrum [(de)oxy-Hb] before and after BP (Mayer wave) and Respiration removement with CAR' ,'FontSize',12,'FontWeight','demi','Interpreter','none')
                    legend('[oxy-Hb]_r_a_w','[deoxy-Hb]_r_a_w','[oxy-Hb]_c_l_e_a_n','[deoxy-Hb]_c_l_e_a_n')