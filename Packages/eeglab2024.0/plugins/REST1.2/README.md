All copyrights of the software are reserved by the Key Laboratory for NeuroInformation of Ministry of Education, School of Life Science and Technology, University of Electronic Science and Technology of China. This software is for non-commercial use only. It is freeware but not in the public domain.

For more see http://www.neuro.uestc.edu.cn/rest/
Reference: Yao D (2001) A method to standardize a reference of scalp EEG recordings to a point at infinity.
                      Physiol Meas 22:693?11. doi: 10.1088/0967-3334/22/4/305

Written by Li Dong (Li_dong729@163.com) and Shiang Hu (hushiang@126.com)
Date: Aug. 19, 2019
% ---------------------------------------------------------------------------------------------
######################Please cite this toolbox as:###############################

Li Dong*, Fali Li, Qiang Liu, Xin Wen, Yongxiu Lai, Peng Xu and Dezhong Yao*. MATLAB Toolboxes for Reference Electrode Standardization Technique (REST) of Scalp EEG. Frontiers in Neuroscience,  2017:11(601).
% --------------------------------------------------------------------------------------------

Log

REST_reference_v1.0_20170411
      [1] Retain unselected channels (e.g. EMG, EOG, band channels etc.), while saving re-reference EEG data.

REST_reference_v1.1_20190819
      [1] Calculate leadfield at once, based on canonical concentric-three-spheres head model.
      [2] remove the leadfield.exe.
      [3] add reference papers on the interface
% -----------------------------------------------------------
To install REST, download the zip file (linked above), unzip and place the folder in the 'plugins' folder of your existing EEGLAB installation (so something like ~/eeglab12_0_6b/plugins/REST_reference_v1.0/eegplugin_rest.m exists).

To run REST, ensure that the correct EEGLAB folder is in your current Matlab path, and run 'eeglab' as a command from the Matlab Command Window.

Then, load data using EEGLAB, click 'REST'--> 'Re-referencing to REST'

% -----------------------------------------------------------
Re-referencing to REST steps:
        '   [1] Select original reference;';...
        '       Average: average reference.';...
        '       A Fixed Electrode: reference is a fixed a electrode, e.g. Cz.';...
        '   [2] Select Channels: Select channels you want to re-reference;';...
        '       Retain remaining channels?: if you want to keep un-selected channels in the data, check the box;';...
        '   [3] If channel location has been imported in EEG.channlocs, Press button ''Run'' directly;';... 
        '       OR you can selcet a user defined lead field file, and then Press button ''Run'';';...
        '    (a) Load channel location in EEGLAB';...
        '        Edit--> channel locations --> OK';...
        '        Default head model is a 3-concentric sphere head model,more details can be seen in Dong et al., 2017';...
        '        headmodel.cond = [1,0.0125,1]';...
        '        headmodel.tissue = { 'brain','skull', 'scalp' }';...
        '        headmodel.r = [0.87,0.92,1]';...
                 It calculates the leadfield matrix from the 3000 cortical dipoles (spherical equivalent dipoles, see 'corti869-3000dipoles.dat') and the newly given electrode array for the canonical concentric-three-spheres head model.The radii of the three concentri spheres are 0.87(inner radius of the skull), 0.92(outer radius of the skull) and 1.0(radius of the head), while the conductivities are 1.0(brain and scalp) and 0.0125 (skull). The electorde array shoved as *.txt ASCII files with their Cartesian x (the left ear is defined as -x axis),y (the nasion is the +y axis),z coordinates in three columns.
                 Noting that because the default sources/dipoles (concentric-three-spheres head model) are symmetric,if the left ear is defined as +x axis, it will generate same re-referencing results while the left ear is defined as -x axis.

        '    (b) Select a lead field file (has been calculated and saved as *.txt/*.xls/*.xlsx/*.dat):';...
        '        The size of lead field matrix should be No. of sources X No. channels';...
        '  [4] Press button ''OK'' to save the re-referencing data to workspace (ALLEEG).In EEGLAB, click 'Datasets'-->'*_REST'';...
% --------------------------------------------------------------

Version history

v1.2 - fixed by SCCN/UCSD


