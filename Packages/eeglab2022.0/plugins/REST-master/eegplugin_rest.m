% All copyrights of the software are reserved by the Key Laboratory for 
% NeuroInformation of Ministry of Education, School of Life Science and 
% Technology, University of Electronic Science and Technology of China. 
% This software is for non-commercial use only. It is freeware but not 
% in the public domain.
% 
% For more see http://www.neuro.uestc.edu.cn/rest/
% Reference: Yao D (2001) A method to standardize a reference of scalp EEG 
%            recordings to a point at infinity.
% Physiol Meas 22:693?11. doi: 10.1088/0967-3334/22/4/305
% 
% Written by Li Dong (Li_dong729@163.com)
% Date: Oct. 12, 2016
% update by Li Dong Aug. 18 2020

function currvers = eegplugin_rest(fig, trystrs, catchstrs)

currvers  = ['REST_v1.2']; % version

if nargin < 3
   error('eegplugin_rest requires 3 arguments');
end

%
% ADD FOLDER TO PATH
%
p = which('eegplugin_rest','-all');
if length(p)>1
        fprintf('\nREST WARNING: More than one REST folder was found.\n\n');
end
p = p{1};
p = p(1:findstr(p,'eegplugin_rest.m')-1);
% add all ERPLAB subfolders
addpath(genpath(p))

%
% REST NEST-MENU  (REST at the EEGLAB's Main Menu)
%
if ispc      % windows
        wfactor1 = 1.20;
        wfactor2 = 1.21;
elseif ismac % Mac OSX
        wfactor1 = 1.45;
        wfactor2 = 1.46;
else
        wfactor1 = 1.30;
        wfactor2 = 1.31;
end
% posmainfig = get(gcf,'Position');
% hframe     = findobj('parent', gcf,'tag','Frame1');
% posframe   = get(hframe,'position');
% set(gcf,'position', [posmainfig(1:2) posmainfig(3)*wfactor1 posmainfig(4)]);
% set(hframe,'position', [posframe(1:2) posframe(3)*wfactor2 posframe(4)]);
% 
% menuREST = findobj(fig,'tag','EEGLAB');   % At EEGLAB Main Menu


%% ************************************************************************
%  **********************|        MENU      |******************************
%  **********************|      CALLBACKS   |******************************
%  ************************************************************************

%
% REST callback
%
comRESTman   = [trystrs.no_check,'pop_REST_reref;',catchstrs.add_to_hist] ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%        MAIN      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%        MENU      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create REST menu
%
% submenu = uimenu(menuREST,'Label','REST','separator','on','tag','REST');
% set(submenu,'position', 4); % menu position
% 
% restMenu = uimenu( submenu,'Label',['Re-referencing to REST'],...
%     'CallBack' , comRESTman,...
%     'separator', 'off');

menu = findobj(fig, 'tag', 'tools');

submenu = uimenu(menu,'Label','Re-referencing to REST','CallBack' , comRESTman,'separator', 'on');
