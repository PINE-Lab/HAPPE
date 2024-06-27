% *** This function is part of ERPLAB Toolbox ***
% Author: Javier Lopez-Calderon
% Center for Mind and Brain
% University of California, Davis,
% Davis, CA
% 2014

%b8d3721ed219e65100184c6b95db209bb8d3721ed219e65100184c6b95db209b
%
% ERPLAB Toolbox
% Copyright � 2007 The Regents of the University of California
% Created by Javier Lopez-Calderon and Steven Luck
% Center for Mind and Brain, University of California, Davis,
% javlopez@ucdavis.edu, sjluck@ucdavis.edu
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [ERP, erpcom] = pop_getFFTfromERP(ERP, varargin)
erpcom = '';
if nargin < 1
      help pop_getFFTfromERP
      return
end
if isfield(ERP(1), 'datatype')
      datatype = ERP.datatype;
else
      datatype = 'ERP';
end
if nargin==1
      title_msg  = 'ERPLAB: pop_getFFTfromERP() error:';
      if isempty(ERP)
            ERP = preloadERP;
            if isempty(ERP)
                  msgboxText =  'No ERPset was found!';
                  
                  errorfound(msgboxText, title_msg);
                  return
            end
      end
      if isempty(ERP.bindata)
            msgboxText = 'cannot work with an empty ERP erpset';
            errorfound(msgboxText, title_msg);
            return
      end
      if ~strcmpi(datatype, 'ERP')
            msgboxText =  'This ERPset is already converted into frequency domain!';
            errorfound(msgboxText, title_msg);
            return
      end     
      
      %
      % FFT points will be as much as needed
      % to get 1 point each 0.25 Hz, at least.
      % Users can change this value using scripting. Jav
      %
      fnyqx = round(ERP.srate/2);
      np    = 2.^nextpow2(4*fnyqx); % NFFT
      
      erpcom = pop_getFFTfromERP(ERP, 'NFFT', np, 'Saveas', 'on','History','gui');
      return
end

%
% Parsing inputs
%
p = inputParser;
p.FunctionName  = mfilename;
p.CaseSensitive = false;
p.addRequired('ERP');
% option(s)
p.addParamValue('TaperWindow', 'on');       % 'ERP': compute the average of epochs per bin;
p.addParamValue('NFFT', []);                % number of points for the FFT
p.addParamValue('Saveas', 'off', @ischar);  % 'on', 'off'
p.addParamValue('Warning', 'off', @ischar);
p.addParamValue('History', 'script', @ischar); % history from scripting

p.parse(ERP, varargin{:});

if iseegstruct(ERP)
      if length(ERP)>1
            msgboxText =  'ERPLAB says: Unfortunately, this function does not work with multiple ERPsets';
            error(msgboxText);
      end
end


% NFFT, iswindowed

if strcmpi(p.Results.TaperWindow,'off')
      iswindowed = 0;
elseif strcmpi(p.Results.TaperWindow,'on')
      iswindowed = 1;
else
      if ~isempty(p.Results.TaperWindow) && ischar(p.Results.TaperWindow)
            iswindowed = p.Results.TaperWindow;
      else
            error('Unknow value for "TaperWindow"')
      end
end

NFFT = p.Results.NFFT;
if ismember_bc2({p.Results.Saveas}, {'on','yes'})
      issaveas  = 1;
else
      issaveas  = 0;
end

% if strcmpi(p.Results.Warning, 'on')
%         rwwarn = 1;
% else
%         rwwarn = 0;
% end
if strcmpi(p.Results.History,'implicit')
      shist = 3; % implicit
elseif strcmpi(p.Results.History,'script')
      shist = 2; % script
elseif strcmpi(p.Results.History,'gui')
      shist = 1; % gui
else
      shist = 0; % off
end

%
% subroutine
%
ERP = getFFTfromERP(ERP, NFFT, iswindowed);

%
% History
%
skipfields = {'ERP', 'Saveas','History'};
fn     = fieldnames(p.Results);
erpcom = sprintf('%s = pop_getFFTfromERP( %s ', inputname(1), inputname(1));

for q=1:length(fn)
      fn2com = fn{q};
      if ~ismember_bc2(fn2com, skipfields)
            fn2res = p.Results.(fn2com);
            if ~isempty(fn2res)
                  if ischar(fn2res)
                        if ~strcmpi(fn2res,'off')
                              erpcom = sprintf( '%s, ''%s'', ''%s''', erpcom, fn2com, fn2res);
                        end
                  else
                        if iscell(fn2res)
                              if ischar([fn2res{:}])
                                    fn2resstr = sprintf('''%s'' ', fn2res{:});
                              else
                                    fn2resstr = vect2colon(cell2mat(fn2res), 'Sort','on');
                              end
                              fnformat = '{%s}';
                        else
                              fn2resstr = vect2colon(fn2res, 'Sort','on');
                              fnformat = '%s';
                        end
                        if strcmpi(fn2com,'Criterion')
                              if p.Results.Criterion<100
                                    erpcom = sprintf( ['%s, ''%s'', ' fnformat], erpcom, fn2com, fn2resstr);
                              end
                        else
                              erpcom = sprintf( ['%s, ''%s'', ' fnformat], erpcom, fn2com, fn2resstr);
                        end
                  end
            end
      end
end
erpcom = sprintf( '%s );', erpcom);

%
% Save ERPset
%
if issaveas
      [ERP, issave, erpcom_save] = pop_savemyerp(ERP,'gui','erplab', 'History', 'implicit');
      if issave>0
            if issave==2
                  erpcom = sprintf('%s\n%s', erpcom, erpcom_save);
                  msgwrng = '*** Your ERPset was saved on your hard drive.***';
            else
                  msgwrng = '*** Warning: Your ERPset was only saved on the workspace.***';
            end
      else
            msgwrng = 'ERPLAB Warning: Your changes were not saved';
      end
      try cprintf([1 0.52 0.2], '%s\n\n', msgwrng); catch,fprintf('%s\n\n', msgwrng);end ;
end
% get history from script. ERP
switch shist
      case 1 % from GUI
            displayEquiComERP(erpcom);
      case 2 % from script
            ERP = erphistory(ERP, [], erpcom, 1);
      case 3
            % implicit
      otherwise %off or none
            erpcom = '';
            return
end

%
% Completion statement
%
msg2end
return