function [traindata,testdata,metadata]=load_cifar100()
% LOAD_CIFAR100 loads cifar100 data [1,2]. 
%    the data files will be downloaded and uncompressed into 
%    directory $HOME/.cache/ 
% 
% Usage: 
%    [traindata, testdata, meta]=load_cifar100();
% 
% References: 
% [1] Alex Krizhevsky, CIFAR-100 dataset
%     https://www.cs.toronto.edu/~kriz/cifar.html
% [2] https://www.cs.toronto.edu/~kriz/cifar-100-matlab.tar.gz


%       Copyright (C) 2019 Alois Schlögl
%    	This is part of the NaN-toolbox 
%	https://octave.sourceforge.io/nan/index.html
%	https://pub.ist.ac.at/~schloegl/matlab/NaN/
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.



URL="https://www.cs.toronto.edu/~kriz/cifar-100-matlab.tar.gz";
DOWNLOAD_DIRECTORY = fullfile(getenv('HOME'),'.cache');
if ~exist(DOWNLOAD_DIRECTORY,'dir'),
    mkdir(DOWNLOAD_DIRECTORY); 
end;
DOWNLOAD = fullfile(DOWNLOAD_DIRECTORY,'cifar-100-matlab.tar.gz');
DATAFILE = fullfile(DOWNLOAD_DIRECTORY,'cifar-100-matlab','train.mat');

if ~exist(sprintf(DATAFILE,1))
    if ~exist(DOWNLOAD,'file')
        fprintf(1,'Downloading cifar-100 database (~170 MB) to %s/\n', DOWNLOAD_DIRECTORY);
        system(sprintf('wget %s -O %s',URL, DOWNLOAD));
    end
    untar(DOWNLOAD, DOWNLOAD_DIRECTORY)
end;

traindata = load(fullfile(fileparts(DATAFILE),'train.mat'));
testdata  = load(fullfile(fileparts(DATAFILE),'test.mat'));
metadata  = load(fullfile(fileparts(DATAFILE),'meta.mat'));



