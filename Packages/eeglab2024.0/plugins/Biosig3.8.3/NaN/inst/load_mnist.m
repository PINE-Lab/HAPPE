function [train_data, train_labels, test_data, test_labels] = load_mnist(f)
% LOAD_MNIST load MNIST database [1] 
% 
% Usage:  
%     [train_data, train_labels, test_data, test_labels] = load_mnist();
%     
% 
% References: 
% [1] Yann LeCun, Corinna Cortes, Christopher J.C. Burges, 
%     THE MNIST DATABASE of handwritten digits 
%     http://yann.lecun.com/exdb/mnist/


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Download and MNIST data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
URL   = 'http://yann.lecun.com/exdb/mnist/';
files = {'train-images-idx3-ubyte.gz', 't10k-images-idx3-ubyte.gz', 'train-labels-idx1-ubyte.gz', 't10k-labels-idx1-ubyte.gz'};
DOWNLOAD_DIRECTORY = fullfile(getenv('HOME'),'.cache/mnist');
if ~exist(DOWNLOAD_DIRECTORY,'dir'), 
      mkdir(DOWNLOAD_DIRECTORY); 
end;
for k = 1:length(files)
        DOWNLOAD = fullfile(DOWNLOAD_DIRECTORY,files{k});
        if ~exist(DOWNLOAD,'file')
                system(sprintf('wget %s/%s -O %s',URL, files{k},DOWNLOAD));
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load all files 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ((nargin < 1) || (nargout>1)),
        train_data   = load_mnist(fullfile(DOWNLOAD_DIRECTORY,files{1}));
        train_labels = load_mnist(fullfile(DOWNLOAD_DIRECTORY,files{2}));
        test_data    = load_mnist(fullfile(DOWNLOAD_DIRECTORY,files{3}));
        test_labels  = load_mnist(fullfile(DOWNLOAD_DIRECTORY,files{4}));
        return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open and read content of file(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid   = fopen(f, 'rz', 'ieee-be');
if fid<0, error('can not open file'); end

magic = fread(fid, 1, 'int32');
N     = fread(fid, 1, 'int32');
if magic==2051,
    sz = fread(fid, [1,2], 'int32');
    pixel = reshape(fread(fid, [prod(sz),N], 'uint8=>uint8')',[N,sz]);
elseif magic==2049,
    pixel = fread(fid, [N,1], 'uint8=>uint8');
else 
    error('unknown file type');
end
fclose(fid);
train_data = pixel;
return

