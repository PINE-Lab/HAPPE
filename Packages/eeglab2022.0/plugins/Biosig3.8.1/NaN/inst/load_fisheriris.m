% LOAD_FISHERIRIS 
%        loads famous iris data set from Fisher, 1936 [1]. 
%
% References: 
% [1] Fisher,R.A. "The use of multiple measurements in taxonomic problems" 
%        Annual Eugenics, 7, Part II, 179-188 (1936); also in "Contributions to Mathematical Statistics" (John Wiley, NY, 1950).
% [2] Duda,R.O., & Hart,P.E. (1973) Pattern Classification and Scene Analysis. 
%        (Q327.D83) John Wiley & Sons. ISBN 0-471-22361-1. See page 218.

% Copyright (C) 2009,2010,2016,2019,2020 by Alois Schloegl <alois.schloegl@gmail.com>
%   This function is part of the NaN-toolbox
%   http://pub.ist.ac.at/~schloegl/matlab/NaN/

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
% of the  License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street - Fifth Floor, Boston, MA 02110-1301, USA.


if exist('OCTAVE_VERSION','builtin')
	IRIS=fullfile(tempdir(),'iris.dat');

	if ~exist(IRIS,'file')
		if exist('webread','file')     % available in Octave 6 or 7
			species = webread('http://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data');
			fid=fopen(IRIS,'w'); fwrite(fid,species,'char'); fclose(fid);
		elseif strncmp(computer,'PCWIN',5)     % on MXE-OCTAVE for Windows, wget is available
			fprintf(1,'Download http://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data and save in under %s\nPress any key to continue ...\n',IRIS);
		else 
		        system(['wget http://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data -O ',IRIS]);
		end;         
	end;      

	if exist('str2array','file')==3,
		if ~exist('species','var')
		        tmp = fopen(IRIS); species=fread(tmp,[1,inf],'uint8=>char'); fclose(tmp);
		end
        	[meas, tmp, species] = str2array(species,',');        
	        meas = meas(1:150, 1:4);
	        species = species(1:150, 5);

	elseif exist('textread', 'file'),
		[a,b,c,d,species] = textread (IRIS, '%f,%f,%f,%f,%s\n');
		meas = [a,b,c,d];
		clear a b c d;
		if (size(meas,1)==151)
			% remove empty line at the end
			meas(151,:)=[];
			species(151)=[];
		end
	else 
		tmp  = fopen(IRIS); species=fread(tmp,[1,inf],'uint8=>char'); fclose(tmp);
		[meas,tmp,species]=str2double(species,',');
	        meas = meas(:,1:4);
	        species = species(:,5);        
	end
else 
        load fisheriris; 
end; 

%!xtest
%! load_fisheriris
%! assert(all(size(meas)==[150,4]))
%! assert(all(size(species)==[150,1]))

