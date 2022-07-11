function k=convert2single_sweep_atf(filename)
% Demostration for convert multisweep data into single sweep ATF files
% 

%    Copyright (C) 2016 by Alois Schloegl <alois.schloegl@ist.ac.at>	
%    This is part of the BIOSIG-toolbox http://biosig.sf.net/
%
%    BioSig is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    BioSig is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with BioSig.  If not, see <http://www.gnu.org/licenses/>.


% load file 
[s,HDR]=sload(filename);
if ~isfield(HDR,'PhysDim')
	HDR.PhysDim=physicalunits(HDR.PhysDimCode);
end; 

H=sprintf('ATF\t1.0\r\n0\t%i\r\n"Time (ms)"',HDR.NS+1);
FMTSTR='\r\n%f';
for n=1:HDR.NS
	H=[H,sprintf('\t"%s (%s)"',HDR.Label{n},HDR.PhysDim{n})];
	FMTSTR=[FMTSTR,'\t%f'];
end; 

ix = [1;HDR.EVENT.POS(HDR.EVENT.TYP==hex2dec('7ffe'));size(s,1)+1];
for k=1:length(ix)-1,
	F = fullfile(HDR.FILE.Path,[HDR.FILE.Name,sprintf('_%04i.atf',k)]);
	fid=fopen(F,'wt');
	fprintf(fid,'%s',H);	

	t = [0:ix(k+1)-1-ix(k)]/HDR.SampleRate; 
	fprintf(fid,FMTSTR,[t;s(ix(k):ix(k+1)-1,:)']);

	fprintf(fid,'\r\n');
	fclose(fid); 
end; 


