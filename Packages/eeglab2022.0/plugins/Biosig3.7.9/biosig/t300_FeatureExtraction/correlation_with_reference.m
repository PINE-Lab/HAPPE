function R = correlation_with_reference(fn)
% CORRELATION_WITH_REFERENCE estimates the the correlation of each 
%    channel with the (common, global) activity at the references
%    electrode.
%
%  R = CORRELATION_WITH_REFERENCE(filename)
%  R = CORRELATION_WITH_REFERENCE(data)
%

%	$Id$
%	Copyright (C) 2005 by Alois Schloegl <alois.schloegl@gmail.com>
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% Version 2 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Library General Public License for more details.
%
% You should have received a copy of the GNU Library General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc., 59 Temple Place - Suite 330,
% Boston, MA  02111-1307, USA.


if isnumeric(fn)
        Y = fn;
        HDR.NS = size(S,2);
else
        [Y,HDR]=sload(fn);
        R.Label = HDR.Label; 
        if isfield(HDR,'ELPOS')
                R.ELPOS = HDR.ELPOS; 
        end;
        if isfield(HDR,'ELEC')
                R.ELPOS = HDR.ELEC.XYZ; 
        end;
end;

Y(any(isnan(Y),2),:)=NaN;
m = mean(Y);
for k1 = 1:HDR.NS,
        Y(:,k1) = Y(:,k1)-m(k1);
end;
w = toeplitz([0,ones(1,HDR.NS-1)/(HDR.NS-1)]);
v = eye(HDR.NS)-w;

R.ssk=meansq(Y*v); % sk
R.ss0=meansq(Y*w); % s0
R.sk0=mean((Y*v).*(Y*w));

R.corr = -R.sk0./sqrt(R.ssk.*R.ss0);
R.datatype = 'CORRELATION_WITH_REFERENCE';

return; 

if 0, % obsolote
        % calculation of angle alpha
        R.sy=meansq(Y);
        a = acos((R.ss0+R.ssk-R.sy)./(sqrt(R.ss0.*R.ssk)*2))*180/pi;
        a(imag(a)~=0)=NaN;
        R.alpha = 90-real(a);
end;



