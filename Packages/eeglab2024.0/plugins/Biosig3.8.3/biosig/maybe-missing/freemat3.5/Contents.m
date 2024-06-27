% BIOSIG/MAYBE-MISSING/FREEMAT3.5  contains functions that 
% are missing in FreeMat v3.5. These improve the compatibility 
% of BioSig with FreeMat

%	$Id: Contents.m,v 1.2 2008-01-23 11:37:29 schloegl Exp $
%	Copyright (C) 2008 by Alois Schloegl <alois.schloegl@gmail.com>
%	This is part of the BIOSIG project http://biosig.sf.net/


see the list of files in this directory for missing files. 

Other known incompatibilities are: 
- FSEEK: no return argument
- FREAD: fails if the number of avaiable bytes is insufficient.
- NDIMS: trailing singleton dimensions are not ignored.



