% BIOSIG/T330_StimFit contains Matlab/Octave implementation of 
%  StimFit [1] and related functions. 
%
% The following functions are included:
%    SIMUL001		generate simulated data (file: test01.gdf)
%    SIMUL002		generate simulated data (file: test02.gdf)
%    MICROSTIMFIT  	stimfit event analysis 
%    MINIDET 		miniature EPSP detection based on template matching		
%    DEMO
%
% REFERENCES: 
% [1] Jose Guzman, Alois Schlögl, Christoph Schmidt-Hieber
%     Stimfit: quantifying electrophysiological data with Python.
%     Front. Neuroinform. 8:16, 2014
%     available online: doi: http://dx.doi.org/10.3389/fninf.2014.00016
%     https://pub.ist.ac.at/~schloegl/publications/GuzmanEtAl2014.fninf-08-00016.pdf
% [2] http://pub.ist.ac.at/~schloegl/biosig/
% [3] http://biosig.sf.net/
% [4] https://github.com/neurodroid/stimfit
% [5] Jonas P, Major G, Sakmann B. Quantal components of unitary EPSCs at the mossy
%     fibre synapse on CA3 pyramidal cells of rat hippocampus. J Physiol. 1993
%     Dec;472:615-63.


%    Copyright (C) 2013, 2019 by Alois Schloegl <alois.schloegl@ist.ac.at>
%    This is part of the BIOSIG project http://biosig.sf.net/
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


