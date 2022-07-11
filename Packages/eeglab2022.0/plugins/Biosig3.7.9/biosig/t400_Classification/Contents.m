% BIOSIG/T400 contains classifiers.  
%
%
% TRAIN_SC      train classifier
% TEST_SC       test  classifier
% CLASSIFY	to classify samples into categories
% XVAL		cross-validation procedure
% UNTRAIN_SC	decremental training of classifier (depre
% FC0   generates classifier from "asychronous" data
% FINDCLASSIFIER generates BCI classifier [1-3]
%
% obsolete: 
% FINDCLASSIFIER1 generate classifier for a BCI [1,2]
% FINDCLASSIFIER2 2nd generation classifier for BCI [2,3]
%
%  helper function
%
% PERM	permutations of indices in trials 
% TRAIN_LDA_SPARSE sparse LDA classifier
%
% Several core functions (like train_sc, test_sc, xval) have been included in the NaN-toolbox
% and will primarily be maintained there. For more details see: 
%   http://pub.ist.ac.at/~schloegl/matlab/NaN 
%
% REFERENCES: 
% [1] Schloegl A, Neuper C, Pfurtscheller G
% 	Estimating the mutual information of an EEG-based Brain-Computer-Interface
%  	Biomedizinische Technik 47(1-2): 3-8, 2002.
% [2] Schloegl A, Keinrath C, Scherer R, Pfurtscheller G,
%	Information transfer of an EEG-based Bran-computer interface.
%	Proceedings of the 1st International IEEE EMBS Conference on Neural Engineering, Capri, Italy, Mar 20-22, 2003 
% [3] Schloegl A, Lee FY, Bischof H, Pfurtscheller G
%	Characterization of Four-Class Motor Imagery EEG Data for the BCI-Competition 2005.
%	Journal of neural engineering 2 (2005) 4, S. L14-L22
% [4] Schloegl A, Kronegg J, Huggins JE, Mason SG;
%	Evaluation criteria in BCI research.
%	(Eds.) G. Dornhege, J.R. Millan, T. Hinterberger, D.J. McFarland, K.-R.Mueller;
%	Towards Brain-Computer Interfacing, MIT press, 2007, p. 327-342.
% [5] Schloegl A, Brunner C, Scherer R, Glatz A
%       BioSig - an open source software library for BCI research.
%       (Eds.) G. Dornhege, J.R. Millan, T. Hinterberger, D.J. McFarland, K.-R. Mueller;
%       Towards Brain-Computer Interfacing, MIT Press, 2007, p.347-358.
% [6] Schloegl A., Brunner C.
%	BioSig: A Free and Open Source Software Library for BCI Research,
%	Computer, vol. 41, no. 10, pp. 44-50, October, 2008.
%
%	$Id$
%	Copyright (c) 1997-2005,2006,2008,2009 by Alois Schloegl
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

