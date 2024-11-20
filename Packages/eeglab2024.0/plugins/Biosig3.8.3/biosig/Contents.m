% BIOSIG toolbox contains many useful functions for biomedical signal processing
%     http://biosig.sf.net/
%
% Copyright (C) 2003,2004,2007,2008 by Alois Schloegl <alois.schloegl@gmail.com>
% WWW: http://biosig.sf.net/
% $Id: Contents.m,v 1.8 2008/01/23 22:04:41 schloegl Exp $
% This is part of the BIOSIG-toolbox http://biosig.sf.net/
%
% LICENSE:
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program; if not, write to the Free Software
%     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA 
% 
% 
% 
% 
% === INDEX  BIOSIG ===
% 
% BIOSIG/DOC: Documentation
% ---------------------------------------------
% 	header.txt	specification of header (Dokumentation) 
%	eventcodes.txt	Codes for Events, Markers, Annotations as used in GDF and BioSig
% 	DecimalFactors.txt   
%	units.csv	codes for phyiscal dimensions	
% 
% 
% BIOSIG/demo: [demo's, examples] 
% ---------------------------------------------
% 	demo1	QRS-detection
% 	demo2	estimates and validates BCI classifier
%	demo3	demonstrates how the generate an EDF/GDF/BDF file
%	demo4	Demonstrates how the generate an BKR file
%	demo5	Demonstrates how the generate an WAV file
%       demo6   transfer functions of lumped circuit model
%       demo7   simulations for MVAR estimates        
%       scptest tests loading routine of SCP-ECG data
% 	bench_biosig    benchmark test based on BioSig4OctMat
% 
%
% T100: [Data Acquistion] 
% ---------------------------------------------
% 
% 
% 
% BIOSIG/T200: Data Formats
% ---------------------------------------------
% 	SOPEN		opens biosig data and reads header information
% 	SREAD		reads biosig data
% 	SCLOSE	closes biosig file
% 	SWRITE	writes biosig data (currently only BKR, EDF, BDF implemented)
% 	SSEEK		set file positon indicator
% 	STELL		returns file position indicator
% 	SEOF		checks for end-of-file
%	SREWIND		sets file pointer to the start
% 
% 	SLOAD 	Opens, reads all data and closes file Biosig files. 
% 	TLOAD 	triggered loading of data
%  	and some utility functions
% 
% 
% BIOSIG/T250: Quality Control and Artifact Processing
% ---------------------------------------------
%	ARTIFACT_SELECTION	converts artifact scorings into trial selections
% 	EEG2HIST	calculates histogram
% 	GETTRIGGER	gets trigger points
% 	TRIGG		extract fixed-length trials around trigger points	
%	DETECT_MUSCLE	detection of muscle artefacts using an inverse filter
%       REGRESS_EOG     reduce EOG artifacts with regression analysis 
%       REMOVE5060HZ    methods for removing line interference
% 
% 
% BIOSIG/T300: Signal Processing and Feature extraction
% ---------------------------------------------
%    general: 
% 	processing	general framework for blockwise-dataprocessing
%    EEG: 
%       BARLOW          Barlow parameters 
%       HJORTH          Hjorth parameters 
% 	WACKERMANN      Wackermann parameters
%       EVOKED_POTENTIAL        EP estimation
%	TFMVAR		time-frequency multivariate autoregressive modelling
% 	LUMPED          Lumped Circuit model for the EEG alpha rhythm
%    EMG:
%       PAYNTER         Paynter filter for Amplitude demodulation of EMG
%    ECG: 
%       QRSDETECT       QRS-Detection methods 
%       BERGER          resampling of HRV data
%       HEARTRATEVARIABILITY    estimates various HRV parameters
%	QRScorr		correctiong of QRS-detection
%	ECTBcorr	correction of Ectopic beat effect
%       TVAAR           Time-varying AAR estimation 
%
%
% BIOSIG/T310: Wavelet and Fourier Analysis,
% ---------------------------------------------
% 	ERD/ERS
%	PLV
%	HRV
%
%
% BIOSIG/T400: Classification, Single Trial Analysis
% ---------------------------------------------
% 	DECOVM		decomposes an "extended" covariance matrix	
% 	TRAIN_SC	estimates classifier from labelled data
% 	TEST_SC         applies classifier to test data 
% 	FINDCLASSIFIER	obtains classifier includeding performance test. 
%	XVAL		cross-validation procedure
%
%
% BIOSIG/T450: Statistical analysis, False Discovery Rate (FDR)
% -------------------------------------------------------------
%	FDR	false discovery rate
%	several utility functions
%
%
% BIOSIG/T490: Evaluation criteria
% ---------------------------------------------
% 	AUC		area under the curve
% 	ROC		receiver-operator-characteristics 	
% 	MUTINFO		mutual information
% 	QCMAHAL		quality check of multiple discriminator
% 	KAPPA		kappa statistics
% 	BCI3EVAL	Evaluation of BCI results (triggered output)
% 	BCI4EVAL	Evaluation of BCI results (continous output)
%	CRITERIA2005IIIb   calculates maximum steepness of mutual information
% 	
% 
% BIOSIG/T500: Presentation, Output
% ---------------------------------------------
% 	PLOTA		general plot functions for various data structures
%	SVIEW		simple signal viewer 
%       ELPOS           2-D electrode positions
%       ELPOS3          3-D electrode positions
% 
% 
% BIOSIG/T510: Visualization II
% ---------------------------------------------
%      various functions to display results of T310	
%
% BIOSIG/T501: Visualization of EEG Cooupling 
% ---------------------------------------------
%	PLOT_COUPLING	displays EEG coupling
% 
% 
% BIOSIG/T600: Interactive Viewer and Scoring  
% ---------------------------------------------
%	SVIEWER		Interactive Viewer and Scoring
% 	VIEWEDF		EDF-Viewer 
%
% 
% 
% TSA: Time Series Analysis 
% ---------------------------------------------
%
%
% NaN: Statistical Toolbox for data with missing Samples 
% ---------------------------------------------
%
%
%
%
%
%


