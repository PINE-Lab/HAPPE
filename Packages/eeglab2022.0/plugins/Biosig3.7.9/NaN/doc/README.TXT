NaN-Tb: A statistics toolbox 
------------------------------------------------------------
Copyright (C) 2000-2005,2009,2010,2011,2014 Alois Schloegl <alois.schloegl@gmail.com>


FEATURES of the NaN-tb:
-----------------------
 - statistical toolbox
 - machine learning and classification toolbox
 - NaN's are treated as missing values 
 - supports weightening of data
 - usage of multiple CPU cores 

 - supports DIM argument
 - less round-off errors using extended double
 - less but more powerful functions (no nan-FUN needed)
 - supports unbiased estimation
 - fixes known bugs
 - compatible with Matlab and Octave
 - easy to use
 - The toolbox is tested with Octave 3.x and Matlab 7.x


Currently are implemented:
-------------------------- 
level 1: basic functions (not derived)
	SUMSKIPNAN 	SUM is a built-in function and cannot not be replaced,
			For this reason, a different name (than SUM) had to be chosen. 
			SUMSKIPNAN is central, it implements skipping NaN's, the 
			DIM-argument and returns the number of valid elements, too.
	COVM		covariance estimation (several modes)
			Round-off errors avoided by using internally extended accuracy
	CUMSUMSKIPNAN	Cumulative sum, skipping NaN's
	DECOVM		decomposes the extended covarianced matrix into mean and cov 
	XCOVF		cross-correlation function
	FLAG_NANS_OCCURED	returns 0 if no NaN's appeared in the input data 
			of the last call to one of the following functions, and 1 otherwise:  
			sumskipnan, covm, center, cor, coefficient of variation, corrcoef, geomean, harmmean, 
			kurtosis, mad, mean, meandev, meansq, moment, nanmean, nanstd, nansum, 
			rms, sem, skewness, statistic, std, var	
	FLAG_IMPLICIT_SKIP_NAN  can be used to turn off and on the NaN-skipping behaviour. This can 
	                be useful for debugging or for compatibility reasons. 
	FLAG_ACCURACY_LEVEL can be used to increase the accuracy of summations (sumskipnan and covm) 
	                at the cost of speed.
	LOAD_FISHERIRIS	loads famous fisher iris data set 
	STR2ARRAY	convert string to array - useful to extract numeric data from 
			delimiter files
	XPTOPEN		read and write SAS Transport Format (XPT); reads ARFF and STATA files

		The following functions are experimental, not all effects of missing values are fully understood. 
		E.g. Missing values can cause aliasing, also effects on bandpass und highpass filters need to be investigated. 
	NANCONV         convolution
	NANCONV2        2-dimensional convolution
	NANFILTER	filter function
	NANFFT		Fourier transform

level 2a: derived functions
	MEAN	 	mean (options: arithmetic, geometric, harmonic)
	VAR 		variance
	STD 		standard deviation
	MEDIAN		median (currently only for 2-dim matrices)
	SEM		standard error of the mean (does not depend on distribution)
	TRIMMEAN 	trimmed mean 
	medAbsDev	median absolute deviation 

	MEANSQ		mean square
	RMS		root mean square

	STATISTIC	estimates various statistics at once
	MOMENT		moment 	
	SKEWNESS	skewness
	KURTOSIS	excess 

*	IQR		interquartile range
	MAD		mean absolute deviation
*	RANGE		range (max-min)

	CENTER		removes mean
	ZSCORE		normalizes x to zero mean and variance 1 (z = (x-mean)/std) 
	zScoreMedian	non-parametric z-score, normalizes is to zero median and 1/(1.483*median absolute deviation)  

	HARMMEAN	harmonic mean
	GEOMEAN		geometric mean

	NANTEST		checks whether all functions have been replaced
	DETREND		detrending of data with missing values and non-equidistant sampled data

	COR		correlation matrix
	COV		covariance matrix
	CORRCOEF	correlation coefficient, including rank correlation, 
				significance test and confidence intervals 
	SPEARMAN, RANKCORR spearman's rank correlation coefficient. They might be replaced by CORRCOEF. 
	PARTCORRCOEF	partial correlation coefficient			
	RANKS		calculates ranks for non-parametric statistics 
	TIEDRANK	similar to RANKS, used for compatibility reasons 

	QUANTILE	q-th quantile 
	PRCTILE,PERCENTILE	p-th percentile
	TRIMEAN		trimean

	BLAND_ALTMANN	Bland-Altmann plot 
	ECDF		empirical cumulative distribution function 
	CDFPLOT		plot empirical cumulative distribution function 
	GSCATTER	scatter plot of grouped data 
	NORMPDF		normal probability distribution
	NORMCDF		normal cumulative distribution 
	NORMINV		inverse of the normal cumulative distribution
	TPDF		student probability distribution
	TCDF		student cumulative distribution 
	TINV		inverse of the student cumulative distribution
	NANSUM, NANSTD	fixes for buggy versions included
	TTEST		paired t-test
	TTEST2		(unpaired) t-test
	SIGNRANK	Wilcoxon's signed-rank test

level 2b: classification, cross-validation 
	TRAIN_SC	train classifier
	TEST_SC		test classifier
	CLASSIFY	classify data (no cross validation)
	XVAL		classify data with cross validation
	KAPPA		performance evaluation 
	TRAIN_LDA_SPARSE	utility function 
	FSS		feature subset selection and feature ranking
	CAT2BIN         converts categorial to binary data
	SVMTRAIN_MEX	libSVM-training algorithm 
	ROW_COL_DELETION heuristic to select rows and columns to remove missing values


REFERENCE(S):
----------------------------------
[1] http://www.itl.nist.gov/
[2] http://mathworld.wolfram.com/


What is the difference to previous implementations?
===================================================
1) The default behavior of previous implementations is that NaNs in the input 
data results in NaNs in the output data. In many applications this behavior 
is not what you want. In this implementation, NaNs are handled as missing values and 
are skipped. 

2) In previous implementations the workaround was using different functions 
like NANSUM, NANMEAN etc. In this toolbox, the same routines can be applied to
data with and without NaNs. This enables more natural (better read- and 
understandable) applications. 

3) SUMSKIPNAN is central to the other functions. 
It implements 
- the DIMENSION-argument, 
- handles NaNs as missing values or as exception signal (depending on a 
  hidden FLAG), 
- and returns the number of valid elements (which are not NaNs) in the 
  second output argument.
(Note, NANSUM from Matlab does not support the DIM-argument, and NANSUM(NaN)
gives NaN instead of 0);

4) [obsolete] 

5) The DIMENSION argument is implemented in most routines. 
These should work in all Matlab and Octave versions. A workaround for a bug in 
Octave versions <=2.1.35 is implemented. Also several functions from Matlab 
have no support for the DIM argument (e.g. SKEWNESS, KURTOSIS, VAR)

6) Compatible to previous Octave implementation
MEAN implements also the GEOMETRIC and HARMONIC mean. Handling of some special 
cases has been removed because its not necessary, anymore. 
MOMENT implements Mode 'ac' (absolute and/or central) moment as implemented
in Octave. 

7) Performance increase
In most numerical applications, NaN's should be simply skipped. Therefore, 
it is efficient to skip NaN's in the default case. 
In case an explicit check for NaN's is necessary, implicit exception 
handling could be avoided. Eventually the overall performance could increase. 

8) More readable code
An explicit check for NaN's display the importance of this special case. 
Therefore, the application program might be more readable.

9) ZSCORE, MAD, HARMMEAN and GEOMEAN
DIM-argument and skipping of NaN's implemented. None of these features is
implemented in the Matlab versions.

10a) NANMEAN, NANVAR, NANMEDIAN
These are not necessary anymore. They are implemented in SUMSKIPNAN, MEAN, 
VAR, STD and MEDIAN, respectively. 

10b) NANSUM, NANSTD
These functions are obsolete, too. However, previous implementations 
do not always provide the expected result. Therefore, a correct
version is included for backward compatibility. 

11) GPL license
Permits to implement useful modifications. 

12) NORMPDF, NORMCDF, NORMINV
In the Matlab statistics toolbox V 3.0, NORMPDF, NORMCDF and NORMINV gave 
incorrect results for SIGMA=0; A similar problem was observed in Octave 
with NORMAL_INV, NORMAL_PDF, and NORMALCDF. 

The problem is fixed with this version. Furthermore, the check of the input 
arguments is implemented simpler and easier in this versions. 

13) TPDF, TCDF, TINV
In the Matlab statistics toolbox V3.0(12.1) and V4.0(13), TCDF and TINV do not handle NaNs 
correctly. TINV returns 0 instead of NaN, TCDF stops with an error message. 
In Stats-tb V2.2(R11) TINV has also the same problem.

For these reasons, the NaN-tb is a bug fix. Furthermore, the check of the input 
arguments is implemented simpler. Overall, the code becomes cleaner and leaner. 

14) NANCONV, NANCONV2, NANFFT, NANFILTER, NANFILTER1UC
are signal processing functions for graceful handling of data with
missing values. These functions are very experimental, because the behavior in 
case of data with missing values is not fully investigated. 
E.g. missing values can cause aliasing, and also the behavior of bandpass and highpass
filters is not sufficiently investigated. Therefore, these functions should be
used with care. 


Q: WHY SKIPPING NaN's?:
------------------------
A: Usually, NaN means that the value is not available. This meaning is most 
common, even many different reasons might cause NaN's. In statistics, NaN's 
represent missing values, in biosignal processing such missing values might 
have been caused by some recording error. Other reasons for NaN's are, 
undetermined expressions like e.g. 0/0, inf-inf, data not available, unknown value, 
not a numeric value, etc. 

If NaN has the meaning of a missing value, it is only consequent to say, the 
sum of NaN's should be zero. Similar arguments hold for the other functions. 
The mean of X is undefined if and only if X contains no numbers. The 
implementation sum(X)/sum(~isnan(X)) gives 0/0=NaN, which is the desired 
result. The variance of X is undefined if and only if X contains less than 
2 numbers.  

In most numerical applications, NaN's should be simply skipped. Therefore, 
it is efficient to skip NaN's in the default case. In the other cases, the 
NaN's can still be checked explicitly. This could eventually result in a 
more readable code and in improved performance, too.  


Q: What if I need to check for NaN's:
-------------------------------------
A: You can always check whether there were some skipped NaN's in your 
data with the command FLAG_NANS_OCCURED().

m = mean(x);
if flag_nans_occured()
	% do your error handling, e.g. 
	error('there were NaN's in x, ignore m'); 
end; 

Its also easy to control the granularity of the checks

flag_nans_occured(); 	% reset flag
 % do any statistical analysis you want 
if flag_nans_occured()
	% check, whether some NaN's occured. 
end; 


Installing the NaN-tb for Octave and Matlab:
--------------------------------------------
a) Extract files and save them in /your/directory/structure/to/NaN/

b) Include the path with one of the following commands: 
	addpath('/your/directory/structure/to/NaN/')
	path('/your/directory/structure/to/NaN/',path)
   Make sure the functions in the NaN-toolbox are found before the default functions.  
   	
c) run NANINSTTEST 
This checks whether the installation was successful. 

d) Compile mex files:  
  This is useful to improve speed, and is required if you used weighted samples. 
  Check if precompiled binaries are provided. If your platform is not supported, 
  compile the C-Mex-functions using "make". 
	
  Run NANINSTTEST again to check the stability of the compiled SUMSKIPNAN.  

	$Id$
	Copyright (C) 2000-2005,2009,2010,2011,2014 by Alois Schloegl <alois.schloegl@gmail.com>	
        http://pub.ist.ac.at/~schloegl/matlab/NaN/


LICENSE:
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, see <http://www.gnu.org/licenses/>.

