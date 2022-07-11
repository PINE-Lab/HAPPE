# MOD - A Machine-learning optimal detection method. 
These functions show the essential steps of using MOD for
detecting EPSP/EPSCs as described in [1] and applied in [2]. 
 
In order to run this demo, the following software packages are needed: 

+ Octave (or Matlab) and the Signal Processing toolbox or package
+ Biosig [3]
+ NaN-toolbox [4] 
+ SigViewer [5] - for scoring the data and visualization of the resulting detections
	(either v0.5.1, or the most recent version from git repo)

Copyright (C) 2016-2021 Alois Schlögl, IST Austria


# References:
[1] Zhang X, Schlögl A, Vandael D, Jonas P, MOD: A novel machine-learning optimal-filtering method for accurate and efficient detection of subthreshold synaptic events in vivo. Journal of Neuroscience Methods, 2021.  https://doi.org/10.1016/j.jneumeth.2021.109125

[2] Xiaomin Zhang, Alois Schlögl, Peter Jonas
     Selective Routing of Spatial Information Flow from Input to Output in Hippocampal Granule Cells, Neuron, 2020. https://doi.org/10.1016/j.neuron.2020.07.006 http://www.sciencedirect.com/science/article/pii/S0896627320305237 

[3] Biosig - an Free and Ppen source software library for biomedical signal processing
     http://biosig.sourceforge.net/index.html

[4] The NaN-toolbox: A statistics and machine learning toolbox for Octave and Matlab®
     for data with and w/o MISSING VALUES encoded as NaN's. https://pub.ist.ac.at/~schloegl/matlab/NaN/
     
[5] SigViewer (v0.5.1 or v0.65+)
       	https://github.com/cbrnr/sigviewer/commit/f62f8d9d16db8c5c89c0fc12563b20bf4ad2ab7f (Nov 18, 2020)
 	or later. Binaries for windows are available here:
	https://pub.ist.ac.at/~schloegl/software/sigviewer/sigviewer-0.6.4-win64.exe


# Content:

	% demo - applying the method to the example data available from
	%    https://pub.ist.ac.at/~schloegl/software/mod/k	
	% results will be written into output/*
	demo_mod.m

	% this is the core function for obtaining the parameters of 
        % from the MOD method  (filter coefficients, threshold, and delay)
	mod_optimal_detection_filter.m

	% storing the results in a GDF file for visualization with sigviewer
	minidet2gdf.m


# FAQ

## How to apply it to your own data

Check out demo001.m and address these items: 

1) you need your raw data, and your scorings and adapt loading the data (lines 67-87)
Please note, the scorings consists of two periods (one segment in the beginning, the other from the end)
And you need to tell which periods are scored, when constructing the scoring trace
lines 124-136.

2) In case you want to apply blanking of AP's, you need to detect them and
replace them with NaN's (line 88-90, 115-123)

3) if you want to downsample your data, you can check the lines 92-114,

4) If you want to apply cross-validation, you should enable lines 140 - 194

That should do it. The main data processing of training the detection method happens in lines  199-213, the remainder storing the results,
and in case you still have some HF noise in the raw detection trace, you might want to smooth it before doing the event detection.

## Is the output in the text file simply the event timings?

Yes, the event timing can be found in multiple outputs.

- on the workspace its in the variable "pos",

- and its also stored in all three output files, in one way or another.

   output/*.mat    contains various results, including the variable 'pos'
   
   output/*.txt    is the event times in seconds
   
   output/*.gdf    is best viewed with sigviewer, but can be also loaded into matlab



