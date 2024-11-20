//-------------------------------------------------------------------
//   C-MEX implementation of kth element - this function is part of the NaN-toolbox. 
//
//   usage: x = kth_element(X,k [,flag])
//          returns sort(X)(k)
//
//   References: 
//   [1] https://secure.wikimedia.org/wikipedia/en/wiki/Selection_algorithm
//
//   This program is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 3 of the License, or
//   (at your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program; if not, see <http://www.gnu.org/licenses/>.
//
//
// Input:
//   X    data vector, must be double/real 
//   k    which element should be selected
//   flag [optional]: 
//	 0: data in X might be reorded (partially sorted) in-place and 
//	    is slightly faster because no local copy is generated
//	    data with NaN is not correctly handled. 	
//       1: data in X is never modified in-place, but a local copy is used.    
//	    data with NaN is not correctly handled. 	
//	 2: copies data and excludes all NaN's, the copying might be slower 
//	    than 1, but it enables a faster selection algorithm. 
//	    This is the save but slowest option	
//
// Output:
//      x = sort(X)(k)  
//
//    $Id$
//    Copyright (C) 2010,2011 Alois Schloegl <alois.schloegl@gmail.com>
//    This function is part of the NaN-toolbox
//    http://pub.ist.ac.at/~schloegl/matlab/NaN/
//
//-------------------------------------------------------------------


#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "mex.h"


#ifdef tmwtypes_h
  #if (MX_API_VER<=0x07020000)
    typedef int mwSize;
    typedef int mwIndex;
  #endif 
#endif 

/* 
   math.h has isnan() defined for all sizes of floating point numbers, 
   but c++ assumes isnan(double), causing possible conversions for float and long double
*/
#define ISNAN(a) (a!=a)


#define SWAP(a,b) {temp = a; a=b; b=temp;} 
 
static void findFirstK(double *array, size_t left, size_t right, size_t k)
{
        while (right > left) {
                size_t pivotIndex = (left + right) / 2;

		/* partition */
	        double temp;
	        double pivotValue = array[pivotIndex];
        	SWAP(array[pivotIndex], array[right]);
        	pivotIndex = left;
	        for (size_t i = left; i <= right - 1; ++i ) {
        	        // if (array[i] <= pivotValue || isnan(pivotValue)) // needed if data contains NaN's 
        	        if (array[i] <= pivotValue) 
        	        {
        	        	SWAP(array[i], array[pivotIndex]);
        	                ++pivotIndex;
                	}
        	}
        	SWAP(array[pivotIndex], array[right]);

                if (pivotIndex > k)
                	right = pivotIndex - 1;
                else if (pivotIndex < k)
                        left = pivotIndex + 1;
                else break;        
        }
}
 

void mexFunction(int POutputCount,  mxArray* POutput[], int PInputCount, const mxArray *PInputs[]) 
{
    	size_t k, n;	// running indices 
    	size_t szK, szX;
    	double 	*T,*X,*Y,*K; 
    	char flag = 0; // default value 

	// check for proper number of input and output arguments
	if ( PInputCount < 2 || PInputCount > 3 ) {
		mexPrintf("KTH_ELEMENT returns the K-th smallest element of vector X\n");
		mexPrintf("\nusage:\tx = kth_element(X,k)\n");
		mexPrintf("\nusage:\tx = kth_element(X,k,flag)\n");
		mexPrintf("\nflag=0: the elements in X can be modified in-place, and data with NaN's is not correctly handled. This can be useful for performance reasons, but it might modify data in-place and is not save for data with NaN's. You are warned.\n");
		mexPrintf("flag=1: prevents in-place modification of X using a local copy of the data, but does not handle data with NaN in the correct way.\n");
		mexPrintf("flag=2: prevents in-place modification of X using a local copy of the data and handles NaN's correctly. This is the save but slowest option.\n");
		
	    	mexPrintf("\nsee also: median, quantile\n\n");
	        mexErrMsgTxt("KTH_ELEMENT requires two or three input arguments\n");
	}  
	else if (PInputCount == 3) {
		// check value of flag
		size_t N = mxGetNumberOfElements(PInputs[2]); 
		if (N>1)
		        mexErrMsgTxt("KTH_ELEMENT: flag argument must be scalar\n");
		else if (N==1) {
	    		switch (mxGetClassID(PInputs[2])) {
    			case mxLOGICAL_CLASS:
    			case mxCHAR_CLASS:
	    		case mxINT8_CLASS:
    			case mxUINT8_CLASS:
    				flag = (char)*(uint8_t*)mxGetData(PInputs[2]);
    				break; 
	    		case mxDOUBLE_CLASS:
    				flag = (char)*(double*)mxGetData(PInputs[2]);
    				break; 
	    		case mxSINGLE_CLASS:
    				flag = (char)*(float*)mxGetData(PInputs[2]);
    				break; 
	    		case mxINT16_CLASS:
    			case mxUINT16_CLASS:
    				flag = (char)*(uint16_t*)mxGetData(PInputs[2]);
    				break; 
	    		case mxINT32_CLASS:
    			case mxUINT32_CLASS:
    				flag = (char)*(uint32_t*)mxGetData(PInputs[2]);
    				break; 
	    		case mxINT64_CLASS:
    			case mxUINT64_CLASS:
    				flag = (char)*(uint64_t*)mxGetData(PInputs[2]);
    				break; 
	    		case mxFUNCTION_CLASS:
    			case mxUNKNOWN_CLASS:
    			case mxCELL_CLASS:
	    		case mxSTRUCT_CLASS:
    			default: 
    				mexErrMsgTxt("KTH_ELEMENT: Type of 3rd input argument not supported.");
			}
		}
		// else flag = default value 
	}      
	// else flag = default value 

	if (POutputCount > 2)
	        mexErrMsgTxt("KTH_ELEMENT has only one output arguments.");

	// get 1st argument
	if (mxIsComplex(PInputs[0]) || mxIsComplex(PInputs[1]))
		mexErrMsgTxt("complex argument not supported (yet). ");
	if (!mxIsDouble(PInputs[0]) || !mxIsDouble(PInputs[1]))
		mexErrMsgTxt("input arguments must be of type double . ");
	// TODO: support of complex, and integer data	
		

	szK = mxGetNumberOfElements(PInputs[1]);
	K = (double*)mxGetData(PInputs[1]);

	szX = mxGetNumberOfElements(PInputs[0]);
	X = (double*)mxGetData(PInputs[0]);

	if (flag==0)
		T = X;
	else {
		//***** create temporary copy for avoiding unintended side effects (in-place sort of input data) */ 
		T = (double*)mxMalloc(szX*sizeof(double)); 
		if (flag==1)
			memcpy(T,X,szX*sizeof(double));
		else  {
			/* do not copy NaN's */
			for (k=0,n=0; k < szX; k++) {
				if (!ISNAN(X[k])) T[n++]=X[k]; 
			}
			szX = n; 
		}	
	}
	
	/*********** create output arguments *****************/
	POutput[0] = mxCreateDoubleMatrix(mxGetM(PInputs[1]),mxGetN(PInputs[1]),mxREAL);
	Y = (double*) mxGetData(POutput[0]);
	for (k=0; k < szK; k++) {
		if (K[k] > szX || K[k] < 1)
			Y[k] = 0.0/0.0;	// NaN: result undefined
		else {
        		n = (size_t)(K[k]-1);       // convert to zero-based indexing, round towards 0 
        		findFirstK(T, 0, szX-1, n);
        		Y[k] = T[n];
		}	
	}

	if (flag) mxFree(T);
	
	return; 
}

