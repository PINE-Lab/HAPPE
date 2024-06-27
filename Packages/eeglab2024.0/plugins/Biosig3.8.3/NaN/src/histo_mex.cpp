//-------------------------------------------------------------------
//   C-MEX implementation of Histogram - this function is part of the NaN-toolbox. 
//
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
// histo_mex: computes histogram 
//
// Input:
// - data matrix 
// - flag for row-wise histogram
//
// Output:
// - histogram 
//     HIS.X
//     HIS.H
//
//    $Id$
//    Copyright (C) 2009,2010,2011 Alois Schloegl <alois.schloegl@gmail.com>
//    This function is part of the NaN-toolbox
//    http://pub.ist.ac.at/~schloegl/matlab/NaN/
//
//-------------------------------------------------------------------

/* TODO: 
	speed: its slower than the m-functions histo2/3/4
		|-> use a more efficient sorting function 
	resembling of histo3 for multicolumn data. 
	support of complex data and char-strings 
*/

#include <math.h>
#include <stdint.h>
#include <string.h>
#include "mex.h"

/* 
   math.h has isnan() defined for all sizes of floating point numbers, 
   but c++ assumes isnan(double), causing possible conversions for float and long double
*/
#define ISNAN(a) (a!=a)


#ifdef tmwtypes_h
  #if (MX_API_VER<=0x07020000)
    typedef int mwSize;
  #endif 
#endif 

struct sort_t {
	uint8_t *Table;	// data table
	size_t Size;	// sizeof elements e.g. 4 for single
	size_t Stride; 	// for multicolumn data 
	size_t N; 	// number of rows
	mxClassID Type;	// data type  
} Sort; 	

//inline int compare(const sqize_t *a, const size_t *b) {
int compare(const void *a, const void *b) {
	int z = 0; 
	size_t i = 0;
	size_t ix1 = *(size_t*)a;
	size_t ix2 = *(size_t*)b;
	
	while ((i<Sort.N) && !z) {
		switch (Sort.Type) {
		case mxCHAR_CLASS:
			z = memcmp(Sort.Table+ix1*Sort.Size,Sort.Table+ix2*Sort.Size,Sort.Size); 
			break;
		case mxINT32_CLASS: {
			int32_t f1,f2;
			f1 = ((int32_t*)Sort.Table)[ix1];
			f2 = ((int32_t*)Sort.Table)[ix2];
			if (f1<f2) z = -1; 
			else if (f1>f2) z = 1; 
			break;
			}
		case mxUINT32_CLASS: {
			uint32_t f1,f2;
			f1 = ((uint32_t*)Sort.Table)[ix1];
			f2 = ((uint32_t*)Sort.Table)[ix2];
			if (f1<f2) z = -1; 
			else if (f1>f2) z = 1; 
			break;
			}
		case mxINT64_CLASS: {
			int64_t f1,f2;
			f1 = ((int64_t*)Sort.Table)[ix1];
			f2 = ((int64_t*)Sort.Table)[ix2];
			if (f1<f2) z = -1; 
			else if (f1>f2) z = 1; 
			break;
			}
		case mxUINT64_CLASS: {
			uint64_t f1,f2; 
			f1 = ((uint64_t*)Sort.Table)[ix1];
			f2 = ((uint64_t*)Sort.Table)[ix2];
			if (f1<f2) z = -1; 
			else if (f1>f2) z = 1; 
			break;
			}
		case mxSINGLE_CLASS: {
			float f1,f2;
			f1 = ((float*)Sort.Table)[ix1];
			f2 = ((float*)Sort.Table)[ix2];
			z = ISNAN(f1) - ISNAN(f2);
			if (z) break;
			
			if (f1<f2) z = -1; 
			else if (f1>f2) z = 1; 
			// else f1==f2 || (isnan(f1) && isnan(f2))	
			break;
			}
		case mxDOUBLE_CLASS: {
			double f1,f2;
			f1 = ((double*)Sort.Table)[ix1];
			f2 = ((double*)Sort.Table)[ix2];
			z = ISNAN(f1) - ISNAN(f2);
			if (z) break;
			
			if (f1<f2) z = -1; 
			else if (f1>f2) z = 1; 
			// else f1==f2 || (isnan(f1) && isnan(f2))	
			break;
			}
		case mxINT16_CLASS: {
			int16_t f1,f2;
			f1 = ((int16_t*)Sort.Table)[ix1];
			f2 = ((int16_t*)Sort.Table)[ix2];
			if (f1<f2) z = -1; 
			else if (f1>f2) z = 1; 
			break;
			}
		case mxUINT16_CLASS: {
			uint16_t f1,f2;
			f1 = ((uint16_t*)Sort.Table)[ix1];
			f2 = ((uint16_t*)Sort.Table)[ix2];
			if (f1<f2) z = -1; 
			else if (f1>f2) z = 1; 
			break;
			}
		case mxINT8_CLASS: {
			int8_t f1,f2;
			f1 = ((int8_t*)Sort.Table)[ix1];
			f2 = ((int8_t*)Sort.Table)[ix2];
			if (f1<f2) z = -1; 
			else if (f1>f2) z = 1; 
			break;
			}
		case mxUINT8_CLASS: {
			uint8_t f1,f2;
			f1 = ((uint8_t*)Sort.Table)[ix1];
			f2 = ((uint8_t*)Sort.Table)[ix2];
			if (f1<f2) z = -1; 
			else if (f1>f2) z = 1; 
			break;
			}
		default:
			mexErrMsgTxt("unsupported input type");	
		}
		i++;	
		ix1 += Sort.Stride;
		ix2 += Sort.Stride;
	}	
	return(z);
}


void mexFunction(int POutputCount,  mxArray* POutput[], int PInputCount, const mxArray *PInputs[]) 
{

    	const mwSize	*SZ;	    
	char 		flag_rows = 0; 
	char 		done = 0; 
    	size_t    	j, k, l;	// running indices 
	const mxArray	*W = NULL; 
    	double 		*w = NULL; 

	// check for proper number of input and output arguments
	if ((PInputCount <= 0) || (PInputCount > 3)) {
		mexPrintf("HISTO_MEX computes histogram from vector or column matrices\n\n");
		mexPrintf("usage:\tHIS = histo_mex(Y)\n\t\tComputes histogram from each column\n");
		mexPrintf("\t[HIS,tix] = histo_mex(Y,'rows')\n\t\tComputes row-wise histogram, tix is useful for data compression.\n\t\t Y = HIS.X(tix,:); \n\n");
		
	    	mexPrintf("see also: HISTO2, HISTO3, HISTO4\n\n");
	        mexErrMsgTxt("HISTO_MEX requires 1 or 2 input arguments\n");
	}        
	if (POutputCount > 2)
	        mexErrMsgTxt("histo.MEX has 1 output arguments.");

	// get 1st argument
	if (mxIsComplex(PInputs[0]))
		mexErrMsgTxt("complex argument not supported (yet). ");
	// TODO: support complex argument!	
		
	if (PInputCount==1) 
		; // histo_mex(X)
	else if (mxIsChar(PInputs[1])) {
		// histo_mex(X,'rows')
		char *t = mxArrayToString(PInputs[1]);
		flag_rows = !strcmp(t,"rows");
		mxFree(t); 
		// histo_mex(X,'rows',W)
		if ((PInputCount>2) && mxIsDouble(PInputs[2])) 	W = PInputs[2]; 
	} 
		// histo_mex(X,W)
	else if (mxIsDouble(PInputs[1])) {	
		W = PInputs[1]; 
	}
	else 
		mexErrMsgTxt("Weight vector must be REAL/DOUBLE.");

	if (W != NULL) 	{
	 	if (mxGetM(PInputs[0])==mxGetM(W) ) 
	 		w = (double*)mxGetData(W);
	 	else 
			mexErrMsgTxt("number of rows in X and W do not match.");

		for (k=0; (k<mxGetM(W)) && (w[k]>=0.0); k++); 
		if (k<mxGetM(W)) 
			mexWarnMsgTxt("weight vector contains also non-negative numbers or NaN.");
	} 			
		
	// get size
    	mwSize ND = mxGetNumberOfDimensions(PInputs[0]);	
    	// NN = mxGetNumberOfElements(PInputs[0]);
    	SZ = mxGetDimensions(PInputs[0]);		

	if (ND>2) 
		mexErrMsgTxt("Error HISTO.MEX: input must be vector or matrix (no more than two dimensions)");

	size_t n  = SZ[0];
	
	const char *fnames[] = {"datatype","X","H"};
	mxArray	*HIS = mxCreateStructMatrix(1, 1, 3, fnames);
	mxSetField(HIS,0,"datatype",mxCreateString("HISTOGRAM"));
	
	if (flag_rows || (SZ[1]==1)) { 

		///***** SORT each column: initialize sorting algorithm  
		size_t *idx = NULL;
		idx = (size_t*) mxMalloc(SZ[0]*sizeof(size_t));
		for (n=0; n<SZ[0]; n++) {
			idx[n]=n;			
		}
		Sort.Type = mxGetClassID(PInputs[0]); 
		Sort.Table = (uint8_t*) mxGetData(PInputs[0]);

		switch (mxGetClassID(PInputs[0])) {
		case mxLOGICAL_CLASS:
		        Sort.Size = sizeof(mxLogical);
		        break; 
		case mxINT8_CLASS: 
		case mxUINT8_CLASS: 
			Sort.Size = 1;
			break; 
		case mxINT16_CLASS: 
		case mxUINT16_CLASS: 
			Sort.Size = 2;
			break; 
		case mxINT32_CLASS: 
		case mxUINT32_CLASS: 
		case mxSINGLE_CLASS:
			Sort.Size = 4;
			break; 
		case mxINT64_CLASS: 
		case mxUINT64_CLASS: 
		case mxDOUBLE_CLASS: 
			Sort.Size = 8;
			break;
		case mxCHAR_CLASS: 
		        Sort.Size = sizeof(mxChar);
			break;
		default:
			mexErrMsgTxt("unsupported input type"); 
		}	
		Sort.N = flag_rows ? SZ[1] : 1; 
		qsort(idx,SZ[0],sizeof(*idx),compare);

		// count number of different elements 	
		n = SZ[0] ? 1 : 0; 
		size_t k;
		for (k=1; k<SZ[0]; k++) {
			if (compare(idx+k-1,idx+k)) n++;
		}

		uint64_t *tix = NULL; 
		if (POutputCount>1) {
			POutput[1] = mxCreateNumericMatrix(SZ[0], 1, mxUINT64_CLASS,mxREAL);
			tix = (uint64_t*)mxGetData(POutput[1]);
		}	

		// fill HIS.H and HIS.X 
		mxArray *H = mxCreateNumericMatrix(n, 1, mxDOUBLE_CLASS,mxREAL);
		mxArray *X = mxCreateNumericMatrix(n, SZ[1], mxGetClassID(PInputs[0]),mxREAL); 
		mxSetField(HIS,0,"H",H);
		mxSetField(HIS,0,"X",X);
		double  *h = (double*)mxGetData(H);
		uint8_t *x = (uint8_t*)mxGetData(X);
		
		l = 0;
		if (tix) tix[idx[0]] = 1;
		for (k=0; k<SZ[0]; k++) {
			if ((k==0) || compare(&idx[k-1], &idx[k])) {
				for (j=0; j<SZ[1]; j++) {
					memcpy(x + (l+j*n)*Sort.Size, Sort.Table+(idx[k] + j*Sort.Stride)*Sort.Size, Sort.Size);
				}
				l++;
			}
			if (tix) tix[idx[k]] = l;
			if (w != NULL) 
				h[l-1]+=w[idx[k]];
			else 
				h[l-1]+=1.0;	
		}
		mxFree(idx); 
		done = 1;
	}

	else { 
		switch (mxGetClassID(PInputs[0])) {

		case mxINT8_CLASS: { 	
			mxArray *H = mxCreateNumericMatrix(256, SZ[1], mxDOUBLE_CLASS,mxREAL);
			mxArray *X = mxCreateNumericMatrix(256, 1, mxINT8_CLASS,mxREAL); 
			mxSetField(HIS,0,"H",H);
			mxSetField(HIS,0,"X",X);

			int8_t *x;
			x = (int8_t*)mxGetData(X);
			for (k=0; k<0x0100; k++)
				x[k]=(int8_t)(k-128);

			x = (int8_t*)mxGetData(PInputs[0]);
			double *h = (double*)mxGetData(H);
			for (k=0; k<SZ[0]*SZ[1]; k++)
				h[x[k]+128+(k/SZ[0]<<8)] += (w!=NULL ? w[k%SZ[0]] : 1.0);
		
			done = 1; 	
			break;
			}

		case mxUINT8_CLASS: { 	
			mxArray *H = mxCreateNumericMatrix(256, SZ[1], mxDOUBLE_CLASS,mxREAL);
			mxArray *X = mxCreateNumericMatrix(256, 1, mxUINT8_CLASS,mxREAL); 
			mxSetField(HIS,0,"H",H);
			mxSetField(HIS,0,"X",X);

			uint8_t *x = (uint8_t*)mxGetData(X);
			for (k=0; k<0x0100; k++) x[k]=(uint8_t)k;

			x = (uint8_t*)mxGetData(PInputs[0]);
			double *h = (double*)mxGetData(H);
			for (k=0; k<SZ[0]*SZ[1]; k++)
				h[x[k]+(k/SZ[0]<<8)] += (w!=NULL ? w[k%SZ[0]] : 1.0);
				
			done = 1; 	
			break;
			}

		case mxINT16_CLASS: {	
			mxArray *H = mxCreateNumericMatrix(0x10000, SZ[1], mxDOUBLE_CLASS,mxREAL);
			mxArray *X = mxCreateNumericMatrix(0x10000, 1, mxINT16_CLASS,mxREAL); 
			mxSetField(HIS,0,"H",H);
			mxSetField(HIS,0,"X",X);

			double *h = (double*)mxGetData(H);
			int16_t *x = (int16_t*)mxGetData(X);
			for (k=0; k<0x10000; k++)
				x[k]=(int16_t)(k-0x8000);

			x = (int16_t*)mxGetData(PInputs[0]);
			for (k=0; k<SZ[0]*SZ[1]; k++)
				h[x[k]+0x8000+(k/SZ[0]<<16)] += (w!=NULL ? w[k%SZ[0]] : 1.0);
			
			done = 1; 	
			break;
			}

		case mxUINT16_CLASS: {	
			mxArray *H = mxCreateNumericMatrix(0x10000, SZ[1], mxDOUBLE_CLASS,mxREAL);
			mxArray *X = mxCreateNumericMatrix(0x10000, 1, mxINT16_CLASS,mxREAL); 
			mxSetField(HIS,0,"H",H);
			mxSetField(HIS,0,"X",X);
	
			double *h = (double*)mxGetData(H);
			int16_t *x = (int16_t*)mxGetData(X);
			for (k=0; k<0x10000; k++)
				x[k]=(uint16_t)(k-0x8000);

			uint16_t *x16 = (uint16_t*)mxGetData(PInputs[0]);
			for (k=0; k<SZ[0]*SZ[1]; k++)
				h[x16[k]+(k/SZ[0]<<16)] += (w!=NULL ? w[k%SZ[0]] : 1.0);
			done = 1; 	
			break;	
			}

		default: {
			/* FIXME */
			mexErrMsgTxt("multicolumns with int32 or larger not supported!"); 

			mxArray *H = mxCreateNumericMatrix(0x10000, SZ[1], mxDOUBLE_CLASS,mxREAL);
			mxArray *X = mxCreateNumericMatrix(0x10000, SZ[1], mxGetClassID(PInputs[0]),mxREAL); 
			mxSetField(HIS,0,"H",H);
			mxSetField(HIS,0,"X",X);

			/*
			double *h = (double*)mxGetData(H);
			int16_t *x = (int16_t*)mxGetData(X);

			for (n=0; n<SZ[1]; n++) {
			}	
			*/
			}
		} // end switch 	
	}
	

	/*******  output    *******/
	if (done) POutput[0] = HIS; 
	
	return; 
}
