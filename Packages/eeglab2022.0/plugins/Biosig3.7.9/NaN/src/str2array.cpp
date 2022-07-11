//-------------------------------------------------------------------
//   C-MEX implementation of STR2ARRAY - this function is part of the NaN-toolbox. 
//   Actually, it also fixes a problem in STR2ARRAY.m described here:
//   http://www-old.cae.wisc.edu/pipermail/help-octave/2007-December/007325.html
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
// usage:
//	[...] = STR2ARRAY(s)
//	[...] = STR2ARRAY(sa)
//	[...] = STR2ARRAY(s,cdelim)
//	[...] = STR2ARRAY(s,cdelim,rdelim)
//	[...] = STR2ARRAY(s,cdelim,rdelim,ddelim)
//	[num,status,strarray] = STR2ARRAY(...)
//
// Input:
//  s 	        char string 
//  sa 	        cell array of strings 
//  cdelim	column delimiter
//  rdelim	row delimiter
//  ddelim      decimal delimiter
//
// Output:
//    $Id: STR2ARRAY.cpp 7142 2010-03-30 18:48:06Z schloegl $
//    Copyright (C) 2010,2011 Alois Schloegl <alois.schloegl@gmail.com>
//    This function is part of the NaN-toolbox
//    http://pub.ist.ac.at/~schloegl/matlab/NaN/
//
//-------------------------------------------------------------------



#include <ctype.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include "mex.h"

#ifdef tmwtypes_h
  #if (MX_API_VER<=0x07020000)
    typedef int mwSize;
  #endif 
#endif 


int str2val(char *s, double *r, double *i) 
{
/* 
	str2val converts string into numeric value. real and complex numbers are supported. 
	complex numbers are "3.4 + 5.6i" or "3.4 + i * 5.6" (spaces are optional)

	input:
		s	char string 
	output: 
		*r	real value 
		*i	imaginary value
	return values: 
		0: conversion failed 
		1: real number returned: 
		2: complex number returned 	
	
*/
	char *endptr = NULL;
	double val = strtod(s, &endptr);	// conversion 

	while (isspace(*endptr)) endptr++;
#ifdef DEBUG		
	mexPrintf("123<%s>\t,%f,\t[%s]\n",s,val,endptr);
#endif
	if (!*endptr) {
		// conversion successful 
		*r = val; 
		return(1); 
	}
	else if ((*endptr=='+') || (*endptr=='-')) {
		// imaginary part 
		double sgn = (*endptr=='+') ? 1.0 : -1.0; 
		double ival; 
		while (isspace(*(++endptr))); 

		if (*endptr=='i') {
			// case " a + i * b " 		
			while (isspace(*(++endptr))); 

			if (*endptr=='*') {
				ival = strtod(endptr+1, &endptr);	// conversion 
				if (*endptr && !isspace(*endptr)) {
					return(0); 	// failed 
				}	
				else {
					*r = val; 
					*i = sgn*ival; 
					return(2);	// 
				}	
			}
			else 
				return(0); 	//failed 
		}
		else {
			// case " a + bi " 		
			ival = strtod(endptr, &endptr);	// conversion 
			if (*endptr != 'i') return(0); 
			endptr++;
			while (*endptr) {
				if (!isspace(*endptr)) return(0); 
				endptr++;
			}	
			*r = val; 
			*i = sgn*ival; 
			return(2); 
		}
	}
	else //if (*endptr && !isspace(*endptr)) 
	{
		// conversion failed 
		return(0); 
	}
}



void mexFunction(
    int              nlhs,        /* number of expected outputs */
    mxArray       *plhs[],        /* array of pointers to output arguments */
    int              nrhs,        /* number of inputs */
    const mxArray *prhs[]         /* array of pointers to input arguments */
)

{
	char *s = NULL;
	const char *cdelim = "\x09,";
	const char *rdelim = "\x0a;";
	const char *ddelim = NULL;
	const char *valid_delim = " ()[]{},;:\"|/\x21\x22\x09\0x0a\0x0b\0x0c\0x0d\x00";	// valid delimiter
	uint8_t *u;	
	size_t slen = 0,k;
	size_t maxcol=0, maxrow=0, nr, nc;
	
	if (nrhs<1) {
		mexPrintf("   STR2ARRAY.MEX converts delimiter text files into arrays of numerics and cell-strings\n");
		mexPrintf("   STR2ARRAY.MEX converts delimiter text files into numeric arrays\n");
		mexPrintf("   It fixes a problem of the old STR2DOUBLE discussed here: http://www-old.cae.wisc.edu/pipermail/help-octave/2007-December/007325.html\n");
		mexPrintf("   at avoids using the insecure STR2NUM using EVAL\n");
		mexPrintf("\n   Usage of STR2ARRAY:\n");
		mexPrintf("\t[...] = STR2ARRAY(s)\n");
		mexPrintf("\t[...] = STR2ARRAY(sa)\n");
		mexPrintf("\t[...] = STR2ARRAY(s,cdelim)\n");
		mexPrintf("\t[...] = STR2ARRAY(s,cdelim,rdelim)\n");
		mexPrintf("\t[...] = STR2ARRAY(s,cdelim,rdelim,ddelim)\n");
		mexPrintf("\t[num,status,strarray] = STR2ARRAY(...)\n");
		mexPrintf("   Input:\n\ts\tstring\n\tsa\tcell array of strings\n\tcdelim\tlist of column delimiters (default: \"<Tab>,\"\n\trdelim\tlist of row delimiter (default: \"<LF><CR>;\")");
		mexPrintf("\n\tddelim\tdecimal delimiter (default: \".\"). This is useful if decimal delimiter is a comma (e.g. after Excel export in Europe)\n");
		mexPrintf("   Output:\n\tnum\tnumeric array\n\tstatus\tflag failing conversion\n\tstrarray\tcell array of strings contains strings of failed conversions\n");
		mexPrintf("\nExamples:\n\tSTR2ARRAY('4.12')\n\tSTR2ARRAY('1.2 - 3.4e2i')   complex numbers\n\tSTR2ARRAY('101.01 , 0-i4; 1.2 - i * 3.4, abc')\n\tSTR2ARRAY({'101.01', '0-i4'; '1.2 - i * 3.4', 'abc'})\n\tSTR2ARRAY('1,2;a,b,c;3,4')\n");
		mexPrintf("\tSTR2ARRAY('1;2,3;4',';',',')   exchange row- and column delimiter\n\tSTR2ARRAY('1,200 4;3,400 5',' ',';',',')  replace decimal delimter\n");
		return; 
	}

	/* sanity check of input arguments */

	if ((nrhs==1) && mxIsCell(prhs[0])) {	
		// cell array of strings 
		maxrow = mxGetM(prhs[0]);
		maxcol = mxGetN(prhs[0]);

		/* allocate output memory */
		if (nlhs>2) plhs[2] = mxCreateCellMatrix(maxrow, maxcol);
		uint8_t *v = NULL;
		if (nlhs>1) {
			plhs[1] = mxCreateLogicalMatrix(maxrow, maxcol);
			v = (uint8_t*)mxGetData(plhs[1]);	
			memset(v, 1, maxrow*maxcol);
		}
		plhs[0] = mxCreateDoubleMatrix(maxrow, maxcol, mxREAL);
		double *o = (double*)mxGetData(plhs[0]);
		double *oi= NULL; 
		for (k=0; k<maxrow*maxcol; k++) {
			o[k] = 0.0/0.0;
		}
		for (k = 0; k < maxrow*maxcol; k++) {
			double ival; 
			char *s = mxArrayToString(mxGetCell(prhs[0],k)); 
			if (s==NULL) {
				mxArray *a = mxGetCell(prhs[0],k);
				/*
					this does not work because a must not be modified 
					if (nlhs>2) mxSetCell(plhs[2], k, a);
				*/
			}
			else {
				int typ = str2val(s, o+k, &ival);
				if ((nlhs>2) && (typ==0)) mxSetCell(plhs[2], (mwSize)k, mxCreateString(s));
				if ((nlhs>1) && (typ> 0)) v[k] = 0;
				if (typ==2) {
					if (mxGetPi(plhs[0])==NULL) {
						oi = (double*) mxCalloc(maxrow*maxcol, sizeof(double));
						mxSetPi(plhs[0], oi);
					}
					oi[k] = ival;
				}	
			}	
		}
		// cell-array input is finished 
		return; 
	}

	if (nrhs>0) {
		if (mxIsChar(prhs[0])) {	
			s = mxArrayToString(prhs[0]);
			slen = mxGetNumberOfElements(prhs[0]);
		}
		else 
			mexErrMsgTxt("arg1 is not a char array");
	}		
	if (nrhs>1) {
		if (mxIsChar(prhs[1])) 	
			cdelim = mxArrayToString(prhs[1]);
		else 
			mexErrMsgTxt("arg2 is not a char array");
	}		
	if (nrhs>2) {
		if (mxIsChar(prhs[2])) 	
			rdelim = mxArrayToString(prhs[2]);
		else 
			mexErrMsgTxt("arg3 is not a char array");
	}		
	if (nrhs>3) {
		if (mxIsChar(prhs[3]) && (mxGetNumberOfElements(prhs[3])==1) ) {
			ddelim = mxArrayToString(prhs[3]);
			for (k=0; k<slen; k++) {
				if (s[k]==ddelim[0]) 
					s[k] = '.'; 
			}
		}	
		else 
			mexErrMsgTxt("arg4 is not a single char");
	}		

	/* identify separators */
	u = (uint8_t*) mxCalloc(1,slen+1);
	for (k = 0; k < slen; ) {
		if (strchr(cdelim,s[k]) != NULL) {
			u[k] = 1;      // column delimiter
			while (s[++k]==' ');    // ignore extra space characters
		}	
		else if (strchr(rdelim,s[k]) != NULL)
			u[k++] = 2;    // row delimiter
		else 
			k++; 	       // ordinary character 
	}

	/* count dimensions and set delimiter elements to 0 */
	nc=0, nr=0;
	if (u[slen-1]<2) {
		// when terminating char is not a row delimiter 
		nr = (slen>0);
		u[slen] = 2; 	
	}	
	for (k = 0; k < slen; ) {
		if (u[k]==2) {
			s[k] = 0; 
			nr++;
			if (nc > maxcol) maxcol=nc;
			nc = 0; 
		}
		else if (u[k]==1) {
			s[k] = 0;
			nc++;
		}
		k++; 
	}
	if (nc > maxcol) maxcol=nc;
	maxcol += (slen>0);
	maxrow = nr;

	/* allocate output memory */
	if (nlhs>2) plhs[2] = mxCreateCellMatrix(maxrow, maxcol);
	uint8_t *v = NULL;
	if (nlhs>1) {
		plhs[1] = mxCreateLogicalMatrix(maxrow, maxcol);
		v = (uint8_t*)mxGetData(plhs[1]);	
		memset(v,1,maxrow*maxcol);
	}
	plhs[0] = mxCreateDoubleMatrix(maxrow, maxcol, mxREAL);
	double *o = (double*)mxGetData(plhs[0]);
	double *oi = NULL; 
	for (k=0; k<maxrow*maxcol; k++) {
		o[k] = 0.0/0.0;
	}
	
	nr = 0; nc = 0; 
	size_t last=0; 
	for (k = 0; k <= slen; k++) {
		if (u[k]) {
			// delimiter triggers action 
			size_t idx = nr+nc*maxrow;
			if (last==k) {
				// empty field 
				o[idx] = 0.0/0.0; 
			}
			else {
				double ival;
				int typ = str2val(s+last, o+idx, &ival);

				if ((nlhs>2) && (typ==0)) mxSetCell(plhs[2], idx, mxCreateString(s+last));
				if ((nlhs>1) && (typ> 0)) v[idx] = 0;
				if (typ==2) {
					if (oi==NULL) {
						oi = (double*) mxCalloc(maxrow*maxcol, sizeof(double));
						mxSetPi(plhs[0], oi);
					}
					oi[idx] = ival;
				}	
			}	

			nc++;	// next element 
			if (u[k]==2) {
				nr++;	// next row 
				nc = 0;
			}	
			last = k+1; 
		}
	}	
	mxFree(u); 
};

