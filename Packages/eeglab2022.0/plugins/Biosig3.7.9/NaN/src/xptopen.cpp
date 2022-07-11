//-------------------------------------------------------------------
//   XPTOPEN is C-MEX implementation for reading various  
//   statistical data formats including SAS/XPT, SPSS/PASW, 
//   STATA and ARFF data formats. Basic support for writing 
//   SAS/XPT is also supported. 
//   Endian conversion is done automatically.
//
//   usage: x = xptopen(filename)
//   usage: x = xptopen(filename,'r')
//		read filename and return variables in struct x
//   usage: xptopen(filename,'w',x)
//		save fields of struct x in filename
//   usage: x = xptopen(filename,'a',x)
//		append fields of struct x to filename
//
//   References:
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
//    $Id$
//    Copyright (C) 2010,2011,2012,2013 Alois Schloegl <alois.schloegl@ist.ac.at>
//    This function is part of the NaN-toolbox
//    http://pub.ist.ac.at/~schloegl/matlab/NaN/
//
// References:
// [1]	TS-140 THE RECORD LAYOUT OF A DATA SET IN SAS TRANSPORT (XPORT) FORMAT
//	http://support.sas.com/techsup/technote/ts140.html
// [2] IBM floating point format
//	http://en.wikipedia.org/wiki/IBM_Floating_Point_Architecture
// [3] see http://old.nabble.com/Re%3A-IBM-integer-and-double-formats-p20428979.html
// [4] STATA File Format
//	http://www.stata.com/help.cgi?dta
//	http://www.stata.com/help.cgi?dta_113
//-------------------------------------------------------------------

/*
SPSS file format
// http://cvs.savannah.gnu.org/pspp/doc/data-file-format.texi?root=pspp&content-type=text%2Fplain
*/

#define TEST_CONVERSION 0  // 0: ieee754, 1: SAS converter (big endian bug), 2: experimental
#define DEBUG 1

#include <ctype.h>
#include <math.h>
//#include <sqlite.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <sys/param.h>
#include <time.h>
#include "mex.h"

#ifndef typeof
#define typeof __typeof__
#endif

#ifdef tmwtypes_h
  #if (MX_API_VER<=0x07020000)
    typedef int mwSize;
    typedef int mwIndex;
  #endif
#endif

#define NaN  		(0.0/0.0)
#define fix(m)     	(m<0 ? ceil(m) : floor(m))
#define max(a,b)	(((a) > (b)) ? (a) : (b))
#define min(a,b)	(((a) < (b)) ? (a) : (b))


#if 0

#elif defined(__linux__)
#  include <endian.h>
#  include <byteswap.h>

#elif defined(__CYGWIN__)
#  include <endian.h>
#  include <byteswap.h>

#elif defined(__GLIBC__)	// for Hurd
#  include <endian.h>
#  include <byteswap.h>

#elif defined(__MINGW32__) 
   /* use local version because MINGW does not provide byteswap.h */
#  define __BIG_ENDIAN		4321
#  define __LITTLE_ENDIAN  	1234
#  define __BYTE_ORDER 		__LITTLE_ENDIAN

#elif defined(__NetBSD__)
#  include <sys/bswap.h>
#  define __BIG_ENDIAN _BIG_ENDIAN
#  define __LITTLE_ENDIAN _LITTLE_ENDIAN
#  define __BYTE_ORDER _BYTE_ORDER
#  define bswap_16(x) bswap16(x)
#  define bswap_32(x) bswap32(x)
#  define bswap_64(x) bswap64(x)

#elif defined(_APPLE_) && defined(_MACH_)
#  include <machine/endian.h>
#  define _BYTE_ORDER __DARWIN_BYTE_ORDER
#  define _LITTLE_ENDIAN __DARWIN_LITTLE_ENDIAN
#  define _BIG_ENDIAN __DARWIN_BIG_ENDIAN
#  define bswap_16(x) __bswap16(x)
#  define bswap_32(x) __bswap32(x)
#  define bswap_64(x) __bswap64(x)

#elif defined(__APPLE__)
#  include <CoreFoundation/CFByteOrder.h>
#  define __BIG_ENDIAN       4321
#  define __LITTLE_ENDIAN  1234
#if (defined(__LITTLE_ENDIAN__) && (__LITTLE_ENDIAN__ == 1))
    #define __BYTE_ORDER __LITTLE_ENDIAN
#else
    #define __BYTE_ORDER __BIG_ENDIAN
#endif
#  define bswap_16(x) CFSwapInt16(x)
#  define bswap_32(x) CFSwapInt32(x)
#  define bswap_64(x) CFSwapInt64(x)

#elif (defined(BSD) && (BSD >= 199103)) && !defined(__GLIBC__)
#  include <machine/endian.h>
#  define __BIG_ENDIAN _BIG_ENDIAN
#  define __LITTLE_ENDIAN _LITTLE_ENDIAN
#  define __BYTE_ORDER _BYTE_ORDER
#  define bswap_16(x) __bswap16(x)
#  define bswap_32(x) __bswap32(x)
#  define bswap_64(x) __bswap64(x)

#elif defined(__GNUC__) 
   /* use byteswap macros from the host system, hopefully optimized ones ;-) */
#  include <endian.h>
#  include <byteswap.h>
#  define bswap_16(x) __bswap_16 (x)
#  define bswap_32(x) __bswap_32 (x)
#  define bswap_64(x) __bswap_64 (x)

#elif defined(__sparc__) 
#  define __BIG_ENDIAN  	4321
#  define __LITTLE_ENDIAN  	1234
#  define __BYTE_ORDER 	__BIG_ENDIAN

#else
#  error Unknown platform
#endif 

#if defined(__MINGW32__) || defined(__sparc__) 

# ifndef bswap_16
#  define bswap_16(x)   \
	((((x) & 0xff00) >> 8) | (((x) & 0x00ff) << 8))
# endif

# ifndef bswap_32
#  define bswap_32(x)   \
	 ((((x) & 0xff000000) >> 24) \
        | (((x) & 0x00ff0000) >> 8)  \
	| (((x) & 0x0000ff00) << 8)  \
	| (((x) & 0x000000ff) << 24))

# endif

# ifndef bswap_64
#  define bswap_64(x) \
      	 ((((x) & 0xff00000000000000ull) >> 56)	\
      	| (((x) & 0x00ff000000000000ull) >> 40)	\
      	| (((x) & 0x0000ff0000000000ull) >> 24)	\
      	| (((x) & 0x000000ff00000000ull) >> 8)	\
      	| (((x) & 0x00000000ff000000ull) << 8)	\
      	| (((x) & 0x0000000000ff0000ull) << 24)	\
      	| (((x) & 0x000000000000ff00ull) << 40)	\
      	| (((x) & 0x00000000000000ffull) << 56))
# endif

#endif


#if !defined(__BIG_ENDIAN) && !defined(__LITTLE_ENDIAN) 
#error  ENDIANITY is not known 
#endif 


#if __BYTE_ORDER == __BIG_ENDIAN
#define l_endian_u16(x) ((uint16_t)bswap_16((uint16_t)(x)))
#define l_endian_u32(x) ((uint32_t)bswap_32((uint32_t)(x)))
#define l_endian_u64(x) ((uint64_t)bswap_64((uint64_t)(x)))
#define l_endian_i16(x) ((int16_t)bswap_16((int16_t)(x)))
#define l_endian_i32(x) ((int32_t)bswap_32((int32_t)(x)))
#define l_endian_i64(x) ((int64_t)bswap_64((int64_t)(x)))

#define b_endian_u16(x) ((uint16_t)(x))
#define b_endian_u32(x) ((uint32_t)(x))
#define b_endian_u64(x) ((uint64_t)(x))
#define b_endian_i16(x) ((int16_t)(x))
#define b_endian_i32(x) ((int32_t)(x))
#define b_endian_i64(x) ((int64_t)(x))

#elif __BYTE_ORDER==__LITTLE_ENDIAN
#define l_endian_u16(x) ((uint16_t)(x))
#define l_endian_u32(x) ((uint32_t)(x))
#define l_endian_u64(x) ((uint64_t)(x))
#define l_endian_i16(x) ((int16_t)(x))
#define l_endian_i32(x) ((int32_t)(x))
#define l_endian_i64(x) ((int64_t)(x))

#define b_endian_u16(x) ((uint16_t)bswap_16((uint16_t)(x)))
#define b_endian_u32(x) ((uint32_t)bswap_32((uint32_t)(x)))
#define b_endian_u64(x) ((uint64_t)bswap_64((uint64_t)(x)))
#define b_endian_i16(x) ((int16_t)bswap_16((int16_t)(x)))
#define b_endian_i32(x) ((int32_t)bswap_32((int32_t)(x)))
#define b_endian_i64(x) ((int64_t)bswap_64((int64_t)(x)))

#endif /* __BYTE_ORDER */


/* 
	Including ZLIB enables reading gzipped files (they are decompressed on-the-fly)  
	The output files can be zipped, too. 
 */

#ifdef WITH_ZLIB
#include <zlib.h>
#endif


double xpt2d(uint64_t x);
uint64_t d2xpt(double x);
double tm_time2gdf_time(struct tm *t);

void mexFunction(int POutputCount,  mxArray* POutput[], int PInputCount, const mxArray *PInputs[])
{
	const char L1[] = "HEADER RECORD*******LIBRARY HEADER RECORD!!!!!!!000000000000000000000000000000  ";
	const char L2[] = "SAS     SAS     SASLIB 6.06     bsd4.2                          13APR89:10:20:06";
	//const char L3[] = "";
	const char L4[] = "HEADER RECORD*******MEMBER  HEADER RECORD!!!!!!!000000000000000001600000000140  ";
	const char L5[] = "HEADER RECORD*******DSCRPTR HEADER RECORD!!!!!!!000000000000000000000000000000  ";
	const char L6[] = "SAS     ABC     SASLIB 6.06     bsd4.2                          13APR89:10:20:06";
	//const char L7[] = "";
	const char L8[] = "HEADER RECORD*******NAMESTR HEADER RECORD!!!!!!!000000000200000000000000000000  ";
	const char LO[] = "HEADER RECORD*******OBS     HEADER RECORD!!!!!!!000000000000000000000000000000  ";

	const  char DATEFORMAT[] = "%d%b%y:%H:%M:%S";
	char   *fn = NULL;
	char   Mode[3] = "r";
	size_t count = 0, HeadLen0=80*8, HeadLen2=0, sz2 = 0;
	uint32_t NS = 0;
	char   H0[HeadLen0];
	char   *H2 = NULL;
	char   SWAP = 0;

#ifndef ZLIB_H
	FILE   *fid;
#else
	gzFile  fid;
	#define fopen 		gzopen	
	#define fread(a,b,c,d) 	(gzread(d,a,b*c)/b)	
	#define fwrite(a,b,c,d)	(gzwrite(d,a,b*c)/b)	
	#define feof 		gzeof	
	#define fseek 		gzseek	
	#define fclose 		gzclose	
	#define	rewind(fid) 	(gzseek(fid,0,SEEK_SET)) 
#endif

	// check for proper number of input and output arguments
	if ( PInputCount > 0 && mxGetClassID(PInputs[0])==mxCHAR_CLASS) {
		size_t buflen = (mxGetM(PInputs[0]) * mxGetN(PInputs[0]) * sizeof(mxChar)) + 1;
		fn = (char*)malloc(buflen);
		mxGetString(PInputs[0], fn, buflen);
	}
	else {
		mexPrintf("XPTOPEN read of several file formats and writing of the SAS Transport Format (*.xpt)\n");
		mexPrintf("\n\tX = xptopen(filename)\n");
		mexPrintf("\tX = xptopen(filename,'r')\n");
		mexPrintf("\t\tread filename and return variables in struct X\n");
#ifdef ZLIB_H		
		mexPrintf("\tSupported are ARFF, SAS-XPT and STATA files with or w/o zlib/gzip compression.\n");
#else
		mexPrintf("\tSupported are ARFF, SAS-XPT and STATA files.\n");
#endif
		mexPrintf("\n\tX = xptopen(filename,'w',X)\n");
		mexPrintf("\t\tsave fields of struct X in filename.\n\n");
		mexPrintf("\tThe fields of X must be column vectors of equal length.\n");
		mexPrintf("\tEach vector is either a numeric vector or a cell array of strings.\n");
		mexPrintf("\nThe SAS-XPT format stores Date/Time as numeric value counting the number of days since 1960-01-01.\n\n");
		return;
	}

	if (PInputCount > 1)
	if (mxGetClassID(PInputs[1])==mxCHAR_CLASS && mxGetNumberOfElements(PInputs[1])) {
		mxGetString(PInputs[1],Mode,3);
		Mode[2]=0;
	}

	fid = fopen(fn,Mode);
	if (fid == NULL) {
	        mexErrMsgTxt("Can not open file!\n");
		}

	if (Mode[0]=='r' || Mode[0]=='a' ) {
		count += fread(H0,1,80*8,fid);
		enum FileFormat {
			noFile, unknown, ARFF, SASXPT, SPSS, SQLite, STATA
		};
		enum FileFormat TYPE; 		/* type of file format */
		uint8_t		LittleEndian;   /* 1 if file is LittleEndian data format and 0 for big endian data format*/

		TYPE = unknown;
		if (!memcmp(H0,"$FL2@(#) SPSS DATA FILE",23) || !memcmp(H0,"$FL2@(#) PASW STATISTICS DATA FILE",27)) {
		/*
			SPSS file format
		*/
                        uint32_t M=0; 

		        mexWarnMsgTxt("XPTOPEN: support of for SPSS file format is very experimental (do not use it for production use)\n");

			TYPE = SPSS;
			switch (*(uint32_t*)(H0+64)) {
			case 0x00000002:
			case 0x00000003:
		    		LittleEndian = 1;
		    		SWAP = __BYTE_ORDER==__BIG_ENDIAN;
		    		NS = l_endian_u32(*(uint32_t*)(H0+68));
		    		M  = l_endian_u32(*(uint32_t*)(H0+80));
			    	break;
			case 0x02000000:
			case 0x03000000:
		    		SWAP = __BYTE_ORDER==__LITTLE_ENDIAN;
		    		LittleEndian = 0;
		    		NS = b_endian_u32(*(uint32_t*)(H0+68));
		    		M  = b_endian_u32(*(uint32_t*)(H0+80));
			    	break;
			default:
				TYPE = unknown;
			}
			NS = *(int32_t*)(H0+80);
			M  = *(int32_t*)(H0+80);
			if (SWAP) {
				NS = bswap_32(NS);
				M  = bswap_32(M);
			}
			HeadLen0 = 184;
			char *H2 = (char*)malloc(NS*32);
			size_t c2 = 0;

			/*
				Read Variable SPSS header
			*/
			int ns = 0;
			const char **ListOfVarNames = (const char**)malloc((NS+1) * sizeof(char*));
			char *VarNames   = (char*)malloc((NS+1) * sizeof(char) * 9);
			double *MISSINGS = (double*)malloc((NS+1) * sizeof(double));
			for (uint32_t k=0; k<NS; k++) {
				int32_t rec_type, type, FlagHasLabel, FlagMissing;
				c2 += fread(&rec_type,1,4,fid);
				c2 += fread(&type,1,4,fid);
				c2 += fread(&FlagHasLabel,1,4,fid);
				c2 += fread(&FlagMissing,1,4,fid);
				fseek(fid,4,SEEK_CUR);
				if (SWAP) {
					rec_type     = bswap_32(rec_type);
					type         = bswap_32(type);
					FlagHasLabel = bswap_32(FlagHasLabel);
				}

				if (rec_type != 2) 
					mexErrMsgTxt("invalid SPSS file\n");

				c2 += fread(VarNames+9*ns,1,8,fid);
				VarNames[9*ns+8] = 0;
				ListOfVarNames[ns] = VarNames+9*ns;
				if (FlagHasLabel==1) {
					int32_t LenLabel;
					c2 += fread(&LenLabel,1,4,fid);
					if (SWAP) LenLabel = bswap_32(LenLabel);
					if (LenLabel%4) LenLabel += 4 - LenLabel % 4;
					fseek(fid,LenLabel,SEEK_CUR);
				}
				if (FlagMissing)
					c2 += fread(MISSINGS+ns,1,8,fid);

				if (type != -1) ns++;
			}


			NS = ns;
			mxArray **R = (mxArray**) mxMalloc(NS * sizeof(mxArray*));
			/* ToDo:
				EXTRACT data
			*/

			/* convert into output */
			POutput[0] = mxCreateStructMatrix(1, 1, NS, ListOfVarNames);
			for (uint32_t k = 0; k < NS; k++) {
				 mxSetField(POutput[0], 0, ListOfVarNames[k], R[k]);
			}

			if (MISSINGS) 	    free(MISSINGS);
			if (VarNames) 	    free(VarNames);
			if (ListOfVarNames) free(ListOfVarNames);
			if (H2) 	    free(H2);
		}

		if (TYPE == SPSS) {
/*
The records must appear in the following order:
- File header record.
- Variable records.
- All pairs of value labels records and value label variables records,
if present.
- Document record, if present.
- Any of the following records, if present, in any order:
	Machine integer info record.
	Machine floating-point info record.
	Variable display parameter record.
	Long variable names record.
	Miscellaneous informational records.
- Dictionary termination record.
- Data record.

*/			;
		}

		else if (!memcmp(H0,"SQLite format 3\000",16) && H0[21]==64 && H0[22]==32 && H0[23]==32 ) {
			TYPE = SQLite;

			fclose(fid);
			mexErrMsgTxt("SQLite format not supported yet");
			return;
		}

		else if ((H0[0]>=0x6e || H0[0]<=114) && (H0[1]==1 || H0[1]==2) && H0[2]==1 && H0[3]==0) {
		/*
			STATA File Format
			http://www.stata.com/help.cgi?dta
			http://www.stata.com/help.cgi?dta_113
			Stata files written by R start with 0x6e
		*/
                        uint32_t M=0; 

			TYPE = STATA;
			// Header 119 bytes
	    		LittleEndian = H0[1]==2;
	    		if (LittleEndian) {
	    			NS = l_endian_u16(*(uint16_t*)(H0+4));
	    			M  = l_endian_u32(*(uint32_t*)(H0+6));
	    		}
	    		else {
	    			NS = b_endian_u16(*(uint16_t*)(H0+4));
	    			M  = b_endian_u32(*(uint32_t*)(H0+6));
	    		}

	    		// Descriptors
	    		int fmtlen = (H0[0]==113) ? 12 : 49;
	    		fseek(fid,109,SEEK_SET);
	    		size_t HeadLen2 = 2+NS*(1+33+2+fmtlen+33+81);
			char *H1 = (char*)malloc(HeadLen2);
			HeadLen2 = fread(H1,1,HeadLen2,fid);

			// expansion fields
			char typ; int32_t len;
			char flagSWAP = (((__BYTE_ORDER == __BIG_ENDIAN) && LittleEndian) || ((__BYTE_ORDER == __LITTLE_ENDIAN) && !LittleEndian));
			do {
				if (!fread(&typ,1,1,fid)) break;
				if (!fread(&len,4,1,fid)) break;
				if (flagSWAP) bswap_32(len);
				fseek(fid,len,SEEK_CUR);
			} while (len);
			uint8_t *typlist = (uint8_t*)H1;

/*
			char *varlist = H1+NS;
			char *srtlist;
			char *fmtlist = H1+NS*36+2;
			char *lbllist = H1+NS*(36+fmtlen)+2;
*/

			mxArray **R = (mxArray**) mxMalloc(NS*sizeof(mxArray*));
			size_t *bi = (size_t*) malloc((NS+1)*sizeof(size_t*));
			const char **ListOfVarNames = (const char**)malloc(NS * sizeof(char*));
			bi[0] = 0;
			for (size_t k = 0; k < NS; k++) {
				size_t sz;
				ListOfVarNames[k] = H1+NS+33*k;
				switch (typlist[k]) {
				case 0xfb: sz = 1; break;
				case 0xfc: sz = 2; break;
				case 0xfd: sz = 4; break;
				case 0xfe: sz = 4; break;
				case 0xff: sz = 8; break;
				default: sz = typlist[k];
				}
				bi[k+1] = bi[k]+sz;
			}

			// data
			uint8_t *data = (uint8_t *) malloc(bi[NS] * M);
			fread(data, bi[NS], M, fid);

			char *f = (char*)malloc(bi[NS]+1);
			for (size_t k = 0; k < NS; k++) {
				switch (typlist[k]) {
				case 0xfb:
					R[k] = mxCreateDoubleMatrix(M, 1, mxREAL);
					for (typeof(M) m = 0; m < M; m++) {
						int8_t d = *(int8_t*)(data+bi[k]+m*bi[NS]);
						((double*)mxGetData(R[k]))[m] = (d>100) ? NaN : d;
					}
					break;
				case 0xfc:
					R[k] = mxCreateDoubleMatrix(M, 1, mxREAL);
					if (flagSWAP) for (size_t m = 0; m < M; m++) {
						int16_t d = (int16_t) bswap_16(*(uint16_t*)(data+bi[k]+m*bi[NS]));
						((double*)mxGetData(R[k]))[m] = (d>32740) ? NaN : d;
					}
					else for (typeof(M) m = 0; m < M; m++) {
						int16_t d = *(int16_t*)(data+bi[k]+m*bi[NS]);
						((double*)mxGetData(R[k]))[m] = (d>32740) ? NaN : d;
					}
					break;
				case 0xfd:
					R[k] = mxCreateDoubleMatrix(M, 1, mxREAL);
					if (flagSWAP) for (size_t m = 0; m < M; m++) {
						int32_t d = (int32_t)bswap_32(*(uint32_t*)(data+bi[k]+m*bi[NS]));
						((double*)mxGetData(R[k]))[m] = (d>2147483620) ? NaN : d;
					}
					else for (typeof(M) m = 0; m < M; m++) {
						int32_t d = *(int32_t*)(data+bi[k]+m*bi[NS]);
						((double*)mxGetData(R[k]))[m] = (d>2147483620) ? NaN : d;
					}
					break;
				case 0xfe:
					R[k] = mxCreateNumericMatrix(M, 1, mxSINGLE_CLASS, mxREAL);
					if (flagSWAP) for (size_t m = 0; m < M; m++) {
						((uint32_t*)mxGetData(R[k]))[m] = bswap_32(*(uint32_t*)(data+bi[k]+m*bi[NS]));;
					}
					else for (size_t m = 0; m < M; m++) {
						((uint32_t*)mxGetData(R[k]))[m] = *(uint32_t*)(data+bi[k]+m*bi[NS]);
					}
					break;
				case 0xff:
					R[k] = mxCreateDoubleMatrix(M, 1, mxREAL);
					if (flagSWAP) for (typeof(M) m = 0; m < M; m++) {
						((uint64_t*)mxGetData(R[k]))[m] = bswap_64(*(uint64_t*)(data+bi[k]+m*bi[NS]));
					}
					else for (typeof(M) m = 0; m < M; m++) {
						((uint64_t*)mxGetData(R[k]))[m] = *(uint64_t*)(data+bi[k]+m*bi[NS]);
					}
					break;
				default:
					R[k] =	mxCreateCellMatrix(M, 1);
					size_t sz = typlist[k];
					for (typeof(M) m = 0; m < M; m++) {
						memcpy(f, data+bi[k]+m*bi[NS], sz);
						f[sz] = 0;
						mxSetCell(R[k], m, mxCreateString(f));
					}
				}
			}
			if (f)  free(f);
			if (H1) free(H1);
			if (bi) free(bi);

			/* convert into output */
			POutput[0] = mxCreateStructMatrix(1, 1, NS, ListOfVarNames);
			for (size_t k = 0; k < NS; k++) {
				 mxSetField(POutput[0], 0, ListOfVarNames[k], R[k]);
			}

			if (ListOfVarNames) free(ListOfVarNames);
		}

		else if (H0[0]=='%' || H0[0]=='@') {
		/*
			 ARFF
		*/
                        uint32_t M=0; 

			TYPE = ARFF;
			rewind(fid);

			char *H1 = NULL;
			count = 0;
			size_t ns = 0;
			char *vartyp = NULL;
			char **datestr = NULL;
			const char **ListOfVarNames = NULL;
			mxArray **R = NULL;
			size_t m = 0;

			while (!feof(fid)) {
				HeadLen0 = max(1024,HeadLen0*2);
				H1 = (char*)realloc(H1,HeadLen0);
				count += fread(H1+count,1,HeadLen0-count-1,fid);
			}
			H1[count]   = 0;

			switch (H1[count-1]) {
			case 0x0a:
			case 0x0d:
				H1[count]   = 0;
				break;
			default:
				H1[count]   = 0x0a;
			}
			H1[count+1] = 0;

			char *line = strtok(H1,"\x0a\0x0d");

			int status = 0;
			while (line) {

				if (!strncasecmp(line,"@relation",9)) {
					status = 1;
				}

				else if (status == 1 && !strncasecmp(line,"@attribute",10)) {
					if (ns<=NS) {
						ns = max(16, ns*2);
						ListOfVarNames = (const char**)realloc(ListOfVarNames,ns*sizeof(char*));
						vartyp         = (char*)realloc(vartyp,ns*sizeof(char));
						R              = (mxArray**) mxRealloc(R,ns*sizeof(mxArray*));
					}
					size_t k = 10;
					char *p1, *p2;
					while (isspace(line[k])) k++;
					p1 = line+k;
					while (!isspace(line[k])) k++;
					line[k++]=0;
					while (isspace(line[k])) k++;
					p2 = line+k;

					ListOfVarNames[NS] = p1;
					if      (!strncasecmp(p2,"numeric",7)) {
						vartyp[NS] = 1;
					}
					else if (!strncasecmp(p2,"integer",7)) {
						vartyp[NS] = 2;
					}
					else if (!strncasecmp(p2,"real",4)) {
						vartyp[NS] = 3;
					}
					else if (!strncasecmp(p2,"string",6)) {
						vartyp[NS] = 4;
					}
					else if (!strncasecmp(p2,"{",1)) {
						vartyp[NS] = 5;
					}
					else if (!strncasecmp(p2,"date",4)) {
						vartyp[NS] = 6;
						datestr = (char**)realloc(datestr,(NS+1)*sizeof(char*));
						p2+=4;
						while (isspace(*p2)) p2++;
						datestr[NS] = p2;
						if (p2[0]==34) {
							p2++;
							while (p2[0]!=34 && p2[0]) p2++;
							p2[1]=0;
						}
					}
					else if (!strncasecmp(p2,"relational",10)) {
						vartyp[NS] = 7;
					}
					else vartyp[NS] = 99;

					NS++;
				}

				else if (status == 1 && !strncasecmp(line,"@data",5)) {
					status = 2;
					char *p = line;
					while (*p) p++;  // goto end of current line
					p++; 		   // skip \x00
					M = 0;
					while (*p) {
						if (p[0]==0x0a || p[0]==0x0d) {
							// count number of <CR>
							M++;
							// skip next char (deals with <CR><NL>)
							p+=2;
						}
						else p++;
					}
					for (size_t k=0; k<NS; k++) {
						if (vartyp[k]==4 || vartyp[k]==5)
							R[k] = mxCreateCellMatrix(M, 1);
						else
							R[k] = mxCreateDoubleMatrix(M, 1, mxREAL);
					}
				}

				else if (status == 2) {

					size_t p = 0,k;
					for (ns = 0; ns<NS; ns++) {
						// read next token
						while (isspace(line[p])) p++;
						if (line[p]==39) {
							p++; k=p;
							while (line[k]!=39 && line[k]) k++;
							// if (!line[k]) ; // error
							line[k++] = 0;
						}
						else
							k=p;
						while (line[k] != ',' && line[k] != 0) k++;
						line[k] = 0;

						if (vartyp[ns] < 4) {
							double d = atof(line+p);
							*(mxGetPr(R[ns])+m) = d;
						}
						else if (vartyp[ns] < 6) {
							mxSetCell(R[ns], m, mxCreateString(line+p));
						}
						else if (vartyp[ns] == 6) {
							size_t kk[6],n=0, N=strlen(datestr[ns]);
							char T0[6][5];
							int ix = 0;
							struct tm t;

							for (n=0; n < N; n++) {
								switch (datestr[ns][n]) {
								case 'Y':
									ix = 0;
									break;
								case 'M':
									ix = 1;
									break;
								case 'd':
									ix = 2;
									break;
								case 'H':
									ix = 3;
									break;
								case 'm':
									ix = 4;
									break;
								case 's':
									ix = 5;
									break;
								default:
									ix = 99;
								}

								if (ix < 6) {
									T0[ix][kk[ix]++] = line[p+n];
								}
							}
							for (n=0; n<6; n++) {
								T0[n][kk[n]] = 0;
							}
							t.tm_year = atoi(T0[0]);
							t.tm_mon  = atoi(T0[1]);
							t.tm_mday = atoi(T0[2]);
							t.tm_hour = atoi(T0[3]);
							t.tm_min  = atoi(T0[4]);
							t.tm_sec  = atoi(T0[5]);

							*(mxGetPr(R[ns])+m) = tm_time2gdf_time(&t);
						}
						p = k+1;
					}
					m++;
				}
				line = strtok(NULL, "\x0a\x0d");
			}

			/* convert into output */
			POutput[0] = mxCreateStructMatrix(1, 1, NS, ListOfVarNames);
			for (size_t k = 0; k < NS; k++) {
				mxSetField(POutput[0], 0, ListOfVarNames[k], R[k]);
			}

			if (ListOfVarNames) free(ListOfVarNames);
			if (vartyp) free(vartyp);
			if (datestr) free(datestr);
			if (H1) free(H1);
		}

		else if (!memcmp(H0,"HEADER RECORD*******LIBRARY HEADER RECORD!!!!!!!000000000000000000000000000000",78)) {
		/*
			 SAS Transport file format (XPORT)
		*/
                        size_t M=0; 
			TYPE = SASXPT;

			/* TODO: sanity checks */

			char tmp[5];
			memcpy(tmp,H0+7*80+54,4);
			tmp[4] = 0;
			NS = atoi(tmp);

			char *tmp2;
			sz2 = strtoul(H0+4*80-6, &tmp2, 10);

			HeadLen2 = NS*sz2;
			if (HeadLen2 % 80) HeadLen2 = (HeadLen2 / 80 + 1) * 80;

			/* read namestr header, and header line "OBS" */
			H2 = (char*) realloc(H2, HeadLen2+81);
			count  += fread(H2,1,HeadLen2+80,fid);

			/* size of single record */
			size_t pos=0, recsize = 0, POS = 0;
			for (size_t k = 0; k < NS; k++)
				recsize += b_endian_u16(*(int16_t*)(H2+k*sz2+4));

			/* read data section */
			size_t szData = 0;
			uint8_t *Data = NULL;
			while (!feof(fid)) {
				size_t szNew = max(16,szData*2);
				Data         = (uint8_t*)realloc(Data,szNew);
				szData      += fread(Data+szData,1,szNew-szData,fid);
			}

			M = szData/recsize;

			mxArray **R = (mxArray**) mxMalloc(NS*sizeof(mxArray*));
			const char **ListOfVarNames = (const char**)malloc(NS * sizeof(char*));
			char *VarNames = (char*)malloc(NS * 9);

			for (size_t k = 0; k < NS; k++) {
				size_t maxlen = b_endian_u16(*(int16_t*)(H2+k*sz2+4));

				ListOfVarNames[k] = VarNames+pos;
				size_t n = k*sz2+8;
				// int flagDate = (!memcmp(H2+n+48,"DATE    ",8) || !memcmp(H2+n+48,"MONNAME ",8)); // not used
				do {
					VarNames[pos++] = H2[n];
				} while (isalnum(H2[++n]) && (n < k*sz2+16));


				VarNames[pos++] = 0;

				if ((*(int16_t*)(H2+k*sz2)) == b_endian_u16(1) && (*(int16_t*)(H2+k*sz2+4)) == b_endian_u16(8) ) {
					// numerical data
					R[k] = mxCreateDoubleMatrix(M, 1, mxREAL);
					for (size_t m=0; m<M; m++) {
						double d = xpt2d(b_endian_u64(*(uint64_t*)(Data+m*recsize+POS)));

//						if (flagDate) d += 715876;  // add number of days from 0000-Jan-01 to 1960-Jan-01

						*(mxGetPr(R[k])+m) = d;

					}
				}
				else if ((*(int16_t*)(H2+k*sz2)) == b_endian_u16(2)) {
					// character string
					R[k] = mxCreateCellMatrix(M, 1);
					char *f = (char*)malloc(maxlen+1);
					for (size_t m=0; m<M; m++) {
						memcpy(f, Data+m*recsize+POS, maxlen);
						f[maxlen] = 0;
						mxSetCell(R[k], m, mxCreateString(f));
					}
					if (f) free(f);
				}
				POS += maxlen;
			}

			POutput[0] = mxCreateStructMatrix(1, 1, NS, ListOfVarNames);
			for (size_t k = 0; k < NS; k++) {
				 mxSetField(POutput[0], 0, ListOfVarNames[k], R[k]);
			}

			if (VarNames) 	    free(VarNames);
			if (ListOfVarNames) free(ListOfVarNames);
			if (Data)	    free(Data);
			/* end of reading SAS format */
		}

		else {
			fclose(fid);
			mexErrMsgTxt("file format not supported");
			return;
		}



	}

//	if (Mode[0]=='w' || Mode[0]=='a' ) {
	if (Mode[0]=='w') {

		NS += mxGetNumberOfFields(PInputs[2]);

		// generate default (fixed) header
		if (Mode[0]=='w') {
			memset(H0,' ',80*8);
			memcpy(H0,L1,strlen(L1));
			memcpy(H0+80,L2,strlen(L2));

			memcpy(H0+80*3,L4,strlen(L4));
			memcpy(H0+80*4,L5,strlen(L5));
			memcpy(H0+80*5,L6,strlen(L6));

			memcpy(H0+80*7,L8,strlen(L8));
		}

		time_t t;
		time(&t);
		char tt[20];
		strftime(tt, 17, DATEFORMAT, localtime(&t));
		memcpy(H0+80*2-16,tt,16);
		memcpy(H0+80*2,tt,16);
		memcpy(H0+80*5+8,fn,min(8,strcspn(fn,".\x00")));
		memcpy(H0+80*5+32,"XPTOPEN.MEX (OCTAVE/MATLAB)",27);
		memcpy(H0+80*6-16,tt,16);
		memcpy(H0+80*6,tt,16);

		char tmp[17];
		sprintf(tmp,"%04i", NS);	// number of variables
		memcpy(H0+80*7+54, tmp, 4);

		if (sz2==0) sz2 = 140;
		if (sz2 < 136)
			mexErrMsgTxt("error XPTOPEN: incorrect length of namestr field");

		/* generate variable NAMESTR header */
		HeadLen2 = NS*sz2;
		if (HeadLen2 % 80) HeadLen2 = (HeadLen2 / 80 + 1) * 80;
		H2 = (char*) realloc(H2,HeadLen2);
		memset(H2,0,HeadLen2);

		size_t M = 0;
		mxArray **F = (mxArray**) mxMalloc(NS*sizeof(mxArray*));
		char **Fstr = (char**) malloc(NS*sizeof(char*));
		size_t *MAXLEN = (size_t*) malloc(NS*sizeof(size_t*));
		for (uint16_t k = 0; k < NS; k++) {
			Fstr[k] = NULL;
			MAXLEN[k]=0;
			F[k] = mxGetFieldByNumber(PInputs[2],0,k);
			if (k==0) M = mxGetM(F[k]);
			else if (M != mxGetM(F[k])) {
				if (H2) free(H2);
				if (F)  free(F);
				mexErrMsgTxt("Error XPTOPEN: number of elements (rows) do not fit !!!");
			}

			if (mxIsChar(F[k])) {
				*(int16_t*)(H2+k*sz2) = b_endian_u16(2);
				*(int16_t*)(H2+k*sz2+4) = b_endian_u16(mxGetN(F[k]));
			}
			else if (mxIsCell(F[k])) {
				size_t maxlen = 0;
				for (size_t m = 0; m<M; m++) {
					mxArray *f = mxGetCell(F[k],m);
					if (mxIsChar(f) || mxIsEmpty(f)) {
						size_t len = mxGetNumberOfElements(f);
						if (maxlen<len) maxlen = len;
					}
				}
				Fstr[k] = (char*) malloc(maxlen+1);
				*(int16_t*)(H2+k*sz2) = b_endian_u16(2);
				*(int16_t*)(H2+k*sz2+4) = b_endian_u16(maxlen);
				MAXLEN[k] = maxlen;
			}
			else {
				*(int16_t*)(H2+k*sz2) = b_endian_u16(1);
				*(int16_t*)(H2+k*sz2+4) = b_endian_u16(8);
			}
			*(int16_t*)(H2+k*sz2+6) = b_endian_u16(k);
			strncpy(H2+k*sz2+8,mxGetFieldNameByNumber(PInputs[2],k),8);
		}

		count  = fwrite(H0, 1, HeadLen0, fid);
		count += fwrite(H2, 1, HeadLen2, fid);
		/* write OBS header line */
		count += fwrite(LO, 1, strlen(LO), fid);
		for (size_t m = 0; m < M; m++) {
			for (uint16_t k = 0; k < NS; k++) {

				if (*(int16_t*)(H2+k*sz2) == b_endian_u16(1)) {
					// numeric
					uint64_t u64 = b_endian_u64(d2xpt(*(mxGetPr(F[k])+m)));
					count += fwrite(&u64, 1, 8, fid);
				}

/*				else if (mxIsChar(F[k])) {
					*(int16_t*)(H2+k*sz2) = b_endian_u16(2);
					*(int16_t*)(H2+k*sz2+4) = b_endian_u16(mxGetN(F[k]));
				}
*/
				else if (mxIsCell(F[k])) {
					size_t maxlen = MAXLEN[k];
					mxArray *f = mxGetCell(F[k],m);
					mxGetString(f, Fstr[k], maxlen+1);
					count += fwrite(Fstr[k], 1, maxlen, fid);
				}
			}
		}
/*
		// padding to full multiple of 80 byte:
		// FIXME: this might introduce spurios sample values
		char d = count%80;
		while (d--) fwrite("\x20",1,1,fid);
*/
		// free memory
		for (size_t k=0; k<NS; k++) if (Fstr[k]) free(Fstr[k]);
		if (Fstr)   free(Fstr);
		if (MAXLEN) free(MAXLEN);
		Fstr = NULL;
	}
	fclose(fid);

	if (H2) free(H2);
	H2 = NULL;

	return;
}


/*
	XPT2D converts from little-endian IBM to little-endian IEEE format
*/

double xpt2d(uint64_t x) {
	// x is little-endian 64bit IBM floating point format
#if __BYTE_ORDER == __BIG_ENDIAN

	mexErrMsgTxt("IEEE-to-IBM conversion on big-endian platform not supported, yet");
	return(NaN);

#elif __BYTE_ORDER==__LITTLE_ENDIAN
	char c = *((char*)&x+7) & 0x7f;
	uint64_t u = x;
	*((char*)&u+7)=0;


#if DEBUG
	mexPrintf("xpt2d(%016Lx): [0x%x]\n",x,c);
#endif

	// missing values
	if ((c==0x2e || c==0x5f || (c>0x40 && c<0x5b)) && !u )
		return(NaN);

	int s,e;
	s = *(((char*)&x) + 7) & 0x80;		// sign
	e = (*(((char*)&x) + 7) & 0x7f) - 64;	// exponent
	*(((char*)&x) + 7) = 0; 		// mantisse x

#if DEBUG
	mexPrintf("%x %x %016Lx\n",s,e,x);
#endif

	double y = ldexp((double)x, e*4-56);
	if (s) return(-y);
	else   return( y);

#endif
}


/*
	D2XPT converts from little-endian IEEE to little-endian IBM format
*/

uint64_t d2xpt(double x) {
#if __BYTE_ORDER == __BIG_ENDIAN

	mexErrMsgTxt("IEEE-to-IBM conversion on big-endian platform not supported, yet");
	return(0x2eLL << 56);	// NaN - not a number

#elif __BYTE_ORDER==__LITTLE_ENDIAN
	uint64_t s,m;
	int e;


	if (x != x) return(0x2eLL << 56);	// NaN - not a number

	if (fabs(x) == 1.0/0.0) return(0x5fLL << 56); 	// +-infinity

	if (x == 0.0) return(0);

	if (x > 0.0) s=0;
	else s=1;

	x = frexp(x,&e);

#if DEBUG
	mexPrintf("d2xpt(%f)\n",x);
#endif
	// see http://old.nabble.com/Re%3A-IBM-integer-and-double-formats-p20428979.html
	memcpy(&m, &x, 8);
	*(((char*)&m) + 6) &= 0x0f; //
	if (e) *(((char*)&m) + 6) |= 0x10; // reconstruct implicit leading '1' for normalized numbers
	m <<= (3-(-e & 3));
	*(((uint8_t*)&m) + 7)  = s ? 0x80 : 0;
	e = (e + (-e & 3)) / 4 + 64;

	if (e >= 128) return(0x5f); // overflow
	if (e < 0) {
		uint64_t h = 1<<(4*-e - 1);
		m = m / (2*h) + (m & h && m & (3*h-1) ? 1 : 0);
		e = 0;
	}
	return (((uint64_t)e)<<56 | m);

#endif

}


double tm_time2gdf_time(struct tm *t) {
	/* based Octave's datevec.m
	it referes Peter Baum's algorithm at http://vsg.cape.com/~pbaum/date/date0.htm
	but the link is not working anymore as of 2008-12-03.

	Other links to Peter Baum's algorithm are
	http://www.rexswain.com/b2mmddyy.rex
	http://www.dpwr.net/forums/index.php?s=ecfa72e38be61327403126e23aeea7e5&showtopic=4309
	*/

	int Y,M,s; //h,m,
	double D;

	D = t->tm_mday;
	M = t->tm_mon+1;
	Y = t->tm_year+1900;

	// Set start of year to March by moving Jan. and Feb. to previous year.
  	// Correct for months > 12 by moving to subsequent years.
  	Y += (int)fix ((M-14.0)/12);

  	const int monthstart[] = {306, 337, 0, 31, 61, 92, 122, 153, 184, 214, 245, 275};
	// Lookup number of days since start of the current year.
  	D += monthstart[t->tm_mon % 12] + 60;

	// Add number of days to the start of the current year. Correct
  	// for leap year every 4 years except centuries not divisible by 400.
  	D += 365*Y + floor (Y/4) - floor (Y/100) + floor (Y/400);

  	// Add fraction representing current second of the day.
  	s = t->tm_hour*3600 + t->tm_min*60 + t->tm_sec;

	// s -= timezone;
	return(D + s/86400.0);
}

