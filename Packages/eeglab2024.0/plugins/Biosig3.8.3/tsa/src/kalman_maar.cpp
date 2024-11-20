// Copyright by Martin Billinger

// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with
// this program; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place - Suite 330, Boston, MA  02111-1307, USA.

#include <mex.h>
#include <memory.h>

#include <blas.h>
#include <lapack.h>

#include "kalman_maar.h"

// input to the matlab function. for convenience elements
// can be accessed by name as well as array element.
union mexin
{
	struct {
		const mxArray *y, *p, *UC, *Kalman, *mode;
	};
	const mxArray *list[5];
};

// output of the matlab function. for convenience elements
// can be accessed by name as well as array element.
union mexout
{
	struct {
		mxArray *x, *e, *Kalman;
	};
	mxArray *list[3];
};

struct mxKalman
{
	mxArray *G, *H, *Kp, *Q1, *Q2, *errCov, *x, *yr;
	struct Kalman pointers;
};

// true when persistent variables are initialized
static bool initialized = false;

// holds persistent variables used by the MVAAR algorithm
static Kalman_Helper *kh = NULL;

// parameters
int M = 0, L = 0, p = 0;

void initialize( ptrdiff_t M, ptrdiff_t p )
{
    mexPrintf( "Initializing Persistent variables.\n" );
	kh = new Kalman_Helper( M, p );
    initialized = true;
}

void cleanup( )
{
    mexPrintf( "Destroying Persistent variables.\n" );
	delete kh;
	kh = NULL;
    initialized = false;
}

void mexFunction( 
	int nlhs, mxArray *plhs[],			// output
	int nrhs, const mxArray *prhs[] )	// input
{
	union mexin in;
	union mexout out;
	struct mxKalman inKalman, outKalman;

	// some pointers for various uses
	double *tmp1, *tmp2, *tmp3;

	int n;
	double UC;
	double *y, *err;
	double mode;

	if( nrhs < 3 )
		mexErrMsgTxt( "Not enough input arguments" );
		
	if( nrhs > 5 )
		mexErrMsgTxt( "Too many input arguments" );

	for( n=0; n<nrhs; n++ )
		in.list[n] = prhs[n];
		
	if( nrhs < 5 )
		mode = 7;
	else
		mode = *mxGetPr( in.mode );
		 
	y = mxGetPr( in.y );

	// if any parameters is changed the persistent variables have to be
	// initialized anew to accomodate for eventually changed memory requirements
	if( ( p!=(int)mxGetPr(in.p)[0] ) || ( M!=mxGetN(in.y) ) )
		if( initialized )
			cleanup( );
    
	if( !initialized )
    {
       	p = (int)mxGetPr(in.p)[0];
		M = mxGetN( in.y );
		L = M*M*p;
		initialize( M, p );
		mexAtExit( cleanup );
	}

	UC = mxGetPr( in.UC )[0];
    
    // create output variables	
	{	// (in C variables can only be declared et the beginning of a block)
		const char *fieldnames[8] = { "G", "H", "Kp", "Q1", "Q2", "errCov", "x", "yr" };
		plhs[0] = mxCreateDoubleMatrix( 1, L, mxREAL );			// x
		plhs[1] = mxCreateDoubleMatrix( 1, M, mxREAL );			// err
		plhs[2] = mxCreateStructMatrix( 1, 1, 8, fieldnames );	// Kalman	
	}
	
	for( n=0; n<3; n++ )
		out.list[n] = plhs[n];
		
	err = mxGetPr( out.e );

	if( nrhs == 3 ) // no filter-state argument given. initialize filter
	{		
		printf( "Initialising Kalman Filter\n" );
		outKalman.G = mxCreateDoubleMatrix( L, M, mxREAL );
		outKalman.H = mxCreateDoubleMatrix( M, L, mxREAL );
		outKalman.Kp = eye( L, 1.0 );
		outKalman.Q1 = eye( L, UC );
		outKalman.Q2 = eye( M, 1.0 );
        outKalman.errCov = eye( M, 1.0 );
		outKalman.x = mxCreateDoubleMatrix(L,1,mxREAL);
		outKalman.yr = mxCreateDoubleMatrix(1,M*p,mxREAL);
        
		// Kalman.yr = [y(:)' zeros(1,M*(p-1))];
		// err = y;
		/*tmp1 = mxGetPr( in.y );
		tmp2 = mxGetPr( outKalman.yr );
		tmp3 = mxGetPr( out.e );
		for( n=0; n<M; n++ )
		{
			*(tmp2++) = *(tmp1);
			*(tmp3++) = *(tmp1++);
		}*/
	}
	else
	{
		// extract filter-state from input data
		inKalman.G = mxGetField( in.Kalman, 0, "G" );
		inKalman.H = mxGetField( in.Kalman, 0, "H" );
		inKalman.Kp = mxGetField( in.Kalman, 0, "Kp" );
		inKalman.Q1 = mxGetField( in.Kalman, 0, "Q1" );
		inKalman.Q2 = mxGetField( in.Kalman, 0, "Q2" );
        inKalman.errCov = mxGetField( in.Kalman, 0, "errCov" );
		inKalman.x = mxGetField( in.Kalman, 0, "x" );
		inKalman.yr = mxGetField( in.Kalman, 0, "yr" );
        
		// set pointers
		inKalman.pointers.G = mxGetPr( inKalman.G );
		inKalman.pointers.H = mxGetPr( inKalman.H );
		inKalman.pointers.Kp = mxGetPr( inKalman.Kp );
		inKalman.pointers.Q1 = mxGetPr( inKalman.Q1 );
		inKalman.pointers.Q2 = mxGetPr( inKalman.Q2 );
		inKalman.pointers.errCov = mxGetPr( inKalman.errCov );
		inKalman.pointers.x = mxGetPr( inKalman.x );
		inKalman.pointers.yr = mxGetPr( inKalman.yr );

		// initialize filter-state output
		outKalman.G = mxCreateDoubleMatrix( L, M, mxREAL );
		outKalman.H = mxCreateDoubleMatrix( M, L, mxREAL );
		outKalman.Kp = mxCreateDoubleMatrix(L, L, mxREAL );
		outKalman.Q1 = mxCreateDoubleMatrix(L, L, mxREAL );
		outKalman.Q2 = mxCreateDoubleMatrix(M, M, mxREAL );
		outKalman.errCov = mxCreateDoubleMatrix(M, M, mxREAL );
		outKalman.x = mxCreateDoubleMatrix( L, 1, mxREAL );
		outKalman.yr = mxCreateDoubleMatrix(1,M*p,mxREAL );
		
		// set pointers
		outKalman.pointers.G = mxGetPr( outKalman.G );
		outKalman.pointers.H = mxGetPr( outKalman.H );
		outKalman.pointers.Kp = mxGetPr( outKalman.Kp );
		outKalman.pointers.Q1 = mxGetPr( outKalman.Q1 );
		outKalman.pointers.Q2 = mxGetPr( outKalman.Q2 );
		outKalman.pointers.errCov = mxGetPr( outKalman.errCov );
		outKalman.pointers.x = mxGetPr( outKalman.x );
		outKalman.pointers.yr = mxGetPr( outKalman.yr );

		// filtering step
		kalman_maar( err, y, UC, M, L, p, 
			&inKalman.pointers, &outKalman.pointers, kh, mode );
	}
	
	// x = Kalman.x'
	memcpy( mxGetPr( out.x ), mxGetPr( outKalman.x ), L*sizeof(double) );
    
	mxSetField( out.Kalman, 0, "G", outKalman.G );
	mxSetField( out.Kalman, 0, "H", outKalman.H );
	mxSetField( out.Kalman, 0, "Kp", outKalman.Kp );
	mxSetField( out.Kalman, 0, "Q1", outKalman.Q1 );
	mxSetField( out.Kalman, 0, "Q2", outKalman.Q2 );
	mxSetField( out.Kalman, 0, "errCov", outKalman.errCov );
	mxSetField( out.Kalman, 0, "x", outKalman.x );
	mxSetField( out.Kalman, 0, "yr", outKalman.yr );
}
