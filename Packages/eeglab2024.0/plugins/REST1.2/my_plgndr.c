/*
 *  The original plgndr function was implemented using some code from
 *  Numerical Recipes in C. However, it was not allowed to release that code
 *  as open source. The new implementation below is using some code from
 *  the GNU Scientific Library (http://www.gnu.org/software/gsl).
 *
 *  Copyright (C) 2002-2006 Robert Oostenveld
 *  Copyright (C) 2006, Thomas Hartmann
 *  Copyright (C) 2016, Ricardo Bruna
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/*  Based on FieldTrip 20160222 functions:
 *  * plgndr
 */


#include <string.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"

/* Solution for the Legendre polinomial at l=m. */
double legendre_Pmm ( double m, double x ) {
    
    if ( m == 0 )
        return 1.0;
    
    double p_mm = 1.0;
    double root_factor = sqrt ( 1.0 - x * x );
    double fact_coeff = 1.0;
    int i;
    for ( i = 1; i <= m; i ++ ) {
        p_mm *= -fact_coeff * root_factor;
        fact_coeff += 2.0;
    }
    
    return p_mm;
}

/* Iterative solution for the Legendre polinomial. */
void plgndr ( int l, int m, double *x, int nx, double *p ) {
    int index;
    int iter;
    
    /* Uses the analytical formula to construct the first iteration. */
    for ( index = 0; index < nx; index ++ )
        p [ index + m * nx ] = legendre_Pmm ( m, x [ index ] );
    
    if ( l == m )
        return;
    
    /* Uses a simplified version of the Legendre formula for the second iteration. */
    for ( index = 0; index < nx; index ++ )
        p [ index + ( m + 1 ) * nx ] = x [ index ] * ( 2 * ( m + 1 ) - 1 ) * p [ index + m * nx ];
    
    /* For values greater than 2 uses the iterative construction. */
    for ( iter = 2; iter <= l - m; iter ++ )
        for ( index = 0; index < nx; index ++ )
            p [ index + ( m + iter ) * nx ] = ( x [ index ] * ( 2 * ( m + iter ) - 1 ) * p [ index + ( m + iter - 1 ) * nx ] - ( iter + 2 * m - 1 ) * p [ index + ( m + iter - 2 ) * nx ] ) / iter;
}

void mexFunction ( int nlhs, mxArray * plhs [], int nrhs, const mxArray * prhs [] ) {
    int l, m;
    double *x;
    double *pd;
    
    int n, index;
    
    /* Checks the inputs. */
    if ( nrhs != 3 ) mexErrMsgTxt ( "Invalid number of arguments for PLGNDR." );
    if ( !mxIsDouble ( prhs [2] ) ) mexErrMsgTxt ( "This function requires double data as input." );
    
    /* Gets the input variables */
    l = mxGetScalar ( prhs [0] );
    m = mxGetScalar ( prhs [1] );
    x = mxGetData   ( prhs [2] );
    n = mxGetNumberOfElements ( prhs [2] );
    
    /* Checks the inputs. */
    if ( m < 0 ) mexErrMsgTxt ( "m must be a positive integer." );
    if ( m > l ) mexErrMsgTxt ( "l must be a positive integer equal or greater than m." );
    for ( index = 0; index < n; index ++ )
        if ( fabs ( x [ index ] ) > 1.0 ) mexErrMsgTxt ( "abs(x) cannot be greater than 1." );
    
    
    /* Reserves memory for the output array. */
    plhs [0] = mxCreateDoubleMatrix ( n, l + 1, mxREAL );
    pd = mxGetData ( plhs [0] );
    
    /* Calculates the Legendre polinomial. */
    plgndr ( l, m, x, n, pd );
}
