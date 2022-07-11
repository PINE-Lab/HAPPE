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

// struct Kalman
//  holds all information to describe a 
//  kalman filter and it's current state
struct Kalman
{
	// C-style pointers to matrices and vectors that
	// describe a kalman filter and it's current state	
	double *G, *H, *Kp, *Q1, *Q2, *errCov, *x, *yr;
};

// struct Kalman_Helper
//  holds temporary variables that are
//  used by the algorithm
struct Kalman_Helper
{
	Kalman_Helper( ptrdiff_t M, ptrdiff_t p )
	{
		ptrdiff_t L = M*M*p;
		ye = new double[M];
		ipiv = new ptrdiff_t[M];
		tmp = new double[M*M];
		KpH = new double[L*M];
		HKp = new double[M*L];
		HKpH = new double[M*M];
		K = new double[L*L];
	}
	
	~Kalman_Helper( )
	{
		delete[] ye;
		delete[] ipiv;
		delete[] tmp;
		delete[] KpH;
		delete[] HKp;
		delete[] HKpH;
		delete[] K;
	}
	
	double *ye;
	ptrdiff_t *ipiv;
	double *tmp;
	double *KpH;
	double *HKp;
	double *HKpH;
	double *K;
};

// BLAS' and LAPACK's fortran functions need
// all their arguments to be passed by reference.
// Therefore we need constants that can be referenced.
double one = 1.0, zero = 0.0; 
double minusone = -1.0;
ptrdiff_t inc = 1;
char *chn="N"; char *cht="T";

// Helper functions to create identity matrices
void deye( int N, double *i, double value );
mxArray *eye( int N, double value );


// Implementation of the mvaar algorithm
void kalman_maar( double *err, const double *y,
	//double UC, int M, int L, int p,
	double UC, ptrdiff_t M, ptrdiff_t L, ptrdiff_t p,
	struct Kalman *inKalman, struct Kalman *outKalman,
	struct Kalman_Helper *kh, int mode )
{
	double traceK;
	double *tmp1, *tmp2, *tmp3;
	ptrdiff_t n,i,j, info, index;
    
	// % update Measurement Matrix
	// Kalman.H = kron( eye(M), Kalman.yr );
	index = 0;
	for( j=0; j<M; j++ ) {
		for( i=0; i<M*p; i++ ) {
			outKalman->H[index] = inKalman->yr[i];
			index += M;
		}
		index += 1;
	}

	// % calculate prediction error
	// ye = (Kalman.H*Kalman.x)';
	// err = y - ye
	index = 0;
	for( j=0; j<M; j++ ) {
		kh->ye[j] = 0;
		for( i=0; i<M*p; i++ )
			kh->ye[j] += inKalman->yr[i] * inKalman->x[index++];
		err[j] = y[j] - kh->ye[j];
	}
    
    
	// % Adaptive Error covariance
	// Kalman.errCov = (1-UC)*Kalman.errCov + UC*(err'*err);
	index = 0;
	for( i=0; i<M; i++ )
        	for( j=0; j<M; j++ ) {
				outKalman->errCov[index] = (1-UC)*inKalman->errCov[index] + UC * err[i] * err[j];
				index++;
        	}
    
    if( mode == 0 )
    {
		// Kalman.Q2 = Kalman.Q2;	% No update
        memcpy( outKalman->Q2, inKalman->Q2, sizeof(double)*M*M );
    }
    else
    {
		// % update Q2 using the prediction error
		// % of the previous step (corresponds to BioSigs's eMode==4)
		// Kalman.Q2 = Kalman.errCov;
        memcpy( outKalman->Q2, inKalman->errCov, sizeof(double)*M*M );
    }
    
    // KpH = Kalman.Kp * Kalman.H';
    dgemm( chn, cht, &L, &M, &L, &one, inKalman->Kp, &L,
    	outKalman->H, &M, &zero, kh->KpH, &L );
    
    // HKp = Kalman.H * Kalman.Kp;
    dgemm( chn, chn, &M, &L, &L, &one, outKalman->H, &M,
    	inKalman->Kp, &L, &zero, kh->HKp, &M );
    
    // HKpH = Kalman.H * KpH
    dgemm( chn, chn, &M, &M, &L, &one, outKalman->H, &M, 
    	kh->KpH, &L, &zero, kh->HKpH, &M );
					
    
    // % Kalman Gain
    // Kalman.G = KpH * inv( HKpH + Kalman.Q2 );        
    tmp1 = kh->HKpH;
    tmp2 = outKalman->Q2;
    for( n=0; n<M*M; n++ )
        *(tmp1++) += *(tmp2++);
    memset( kh->tmp, 0, M*M*sizeof( double ) );
    deye( M, kh->tmp, 1.0 );
    dgesv( &M, &M, kh->HKpH, &M, kh->ipiv,
    	kh->tmp, &M, &info );
    dgemm( chn, chn, &L, &M, &M, &one, kh->KpH, &L,
    	kh->tmp, &M, &zero, outKalman->G, &L );
    
    // % calculation of the a-posteriori state 
    // % error covariance matrix
    // K = Kalman.Kp-Kalman.G*HKp;
    tmp1 = kh->K;
    tmp2 = inKalman->Kp;
    for( n=0; n<L*L; n++ )
        *(tmp1++) = *(tmp2++);
    dgemm( chn, chn, &L, &L, &M, &minusone, outKalman->G,
    	&L, kh->HKp, &M, &one, kh->K, &L );
		
    
    // % updating Q1
	if( mode == 0 )
    {
		// Kalman.Q1 = Kalman.Q1;	% no update
        memcpy( outKalman->Q1, inKalman->Q1, sizeof(double)*L*L );
    }
	else
    {
		// % corresponds to aMode==17
		// K = 0.5 * (K+K');
		// Kalman.Q1 = UC * Kalman.Kp;
        
		for( j=0; j<L; j++ )
			for( i=0; i<=j; i++ )
			{
				size_t ij = i * L + j;
				size_t ji = j * L + i;
				kh->K[ij] = 0.5 * ( kh->K[ij] + kh->K[ji] );
				kh->K[ji] = kh->K[ij];
			}
        
		for( n=0; n<L*L; n++ )
			outKalman->Q1[n] = UC * inKalman->Kp[n];
    }    
    
    // % a-priori state error covariance matrix
    // % for the next time step
    // Kalman.Kp=K+Kalman.Q1;
    tmp1 = outKalman->Kp;
    tmp2 = kh->K;
    tmp3 = outKalman->Q1;
    for( n=0; n<L*L; n++ )
        *(tmp1++) = *(tmp2++) + *(tmp3++);
    
    // % current estimation of state x
	// Kalman.x = Kalman.x + Kalman.G*err';
    tmp1 = outKalman->x;
    tmp2 = inKalman->x;
    for( n=0; n<L; n++ )
        *(tmp1++) = *(tmp2++);
    dgemv( chn, &L, &M, &one, outKalman->G, &L, err,
    	&inc, &one, outKalman->x, &inc );
    
    // % add new observation, drop oldest
    // Kalman.yr = [y(:)' Kalman.yr(1:M*(p-1))];
    tmp1 = outKalman->yr;
	const double *tmpc = y;
    for( n=0; n<M; n++ )
        *(tmp1++) = *(tmpc++);
    tmp2 = inKalman->yr;
    for( n=0; n<M*(p-1); n++ )
        *(tmp1++) = *(tmp2++);
}


// Implementation of the mvaar algorithm, that uses only one state variable
void kalman_maar( double *err, const double *y,
	//double UC, int M, int L, int p,
	double UC, ptrdiff_t M, ptrdiff_t L, ptrdiff_t p,
	struct Kalman *state,
	struct Kalman_Helper *kh, int mode )
{
	double traceK;
	double *tmp1, *tmp2, *tmp3;
	ptrdiff_t n,i,j, info, index;
    
	// % update Measurement Matrix
	// Kalman.H = kron( eye(M), Kalman.yr );
	// this loop should be faster than calculating the kronecker tensor product.
	index = 0;
	for( j=0; j<M; j++ ) {
		for( i=0; i<M*p; i++ ) {
			state->H[index] = state->yr[i];
			index += M;
		}
		index += 1;
	}
	// % calculate prediction error
	// ye = (Kalman.H*Kalman.x)';
	// err = y - ye
	index = 0;
	for( j=0; j<M; j++ ) {
		kh->ye[j] = 0;
		for( i=0; i<M*p; i++ )
			kh->ye[j] += state->yr[i] * state->x[index++];
		err[j] = y[j] - kh->ye[j];
	}
    
	// % update of Q2 using the prediction error 
	// of the previous step
	// Kalman.Q2 = (1-UC)*Kalman.Q2 + UC*err'*err;
	index = 0;
	for( i=0; i<M; i++ )
        	for( j=0; j<M; j++ ) {
				state->Q2[index] = (1-UC)*state->Q2[index] + UC * err[i] * err[j];
				index++;
        	}
    
    // KpH = Kalman.Kp * Kalman.H';
    dgemm( chn, cht, &L, &M, &L, &one, state->Kp, &L,
    	state->H, &M, &zero, kh->KpH, &L );
    
    // HKp = Kalman.H * Kalman.Kp;
    dgemm( chn, chn, &M, &L, &L, &one, state->H, &M,
    	state->Kp, &L, &zero, kh->HKp, &M );
    
    // HKpH = Kalman.H * KpH
    dgemm( chn, chn, &M, &M, &L, &one, state->H, &M, 
    	kh->KpH, &L, &zero, kh->HKpH, &M );
					
    
    // % Kalman Gain
    // Kalman.G = KpH * inv( HKpH + Kalman.Q2 );        
    tmp1 = kh->HKpH;
    tmp2 = state->Q2;
    for( n=0; n<M*M; n++ )
        *(tmp1++) += *(tmp2++);
    memset( kh->tmp, 0, M*M*sizeof( double ) );
    deye( M, kh->tmp, 1.0 );
    dgesv( &M, &M, kh->HKpH, &M, kh->ipiv,
    	kh->tmp, &M, &info );
    dgemm( chn, chn, &L, &M, &M, &one, kh->KpH, &L,
    	kh->tmp, &M, &zero, state->G, &L );
    
    // % calculation of the a-posteriori state 
    // % error covariance matrix
    // K = Kalman.Kp-Kalman.G*HKp;
    tmp1 = kh->K;
    tmp2 = state->Kp;
    for( n=0; n<L*L; n++ )
        *(tmp1++) = *(tmp2++);
    dgemm( chn, chn, &L, &L, &M, &minusone, state->G,
    	&L, kh->HKp, &M, &one, kh->K, &L );
		
		
	// % updating Q1
	tmp1 = state->Q1;
	if( mode == 0 )
	{
		// no update
	}
	else if( mode == 1 )
	{
		// Kalman.Q1=diag(diag(K)).*UC;
		tmp2 = kh->K;
		for( n=0; n<L*L; n++ )
		{
			if( n%(L+1)==0 )
			{
				*tmp1 = *tmp2 * UC;
				tmp2 += L+1;
			}
			else
				*tmp1 = 0;
			tmp1++;
		}
	}
	else if( mode == 2 )
	{
			// %diagonal matrix containing UC
			// upd = eye(L)/L*UC;
            // Kalman.Q1=upd*trace(K);
			tmp2 = kh->K;
			traceK = *tmp2;
			for( n=1; n<L; n++ )
			{
				tmp2 += L+1;
				traceK += *tmp2;
			}
			
			for( n=0; n<L*L; n++ )
			{
				if( n%(L+1)==0 )
					*tmp1 = traceK * UC / L;
				else
					*tmp1 = 0;
				tmp1++;
			}
			
	}
    
    
    // % a-priori state error covariance matrix
    // % for the next time step
    // Kalman.Kp=K+Kalman.Q1;
    tmp1 = state->Kp;
    tmp2 = kh->K;
    tmp3 = state->Q1;
    for( n=0; n<L*L; n++ )
        *(tmp1++) = *(tmp2++) + *(tmp3++);
    
    // % current estimation of state x
	// Kalman.x = Kalman.x + Kalman.G*err';
    dgemv( chn, &L, &M, &one, state->G, &L, err,
    	&inc, &one, state->x, &inc );
    
    // % add new observation, drop oldest
    // Kalman.yr = [y(:)' Kalman.yr(1:M*(p-1))];    
    
	tmp1 = state->yr + M*p-1;
    tmp2 = tmp1 - M;
    for( n=0; n<M*(p-1); n++ )
        *(tmp1--) = *(tmp2--);
    
    tmp1 = state->yr;
    const double *tmpc = y;
    for( n=0; n<M; n++ )
        *(tmp1++) = *(tmpc++);
}



// ============================
//     Implementation of helper functions
// ============================

mxArray *eye( int N, double value )
{
	mxArray *I;
	I = mxCreateDoubleMatrix( N, N, mxREAL );
	deye( N, mxGetPr( I ), value );
	return I;
}

void deye( int N, double *i, double value )
{
	int n;
	double *p = i;
	for( n=0; n<N; n++ )
	{
		*p = value;
		p = p + (N+1);
	}
}

