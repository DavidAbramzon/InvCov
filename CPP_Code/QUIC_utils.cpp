/*---------------------------------------------------------------------------
Copyright (2014): Eran Treister, Aviva Herman and Irad Yavneh. 
This code is distributed under the terms of the GNU General Public
License 2.0.

Permission to use, copy, modify, and distribute this software
for any purpose without fee is hereby granted, provided that
this entire notice is included in all copies of any software
which is or includes a copy or modification of this software
and in all copies of the supporting documentation for such
software. This software is being provided "as is", without any
express or implied warranty. In particular, the authors do not
make any representation or warranty of any kind concerning the
merchantability of this software or its fitness for any
particular purpose."
---------------------------------------------------------------------------*/

// This is the MEX wrapper for ML_QUIC.  The algorithm is in QUIC_utils.C.

// Invocation form within Matlab or Octave:
// [X W opt time iter] = ML_QUIC(mode, ...)
// [X W opt time iter] = ML_QUIC("default", S, L, tol, msg, maxIter,
//                            X0, W0)
// [X W opt time iter] = ML_QUIC("trace", S, L, tol, msg, maxIter,
//                                 X0, W0)
// See the README file and ML_QUIC.m for more information.

#include "QUIC_utils.h"
#include <stdio.h>
#include <cstring>

											
static inline int calcCholeskyU(double* ansU, uint32_t n)
{
	// this function gets the matrix ansU and perform the Cholesky factorization on it with LAPACK.
	// It runs over ansU. ansU will be upper triangular.
	// returns 1 if succeeds. return 0 if matrix is indefinite.
	ptrdiff_t info = 0;
    ptrdiff_t p0 = n;
    dpotrf_((char*) "U", &p0, ansU, &p0, &info);

    if (info != 0) {
		return 0; // cholesky factorization failed...indefinite matrix.
	}
	return 1;
}


static inline void calcInversionGivenCholeskyFactorU(double* ansU, uint32_t n)
{
	// this function gets the matrix ansU and perform the Cholesky factorization on it with LAPACK.
	// It runs over ansU. ansU will be upper triangular.
	// returns 1 if succeeds. return 0 if matrix is indefinite.
	ptrdiff_t info;
    ptrdiff_t p0 = n;
    dpotri_((char*) "U", &p0, ansU, &p0, &info);
}


static inline void symmetrizeMatrix(double* X, uint32_t n)
{
	for (uint32_t i = 0, k = 0; i < n; i++, k += n) // this is a trick so that k = i*n
		for (uint32_t j = i + 1; j < n; j++)
			X[k + j] = X[j * n + i];
}


static inline double Shrinkage(double q,double val)
{
	double epsilon = 0.0;
	if (val > 0) {
		return (val - q > epsilon) * (val - q);
	} else {
		return (val + q < -epsilon) * (val + q);
	}
}


static inline void CoordinateDescentUpdate(fX_Info& fx_info, Relax_Info& relax_info, uint32_t idx_in_active, double& normD, double& diffD)
{
	uint32_t p			 = fx_info.p;
	double* U 			 = fx_info.X; // We just use the allocation for X as U so we don't need to allocate U. 
	double* W 			 = fx_info.W;
	const double* S 	 = fx_info.S;
	int diagW = 0;

	double *D = relax_info.Dsmall;
	
	uint32_t i = relax_info.activeSet[idx_in_active].i;
	uint32_t j = relax_info.activeSet[idx_in_active].j;
    uint32_t ip = i * p;
    uint32_t jp = j * p;
    uint32_t ij = ip + j;
	
	if (relax_info.supp_size == p) { // diag_newton
		diagW = 1;
	}
	
    
    double a = W[ij] * W[ij];
    if (i != j)
        a += W[ip + i] * W[jp + j];
    double ainv = 1.0 / a;  // multiplication is cheaper than division
    
    double b = S[ij] - W[ij];
	
	if (diagW == 0) {
		for (uint32_t k = 0; k < p ; k++)
			b += W[ip + k] * U[k * p + j];
	} else {
		b +=  W[ip + i] * U[ip + j];
	}
	
    
    double l = fx_info.getLambda(ij) * ainv;
	
	double c = relax_info.Xsmall[idx_in_active] + D[idx_in_active];
    double f = b * ainv;
    double mu;
    normD -= fabs(D[idx_in_active]);
	
	mu = -c + Shrinkage(l, c - f);
	D[idx_in_active] += mu;
 
    diffD += fabs(mu);
    normD += fabs(D[idx_in_active]);
    if (mu != 0.0) {
		if (diagW == 0) {
			for (uint32_t k = 0; k < p; k++)
				U[ip + k] += mu * W[jp + k];
			if (i != j) {
				for (uint32_t k = 0; k < p; k++)
					U[jp + k] += mu * W[ip + k];
			}
		} else {
			U[ip + j] += mu * W[jp + j];
			if (i != j) {
				U[jp + i] += mu * W[ip + i];
			}
		
		}
    }
}


double my_linesearch(fX_Info& fx_info, Relax_Info& relax_info, Trace_Info& trace, double normD)
{
    double& l1normX 	 = fx_info.l1normX;
	double& trSX 		 = fx_info.trSX;
	double& logdetX 	 = fx_info.logdetX;
	double& fX 			 = fx_info.fX;
	double* Xsmall 		 = relax_info.Xsmall;
	double* W 			 = fx_info.W;
	const double* S 	 = fx_info.S;
	
	uint32_t msg = trace.msg;
	
	double* Dsmall = relax_info.Dsmall;
	uint_pair_t* activeSet = relax_info.activeSet;
	
	// This procedure runs over W, and returns the cholesky factor of X+alpha*D.
	uint32_t p = fx_info.p;
	
	unsigned int max_lineiter = 100;
	
    double alpha = 1;
	double beta = 0.5;
    double l1normXD = 0.0;
	double fX1;
    double fX1prev;
	unsigned int lineiter;
	int stop = 0;
	
    for (lineiter = 0; lineiter < max_lineiter; lineiter++) {
        
		double l1normX1 = 0.0;
        double trSX1 = 0.0;
		double tmp;
		
		memset(W, 0, p * p * sizeof(double));
		for (uint32_t i = 0; i < relax_info.numActive; i++ ) {
			uint32_t ij = activeSet[i].i * p + activeSet[i].j;
			tmp = Xsmall[i] + Dsmall[i] * alpha;
			W[ij] = tmp;
			tmp *= (1 + (activeSet[i].i != activeSet[i].j));
			l1normX1 += fabs(tmp) * fx_info.getLambda(ij);
			trSX1 	 += tmp * S[ij];
		}
        // Note that the upper triangular part is the lower
        // triangular part for the C arrays.
		
		// We are using W to hold the Cholesky decomposition; W
        // will hold the inverse later on.
		if (calcCholeskyU(W, p)==0) {
			if (msg >= QUIC_MSG_LINE)
                printf("    Line search step size %e.  Lack of positive definiteness.\n", alpha);
            alpha *= beta;
            continue;
		}
		
        double logdetX1 = 0.0;
		/*The following loop is going on the diagonal*/
        for (unsigned int i = 0, k = 0; i < p; i++, k += (p + 1))
            logdetX1 += log(W[k]);
        logdetX1 *= 2.0;
        fX1 = (trSX1 + l1normX1) - logdetX1;
		if (lineiter==0){
			fX1prev = fX1 + 1;
		}
		//////////////////////////////////////////////// ERAN
		if (stop){
            fX = fX1;
            l1normX = l1normX1;
            logdetX = logdetX1;
            trSX = trSX1;
			break;
		}
		if (fX1prev <= fX1 ){
			alpha /= beta;
			stop = 1;
			continue;
		}
        fX1prev = fX1;
        alpha *= beta;
    }
	if (msg >= QUIC_MSG_LINE)
		printf("    Line search ended in %d alpha decreases, with alpha = %lf\n", lineiter - 1, alpha);
	return alpha;
}



static double linesearch(fX_Info& fx_info, Relax_Info& relax_info, Trace_Info& trace, double normD)
{
    double& l1normX 	 = fx_info.l1normX;
	double& trSX 		 = fx_info.trSX;
	double& logdetX 	 = fx_info.logdetX;
	double& fX 			 = fx_info.fX;
	double* Xsmall 		 = relax_info.Xsmall;
	double* W 			 = fx_info.W;
	const double* S 	 = fx_info.S;
	
	uint32_t msg = trace.msg;
	
	double* Dsmall = relax_info.Dsmall;
	uint_pair_t* activeSet = relax_info.activeSet;
	
	// This procedure runs over W, and returns the cholesky factor of X+alpha*D.
	uint32_t p = fx_info.p;
	uint32_t max_lineiter = 20;
	double sigma = 0.001;
	
	double trgradgD = 0.0;
	
	for (uint32_t i = 0; i < relax_info.numActive; i++ ) {		
		uint32_t ij = activeSet[i].i * p + activeSet[i].j;
		trgradgD += (S[ij] - W[ij]) * Dsmall[i] * (1 + (activeSet[i].i != activeSet[i].j));
	}
	
    double alpha 	= 1.0;
	double beta 	= 0.5;
    double l1normXD = 0.0;
    double fX1prev 	= fX + 1;
	double fX1;
	uint32_t lineiter;
	int stop 		= 0;
	
    for (lineiter = 0; lineiter < max_lineiter; lineiter++) {
		double l1normX1 = 0.0;
        double trSX1 	= 0.0;
		double tmp;
		memset(W, 0, p * p * sizeof(double));
		for (uint32_t i = 0; i < relax_info.numActive; i++ ) {
			uint32_t ij = activeSet[i].i * p + activeSet[i].j;
			tmp = Xsmall[i] + Dsmall[i] * alpha;
			W[ij] = tmp;
			tmp *= (1 + (activeSet[i].i != activeSet[i].j));
			l1normX1 += fabs(tmp) * fx_info.getLambda(ij);
			trSX1 	 += tmp * S[ij];
		}
		
        // Note that the upper triangular part is the lower
        // triangular part for the C arrays.
		
		// We are using W to hold the Cholesky decomposition; W
        // will hold the inverse later on.
		if (calcCholeskyU(W, p) == 0) {
			if (msg >= QUIC_MSG_LINE)
                printf("    Line search step size %e.  Lack of positive definiteness.\n", alpha);
            alpha *= beta;
            continue;
		}
		
        double logdetX1 = 0.0;
        for (uint32_t i = 0, k = 0; i < p; i++, k += (p + 1))
            logdetX1 += log(W[k]);
        logdetX1 *= 2.0;
        fX1 = (trSX1 + l1normX1) - logdetX1;		
		
		if (alpha == 1.0)
            l1normXD = l1normX1;
			
        if (fX1 <= fX + alpha * sigma * (trgradgD + l1normXD - l1normX) || normD == 0) {
            if (msg >= QUIC_MSG_LINE)
                printf("    Line search step size chosen: %e.\n", alpha);
            fX 		= fX1;
            l1normX = l1normX1;
            logdetX = logdetX1;
            trSX 	= trSX1;
            break;
        }
		
        if (msg >= QUIC_MSG_LINE) {
            printf("    Line search step size %e.\n", alpha);
            printf("      Objective value would not decrease sufficiently: "
                    "%e.\n", fX1 - fX);
        }
		
        if (fX1prev < fX1) {
            l1normX = l1normX1;
            logdetX = logdetX1;
            trSX 	= trSX1;
            break;
        }
		
        fX1prev = fX1;
        alpha  *= beta;
    }
	
	if (msg >= QUIC_MSG_LINE)
		printf("    Line search ended in %d iterations \n", lineiter + 1);

	return alpha;
}


static uint32_t computeActiveSetAndSubgrad_debias(fX_Info& fx_info, Relax_Info& relax_info, Trace_Info& trace, double& subgrad)
{
	uint32_t p 			 = fx_info.p;
	double* X  			 = fx_info.X;
	double* W  			 = fx_info.W;
	const double* S  	 = fx_info.S;
	
	uint32_t numActive = 0;
	subgrad = 0.0;
	relax_info.supp_size = 0;
	for (uint32_t k = 0, i = 0; i < p; i++, k += p) { // this is a trick so that k = i*p
		for (uint32_t j = 0; j <= i; j++) {
			double g = S[k + j] - W[k + j];
			
			// here we save the value of |g| into G_k
			int nz = fabs(X[k + j]) > 1e-12; 
			if (nz) {
				relax_info.activeSetGrad[numActive] = fabs(g) + BIG_NUM * (nz);
				
				relax_info.activeSet[numActive].i = i;
				relax_info.activeSet[numActive].j = j;
				
				numActive++;
				
				if (X[k + j] > 0)
					g += fx_info.getLambda(k + j);
				else if (X[k + j] < 0)
					g -= fx_info.getLambda(k + j);
				else
					g = fabs(g) - fx_info.getLambda(k + j);
				
				subgrad += fabs(g);
				relax_info.supp_size += (nz);
			}else{
				X[k + j] = 0.0;
			}			
		}
	}
	if (trace.msg >= QUIC_MSG_NEWTON)
		printf("Supp size: %d\n",relax_info.supp_size);

	return numActive;
}



static uint32_t computeActiveSetAndSubgrad(fX_Info& fx_info, Relax_Info& relax_info, Trace_Info& trace, double& subgrad)
{
	uint32_t p 			 = fx_info.p;
	double* X  			 = fx_info.X;
	double* W  			 = fx_info.W;
	const double* S  	 = fx_info.S;
	
	uint32_t numActive = 0;
	subgrad = 0.0;
	relax_info.supp_size = 0;
	for (uint32_t k = 0, i = 0; i < p; i++, k += p) { // this is a trick so that k = i*p
		for (uint32_t j = 0; j <= i; j++) {
			double g = S[k + j] - W[k + j];
			
			// here we save the value of |g| into G_k
			int nz = fabs(X[k + j]) > 1e-12; 
			if ( nz || (fabs(g) > fx_info.getLambda(k + j))) {
				relax_info.activeSetGrad[numActive] = fabs(g) + BIG_NUM * (nz);
				
				relax_info.activeSet[numActive].i = i;
				relax_info.activeSet[numActive].j = j;
				
				numActive++;
				
				if (X[k + j] > 0)
					g += fx_info.getLambda(k + j);
				else if (X[k + j] < 0)
					g -= fx_info.getLambda(k + j);
				else
					g = fabs(g) - fx_info.getLambda(k + j);
				
				subgrad += fabs(g);
				relax_info.supp_size += (nz);
			}else{
				X[k + j] = 0.0;
			}			
		}
	}
	if (trace.msg >= QUIC_MSG_NEWTON)
		printf("Supp size: %d\n",relax_info.supp_size);

	return numActive;
}


static uint32_t computeActiveSetFromC(fX_Info& fx_info, Relax_Info& relax_info, uint32_t level_size)
{
	uint32_t p 			 = fx_info.p;
	double* X  			 = fx_info.X;
	double* W  			 = fx_info.W;
	const double* S  	 = fx_info.S;
	
	uint32_t numActive = 0;
	relax_info.supp_size = 0;

	for (uint32_t i = 0; i < level_size; i++) { 
		uint32_t ij = relax_info.C[i].i * p + relax_info.C[i].j;
		double g = fabs(S[ij] - W[ij]);
		
		// here we save the value of |g| into G_k
		if (X[ij] != 0.0 || g > fx_info.getLambda(ij)) {
			relax_info.activeSet[numActive].i = relax_info.C[i].i;
			relax_info.activeSet[numActive].j = relax_info.C[i].j;
			numActive++;
			relax_info.supp_size += (X[ij] != 0.0);
		}
	}
	return numActive;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
//sort
///////////////////////////////////////////////////////////////////////////////////////////////////////
static void q_sort_descend(double* numbers, uint_pair_t* indices, size_t left, size_t right)
{
	uint_pair_t pivot_idx;
	size_t l_hold, r_hold, prev_left;
	double pivot;

	l_hold = left;
	r_hold = right;
	pivot = numbers[left];
	pivot_idx.i = indices[left].i;
	pivot_idx.j = indices[left].j;

	while (left < right) {
		while ((numbers[right] <= pivot) && (left < right))
		  right--;

		if (left != right) {
		  numbers[left] = numbers[right];
		  indices[left].i = indices[right].i;
		  indices[left].j = indices[right].j;
		  left++;
		}

		while ((numbers[left] >= pivot) && (left < right))
			left++;

		if (left != right) {
		  numbers[right] = numbers[left];
		  indices[right].i = indices[left].i;
		  indices[right].j = indices[left].j;	  
		  right--;
		}
	}

	numbers[left] = pivot;
	indices[left].i = pivot_idx.i;
	indices[left].j = pivot_idx.j;
	prev_left = left;
	left = l_hold;
	right = r_hold;

	if (left < prev_left)
		q_sort_descend(numbers, indices, left, prev_left - 1);

	if (right > prev_left)
		q_sort_descend(numbers, indices, prev_left + 1, right);
}

static void quickSort_descend(double* numbers, uint_pair_t* indices, size_t array_size)
{
    q_sort_descend(numbers, indices, 0, array_size - 1);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////


static void calc_TraceSX_and_Lambdal1normX(fX_Info& fx_info, Relax_Info& relax_info)
{
	uint32_t p 			 = fx_info.p;
	const double* S  	 = fx_info.S;
	double& l1normX 	 = fx_info.l1normX;
	double& trSX 		 = fx_info.trSX;
	
	l1normX = 0.0;
	trSX 	= 0.0;
	
	for (uint32_t i = 0; i < relax_info.numActive ; i++) {
		uint32_t ij = relax_info.activeSet[i].i * p + relax_info.activeSet[i].j;
		l1normX += fx_info.getLambda(ij) * fabs(relax_info.Xsmall[i]);
		trSX 	+= relax_info.Xsmall[i] * S[ij];
	}

	return;
}



static double ApplyCoordinateDescentIterationsForLASSO(fX_Info& fx_info, Relax_Info& relax_info, Trace_Info& trace)
{
	uint32_t p = fx_info.p;
	
	uint_pair_t* activeSet = relax_info.activeSet;
	
	uint32_t msg = trace.msg;
	
	double diffD, normD = 0.0, tmp;
	double time;
	uint32_t cdSweep = 1;
	srand(1); // this is to set the seed of the random permutation - so we will have a "deterministic code"
	
	memset(fx_info.X, 0, p * p * sizeof(double));
	diffD = 0.0;
	time  = clock();
	if (relax_info.supp_size == p) {
		for (uint32_t l = 0; l < relax_info.numActive; l++) {
				// Here we have 3*p cost.
			CoordinateDescentUpdate(fx_info, relax_info, l, normD, diffD);
		}
	} else {	
		for (cdSweep = 1; cdSweep <= 1 + relax_info.maxNewtonIter / 3; cdSweep++) { // (eran) number of itertions is just a heuristic
			for (uint32_t i = 0; i < relax_info.numActive; i++ ) {
				uint32_t j = i + rand() % (relax_info.numActive - i); 
				
				uint32_t k1 = activeSet[i].i;
				uint32_t k2 = activeSet[i].j;
				activeSet[i].i = activeSet[j].i;
				activeSet[i].j = activeSet[j].j;
				activeSet[j].i = k1;
				activeSet[j].j = k2;
				
				diffD = relax_info.Dsmall[i];
				relax_info.Dsmall[i] = relax_info.Dsmall[j];
				relax_info.Dsmall[j] = diffD;
				
				tmp = relax_info.Xsmall[i];
				relax_info.Xsmall[i] = relax_info.Xsmall[j];
				relax_info.Xsmall[j] = tmp;
			}
			
			for (uint32_t l = 0; l < relax_info.numActive; l++) {
				// Here we have 3*p cost.
				CoordinateDescentUpdate(fx_info, relax_info, l, normD, diffD);
			}
			
			
			if (diffD <= normD * relax_info.cdSweepTol) {
				if (msg >= QUIC_MSG_CD)
					printf("  Coordinate descent stopped through tol\n");
				break;
			}
		}
	}
	if (msg >= QUIC_MSG_CD) {
		printf("CPU time of CD update %d: %.3f seconds\n", cdSweep, (clock() - time) / CLOCKS_PER_SEC);
		printf("  Coordinate descent sweep %ld. norm of D = %e, change in D = %e.\n", cdSweep, normD, diffD);
	}
	return normD;
}

static double Relaxation(fX_Info& fx_info, Relax_Info& relax_info, Trace_Info& trace)
{
	////Applying CD on active set.
	uint32_t p = fx_info.p;
	
	double time; 
	double normD = 0.0;
	memset(relax_info.Dsmall, 0, relax_info.numActive * sizeof(double));
	
	for (uint32_t i = 0 ; i < relax_info.numActive; ++i) {
		relax_info.Xsmall[i] = fx_info.X[relax_info.activeSet[i].i * p + relax_info.activeSet[i].j];
	}
	
	normD = ApplyCoordinateDescentIterationsForLASSO(fx_info, relax_info, trace);		

	time = clock();

	double alpha = linesearch(fx_info, relax_info, trace, normD);
 	
	if (trace.msg >= QUIC_MSG_LINE)
		printf("CPU time of line_search procedure: %.3f seconds\n", (clock() - time) / CLOCKS_PER_SEC);
	
    // compute W = inv(X): W will hold the inverse from now on.
    time = clock();
	calcInversionGivenCholeskyFactorU(fx_info.W, p); // the chol in W was done part of linesearch.
	symmetrizeMatrix(fx_info.W, p);
	if (trace.msg >= QUIC_MSG_INV)
		printf("CPU time of inversion: %.3f seconds\n", (clock() - time) / CLOCKS_PER_SEC);
    
	
	memset(fx_info.X, 0, p * p * sizeof(double));
	for (uint32_t i = 0; i < relax_info.numActive; i++){
		relax_info.Xsmall[i] += alpha * relax_info.Dsmall[i];
		fx_info.X[relax_info.activeSet[i].i * p + relax_info.activeSet[i].j] = relax_info.Xsmall[i];
	}

	if (trace.mode == 'T') {
		if (trace.suppArr)
			trace.suppArr[trace.traceIdx] = relax_info.supp_size;
		if (trace.activeNumArr)
			trace.activeNumArr[trace.traceIdx] = relax_info.numActive;
		if (trace.opt)
			trace.opt[trace.traceIdx] = fx_info.fX;
		if (trace.cputime)
			trace.cputime[trace.traceIdx] = (clock() - trace.timeBegin) / CLOCKS_PER_SEC;
		trace.traceIdx++;
	}
	return alpha;
}


static double Vcycle(fX_Info& fx_info, Relax_Info& relax_info, Trace_Info& trace)
{
	uint32_t p = fx_info.p;	
	double alpha;
	
	uint32_t *level_sizes = (uint32_t *) malloc(VSYCLE_LEVELS * sizeof(uint32_t));
	
	if (!level_sizes) {
		if (trace.msg >= QUIC_MSG_FAILURE)
			printf("Failure: malloc level_sizes\n");
		return -1;
	}
	
	uint32_t minCoarseVars = p;
	
	quickSort_descend(relax_info.activeSetGrad, relax_info.activeSet, relax_info.numActive);
	
	for (int k = 0; k < relax_info.numActive ; k++) {
		relax_info.C[k].i = relax_info.activeSet[k].i;
		relax_info.C[k].j = relax_info.activeSet[k].j;
	}
	
	
	int lev = 0;
	level_sizes[lev] = relax_info.numActive;
    uint32_t wanted_coarse = relax_info.numActive / 2;
	
	if (wanted_coarse > relax_info.supp_size) {
		do {
			if (wanted_coarse <= relax_info.supp_size) {
				level_sizes[++lev] = relax_info.supp_size;
				break;
			} else {
				level_sizes[++lev] = wanted_coarse;
			}
			wanted_coarse = std::max(wanted_coarse/2, relax_info.supp_size);
		} while ((wanted_coarse > minCoarseVars) && (relax_info.supp_size <= wanted_coarse));
		
		//solution for the lowest level of the vcycle
		if (trace.msg >= QUIC_MSG_VCYCLE)
			printf("$$$$$$$$ LOWEST LEVEL %d for size: %d out of %d, supp: %d $$$$$$$\n", lev, level_sizes[lev], level_sizes[0], relax_info.supp_size);
		relax_info.numActive = level_sizes[lev];
		alpha = Relaxation(fx_info, relax_info, trace);
		//alpha = Relaxation(fx_info, relax_info, trace);
		
		//the way up:
		for (--lev; lev >= 1; --lev) {
			relax_info.numActive = computeActiveSetFromC(fx_info, relax_info, level_sizes[lev]);
			if (trace.msg >= QUIC_MSG_VCYCLE)
				printf("$$$$$$$$ LEVEL %d for size: %d out of %d, supp: %d $$$$$$$\n", lev, relax_info.numActive, level_sizes[0], relax_info.supp_size);
			alpha = Relaxation(fx_info, relax_info, trace);
		}
	} else {
		if (trace.msg >= QUIC_MSG_VCYCLE)
			printf("$$$$$$$$ NO LEVELS: for size: %d out of %d, supp: %d $$$$$$$\n", relax_info.numActive, level_sizes[0], relax_info.supp_size);
		relax_info.numActive = level_sizes[lev];
		alpha = Relaxation(fx_info, relax_info, trace);
	}

	free(level_sizes);
	return alpha;
}


static Result InitRelaxInfo(Relax_Info& relax_info, Trace_Info& trace, double cdSweepTol, uint32_t p, uint32_t activeSetSize)
{
	uint32_t msg = trace.msg;
	
	relax_info.Dsmall = (double *) malloc(activeSetSize * sizeof(double));
	
	if (!relax_info.Dsmall) {
		if (msg >= QUIC_MSG_FAILURE)
			printf("Failure: malloc relax_info.Dsmall\n");
		goto ErrDsmall;
	}
	
	relax_info.Xsmall = (double *) malloc(activeSetSize * sizeof(double));
	
	if (!relax_info.Xsmall) {
		if (msg >= QUIC_MSG_FAILURE)
			printf("Failure: malloc relax_info.Xsmall\n");
		goto ErrXsmall;
	}
	
	relax_info.C = (uint_pair_t *) malloc(activeSetSize * sizeof(uint_pair_t));
	
	if (!relax_info.C) {
		if (msg >= QUIC_MSG_FAILURE)
			printf("Failure: malloc relax_info.C\n");
		goto ErrC;
	}
	
	relax_info.activeSet = (uint_pair_t *) malloc((p * (p + 1) / 2) * sizeof(uint_pair_t));
	
	if (!relax_info.activeSet) {
		if (msg >= QUIC_MSG_FAILURE)
			printf("Failure: malloc relax_info.activeSet\n");
		goto ErrActiveSet;
	}
	
	relax_info.activeSetGrad = (double *) malloc((p * (p + 1) / 2) * sizeof(double));
	
	if (!relax_info.activeSetGrad) {
		if (msg >= QUIC_MSG_FAILURE)
			printf("Failure: malloc relax_info.activeSetGrad\n");
		goto ErrActiveSetGrad;
	}
	
	relax_info.supp_size = 0; 
	relax_info.cdSweepTol = cdSweepTol;
	
	return SUCCESS;
	
ErrActiveSetGrad:
	free(relax_info.activeSet);
ErrActiveSet:
	free(relax_info.C);
ErrC:
	free(relax_info.Xsmall);
ErrXsmall:
	free(relax_info.Dsmall);
ErrDsmall:
	return FAILURE;
}


static void DestroyRelaxInfo(Relax_Info& relax_info)
{
	free(relax_info.activeSetGrad);
	free(relax_info.activeSet);
	free(relax_info.C);
	free(relax_info.Xsmall);
	free(relax_info.Dsmall);
}


static Result resizeRelaxInfo(Relax_Info& relax_info, uint32_t& currentAllocatedMemoryForActive, uint32_t msg)
{
	currentAllocatedMemoryForActive = relax_info.numActive + 50;
	
	free(relax_info.Dsmall);
	
	relax_info.Dsmall = (double *) malloc(currentAllocatedMemoryForActive * sizeof(double));

	if (!relax_info.Dsmall) {
		if (msg >= QUIC_MSG_FAILURE)
			printf("Failure: malloc relax_info.Dsmall while resizing\n");
		goto out;
	}
	
	free(relax_info.Xsmall);
	
	relax_info.Xsmall = (double *) malloc(currentAllocatedMemoryForActive * sizeof(double));
	
	if (!relax_info.Xsmall) {
		if (msg >= QUIC_MSG_FAILURE)
			printf("Failure: malloc relax_info.Xsmall while resizing\n");
		goto out;
	}
	
	free(relax_info.C);
	
	relax_info.C = (uint_pair_t *) malloc(currentAllocatedMemoryForActive * sizeof(uint_pair_t));
	
	if (!relax_info.C) {
		if (msg >= QUIC_MSG_FAILURE)
			printf("Failure: malloc relax_info.C while resizing\n");
		goto out;
	}

	return SUCCESS;

out:
	return FAILURE;
}


void QUIC(fX_Info& fx_info, Trace_Info& trace, uint32_t maxIter, QUIC_method method)
{
#ifdef GDEBUG
    startgdb();
#endif
	uint32_t msg = trace.msg;
	uint32_t need_to_calc = 1;
	if (trace.mode >= 'a')
		trace.mode -= ('a' - 'A');

	if (msg >= QUIC_MSG_MIN) {
		if (trace.mode == 'T')
			printf("Running %s version %s in 'trace' mode.\n", METHOD(method), VERSION);
		else
			printf("Running %s version %s in 'default' mode.\n", METHOD(method), VERSION);
	}
	
	uint32_t p = fx_info.p;
	
	uint32_t currentAllocatedMemoryForActive = std::min(10 * p, (p * (p + 1) / 2));
	
	Relax_Info relax_info;	
	Result res = InitRelaxInfo(relax_info, trace, 0.05, p, currentAllocatedMemoryForActive);
	
	if (res != SUCCESS) {
		if (msg >= QUIC_MSG_FAILURE)
			printf("Failure: InitRelaxInfo\n");
		return;
	}
	
	trace.timeBegin = clock();
	double& timeBegin = trace.timeBegin;
	
	double fXprev = INITIAL_FX;

	double time = clock();
	unsigned int NewtonIter;
	
	uint32_t currIdx = 0, prevTraceIdx;
	
	for (NewtonIter = 1; NewtonIter <= maxIter; NewtonIter++) {
		if (msg >= QUIC_MSG_NEWTON)
			printf("Newton iteration %ld.**********************************************\n", NewtonIter);

		fflush(stdout);	
		
		double subgrad = INITIAL_FX;
		relax_info.maxNewtonIter = NewtonIter;
		
		// Compute the active set and the minimum norm subgradient (the initial active set is saved in C):
		if (method == QUIC_DEBIAS) {
			relax_info.numActive = computeActiveSetAndSubgrad_debias(fx_info, relax_info, trace, subgrad);
		} else {
			relax_info.numActive = computeActiveSetAndSubgrad(fx_info, relax_info, trace, subgrad);
		}
		if (relax_info.numActive > currentAllocatedMemoryForActive) {
			// this code is in order to save space in sparse matrices (size of active set): no need for active set to be of size p*(p-1)/2.

			if (msg >= QUIC_MSG_NEWTON)
				printf("Increasing Relax_Info!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
			
			
			Result res = resizeRelaxInfo(relax_info, currentAllocatedMemoryForActive, msg);
			
			if (res != SUCCESS) {
				if (msg >= QUIC_MSG_FAILURE)
					printf("Failure: resize relax_info\n");
				goto out;
			}
		}
		
		if (msg >= QUIC_MSG_NEWTON) {
			printf("  Active set size = %ld, Support: %d.\n", relax_info.numActive, relax_info.supp_size);
			printf("  sub-gradient = %e, l1-norm of X = %e.\n", subgrad, fx_info.l1normX);
		}
		
		
				
		if (need_to_calc) {
			for (uint32_t i = 0 ; i < relax_info.numActive; ++i) {
				relax_info.Xsmall[i] = fx_info.X[relax_info.activeSet[i].i * p + relax_info.activeSet[i].j];
			}
			calc_TraceSX_and_Lambdal1normX(fx_info, relax_info);
			need_to_calc = 0;
		}
		
		fXprev = fx_info.fX;
		
		double alpha;
		prevTraceIdx = trace.traceIdx;
		switch (method)
		{
			case ML_QUIC: 
				alpha = Vcycle(fx_info, relax_info, trace); 
				
				if (trace.mode == 'T' && trace.vcycleIter) {
					trace.vcycleIter[currIdx] = trace.traceIdx;
					currIdx++;
				}
				
				break;
			case QUIC_RELAX:
			case QUIC_DEBIAS:
				alpha = Relaxation(fx_info, relax_info, trace);
				break;
			default:
				if (msg >= QUIC_MSG_FAILURE)
					printf("Invalid method!\n");
				break;
		}
		
		if (msg >= QUIC_MSG_NEWTON)
			printf("  Objective value decreased by %e.\n", fXprev - fx_info.fX);
			
		
		// Check for convergence.
		if (subgrad * alpha < fx_info.l1normX * fx_info.tol || (fabs((fx_info.fX - fXprev) / fx_info.fX) < EPS)){
			NewtonIter++;
			break;
		}
	}
	
	if (trace.mode == 'D') {
		if (trace.suppArr)
			trace.suppArr[0] = relax_info.supp_size;
		if (trace.activeNumArr)
			trace.activeNumArr[0] = relax_info.numActive;
		if (trace.opt)
			trace.opt[0] = fx_info.fX;
		if (trace.iter)
			trace.iter[0] = NewtonIter - 1;
		if (trace.vcycleIter)
			trace.vcycleIter[0] = trace.traceIdx - prevTraceIdx;
	}
	if (trace.mode == 'T' && trace.iter)
		trace.iter[0] = NewtonIter - 1;

	symmetrizeMatrix(fx_info.X, p);

	double elapsedTime = (clock() - timeBegin) / CLOCKS_PER_SEC;
	if (trace.mode == 'D' && trace.cputime)
		trace.cputime[0] = elapsedTime;

	if (msg >= QUIC_MSG_MIN)
		printf("QUIC CPU time: %.3f seconds\n", elapsedTime);
	
out:
	DestroyRelaxInfo(relax_info);
}
