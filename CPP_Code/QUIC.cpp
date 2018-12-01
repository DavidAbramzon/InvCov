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

// This is the MEX wrapper for QUIC.  The algorithm is in QUIC_utils.C.

// Invocation form within Matlab or Octave:
// [X W opt time iter] = QUIC(mode, ...)
// [X W opt time iter] = QUIC("default", S, L, tol, msg, maxIter,
//                            X0, W0)
// [X W opt time iter] = QUIC("trace", S, L, tol, msg, maxIter,
//                                 X0, W0)
// See the README file and QUIC.m for more information.

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "QUIC_utils.h"


int main(int nlhs, void *plhs[], int nrhs, const void * prhs[])
{
    if (nrhs < 2) {
        printf("QUIC:arguments\n"
                          "Missing arguments, please specify\n"
                                  "             S - the empirical covariance matrix, and\n"
                                  "             L - the regularization parameter.");
    }

    fX_Info *fx_info_ptr;
    fX_Info fx_info;
    fX_InfoFullLambda fx_info_lambda;

    Trace_Info trace;
    long argIdx = 0;

    char mode[8];
    mode[0] = 'd';
    char * inputMode = (char*)(&prhs[argIdx]);
    if (!strcmp(inputMode, "path") &&
        !strcmp(inputMode, "trace") &&
        !strcmp(inputMode, "default")) {
            printf("QUIC:arguments\n"
               "Invalid mode, use: 'default', 'path' or "
                       "'trace'.");
        return -1;
    }
    mode[1] = inputMode;
    argIdx++;

    // The empirical covariance matrix:
    const double* S = (double*)(&prhs[argIdx]);
    argIdx++;

    // Regularization parameter matrix:
    fx_info_ptr = &fx_info;
    fx_info_ptr->Lambda = (double *)(&prhs[argIdx]);

    argIdx++;

    double tol = 1e-6;
    if (nrhs > argIdx) {
        tol = (*(double *)(&prhs[argIdx]));
        if (tol < 0) {
            printf("QUIC:tol\n"
                              "Negative tolerance value.");
        }
        argIdx++;
    }

    uint32_t msg = QUIC_MSG_FAILURE;
    if (nrhs > argIdx) {
        msg = (uint32_t) *(int *)(&prhs[argIdx]);
        argIdx++;
    }

    // Maximum number of Newton ierations (whole matrix update):
    uint32_t maxIter = 1000;
    if (nrhs > argIdx) {
        maxIter = (uint32_t) *(int *)(&prhs[argIdx]);
        argIdx++;
    }

    double* X0 = NULL;
    double* W0 = NULL;
    if (nrhs > argIdx) {

        X0 = (double *)(&prhs[argIdx]);
        argIdx++;
        if (nrhs == argIdx)
            printf("QUIC:initializations\n"
                              "Please specify both the initial estimate\n"
                                      "             and the inverse.\n"
                                      "             Maybe incorrect mode is specified?");
        W0 = (double*)(&prhs[argIdx]);
        argIdx++;
    }

    double* X = NULL;
    double* W = NULL;
    uint32_t p1 = *(uint32_t*)(&prhs[argIdx]);
    argIdx++;
    uint32_t p2 = *(uint32_t*)(&prhs[argIdx]);
    int dims[] = {p1, p2};
    double tmp[p1][p2];
    X = (double *) (tmp);
    if (nlhs > 0)
        plhs[0] = tmp;

    double tmp2[p1][p2];
    W = (double *) (tmp2);
    if (nlhs > 1)
        plhs[1] = tmp2;

    if (X0 != NULL) {
        memcpy(X, X0, sizeof(double) * p1 * p2);
        memcpy(W, W0, sizeof(double) * p1 * p2);
    } else {
        memset(X, 0, sizeof(double) * p1 * p2);
        memset(W, 0, sizeof(double) * p1 * p2);
        for (unsigned long i = 0; i < p1 * p2; i += (p1 + 1)) {
            X[i] = 1.0;
            W[i] = 1.0;
        }
    }

    double* opt = NULL;
    double* cputime = NULL;
    uint32_t* iter = NULL;
    uint32_t* suppArr = NULL;
    uint32_t* activeNumArr = NULL;

    unsigned long traceSize = 1;
    unsigned long iterSize = 1;

    if (mode[0] == 't')
        traceSize = maxIter;

    if (nlhs > 2) {
        double tmp3[traceSize];
        plhs[2] = tmp3;
        opt = (double *)(&plhs[2]);
    }

    if (nlhs > 3) {
        double tmp3[traceSize];
        plhs[3] = tmp3;
        cputime = (double*)(&plhs[3]);
    }

    if (nlhs > 4) {
        int dims[] = {iterSize};
        double tmp3[iterSize];
        plhs[4] = tmp3;
        iter = (uint32_t *) (&plhs[4]);
    }

    if (nlhs > 5) {
        int dims[] = {traceSize};
        int tmp3[traceSize];
        plhs[5] = tmp3;
        suppArr = (uint32_t *) (&plhs[5]);
    }

    if (nlhs > 6) {
        int tmp3[traceSize];
        plhs[6] = tmp3;
        activeNumArr = (uint32_t *) (&plhs[6]);
    }

    fx_info_ptr->p 		 = p1*p2;
    fx_info_ptr->l1normX = 0.0;
    fx_info_ptr->trSX 	 = 0.0;
    fx_info_ptr->logdetX = 0.0;
    fx_info_ptr->fX 	 = INITIAL_FX;
    fx_info_ptr->X 		 = X;
    fx_info_ptr->S 		 = S;
    fx_info_ptr->W		 = W;
    fx_info_ptr->tol	 = tol;

    trace.opt 		   = opt;
    trace.iter 		   = iter;
    trace.vcycleIter   = NULL;
    trace.cputime 	   = cputime;
    trace.suppArr 	   = suppArr;
    trace.activeNumArr = activeNumArr;
    trace.traceIdx 	   = 0;
    trace.mode 		   = mode[0];
    trace.msg 		   = msg;

    QUIC(*fx_info_ptr, trace, maxIter, QUIC_RELAX);
    return 0;
}
