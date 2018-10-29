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
#ifndef QUIC_UTILS
#define QUIC_UTILS

#define VERSION "1.1"

#include <algorithm>    // std::max, std::min
#include <math.h>		// fabs
#include <time.h>		// clock, CLOCKS_PER_SEC
#include <iostream>

#ifdef GDEBUG
#include "startgdb.h"
#endif

typedef int int32_t;
typedef unsigned int uint32_t;

#define BIG_NUM 1e+4
#define MSG mexPrintf
#define INITIAL_FX 1e+15
#define VSYCLE_LEVELS 10
#define EPS (double(2.22E-16))
#define METHOD(m) ((m) ? "ML_QUIC" : "QUIC")

typedef enum {
	SUCCESS = 0,
	FAILURE,
} Result;

typedef enum {
	QUIC_MSG_FAILURE = 0,
	QUIC_NO_MSG 	 = 0,
	QUIC_MSG_MIN 	 = 1,
	QUIC_MSG_NEWTON  = 2,
	QUIC_MSG_RELAX   = 2,
	QUIC_MSG_VCYCLE  = 2,
	QUIC_MSG_CD 	 = 2,
	QUIC_MSG_LINE 	 = 2,
	QUIC_MSG_INV	 = 2,
} MsgUrgency;

typedef enum {
	QUIC_RELAX	 = 0,
	ML_QUIC 	 = 1,
	QUIC_DEBIAS = 2,
} QUIC_method;

typedef struct {
    unsigned int i;
    unsigned int j;
} uint_pair_t;


#if !defined(LANG_R) && defined(_WIN32)
#define dpotrf_ dpotrf
#define dpotri_ dpotri
#endif


extern "C" void dpotrf_(char* uplo, ptrdiff_t* n, double* A, ptrdiff_t* lda, ptrdiff_t* info);
extern "C" void dpotri_(char* uplo, ptrdiff_t* n, double* A, ptrdiff_t* lda, ptrdiff_t* info);


class fX_Info {
	public:
		uint32_t 		p;
		double 			*X;
		double 			*W; 
		const double 	*S; 
		const double 	*Lambda;
		double 			l1normX; 
		double 			logdetX; 
		double 			trSX; 
		double 			fX;
		double			tol;
		
		virtual double getLambda(uint32_t pos)
		{	
			return Lambda[0];
		}
};

class fX_InfoFullLambda : public fX_Info {
	public:
		double getLambda(uint32_t pos)
		{	
			return Lambda[pos];
		}
};

typedef struct {
	uint32_t 	maxNewtonIter;
	double  	*Dsmall;
	double		*Xsmall;
	uint32_t	numActive; 
	uint_pair_t *C;
	uint_pair_t *activeSet;
	double 		*activeSetGrad;
	uint32_t 	supp_size; 
	double 		cdSweepTol;
} Relax_Info;

typedef struct {
	double		*opt;
    double		*cputime;
    uint32_t	*iter;
	uint32_t	*vcycleIter;
	uint32_t	*suppArr;
	uint32_t	*activeNumArr;
	uint32_t 	traceIdx;
	double 		timeBegin;
	char 		mode;
	uint32_t  	msg;
} Trace_Info;


void QUIC(fX_Info& fx_info, Trace_Info& trace, uint32_t maxIter, QUIC_method method);


#endif /* QUIC_UTILS */

