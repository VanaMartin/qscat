#ifndef __BLAS__
	#define __BLAS__
	#ifdef INTEL_MKL
		#include "intel.h"
	#else
		#include "openblas.h"
	#endif
#endif
