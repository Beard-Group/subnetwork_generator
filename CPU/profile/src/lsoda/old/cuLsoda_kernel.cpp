/*
 *  cuLsoda_kernel.cu
 *  cuLsoda
 *
 */
#ifndef _CULSODA_CU_H_
#define _CULSODA_CU_H_
 
#include "cuLsoda.cu.h"
 
template<typename Fex, typename Jex>
void cuLsoda(Fex fex, int *neq, double *y, double *t, double *tout, int *itol, double *rtol, double *atol, int *itask, int *istate, int *iopt, double *rwork, int *lrw, int *iwork, int *liw, Jex jac, int *jt, struct cuLsodaCommonBlock *common)
{
	dlsoda_(fex, neq, y, t, tout, itol, rtol, atol, itask, istate, iopt, rwork, lrw, iwork, liw, jac, jt, common);
}


#endif

