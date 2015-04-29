/*
 *  cuLsoda.h
 *
 *	File Notes:	This file is a conversion of the double precision Livermore Solver for
 *	Ordinary Differential Equations with automatic switching for stiff and non-stiff
 *	problems (DLSODA)
 */

/**  barf  [ba:rf]  2.  "He suggested using FORTRAN, and everybody barfed."
 
 - From The Shogakukan DICTIONARY OF NEW ENGLISH (Second edition) */
#ifndef CULSODA_CU_H
#define CULSODA_CU_H

#include <string.h>


using namespace std;
/* Common Block Declarations */
struct cuLsodaCommonBlock
{
	double /*rowns[209],*/ CM_conit, CM_crate, CM_ccmax, CM_el0, CM_h__, CM_hmin, CM_hmxi, CM_hu, CM_rc, CM_tn, CM_uround, CM_pdest, CM_pdlast, CM_ratio, CM_hold, CM_rmax;
	double  CM_el[13], CM_elco[156]	/* was [13][12] */, CM_tesco[36]	/* was [3][12] */;
	double CM_rls[218];
	double CM_tsw, /*rowns2[20],*/ CM_pdnorm;
	double /*rownd2,*/ CM_cm1[12], CM_cm2[5];
	double CM_rlsa[22];
	double CM_sm1[12];
	int CM_init, CM_mxstep, CM_mxhnil, CM_nhnil, CM_nslast, CM_nyh, /*iowns[6],*/ CM_icf, 
	CM_ierpj, CM_iersl, CM_jcur, CM_jstart, CM_kflag, CM_l, CM_lyh, CM_lewt, CM_lacor, CM_lsavf,
	CM_lwm, CM_liwm, CM_meth, CM_miter, CM_maxord, CM_maxcor, CM_msbp, CM_mxncf, CM_n, CM_nq, 
	CM_nst, CM_nfe, CM_nje, CM_nqu;
	int /*iownd[6],*/ CM_ialth, CM_ipup, CM_lmax, /*meo,*/ CM_nqnyh, CM_nslp;
	int CM_ils[37];
	int CM_insufr, CM_insufi, CM_ixpr, /*iowns2[2],*/ CM_jtyp, CM_mused, CM_mxordn, CM_mxords; 
	int /*iownd2[3],*/ CM_icount, CM_irflag;
	int CM_ilsa[9];
};

/* End Common Block */ 


#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))


#define Fex_and_Jex_definition
// The model class is basically Fex
struct myFex
{
        private:
        double *_Hcoeff;
        double *_Hr0, *_Hd, *_Hea, *_Hka64, *_Hka124;
        double t_t, num, den;
        //static const int 1.0 = 1.0;
        public:
        //void set_coeff ( const double *c0, const double c0_v[], const int size) 
        void set_r0 ( const double r0 ) 
        { 
	     _Hr0 = ( double *)malloc(sizeof(double));
             *_Hr0 = r0;
        };
        void set_d ( const double d ) 
        { 
	     _Hd = ( double *)malloc(sizeof(double));
             *_Hd = d;
        };
        void set_ea ( const double ea ) 
        { 
	     _Hea = ( double *)malloc(sizeof(double));
             *_Hea = ea;
        };
        void set_ka_64 ( const double ka_64 ) 
        { 
	     _Hka64 = ( double *)malloc(sizeof(double));
             *_Hka64 = ka_64;
        };
        void set_ka_124 ( const double ka_124 ) 
        { 
	     _Hka124 = ( double *)malloc(sizeof(double));
             *_Hka124 = ka_124;
        };
        void set_coeff ( const double c0_v[], const int size) 
        {
	     _Hcoeff = (double *)malloc(sizeof(double)*size);
	     //cudaMemcpy(_Hcoeff,c0_v,sizeof(double)*size,cudaMemcpyHostToDevice);
	     memcpy(_Hcoeff,c0_v,sizeof(double)*size);
        };
	  void operator()(int *neq, double *t, double *y, double *ydot/*, void *otherData*/)
	{
             t_t = 0.0;
             if ( *t - 1.0 <= 0.0 )
             { 
                num = _Hcoeff[36] * (*_Hka64) + _Hcoeff[72] * (*_Hka124);
                num *= num;
                den = num + 1.0;
                ydot[0] = *_Hr0*((num/den)+*_Hea)-*_Hd*y[0];
	     }
             else if ( (*t-1.0 > 0.0) && (*t-1.0 <= 7.0) )
             {
                t_t = *t - 1.0;
                num = (_Hcoeff[40]+t_t*(_Hcoeff[41]+t_t*(_Hcoeff[42]+t_t*_Hcoeff[43])))*(*_Hka64); 
                num += (_Hcoeff[76]+t_t*(_Hcoeff[77]+t_t*(_Hcoeff[78]+t_t*_Hcoeff[79])))*(*_Hka124); 
                num *= num;
                den = num + 1.0;
                ydot[0] = *_Hr0*((num/den)+*_Hea)-*_Hd*y[0];
             }
             else if ( (*t-1.0 > 7.0) && (*t-1.0 <= 8.0) )
             {
                t_t = *t - 1.0 - 7.0;
                num = (_Hcoeff[44]+t_t*(_Hcoeff[45]+t_t*(_Hcoeff[46]+t_t*_Hcoeff[47])))*(*_Hka64); 
                num += (_Hcoeff[80]+t_t*(_Hcoeff[81]+t_t*(_Hcoeff[82]+t_t*_Hcoeff[83])))*(*_Hka124); 
                num *= num;
                den = num + 1.0;
                ydot[0] = *_Hr0*((num/den)+*_Hea)-*_Hd*y[0];
             }
             else if ( (*t-1.0 > 8.0) && (*t-1.0 <= 9.0) )
             {
                t_t = *t - 1.0 - 8.0;
                num = (_Hcoeff[48]+t_t*(_Hcoeff[49]+t_t*(_Hcoeff[50]+t_t*_Hcoeff[51])))*(*_Hka64); 
                num += (_Hcoeff[84]+t_t*(_Hcoeff[85]+t_t*(_Hcoeff[86]+t_t*_Hcoeff[87])))*(*_Hka124); 
                num *= num;
                den = num + 1.0;
                ydot[0] = *_Hr0*((num/den)+*_Hea)-*_Hd*y[0];
             }
             else if ( (*t-1.0 > 9.0) && (*t-1.0 <= 12.0) )
             {
                t_t = *t - 1.0 - 9.0;
                num = (_Hcoeff[52]+t_t*(_Hcoeff[53]+t_t*(_Hcoeff[54]+t_t*_Hcoeff[55])))*(*_Hka64); 
                num += (_Hcoeff[88]+t_t*(_Hcoeff[89]+t_t*(_Hcoeff[90]+t_t*_Hcoeff[91])))*(*_Hka124); 
                num *= num;
                den = num + 1.0;
                ydot[0] = *_Hr0*((num/den)+*_Hea)-*_Hd*y[0];
             }
             else if ( (*t-1.0 > 12.0) && (*t-1.0 <= 14.0) )
             {
                t_t = *t - 1.0 - 12.0;
                num = (_Hcoeff[56]+t_t*(_Hcoeff[57]+t_t*(_Hcoeff[58]+t_t*_Hcoeff[59])))*(*_Hka64); 
                num += (_Hcoeff[92]+t_t*(_Hcoeff[93]+t_t*(_Hcoeff[94]+t_t*_Hcoeff[95])))*(*_Hka124); 
                num *= num; 
                den = num + 1.0;
                ydot[0] = *_Hr0*((num/den)+*_Hea)-*_Hd*y[0];
             }
             else if ( (*t-1.0 > 14.0) && (*t-1.0 <= 18.0) )
             {
                t_t = *t - 1.0 - 14.0;
                num = (_Hcoeff[60]+t_t*(_Hcoeff[61]+t_t*(_Hcoeff[62]+t_t*_Hcoeff[63])))*(*_Hka64); 
                num += (_Hcoeff[96]+t_t*(_Hcoeff[97]+t_t*(_Hcoeff[98]+t_t*_Hcoeff[99])))*(*_Hka124); 
                num *= num; 
                den = num + 1.0;
                ydot[0] = *_Hr0*((num/den)+*_Hea)-*_Hd*y[0];
             }
             else if ( (*t-1.0 > 18.0) && (*t-1.0 <= 23.0) )
             {
                t_t = *t - 1.0 - 18.0;
                num = (_Hcoeff[64]+t_t*(_Hcoeff[65]+t_t*(_Hcoeff[66]+t_t*_Hcoeff[67])))*(*_Hka64); 
                num += (_Hcoeff[100]+t_t*(_Hcoeff[101]+t_t*(_Hcoeff[102]+t_t*_Hcoeff[103])))*(*_Hka124); 
                num *= num; 
                den = num + 1.0;
                ydot[0] = *_Hr0*((num/den)+*_Hea)-*_Hd*y[0];
             }
             else if ( *t-1.0 > 23.0 ) 
             {
                t_t = *t - 1.0 - 23.0;
                num = (_Hcoeff[68]+t_t*(_Hcoeff[69]+t_t*(_Hcoeff[70]+t_t*_Hcoeff[71])))*(*_Hka64); 
                num += (_Hcoeff[104]+t_t*(_Hcoeff[105]+t_t*(_Hcoeff[106]+t_t*_Hcoeff[107])))*(*_Hka124); 
                num *= num; 
                den = num + 1.0;
                ydot[0] = *_Hr0*((num/den)+*_Hea)-*_Hd*y[0];
             }
	}
};

struct myJex
{
	  void operator()(int *neq, double *t, double *y, int ml, int mu, double *pd, int nrowpd/*, void *otherData*/)
	{
		return;
	}
};




/* dlsoda.f -- translated by f2c (version 20090411).
 You must link the resulting object file with libf2c:
 on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm
 or, if you install libf2c.a in a standard place, with -lf2c -lm
 -- in that order, at the end of the command line, as in
 cc *.o -lf2c -lm
 Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,
 
 http://www.netlib.org/f2c/libf2c.zip
 */




template<typename Fex, typename Jex>
  int dlsoda_(Fex, int *, double *, double *, double *, int *, double *, double *, int *, int *, int *, double *, int *, int *, int *, Jex, int *, struct cuLsodaCommonBlock *);

template<typename Fex, typename Jex> 
  int dstoda_(int *neq, double *y, double *yh, int *NOT_nyh, double *yh1, double *ewt, double *savf, double *acor, double *wm, int *iwm, Fex f, Jex jac, struct cuLsodaCommonBlock *common);

template<typename Fex, typename Jex> 
  int dprja_(int *neq, double *y, double *yh, int *NOT_nyh, double *ewt, double *ftem, double *savf, double *wm, int *iwm, Fex f, Jex jac, struct cuLsodaCommonBlock *common);

  int dsolsy_(double *wm, int *iwm, double *x, double *tem, struct cuLsodaCommonBlock *common);
  int dintdy_(double *t, int k, double *yh, int *NOT_nyh, double *dky, int *iflag, struct cuLsodaCommonBlock *common);
  int dcfode_(int meth, double *DCFODE_elco, double *DCFODE_tesco, struct cuLsodaCommonBlock *common);
  int dsolsy_(double *wm, int *iwm, double *x, double *tem, struct cuLsodaCommonBlock *common);
  int dewset_(int *PARAM_n, int *itol, double *rtol, double *atol, double *ycur, double *ewt, struct cuLsodaCommonBlock *common);
  double dmnorm_(int *PARAM_n, double *v, double *w, struct cuLsodaCommonBlock *common);
  double dfnorm_(int *PARAM_n, double *a, double *w, struct cuLsodaCommonBlock *common);
  double dbnorm_(int *PARAM_n, double *a, int *nra, int *ml, int *mu, double *w, struct cuLsodaCommonBlock *common);
  int dsrcma_(double *rsav, int *isav, int *job, struct cuLsodaCommonBlock *common);
  int dgefa_(double *a, int *lda, int *PARAM_n, int *ipvt, int *info, struct cuLsodaCommonBlock *common);
  int dgesl_(double *a, int *lda, int *PARAM_n, int *ipvt, double *b, int job, struct cuLsodaCommonBlock *common);
  int dgbfa_(double *abd, int *lda, int *PARAM_n, int *ml, int *mu, int *ipvt, int *info, struct cuLsodaCommonBlock *common);
  int dgbsl_(double *abd, int *lda, int *PARAM_n, int *ml, int *mu, int *ipvt, double *b, int job, struct cuLsodaCommonBlock *common);
  double dumach_( struct cuLsodaCommonBlock *common);
//  int xsetf_(int *mflag, struct CommonBlock *common);
//  int xsetun_(int *lun, struct CommonBlock *common);
  int ixsav_(int ipar, int *ivalue, int iset, struct cuLsodaCommonBlock *common);
  int idamax_(int *PARAM_n, double *dx, int incx, struct cuLsodaCommonBlock *common);
  int daxpy_(int *PARAM_n, double *da, double *dx, int incx, double *dy, int incy, struct cuLsodaCommonBlock *common);
  int dumsum_(double a, double b, double *c__, struct cuLsodaCommonBlock *common);
  int dscal_(int *PARAM_n, double *da, double *dx, int incx, struct cuLsodaCommonBlock *common);
  double ddot_(int *PARAM_n, double *dx, int incx, double *dy, int incy, struct cuLsodaCommonBlock *common);
  double d_sign(double *a, double *b);
    void cuLsodaCommonBlockInit(struct cuLsodaCommonBlock *common);

#ifndef use_export
#include "./cuLsoda.cpp"
#endif

#endif



