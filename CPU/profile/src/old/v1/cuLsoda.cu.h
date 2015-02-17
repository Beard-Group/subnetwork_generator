/*
 *  cuLsoda.h
 *
 *	File Notes:	This file is a conversion of the REAL precision Livermore Solver for
 *	Ordinary Differential Equations with automatic switching for stiff and non-stiff
 *	problems (DLSODA)
 */

/**  barf  [ba:rf]  2.  "He suggested using FORTRAN, and everybody barfed."
 
 - From The Shogakukan DICTIONARY OF NEW ENGLISH (Second edition) */
#ifndef CULSODA_CU_H
#define CULSODA_CU_H

#include <string.h>
#define REAL double

using namespace std;

/* Common Block Declarations */
struct cuLsodaCommonBlock
{
	REAL /*rowns[209],*/ CM_conit, CM_crate, CM_ccmax, CM_el0, CM_h__, CM_hmin, CM_hmxi, CM_hu, CM_rc, CM_tn, CM_uround, CM_pdest, CM_pdlast, CM_ratio, CM_hold, CM_rmax;
	REAL  CM_el[13], CM_elco[156]	/* was [13][12] */, CM_tesco[36]	/* was [3][12] */;
	REAL CM_rls[218];
	REAL CM_tsw, /*rowns2[20],*/ CM_pdnorm;
	REAL /*rownd2,*/ CM_cm1[12], CM_cm2[5];
	REAL CM_rlsa[22];
	REAL CM_sm1[12];
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
struct myFex
{
        private:
        vector <double> _Hcoeff;
        int *_Hindex;
        vector <int> _Hnka, _Hkavec, _Hnkd, _Hkdvec;
        value_type *_Htau, *_Ht0;//, *_Hka64, *_Hka124;
        vector <double> _Hr0, _Hd, _Hea, _Hkaval, _Hkdval;
        double t_t, num, den;
        public:
        //void set_coeff ( const double *c0, const double c0_v[], const int size) 
        void set_index ( const int ind ) 
        { 
	     _Hindex = ( int *)malloc(sizeof(int));
	     *_Hindex = ind;
        };
        void set_tau ( const value_type tau ) 
        { 
	     _Htau = ( value_type *)malloc(sizeof(value_type));
	     *_Htau = tau;
        };
        void set_t0 ( const value_type t0 ) 
        { 
	     _Ht0 = ( value_type *)malloc(sizeof(value_type));
	     *_Ht0 = t0;
        };
        void set_r0 ( const vector <double> &r0 ) { _Hr0 = r0; };
        void get_r0 ( vector <double> &r0 ) { r0 = _Hr0; }; 
        void set_d ( const vector <double> &d ) { _Hd = d; };
        void get_d ( vector <double> &d ) { d = _Hd; }; 
        void set_ea ( const vector <double> &ea ) { _Hea = ea; };
        void get_ea ( vector <double> &ea ) { ea = _Hea; }; 
        void set_n_ka ( const vector <int> nka ) { _Hnka = nka; }; 
        //void get_n_ka ( vector <int> &nka ) { nka = _Hnka; }; 
        void set_n_kd ( const vector <int> nkd ) { _Hnkd = nkd; }; 
        //void get_n_kd ( vector <int> &nkd ) { nkd = _Hnkd; }; 
        void set_ka_vec ( const vector <int> kavec ) { _Hkavec = kavec;}; 
        void set_kd_vec ( const vector <int> kdvec ) { _Hkdvec = kdvec;}; 
        void set_ka_val ( const vector <double> &kaval ) { _Hkaval = kaval;}; 
        void get_ka_val ( vector <double> &kaval ) { kaval = _Hkaval;}; 
        void set_kd_val ( const vector <double> &kdval ) { _Hkdval = kdval;}; 
        void get_kd_val ( vector <double> &kdval ) { kdval = _Hkdval;}; 
        /*void set_r0 ( const value_type r0[] ) 
        { 
	     _Hr0 = ( double *)malloc(sizeof(double)*N_gene);
	     memcpy(_Hr0,r0,sizeof(double)*N_gene);
        };
        void set_d ( const value_type d[] ) 
        { 
	     _Hd = ( double *)malloc(sizeof(double)*N_gene);
	     memcpy(_Hd,d,sizeof(double)*N_gene);
        };
        void set_ea ( const value_type ea[] ) 
        { 
	     _Hea = ( double *)malloc(sizeof(double)*N_gene);
	     memcpy(_Hea,ea,sizeof(double)*N_gene);
        };
        void set_n_ka ( const int nka[] ) 
        { 
	     _Hnka = ( int * )malloc( sizeof(int)*N_gene );
	     memcpy(_Hnka,nka,sizeof(int)*N_gene);
        };
        void set_ka_vec ( const int kavec[], const int size ) 
        { 
	     _Hkavec = ( int * )malloc( sizeof(int)*N_gene*size );
	     memcpy(_Hkavec,kavec,sizeof(int)*N_gene*size);
        };
        void set_ka_val ( const value_type kaval[], const int size ) 
        { 
	     _Hkaval = ( value_type * )malloc( sizeof(value_type)*N_gene*size );
	     memcpy(_Hkaval,kaval,sizeof(value_type)*N_gene*size);
        };
        void set_n_kd ( const int nkd[] ) 
        { 
	     _Hnkd = ( int * )malloc( sizeof(int)*N_gene );
	     memcpy(_Hnkd,nkd,sizeof(int)*N_gene);
        };
        void set_kd_vec ( const int kdvec[], const int size ) 
        { 
	     _Hkdvec = ( int * )malloc( sizeof(int)*N_gene*size );
	     memcpy(_Hkdvec,kdvec,sizeof(int)*N_gene*size);
        };
        void set_kd_val ( const value_type kdval[], const int size ) 
        { 
	     _Hkdval = ( value_type * )malloc( sizeof(value_type)*N_gene*size );
	     memcpy(_Hkdval,kdval,sizeof(value_type)*N_gene*size);
        };
        void set_coeff ( const double c0_v[] ) 
        {
	     _Hcoeff = (double *)malloc(sizeof(double)*N_gene*N_time_points*4);
	     memcpy(_Hcoeff,c0_v,sizeof(double)*N_gene*N_time_points*4);
        };*/
        void set_coeff ( const vector<double> &c0_v ) { _Hcoeff = c0_v; };
        double get_coeff0 (int i, int j) const {return _Hcoeff[i*N_time_points*4+j*4];}; 
        double get_coeff1 (int i, int j) const {return _Hcoeff[i*N_time_points*4+j*4+1];};
        double get_coeff2 (int i, int j) const {return _Hcoeff[i*N_time_points*4+j*4+2];};
        double get_coeff3 (int i, int j) const {return _Hcoeff[i*N_time_points*4+j*4+3];};
	void operator()(int *neq, double *t, double *y, double *ydot/*, void *otherData*/)
	{
             t_t = 0.0;
             int ka_index = 0;
             int kd_index = 0;
             if ( *t - *_Htau <= 0.0 )
             {
                for ( int i = 0 ; i < N_gene; i++ )
                { 
                    num = 0.0;
                    den = 0.0;
                    for ( int j = 0; j < _Hnka[i]; j++)
                    {  
                        num += get_coeff0(_Hkavec[ka_index],*_Hindex)*_Hkaval[ka_index];
                        ka_index++; 
                    }
                    num *= num*num*num;
                    num += _Hea[i]*_Hea[i]*_Hea[i]*_Hea[i];
                    for ( int j = 0; j < _Hnkd[i]; j++)
                    {  
                        den += get_coeff0(_Hkdvec[kd_index],*_Hindex)*_Hkdval[kd_index];
                        kd_index++; 
                    }
                    den *= den*den*den;
                    den += (num + 1.0);
                    ydot[i] = _Hr0[i]*(num/den)-_Hd[i]*y[i];
                }
	     }
             else 
             {
                t_t = *t - *_Htau - *_Ht0;
                for ( int i = 0 ; i < N_gene; i++ )
                { 
                    num = 0.0;
                    den = 0.0;
                    for ( int j = 0; j < _Hnka[i]; j++)
                    {  
                        num += (get_coeff0(_Hkavec[ka_index],*_Hindex)+t_t*
                               (get_coeff1(_Hkavec[ka_index],*_Hindex)+t_t*
                               (get_coeff2(_Hkavec[ka_index],*_Hindex)+t_t*
                                    get_coeff3(_Hkavec[ka_index],*_Hindex))))*_Hkaval[ka_index];
                        ka_index++; 
                    }
                    num = num*num*num*num;
                    num = num + (_Hea[i]*_Hea[i]*_Hea[i]*_Hea[i]);
                    for ( int j = 0; j < _Hnkd[i]; j++)
                    {  
                        den += (get_coeff0(_Hkdvec[kd_index],*_Hindex)+t_t*
                               (get_coeff1(_Hkdvec[kd_index],*_Hindex)+t_t*
                               (get_coeff2(_Hkdvec[kd_index],*_Hindex)+t_t*
                                    get_coeff3(_Hkdvec[kd_index],*_Hindex))))*_Hkdval[kd_index];
                        kd_index++; 
                    }
                    den = den*den*den*den ;
                    den = den + (num + 1.0);
                    ydot[i] = _Hr0[i]*(num/den)-_Hd[i]*y[i];
                }
             }
             /*t_t = 0.0;
             int ka_index = 0;
             int kd_index = 0;
             if ( *t - 1.0 <= 0.0 )
             {
                for ( int i = 0 ; i < N_gene; i++ )
                { 
                    num = 0.0;
                    den = 0.0;
                    for ( int j = 0; j < _Hnka[i]; j++)
                    {  
                        num += _Hcoeff[_Hkavec[ka_index]*N_time_points*4]*_Hkaval[ka_index];
                        //num += get_coeff0(_Hkavec[ka_index],0)*_Hkaval[ka_index];
                        ka_index++; 
                    }
                    //num = num*num*num*num;
                    num *= num*num*num;
                    num += _Hea[i]*_Hea[i]*_Hea[i]*_Hea[i];
                    for ( int j = 0; j < _Hnkd[i]; j++)
                    {  
                        den += _Hcoeff[_Hkdvec[kd_index]*N_time_points*4]*_Hkdval[kd_index];
                        //den += get_coeff0(_Hkdvec[kd_index],0)*_Hkdval[kd_index];
                        kd_index++; 
                    }
                    den = den*den*den*den;
                    den += (num + 1.0);
                    //den = den + (num + 1.0);
                    ydot[i] = _Hr0[i]*(num/den)-_Hd[i]*y[i];
                }
	     }
             else if ( (*t-1.0 > 0.0) && (*t-1.0 <= 1.1818) )
             {
                t_t = *t - 1.0;
                for ( int i = 0 ; i < N_gene; i++ )
                { 
                    num = 0.0;
                    den = 0.0;
                    for ( int j = 0; j < _Hnka[i]; j++)
                    {  
                        num += (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+4]+t_t*
                               (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+4+1]+t_t*
                               (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+4+2]+t_t*
                                    _Hcoeff[_Hkavec[ka_index]*N_time_points*4+4+3])))*_Hkaval[ka_index];
                        ka_index++; 
                    }
                    //num = num*num*num*num;
                    num *= num*num*num;
                    num += _Hea[i]*_Hea[i]*_Hea[i]*_Hea[i];
                    for ( int j = 0; j < _Hnkd[i]; j++)
                    {  
                        den += (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+4]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+4+1]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+4+2]+t_t*
                                    _Hcoeff[_Hkdvec[kd_index]*N_time_points*4+4+3])))*_Hkdval[kd_index];
                        //den += (get_coeff0(_Hkdvec[kd_index],1)+t_t*
                        //       (get_coeff1(_Hkdvec[kd_index],1)+t_t*
                        //       (get_coeff2(_Hkdvec[kd_index],1)+t_t*
                        //            get_coeff3(_Hkdvec[kd_index],1))))*_Hkdval[kd_index];
                        kd_index++; 
                    }
                    den = den*den*den*den ;
                    den += (num + 1.0);
                    //den = den + (num + 1.0);
                    ydot[i] = _Hr0[i]*(num/den)-_Hd[i]*y[i];
                }
             }
             else if ( (*t-1.0 > 1.1818) && (*t-1.0 <= 2.3636) )
             {
                t_t = *t - 1.0 - 1.1818;
                for ( int i = 0 ; i < N_gene; i++ )
                { 
                    num = 0.0;
                    den = 0.0;
                    for ( int j = 0; j < _Hnka[i]; j++)
                    {  
                        num += (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+8]+t_t*
                               (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+8+1]+t_t*
                               (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+8+2]+t_t*
                                    _Hcoeff[_Hkavec[ka_index]*N_time_points*4+8+3])))*_Hkaval[ka_index];
                        //num += (get_coeff0(_Hkavec[ka_index],2)+t_t*
                        //       (get_coeff1(_Hkavec[ka_index],2)+t_t*
                        //       (get_coeff2(_Hkavec[ka_index],2)+t_t*
                        //            get_coeff3(_Hkavec[ka_index],2))))*_Hkaval[ka_index];
                        ka_index++; 
                    }
                    //num = num*num*num*num;
                    num *= num*num*num;
                    num += _Hea[i]*_Hea[i]*_Hea[i]*_Hea[i];
                    for ( int j = 0; j < _Hnkd[i]; j++)
                    {  
                        den += (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+8]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+8+1]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+8+2]+t_t*
                                    _Hcoeff[_Hkdvec[kd_index]*N_time_points*4+8+3])))*_Hkdval[kd_index];
                        //den += (get_coeff0(_Hkdvec[kd_index],2)+t_t*
                        //       (get_coeff1(_Hkdvec[kd_index],2)+t_t*
                        //       (get_coeff2(_Hkdvec[kd_index],2)+t_t*
                        //            get_coeff3(_Hkdvec[kd_index],2))))*_Hkdval[kd_index];
                        kd_index++; 
                    }
                    den = den*den*den*den ;
                    den += (num + 1.0);
                    //den = den + (num + 1.0);
                    ydot[i] = _Hr0[i]*(num/den)-_Hd[i]*y[i];
                }
             }
             else if ( (*t-1.0 > 2.3636) && (*t-1.0 <= 3.5455) )
             {
                t_t = *t - 1.0 - 2.3636;
                for ( int i = 0 ; i < N_gene; i++ )
                { 
                    num = 0.0;
                    den = 0.0;
                    for ( int j = 0; j < _Hnka[i]; j++)
                    {  
                        num += (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+12]+t_t*
                               (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+12+1]+t_t*
                               (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+12+2]+t_t*
                                    _Hcoeff[_Hkavec[ka_index]*N_time_points*4+12+3])))*_Hkaval[ka_index];
                        //num += (get_coeff0(_Hkavec[ka_index],3)+t_t*
                        //       (get_coeff1(_Hkavec[ka_index],3)+t_t*
                        //       (get_coeff2(_Hkavec[ka_index],3)+t_t*
                        //            get_coeff3(_Hkavec[ka_index],3))))*_Hkaval[ka_index];
                        ka_index++; 
                    }
                    //num = num*num*num*num;
                    num *= num*num*num;
                    num += _Hea[i]*_Hea[i]*_Hea[i]*_Hea[i];
                    for ( int j = 0; j < _Hnkd[i]; j++)
                    {  
                        den += (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+12]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+12+1]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+12+2]+t_t*
                                    _Hcoeff[_Hkdvec[kd_index]*N_time_points*4+12+3])))*_Hkdval[kd_index];
                        //den += (get_coeff0(_Hkdvec[kd_index],3)+t_t*
                        //       (get_coeff1(_Hkdvec[kd_index],3)+t_t*
                        //       (get_coeff2(_Hkdvec[kd_index],3)+t_t*
                        //            get_coeff3(_Hkdvec[kd_index],3))))*_Hkdval[kd_index];
                        kd_index++; 
                    }
                    den = den*den*den*den ;
                    den += (num + 1.0);
                    //den = den + (num + 1.0);
                    ydot[i] = _Hr0[i]*(num/den)-_Hd[i]*y[i];
                }
             }
             else if ( (*t-1.0 > 3.5455) && (*t-1.0 <= 4.7273) )
             {
                t_t = *t - 1.0 - 3.5455;
                for ( int i = 0 ; i < N_gene; i++ )
                { 
                    num = 0.0;
                    den = 0.0;
                    for ( int j = 0; j < _Hnka[i]; j++)
                    {  
                        num += (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+16]+t_t*
                               (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+16+1]+t_t*
                               (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+16+2]+t_t*
                                    _Hcoeff[_Hkavec[ka_index]*N_time_points*4+16+3])))*_Hkaval[ka_index];
                        //num += (get_coeff0(_Hkavec[ka_index],4)+t_t*
                        //       (get_coeff1(_Hkavec[ka_index],4)+t_t*
                        //       (get_coeff2(_Hkavec[ka_index],4)+t_t*
                        //            get_coeff3(_Hkavec[ka_index],4))))*_Hkaval[ka_index];
                        ka_index++; 
                    }
                    //num = num*num*num*num;
                    num *= num*num*num;
                    num += _Hea[i]*_Hea[i]*_Hea[i]*_Hea[i];
                    for ( int j = 0; j < _Hnkd[i]; j++)
                    {  
                        den += (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+16]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+16+1]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+16+2]+t_t*
                                    _Hcoeff[_Hkdvec[kd_index]*N_time_points*4+16+3])))*_Hkdval[kd_index];
                        //den += (get_coeff0(_Hkdvec[kd_index],4)+t_t*
                        //       (get_coeff1(_Hkdvec[kd_index],4)+t_t*
                        //       (get_coeff2(_Hkdvec[kd_index],4)+t_t*
                        //            get_coeff3(_Hkdvec[kd_index],4))))*_Hkdval[kd_index];
                        kd_index++; 
                    }
                    den = den*den*den*den ;
                    den += (num + 1.0);
                    //den = den + (num + 1.0);
                    ydot[i] = _Hr0[i]*(num/den)-_Hd[i]*y[i];
                }
             }
             else if ( (*t-1.0 > 4.7273) && (*t-1.0 <= 5.9091) )
             {
                t_t = *t - 1.0 - 4.7273;
                for ( int i = 0 ; i < N_gene; i++ )
                { 
                    num = 0.0;
                    den = 0.0;
                    for ( int j = 0; j < _Hnka[i]; j++)
                    {  
                        num += (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+20]+t_t*
                               (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+20+1]+t_t*
                               (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+20+2]+t_t*
                                    _Hcoeff[_Hkavec[ka_index]*N_time_points*4+20+3])))*_Hkaval[ka_index];
                        //num += (get_coeff0(_Hkavec[ka_index],5)+t_t*
                        //       (get_coeff1(_Hkavec[ka_index],5)+t_t*
                        //       (get_coeff2(_Hkavec[ka_index],5)+t_t*
                        //            get_coeff3(_Hkavec[ka_index],5))))*_Hkaval[ka_index];
                        ka_index++; 
                    }
                    //num = num*num*num*num;
                    num *= num*num*num;
                    num += _Hea[i]*_Hea[i]*_Hea[i]*_Hea[i];
                    for ( int j = 0; j < _Hnkd[i]; j++)
                    {  
                        den += (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+20]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+20+1]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+20+2]+t_t*
                                    _Hcoeff[_Hkdvec[kd_index]*N_time_points*4+20+3])))*_Hkdval[kd_index];
                        //den += (get_coeff0(_Hkdvec[kd_index],5)+t_t*
                        //       (get_coeff1(_Hkdvec[kd_index],5)+t_t*
                        //       (get_coeff2(_Hkdvec[kd_index],5)+t_t*
                        //            get_coeff3(_Hkdvec[kd_index],5))))*_Hkdval[kd_index];
                        kd_index++; 
                    }
                    den = den*den*den*den ;
                    den += (num + 1.0);
                    //den = den + (num + 1.0);
                    ydot[i] = _Hr0[i]*(num/den)-_Hd[i]*y[i];
                }
             }
             else if ( (*t-1.0 > 5.9091) && (*t-1.0 <= 7.0909) )
             {
                t_t = *t - 1.0 - 5.9091;
                for ( int i = 0 ; i < N_gene; i++ )
                { 
                    num = 0.0;
                    den = 0.0;
                    for ( int j = 0; j < _Hnka[i]; j++)
                    {  
                        num += (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+24]+t_t*
                               (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+24+1]+t_t*
                               (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+24+2]+t_t*
                                    _Hcoeff[_Hkavec[ka_index]*N_time_points*4+24+3])))*_Hkaval[ka_index];
                        //num += (get_coeff0(_Hkavec[ka_index],6)+t_t*
                        //       (get_coeff1(_Hkavec[ka_index],6)+t_t*
                        //       (get_coeff2(_Hkavec[ka_index],6)+t_t*
                        //            get_coeff3(_Hkavec[ka_index],6))))*_Hkaval[ka_index];
                        ka_index++; 
                    }
                    //num = num*num*num*num;
                    num *= num*num*num;
                    num += _Hea[i]*_Hea[i]*_Hea[i]*_Hea[i];
                    for ( int j = 0; j < _Hnkd[i]; j++)
                    {  
                        den += (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+24]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+24+1]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+24+2]+t_t*
                                    _Hcoeff[_Hkdvec[kd_index]*N_time_points*4+24+3])))*_Hkdval[kd_index];
                        //den += (get_coeff0(_Hkdvec[kd_index],6)+t_t*
                        //       (get_coeff1(_Hkdvec[kd_index],6)+t_t*
                        //       (get_coeff2(_Hkdvec[kd_index],6)+t_t*
                        //            get_coeff3(_Hkdvec[kd_index],6))))*_Hkdval[kd_index];
                        kd_index++; 
                    }
                    den = den*den*den*den ;
                    den += (num + 1.0);
                    //den = den + (num + 1.0);
                    ydot[i] = _Hr0[i]*(num/den)-_Hd[i]*y[i];
                }
             }
             else if ( (*t-1.0 > 7.0909) && (*t-1.0 <= 8.2727) )
             {
                t_t = *t - 1.0 - 7.0909;
                for ( int i = 0 ; i < N_gene; i++ )
                { 
                    num = 0.0;
                    den = 0.0;
                    for ( int j = 0; j < _Hnka[i]; j++)
                    {  
                        num += (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+28]+t_t*
                               (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+28+1]+t_t*
                               (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+28+2]+t_t*
                                    _Hcoeff[_Hkavec[ka_index]*N_time_points*4+28+3])))*_Hkaval[ka_index];
                       // num += (get_coeff0(_Hkavec[ka_index],7)+t_t*
                       //        (get_coeff1(_Hkavec[ka_index],7)+t_t*
                       //        (get_coeff2(_Hkavec[ka_index],7)+t_t*
                       //             get_coeff3(_Hkavec[ka_index],7))))*_Hkaval[ka_index];
                        ka_index++; 
                    }
                    //num = num*num*num*num;
                    num *= num*num*num;
                    num += _Hea[i]*_Hea[i]*_Hea[i]*_Hea[i];
                    for ( int j = 0; j < _Hnkd[i]; j++)
                    {  
                        den += (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+28]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+28+1]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+28+2]+t_t*
                                    _Hcoeff[_Hkdvec[kd_index]*N_time_points*4+28+3])))*_Hkdval[kd_index];
                        //den += (get_coeff0(_Hkdvec[kd_index],7)+t_t*
                        //       (get_coeff1(_Hkdvec[kd_index],7)+t_t*
                        //       (get_coeff2(_Hkdvec[kd_index],7)+t_t*
                        //            get_coeff3(_Hkdvec[kd_index],7))))*_Hkdval[kd_index];
                        kd_index++; 
                    }
                    den = den*den*den*den ;
                    den += (num + 1.0);
                    //den = den + (num + 1.0);
                    ydot[i] = _Hr0[i]*(num/den)-_Hd[i]*y[i];
                }
             }
             else if ( (*t-1.0 > 8.2727) && (*t-1.0 <= 9.4545) )
             {
                t_t = *t - 1.0 - 8.2727;
                for ( int i = 0 ; i < N_gene; i++ )
                { 
                    num = 0.0;
                    den = 0.0;
                    for ( int j = 0; j < _Hnka[i]; j++)
                    {  
                        num += (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+32]+t_t*
                               (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+32+1]+t_t*
                               (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+32+2]+t_t*
                                    _Hcoeff[_Hkavec[ka_index]*N_time_points*4+32+3])))*_Hkaval[ka_index];
                        //num += (get_coeff0(_Hkavec[ka_index],8)+t_t*
                        //       (get_coeff1(_Hkavec[ka_index],8)+t_t*
                        //      (get_coeff2(_Hkavec[ka_index],8)+t_t*
                        //            get_coeff3(_Hkavec[ka_index],8))))*_Hkaval[ka_index];
                        ka_index++; 
                    }
                    //num = num*num*num*num;
                    num *= num*num*num;
                    num += _Hea[i]*_Hea[i]*_Hea[i]*_Hea[i];
                    for ( int j = 0; j < _Hnkd[i]; j++)
                    {  
                        den += (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+32]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+32+1]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+32+2]+t_t*
                                    _Hcoeff[_Hkdvec[kd_index]*N_time_points*4+32+3])))*_Hkdval[kd_index];
                        //den += (get_coeff0(_Hkdvec[kd_index],8)+t_t*
                        //       (get_coeff1(_Hkdvec[kd_index],8)+t_t*
                        //       (get_coeff2(_Hkdvec[kd_index],8)+t_t*
                        //            get_coeff3(_Hkdvec[kd_index],8))))*_Hkdval[kd_index];
                        kd_index++; 
                    }
                    den = den*den*den*den ;
                    den += (num + 1.0);
                    //den = den + (num + 1.0);
                    ydot[i] = _Hr0[i]*(num/den)-_Hd[i]*y[i];
                }
             }
             else if ( (*t-1.0 > 9.4545) && (*t-1.0 <=10.6364) )
             {
                t_t = *t - 1.0 - 9.4545;
                for ( int i = 0 ; i < N_gene; i++ )
                { 
                    num = 0.0;
                    den = 0.0;
                    for ( int j = 0; j < _Hnka[i]; j++)
                    {  
                        num += (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+36]+t_t*
                               (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+36+1]+t_t*
                               (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+36+2]+t_t*
                                    _Hcoeff[_Hkavec[ka_index]*N_time_points*4+36+3])))*_Hkaval[ka_index];
                        //num += (get_coeff0(_Hkavec[ka_index],8)+t_t*
                        //       (get_coeff1(_Hkavec[ka_index],8)+t_t*
                        //      (get_coeff2(_Hkavec[ka_index],8)+t_t*
                        //            get_coeff3(_Hkavec[ka_index],8))))*_Hkaval[ka_index];
                        ka_index++; 
                    }
                    //num = num*num*num*num;
                    num *= num*num*num;
                    num += _Hea[i]*_Hea[i]*_Hea[i]*_Hea[i];
                    for ( int j = 0; j < _Hnkd[i]; j++)
                    {  
                        den += (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+36]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+36+1]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+36+2]+t_t*
                                    _Hcoeff[_Hkdvec[kd_index]*N_time_points*4+36+3])))*_Hkdval[kd_index];
                        //den += (get_coeff0(_Hkdvec[kd_index],8)+t_t*
                        //       (get_coeff1(_Hkdvec[kd_index],8)+t_t*
                        //       (get_coeff2(_Hkdvec[kd_index],8)+t_t*
                        //            get_coeff3(_Hkdvec[kd_index],8))))*_Hkdval[kd_index];
                        kd_index++; 
                    }
                    den = den*den*den*den ;
                    den += (num + 1.0);
                    //den = den + (num + 1.0);
                    ydot[i] = _Hr0[i]*(num/den)-_Hd[i]*y[i];
                }
             }
             else if ( (*t-1.0 >10.6364) && (*t-1.0 <=11.8182) )
             {
                t_t = *t - 1.0 -10.6364;
                for ( int i = 0 ; i < N_gene; i++ )
                { 
                    num = 0.0;
                    den = 0.0;
                    for ( int j = 0; j < _Hnka[i]; j++)
                    {  
                        num += (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+40]+t_t*
                               (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+40+1]+t_t*
                               (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+40+2]+t_t*
                                    _Hcoeff[_Hkavec[ka_index]*N_time_points*4+40+3])))*_Hkaval[ka_index];
                        //num += (get_coeff0(_Hkavec[ka_index],8)+t_t*
                        //       (get_coeff1(_Hkavec[ka_index],8)+t_t*
                        //      (get_coeff2(_Hkavec[ka_index],8)+t_t*
                        //            get_coeff3(_Hkavec[ka_index],8))))*_Hkaval[ka_index];
                        ka_index++; 
                    }
                    //num = num*num*num*num;
                    num *= num*num*num;
                    num += _Hea[i]*_Hea[i]*_Hea[i]*_Hea[i];
                    for ( int j = 0; j < _Hnkd[i]; j++)
                    {  
                        den += (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+40]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+40+1]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+40+2]+t_t*
                                    _Hcoeff[_Hkdvec[kd_index]*N_time_points*4+40+3])))*_Hkdval[kd_index];
                        //den += (get_coeff0(_Hkdvec[kd_index],8)+t_t*
                        //       (get_coeff1(_Hkdvec[kd_index],8)+t_t*
                        //       (get_coeff2(_Hkdvec[kd_index],8)+t_t*
                        //            get_coeff3(_Hkdvec[kd_index],8))))*_Hkdval[kd_index];
                        kd_index++; 
                    }
                    den = den*den*den*den ;
                    den += (num + 1.0);
                    //den = den + (num + 1.0);
                    ydot[i] = _Hr0[i]*(num/den)-_Hd[i]*y[i];
                }
             }
             else if (*t-1.0 >11.8182) 
             {
                t_t = *t - 1.0 -11.8182;
                for ( int i = 0 ; i < N_gene; i++ )
                { 
                    num = 0.0;
                    den = 0.0;
                    for ( int j = 0; j < _Hnka[i]; j++)
                    {  
                        num += (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+44]+t_t*
                               (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+44+1]+t_t*
                               (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+44+2]+t_t*
                                    _Hcoeff[_Hkavec[ka_index]*N_time_points*4+44+3])))*_Hkaval[ka_index];
                        //num += (get_coeff0(_Hkavec[ka_index],8)+t_t*
                        //       (get_coeff1(_Hkavec[ka_index],8)+t_t*
                        //      (get_coeff2(_Hkavec[ka_index],8)+t_t*
                        //            get_coeff3(_Hkavec[ka_index],8))))*_Hkaval[ka_index];
                        ka_index++; 
                    }
                    //num = num*num*num*num;
                    num *= num*num*num;
                    num += _Hea[i]*_Hea[i]*_Hea[i]*_Hea[i];
                    for ( int j = 0; j < _Hnkd[i]; j++)
                    {  
                        den += (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+44]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+44+1]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+44+2]+t_t*
                                    _Hcoeff[_Hkdvec[kd_index]*N_time_points*4+44+3])))*_Hkdval[kd_index];
                        //den += (get_coeff0(_Hkdvec[kd_index],8)+t_t*
                        //       (get_coeff1(_Hkdvec[kd_index],8)+t_t*
                        //       (get_coeff2(_Hkdvec[kd_index],8)+t_t*
                        //            get_coeff3(_Hkdvec[kd_index],8))))*_Hkdval[kd_index];
                        kd_index++; 
                    }
                    den = den*den*den*den ;
                    den += (num + 1.0);
                    //den = den + (num + 1.0);
                    ydot[i] = _Hr0[i]*(num/den)-_Hd[i]*y[i];
                }
             }*/
             /*t_t = 0.0;
             int ka_index = 0;
             if ( *t - 1.0 <= 0.0 )
             {
                for ( int i = 0 ; i < N_gene; i++ )
                { 
                    num = 0.0;
                    den = 0.0;
                    for ( int j = 0; j < _Hnka[i]; j++)
                    {  
                        num += get_coeff0(_Hkavec[ka_index],0)*_Hkaval[ka_index];
                        ka_index++; 
                    }
                    num *= num;
                    den = num + 1.0;
                    ydot[i] = _Hr0[i]*((num/den)+_Hea[i])-_Hd[i]*y[i];
                }
	     }
             else if ( (*t-1.0 > 0.0) && (*t-1.0 <= 7.0) )
             {
                t_t = *t - 1.0;
                for ( int i = 0 ; i < N_gene; i++ )
                { 
                    num = 0.0;
                    den = 0.0;
                    for ( int j = 0; j < _Hnka[i]; j++)
                    {  
                        num += (get_coeff0(_Hkavec[ka_index],1)+t_t*
                               (get_coeff1(_Hkavec[ka_index],1)+t_t*
                               (get_coeff2(_Hkavec[ka_index],1)+t_t*
                                    get_coeff3(_Hkavec[ka_index],1))))*_Hkaval[ka_index];
                        ka_index++; 
                    }
                    num *= num;
                    den = num + 1.0;
                    ydot[i] = _Hr0[i]*((num/den)+_Hea[i])-_Hd[i]*y[i];
                }
             }
             else if ( (*t-1.0 > 7.0) && (*t-1.0 <= 8.0) )
             {
                t_t = *t - 1.0 - 7.0;
                for ( int i = 0 ; i < N_gene; i++ )
                { 
                    num = 0.0;
                    den = 0.0;
                    for ( int j = 0; j < _Hnka[i]; j++)
                    {  
                        num += (get_coeff0(_Hkavec[ka_index],2)+t_t*
                               (get_coeff1(_Hkavec[ka_index],2)+t_t*
                               (get_coeff2(_Hkavec[ka_index],2)+t_t*
                                    get_coeff3(_Hkavec[ka_index],2))))*_Hkaval[ka_index];
                        ka_index++; 
                    }
                    num *= num;
                    den = num + 1.0;
                    ydot[i] = _Hr0[i]*((num/den)+_Hea[i])-_Hd[i]*y[i];
                }
             }
             else if ( (*t-1.0 > 8.0) && (*t-1.0 <= 9.0) )
             {
                t_t = *t - 1.0 - 8.0;
                for ( int i = 0 ; i < N_gene; i++ )
                { 
                    num = 0.0;
                    den = 0.0;
                    for ( int j = 0; j < _Hnka[i]; j++)
                    {  
                        num += (get_coeff0(_Hkavec[ka_index],3)+t_t*
                               (get_coeff1(_Hkavec[ka_index],3)+t_t*
                               (get_coeff2(_Hkavec[ka_index],3)+t_t*
                                    get_coeff3(_Hkavec[ka_index],3))))*_Hkaval[ka_index];
                        ka_index++; 
                    }
                    num *= num;
                    den = num + 1.0;
                    ydot[i] = _Hr0[i]*((num/den)+_Hea[i])-_Hd[i]*y[i];
                }
             }
             else if ( (*t-1.0 > 9.0) && (*t-1.0 <= 12.0) )
             {
                t_t = *t - 1.0 - 9.0;
                for ( int i = 0 ; i < N_gene; i++ )
                { 
                    num = 0.0;
                    den = 0.0;
                    for ( int j = 0; j < _Hnka[i]; j++)
                    {  
                        num += (get_coeff0(_Hkavec[ka_index],4)+t_t*
                               (get_coeff1(_Hkavec[ka_index],4)+t_t*
                               (get_coeff2(_Hkavec[ka_index],4)+t_t*
                                    get_coeff3(_Hkavec[ka_index],4))))*_Hkaval[ka_index];
                        ka_index++; 
                    }
                    num *= num;
                    den = num + 1.0;
                    ydot[i] = _Hr0[i]*((num/den)+_Hea[i])-_Hd[i]*y[i];
                }
             }
             else if ( (*t-1.0 > 12.0) && (*t-1.0 <= 14.0) )
             {
                t_t = *t - 1.0 - 12.0;
                for ( int i = 0 ; i < N_gene; i++ )
                { 
                    num = 0.0;
                    den = 0.0;
                    for ( int j = 0; j < _Hnka[i]; j++)
                    {  
                        num += (get_coeff0(_Hkavec[ka_index],5)+t_t*
                               (get_coeff1(_Hkavec[ka_index],5)+t_t*
                               (get_coeff2(_Hkavec[ka_index],5)+t_t*
                                    get_coeff3(_Hkavec[ka_index],5))))*_Hkaval[ka_index];
                        ka_index++; 
                    }
                    num *= num;
                    den = num + 1.0;
                    ydot[i] = _Hr0[i]*((num/den)+_Hea[i])-_Hd[i]*y[i];
                }
             }
             else if ( (*t-1.0 > 14.0) && (*t-1.0 <= 18.0) )
             {
                t_t = *t - 1.0 - 14.0;
                for ( int i = 0 ; i < N_gene; i++ )
                { 
                    num = 0.0;
                    den = 0.0;
                    for ( int j = 0; j < _Hnka[i]; j++)
                    {  
                        num += (get_coeff0(_Hkavec[ka_index],6)+t_t*
                               (get_coeff1(_Hkavec[ka_index],6)+t_t*
                               (get_coeff2(_Hkavec[ka_index],6)+t_t*
                                    get_coeff3(_Hkavec[ka_index],6))))*_Hkaval[ka_index];
                        ka_index++; 
                    }
                    num *= num;
                    den = num + 1.0;
                    ydot[i] = _Hr0[i]*((num/den)+_Hea[i])-_Hd[i]*y[i];
                }
             }
             else if ( (*t-1.0 > 18.0) && (*t-1.0 <= 23.0) )
             {
                t_t = *t - 1.0 - 18.0;
                for ( int i = 0 ; i < N_gene; i++ )
                { 
                    num = 0.0;
                    den = 0.0;
                    for ( int j = 0; j < _Hnka[i]; j++)
                    {  
                        num += (get_coeff0(_Hkavec[ka_index],7)+t_t*
                               (get_coeff1(_Hkavec[ka_index],7)+t_t*
                               (get_coeff2(_Hkavec[ka_index],7)+t_t*
                                    get_coeff3(_Hkavec[ka_index],7))))*_Hkaval[ka_index];
                        ka_index++; 
                    }
                    num *= num;
                    den = num + 1.0;
                    ydot[i] = _Hr0[i]*((num/den)+_Hea[i])-_Hd[i]*y[i];
                }
             }
             else if ( *t-1.0 > 23.0 ) 
             {
                t_t = *t - 1.0 - 23.0;
                for ( int i = 0 ; i < N_gene; i++ )
                { 
                    num = 0.0;
                    den = 0.0;
                    for ( int j = 0; j < _Hnka[i]; j++)
                    {  
                        num += (get_coeff0(_Hkavec[ka_index],8)+t_t*
                               (get_coeff1(_Hkavec[ka_index],8)+t_t*
                               (get_coeff2(_Hkavec[ka_index],8)+t_t*
                                    get_coeff3(_Hkavec[ka_index],8))))*_Hkaval[ka_index];
                        ka_index++; 
                    }
                    num *= num;
                    den = num + 1.0;
                    ydot[i] = _Hr0[i]*((num/den)+_Hea[i])-_Hd[i]*y[i];
                }
             }*/
	}
};

struct myJex
{
	void operator()(int *neq, REAL *t, REAL *y, int ml, int mu, REAL *pd, int nrowpd/*, void *otherData*/)
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
int dlsoda_(Fex, int *, REAL *, REAL *, REAL *, int *, REAL *, REAL *, int *, int *, int *, REAL *, int *, int *, int *, Jex, int *, struct cuLsodaCommonBlock *);

template<typename Fex, typename Jex> 
int dstoda_(int *neq, REAL *y, REAL *yh, int *NOT_nyh, REAL *yh1, REAL *ewt, REAL *savf, REAL *acor, REAL *wm, int *iwm, Fex f, Jex jac, struct cuLsodaCommonBlock *common);

template<typename Fex, typename Jex> 
int dprja_(int *neq, REAL *y, REAL *yh, int *NOT_nyh, REAL *ewt, REAL *ftem, REAL *savf, REAL *wm, int *iwm, Fex f, Jex jac, struct cuLsodaCommonBlock *common);

int dsolsy_(REAL *wm, int *iwm, REAL *x, REAL *tem, struct cuLsodaCommonBlock *common);
int dintdy_(REAL *t, int k, REAL *yh, int *NOT_nyh, REAL *dky, int *iflag, struct cuLsodaCommonBlock *common);
int dcfode_(int meth, REAL *DCFODE_elco, REAL *DCFODE_tesco, struct cuLsodaCommonBlock *common);
int dsolsy_(REAL *wm, int *iwm, REAL *x, REAL *tem, struct cuLsodaCommonBlock *common);
int dewset_(int *PARAM_n, int *itol, REAL *rtol, REAL *atol, REAL *ycur, REAL *ewt, struct cuLsodaCommonBlock *common);
REAL dmnorm_(int *PARAM_n, REAL *v, REAL *w, struct cuLsodaCommonBlock *common);
REAL dfnorm_(int *PARAM_n, REAL *a, REAL *w, struct cuLsodaCommonBlock *common);
REAL dbnorm_(int *PARAM_n, REAL *a, int *nra, int *ml, int *mu, REAL *w, struct cuLsodaCommonBlock *common);
int dsrcma_(REAL *rsav, int *isav, int *job, struct cuLsodaCommonBlock *common);
int dgefa_(REAL *a, int *lda, int *PARAM_n, int *ipvt, int *info, struct cuLsodaCommonBlock *common);
int dgesl_(REAL *a, int *lda, int *PARAM_n, int *ipvt, REAL *b, int job, struct cuLsodaCommonBlock *common);
int dgbfa_(REAL *abd, int *lda, int *PARAM_n, int *ml, int *mu, int *ipvt, int *info, struct cuLsodaCommonBlock *common);
int dgbsl_(REAL *abd, int *lda, int *PARAM_n, int *ml, int *mu, int *ipvt, REAL *b, int job, struct cuLsodaCommonBlock *common);
REAL dumach_( struct cuLsodaCommonBlock *common);
//int xsetf_(int *mflag, struct CommonBlock *common);
//int xsetun_(int *lun, struct CommonBlock *common);
int ixsav_(int ipar, int *ivalue, int iset, struct cuLsodaCommonBlock *common);
int idamax_(int *PARAM_n, REAL *dx, int incx, struct cuLsodaCommonBlock *common);
int daxpy_(int *PARAM_n, REAL *da, REAL *dx, int incx, REAL *dy, int incy, struct cuLsodaCommonBlock *common);
int dumsum_(REAL a, REAL b, REAL *c__, struct cuLsodaCommonBlock *common);
int dscal_(int *PARAM_n, REAL *da, REAL *dx, int incx, struct cuLsodaCommonBlock *common);
REAL ddot_(int *PARAM_n, REAL *dx, int incx, REAL *dy, int incy, struct cuLsodaCommonBlock *common);
REAL d_sign(REAL *a, REAL *b);
void cuLsodaCommonBlockInit(struct cuLsodaCommonBlock *common);

#ifndef use_export
#include "cuLsoda.cpp"
#endif

#endif
