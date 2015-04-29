/*
 *  cuLsoda.h
 *
 *	File Notes:	This file is a conversion of the value_type precision Livermore Solver for
 *	Ordinary Differential Equations with automatic switching for stiff and non-stiff
 *	problems (DLSODA)
 */

/**  barf  [ba:rf]  2.  "He suggested using FORTRAN, and everybody barfed."
 
 - From The Shogakukan DICTIONARY OF NEW ENGLISH (Second edition) */
#ifndef CULSODA_FEX_H
#define CULSODA_FEX_H
#include "../common.h"

using namespace std;

#define Fex_and_Jex_definition
/*struct myFex
{
        private:
        vector <double> _Hcoeff;
        //int *_Hindex;
        vector <int> _Hnka, _Hkavec, _Hnkd, _Hkdvec;
        value_type *_Htau;
        vector <double> _Hr0, _Hd, _Hea, _Hkaval, _Hkdval;
        double t_t, num, den;
        public:
        void set_tau ( const value_type tau ) 
        { 
	     _Htau = ( value_type *)malloc(sizeof(value_type));
	     *_Htau = tau;
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
        void set_coeff ( const vector<double> &c0_v ) { _Hcoeff = c0_v; };
        double get_coeff0 (int i, int j) const {return _Hcoeff[i*N_time_points*4+j*4];}; 
        double get_coeff1 (int i, int j) const {return _Hcoeff[i*N_time_points*4+j*4+1];};
        double get_coeff2 (int i, int j) const {return _Hcoeff[i*N_time_points*4+j*4+2];};
        double get_coeff3 (int i, int j) const {return _Hcoeff[i*N_time_points*4+j*4+3];};
	void operator()(int *neq, double *t, double *y, double *ydot)
	{*/
             /*t_t = 0.0;
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
             }*/
             /*t_t = 0.0;
             int ka_index = 0;
             int kd_index = 0;
             if ( *t - *_Htau <= 0.0 )
             {
                for ( int i = 0 ; i < N_gene; i++ )
                { 
                    num = 0.0;
                    den = 0.0;
                    //cout << "------------- Gene " << i << "  START----------------" << endl;
                    //cout << " Value of time is  " << *t << endl;
                    for ( int j = 0; j < _Hnka[i]; j++)
                    {  
                        //cout << "calculating numerator " << endl;
                        //cout << "Accessing element " << _Hkavec[ka_index]*48 << " of coeff vector." << endl;
                        //cout << "Value of kaval is: " << _Hkaval[ka_index] << endl;
                        num += _Hcoeff[_Hkavec[ka_index]*N_time_points*4]*_Hkaval[ka_index];
                        //num += get_coeff0(_Hkavec[ka_index],0)*_Hkaval[ka_index];
                        ka_index++; 
                    }
                    //cout << "Value of num(1) is: " << num << endl;
                    num *= num*num*num;
                    //cout << "Value of num(2) is: " << num << endl;
                    num += _Hea[i]*_Hea[i]*_Hea[i]*_Hea[i];
                    //cout << "Value of num(3) is: " << num << endl;
                    //cout << "Value of ea is: " << _Hea[0] << endl;
                    for ( int j = 0; j < _Hnkd[i]; j++)
                    {  
                        den += _Hcoeff[_Hkdvec[kd_index]*N_time_points*4]*_Hkdval[kd_index];
                        kd_index++; 
                    }
                    //cout << "Value of den(1) is: " << den << endl;
                    den = den*den*den*den;
                    //cout << "Value of den(2) is: " << den << endl;
                    den += (num + 1.0);
                    //cout << "Value of den(3) is: " << den << endl;
                    //cout << "Value of r0 is: " << _Hr0[0] << endl;
                    //cout << "Value of d is: " << _Hd[0] << endl;
                    ydot[i] = _Hr0[i]*(num/den)-_Hd[i]*y[i];
                    //cout << "Value of ydot[0] is: " << ydot[i] << endl;
                    //cout << "Value of y[0] is: " << y[i] << endl;
                    //cout << "------------- Gene " << i << "  END------------------" << endl;
                }
	     }
             else if ( (*t-*_Htau > 0.0) && (*t-*_Htau <= 1.1818) )
             {
                t_t = *t - *_Htau;
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
                    ydot[i] = _Hr0[i]*(num/den)-_Hd[i]*y[i];
                }
             }
             else if ( (*t-*_Htau > 1.1818) && (*t-*_Htau <= 2.3636) )
             {
                t_t = *t - *_Htau - 1.1818;
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
                        ka_index++; 
                    }
                    num *= num*num*num;
                    num += _Hea[i]*_Hea[i]*_Hea[i]*_Hea[i];
                    for ( int j = 0; j < _Hnkd[i]; j++)
                    {  
                        den += (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+8]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+8+1]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+8+2]+t_t*
                                    _Hcoeff[_Hkdvec[kd_index]*N_time_points*4+8+3])))*_Hkdval[kd_index];
                        kd_index++; 
                    }
                    den = den*den*den*den ;
                    den += (num + 1.0);
                    ydot[i] = _Hr0[i]*(num/den)-_Hd[i]*y[i];
                }
             }
             else if ( (*t-*_Htau > 2.3636) && (*t-*_Htau <= 3.5455) )
             {
                t_t = *t - *_Htau - 2.3636;
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
                        ka_index++; 
                    }
                    num *= num*num*num;
                    num += _Hea[i]*_Hea[i]*_Hea[i]*_Hea[i];
                    for ( int j = 0; j < _Hnkd[i]; j++)
                    {  
                        den += (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+12]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+12+1]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+12+2]+t_t*
                                    _Hcoeff[_Hkdvec[kd_index]*N_time_points*4+12+3])))*_Hkdval[kd_index];
                        kd_index++; 
                    }
                    den = den*den*den*den ;
                    den += (num + 1.0);
                    ydot[i] = _Hr0[i]*(num/den)-_Hd[i]*y[i];
                }
             }
             else if ( (*t-*_Htau > 3.5455) && (*t-*_Htau <= 4.7273) )
             {
                t_t = *t - *_Htau - 3.5455;
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
                        ka_index++; 
                    }
                    num *= num*num*num;
                    num += _Hea[i]*_Hea[i]*_Hea[i]*_Hea[i];
                    for ( int j = 0; j < _Hnkd[i]; j++)
                    {  
                        den += (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+16]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+16+1]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+16+2]+t_t*
                                    _Hcoeff[_Hkdvec[kd_index]*N_time_points*4+16+3])))*_Hkdval[kd_index];
                        kd_index++; 
                    }
                    den = den*den*den*den ;
                    den += (num + 1.0);
                    ydot[i] = _Hr0[i]*(num/den)-_Hd[i]*y[i];
                }
             }
             else if ( (*t-*_Htau > 4.7273) && (*t-*_Htau <= 5.9091) )
             {
                t_t = *t - *_Htau - 4.7273;
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
                        ka_index++; 
                    }
                    num *= num*num*num;
                    num += _Hea[i]*_Hea[i]*_Hea[i]*_Hea[i];
                    for ( int j = 0; j < _Hnkd[i]; j++)
                    {  
                        den += (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+20]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+20+1]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+20+2]+t_t*
                                    _Hcoeff[_Hkdvec[kd_index]*N_time_points*4+20+3])))*_Hkdval[kd_index];
                        kd_index++; 
                    }
                    den = den*den*den*den ;
                    den += (num + 1.0);
                    ydot[i] = _Hr0[i]*(num/den)-_Hd[i]*y[i];
                }
             }
             else if ( (*t-*_Htau > 5.9091) && (*t-*_Htau <= 7.0909) )
             {
                t_t = *t - *_Htau - 5.9091;
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
                        ka_index++; 
                    }
                    num *= num*num*num;
                    num += _Hea[i]*_Hea[i]*_Hea[i]*_Hea[i];
                    for ( int j = 0; j < _Hnkd[i]; j++)
                    {  
                        den += (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+24]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+24+1]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+24+2]+t_t*
                                    _Hcoeff[_Hkdvec[kd_index]*N_time_points*4+24+3])))*_Hkdval[kd_index];
                        kd_index++; 
                    }
                    den = den*den*den*den ;
                    den += (num + 1.0);
                    ydot[i] = _Hr0[i]*(num/den)-_Hd[i]*y[i];
                }
             }
             else if ( (*t-*_Htau > 7.0909) && (*t-*_Htau <= 8.2727) )
             {
                t_t = *t - *_Htau - 7.0909;
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
                        ka_index++; 
                    }
                    num *= num*num*num;
                    num += _Hea[i]*_Hea[i]*_Hea[i]*_Hea[i];
                    for ( int j = 0; j < _Hnkd[i]; j++)
                    {  
                        den += (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+28]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+28+1]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+28+2]+t_t*
                                    _Hcoeff[_Hkdvec[kd_index]*N_time_points*4+28+3])))*_Hkdval[kd_index];
                        kd_index++; 
                    }
                    den = den*den*den*den ;
                    den += (num + 1.0);
                    ydot[i] = _Hr0[i]*(num/den)-_Hd[i]*y[i];
                }
             }
             else if ( (*t-*_Htau > 8.2727) && (*t-*_Htau <= 9.4545) )
             {
                t_t = *t - *_Htau - 8.2727;
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
                        ka_index++; 
                    }
                    num *= num*num*num;
                    num += _Hea[i]*_Hea[i]*_Hea[i]*_Hea[i];
                    for ( int j = 0; j < _Hnkd[i]; j++)
                    {  
                        den += (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+32]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+32+1]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+32+2]+t_t*
                                    _Hcoeff[_Hkdvec[kd_index]*N_time_points*4+32+3])))*_Hkdval[kd_index];
                        kd_index++; 
                    }
                    den = den*den*den*den ;
                    den += (num + 1.0);
                    ydot[i] = _Hr0[i]*(num/den)-_Hd[i]*y[i];
                }
             }
             else if ( (*t-*_Htau > 9.4545) && (*t-*_Htau <=10.6364) )
             {
                t_t = *t - *_Htau - 9.4545;
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
                        ka_index++; 
                    }
                    num *= num*num*num;
                    num += _Hea[i]*_Hea[i]*_Hea[i]*_Hea[i];
                    for ( int j = 0; j < _Hnkd[i]; j++)
                    {  
                        den += (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+36]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+36+1]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+36+2]+t_t*
                                    _Hcoeff[_Hkdvec[kd_index]*N_time_points*4+36+3])))*_Hkdval[kd_index];
                        kd_index++; 
                    }
                    den = den*den*den*den ;
                    den += (num + 1.0);
                    ydot[i] = _Hr0[i]*(num/den)-_Hd[i]*y[i];
                }
             }
             else if ( (*t-*_Htau >10.6364) && (*t-*_Htau <=11.8182) )
             {
                t_t = *t - *_Htau -10.6364;
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
                        ka_index++; 
                    }
                    num *= num*num*num;
                    num += _Hea[i]*_Hea[i]*_Hea[i]*_Hea[i];
                    for ( int j = 0; j < _Hnkd[i]; j++)
                    {  
                        den += (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+40]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+40+1]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+40+2]+t_t*
                                    _Hcoeff[_Hkdvec[kd_index]*N_time_points*4+40+3])))*_Hkdval[kd_index];
                        kd_index++; 
                    }
                    den = den*den*den*den ;
                    den += (num + 1.0);
                    ydot[i] = _Hr0[i]*(num/den)-_Hd[i]*y[i];
                }
             }
             else if ( (*t-*_Htau >11.8182) && (*t-*_Htau <= 13.0  ) )
             {
                t_t = *t - *_Htau -11.8182;
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
                        ka_index++; 
                    }
                    num *= num*num*num;
                    num += _Hea[i]*_Hea[i]*_Hea[i]*_Hea[i];
                    for ( int j = 0; j < _Hnkd[i]; j++)
                    {  
                        den += (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+44]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+44+1]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+44+2]+t_t*
                                    _Hcoeff[_Hkdvec[kd_index]*N_time_points*4+44+3])))*_Hkdval[kd_index];
                        kd_index++; 
                    }
                    den = den*den*den*den ;
                    den += (num + 1.0);
                    ydot[i] = _Hr0[i]*(num/den)-_Hd[i]*y[i];
                }
             }
             else if (*t-*_Htau > 13.0  )
             {
                t_t = *t - *_Htau - 13.0;
                for ( int i = 0 ; i < N_gene; i++ )
                { 
                    num = 0.0;
                    den = 0.0;
                    for ( int j = 0; j < _Hnka[i]; j++)
                    {  
                        num += (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+48]+t_t*
                               (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+48+1]+t_t*
                               (_Hcoeff[_Hkavec[ka_index]*N_time_points*4+48+2]+t_t*
                                    _Hcoeff[_Hkavec[ka_index]*N_time_points*4+48+3])))*_Hkaval[ka_index];
                        ka_index++; 
                    }
                    num *= num*num*num;
                    num += _Hea[i]*_Hea[i]*_Hea[i]*_Hea[i];
                    for ( int j = 0; j < _Hnkd[i]; j++)
                    {  
                        den += (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+48]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+48+1]+t_t*
                               (_Hcoeff[_Hkdvec[kd_index]*N_time_points*4+48+2]+t_t*
                                    _Hcoeff[_Hkdvec[kd_index]*N_time_points*4+48+3])))*_Hkdval[kd_index];
                        kd_index++; 
                    }
                    den = den*den*den*den ;
                    den += (num + 1.0);
                    ydot[i] = _Hr0[i]*(num/den)-_Hd[i]*y[i];
                }
             }
	}
};*/

/*struct myJex
{
	void operator()(int *neq, value_type *t, value_type *y, int ml, int mu, value_type *pd, int nrowpd)
	{
		return;
	}
};*/
// Definition for a single gene base
struct myFex_single
{
        private:
        vector <double> _Hcoeff;
        //int *_Hindex;
        int *_Hnka, *_Hnkd;
        vector <int>  _Hkavec, _Hkdvec;
        value_type *_Htau;
        double *_Hr0, *_Hd, *_Hea;
        state_type _Hkaval, _Hkdval;
        double t_t, num, den, r0_i, d_i, ea_i;
        public:
        void set_tau ( const value_type tau ) 
        { 
	     _Htau = ( value_type *)malloc(sizeof(value_type));
	     *_Htau = tau;
        };
        void set_r0 ( const double &r0 ) 
        { 
	     _Hr0 = ( double *)malloc(sizeof(double));
             r0_i = r0;
             *_Hr0 = r0;
        };
        void get_r0 ( double &r0 ) { r0 = r0_i; };
        void set_d ( const double &d ) 
        { 
	     _Hd = ( double *)malloc(sizeof(double));
             d_i = d;
             *_Hd = d;
        };
        void get_d ( double &d ) { d = d_i; };
        void set_ea ( const double &ea ) 
        { 
	     _Hea = ( double *)malloc(sizeof(double));
             ea_i = ea;
             *_Hea = ea;
        };
        void get_ea ( double &ea ) { ea = ea_i; };
        void set_n_ka ( const int &nka ) 
        { 
	     _Hnka = ( int *)malloc(sizeof(int));
             *_Hnka = nka; 
        }; 
        void set_n_kd ( const int &nkd ) 
        { 
	     _Hnkd = ( int *)malloc(sizeof(int));
             *_Hnkd = nkd; 
        }; 
        void set_ka_vec ( const vector <int> kavec ) { _Hkavec = kavec;}; 
        void set_kd_vec ( const vector <int> kdvec ) { _Hkdvec = kdvec;}; 
        void set_ka_val ( const vector <double> &kaval ) { _Hkaval = kaval;}; 
        void get_ka_val ( vector <double> &kaval ) { kaval = _Hkaval;}; 
        void set_kd_val ( const vector <double> &kdval ) { _Hkdval = kdval;}; 
        void get_kd_val ( vector <double> &kdval ) { kdval = _Hkdval;}; 
        void set_coeff ( const vector<double> &c0_v ) { _Hcoeff = c0_v; };
        void clear_mem()
        {
             free(_Htau);
             free(_Hr0);
             free(_Hea);
             free(_Hd);
             free(_Hnka);
             free(_Hnkd);
             vector<int>().swap(_Hkavec);
             vector<int>().swap(_Hkdvec);
             vector<double>().swap(_Hkaval);
             vector<double>().swap(_Hkdval);
             vector<double>().swap(_Hcoeff);
        };
	void operator()(int *neq, double *t, double *y, double *ydot)
	{
             /**t_t = 0.0;
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
             }*/
             t_t = 0.0;
             if ( *t - *_Htau <= 0.0 )
             {
                num = 0.0;
                den = 0.0;
                for ( int j = 0; j < *_Hnka; j++)
                {
                    //cout << "calculating numerator " << endl;
                    //cout << "Accessing element " << _Hkavec[j]*48 << " of coeff vector." << endl;
                    //cout << "Value of kaval is: " << _Hkaval[j] << endl;
                    num += _Hcoeff[_Hkavec[j]*N_time_points*4]*_Hkaval[j];
                } 
                //cout << " Value of time is  " << *t << endl;
                //cout << "Value of num(1) is: " << num << endl;
                num *= num*num*num;
                //cout << "Value of num(2) is: " << num << endl;
                num += ((*_Hea)*(*_Hea)*(*_Hea)*(*_Hea));
                //cout << "Value of num(3) is: " << num << endl;
                //cout << "Value of ea is: " << *_Hea << endl;
                for ( int j = 0; j < *_Hnkd; j++)
                    den += _Hcoeff[_Hkdvec[j]*N_time_points*4]*_Hkdval[j];
                //cout << "Value of den(1) is: " << den << endl;
                den = den*den*den*den;
                //cout << "Value of den(2) is: " << den << endl;
                den += (num + 1.0);
                //cout << "Value of den(3) is: " << den << endl;
                ydot[0] = ((*_Hr0)*(num/den))-((*_Hd)*y[0]);
                //cout << "Value of r0 is: " << *_Hr0 << endl;
                //cout << "Value of d is: " << *_Hd << endl;
                //cout << "Value of ydot[0] is: " << ydot[0] << endl;
                //cout << "Value of y[0] is: " << y[0] << endl;
	     }
             else if ( (*t-*_Htau > 0.0) && (*t-*_Htau <= 1.0) )
             {
                t_t = *t - *_Htau;
                num = 0.0;
                den = 0.0;
                for ( int j = 0; j < *_Hnka; j++)
                    num += (_Hcoeff[_Hkavec[j]*N_time_points*4+4]+t_t*
                           (_Hcoeff[_Hkavec[j]*N_time_points*4+4+1]+t_t*
                           (_Hcoeff[_Hkavec[j]*N_time_points*4+4+2]+t_t*
                                _Hcoeff[_Hkavec[j]*N_time_points*4+4+3])))*_Hkaval[j];
                num *= num*num*num;
                num += (*_Hea)*(*_Hea)*(*_Hea)*(*_Hea);
                for ( int j = 0; j < *_Hnkd; j++)
                    den += (_Hcoeff[_Hkdvec[j]*N_time_points*4+4]+t_t*
                           (_Hcoeff[_Hkdvec[j]*N_time_points*4+4+1]+t_t*
                           (_Hcoeff[_Hkdvec[j]*N_time_points*4+4+2]+t_t*
                                _Hcoeff[_Hkdvec[j]*N_time_points*4+4+3])))*_Hkdval[j];
                    //den += (get_coeff0(_Hkdvec[j],1)+t_t*
                    //       (get_coeff1(_Hkdvec[j],1)+t_t*
                    //       (get_coeff2(_Hkdvec[j],1)+t_t*
                    //            get_coeff3(_Hkdvec[j],1))))*_Hkdval[j];
                den = den*den*den*den ;
                den += (num + 1.0);
                ydot[0] = (*_Hr0)*(num/den)-(*_Hd)*y[0];
             }
             else if ( (*t-*_Htau > 1.0) && (*t-*_Htau <= 2.0) )
             {
                t_t = *t - *_Htau - 1.0;
                num = 0.0;
                den = 0.0;
                for ( int j = 0; j < *_Hnka; j++)
                    num += (_Hcoeff[_Hkavec[j]*N_time_points*4+8]+t_t*
                           (_Hcoeff[_Hkavec[j]*N_time_points*4+8+1]+t_t*
                           (_Hcoeff[_Hkavec[j]*N_time_points*4+8+2]+t_t*
                                _Hcoeff[_Hkavec[j]*N_time_points*4+8+3])))*_Hkaval[j];
                num *= num*num*num;
                num += (*_Hea)*(*_Hea)*(*_Hea)*(*_Hea);
                for ( int j = 0; j < *_Hnkd; j++)
                    den += (_Hcoeff[_Hkdvec[j]*N_time_points*4+8]+t_t*
                           (_Hcoeff[_Hkdvec[j]*N_time_points*4+8+1]+t_t*
                           (_Hcoeff[_Hkdvec[j]*N_time_points*4+8+2]+t_t*
                                _Hcoeff[_Hkdvec[j]*N_time_points*4+8+3])))*_Hkdval[j];
                den = den*den*den*den ;
                den += (num + 1.0);
                ydot[0] = (*_Hr0)*(num/den)-(*_Hd)*y[0];
             }
             else if ( (*t-*_Htau > 2.0) && (*t-*_Htau <= 3.0) )
             {
                t_t = *t - *_Htau - 2.0;
                num = 0.0;
                den = 0.0;
                for ( int j = 0; j < *_Hnka; j++)
                    num += (_Hcoeff[_Hkavec[j]*N_time_points*4+12]+t_t*
                           (_Hcoeff[_Hkavec[j]*N_time_points*4+12+1]+t_t*
                           (_Hcoeff[_Hkavec[j]*N_time_points*4+12+2]+t_t*
                                _Hcoeff[_Hkavec[j]*N_time_points*4+12+3])))*_Hkaval[j];
                num *= num*num*num;
                num += (*_Hea)*(*_Hea)*(*_Hea)*(*_Hea);
                for ( int j = 0; j < *_Hnkd; j++)
                    den += (_Hcoeff[_Hkdvec[j]*N_time_points*4+12]+t_t*
                           (_Hcoeff[_Hkdvec[j]*N_time_points*4+12+1]+t_t*
                           (_Hcoeff[_Hkdvec[j]*N_time_points*4+12+2]+t_t*
                                _Hcoeff[_Hkdvec[j]*N_time_points*4+12+3])))*_Hkdval[j];
                den = den*den*den*den ;
                den += (num + 1.0);
                ydot[0] = (*_Hr0)*(num/den)-(*_Hd)*y[0];
             }
             else if ( (*t-*_Htau > 3.0) && (*t-*_Htau <= 4.0) )
             {
                t_t = *t - *_Htau - 3.0;
                num = 0.0;
                den = 0.0;
                for ( int j = 0; j < *_Hnka; j++)
                    num += (_Hcoeff[_Hkavec[j]*N_time_points*4+16]+t_t*
                           (_Hcoeff[_Hkavec[j]*N_time_points*4+16+1]+t_t*
                           (_Hcoeff[_Hkavec[j]*N_time_points*4+16+2]+t_t*
                                _Hcoeff[_Hkavec[j]*N_time_points*4+16+3])))*_Hkaval[j];
                num *= num*num*num;
                num += (*_Hea)*(*_Hea)*(*_Hea)*(*_Hea);
                for ( int j = 0; j < *_Hnkd; j++)
                    den += (_Hcoeff[_Hkdvec[j]*N_time_points*4+16]+t_t*
                           (_Hcoeff[_Hkdvec[j]*N_time_points*4+16+1]+t_t*
                           (_Hcoeff[_Hkdvec[j]*N_time_points*4+16+2]+t_t*
                                _Hcoeff[_Hkdvec[j]*N_time_points*4+16+3])))*_Hkdval[j];
                den = den*den*den*den ;
                den += (num + 1.0);
                ydot[0] = (*_Hr0)*(num/den)-(*_Hd)*y[0];
             }
             else if ( (*t-*_Htau > 4.0) && (*t-*_Htau <= 5.0) )
             {
                t_t = *t - *_Htau - 4.0;
                num = 0.0;
                den = 0.0;
                for ( int j = 0; j < *_Hnka; j++)
                    num += (_Hcoeff[_Hkavec[j]*N_time_points*4+20]+t_t*
                           (_Hcoeff[_Hkavec[j]*N_time_points*4+20+1]+t_t*
                           (_Hcoeff[_Hkavec[j]*N_time_points*4+20+2]+t_t*
                                _Hcoeff[_Hkavec[j]*N_time_points*4+20+3])))*_Hkaval[j];
                num *= num*num*num;
                num += (*_Hea)*(*_Hea)*(*_Hea)*(*_Hea);
                for ( int j = 0; j < *_Hnkd; j++)
                    den += (_Hcoeff[_Hkdvec[j]*N_time_points*4+20]+t_t*
                           (_Hcoeff[_Hkdvec[j]*N_time_points*4+20+1]+t_t*
                           (_Hcoeff[_Hkdvec[j]*N_time_points*4+20+2]+t_t*
                                _Hcoeff[_Hkdvec[j]*N_time_points*4+20+3])))*_Hkdval[j];
                den = den*den*den*den ;
                den += (num + 1.0);
                ydot[0] = (*_Hr0)*(num/den)-(*_Hd)*y[0];
             }
             else if ( (*t-*_Htau > 5.0) && (*t-*_Htau <= 6.0) )
             {
                t_t = *t - *_Htau - 5.0;
                num = 0.0;
                den = 0.0;
                for ( int j = 0; j < *_Hnka; j++)
                    num += (_Hcoeff[_Hkavec[j]*N_time_points*4+24]+t_t*
                           (_Hcoeff[_Hkavec[j]*N_time_points*4+24+1]+t_t*
                           (_Hcoeff[_Hkavec[j]*N_time_points*4+24+2]+t_t*
                                _Hcoeff[_Hkavec[j]*N_time_points*4+24+3])))*_Hkaval[j];
                num *= num*num*num;
                num += (*_Hea)*(*_Hea)*(*_Hea)*(*_Hea);
                for ( int j = 0; j < *_Hnkd; j++)
                    den += (_Hcoeff[_Hkdvec[j]*N_time_points*4+24]+t_t*
                           (_Hcoeff[_Hkdvec[j]*N_time_points*4+24+1]+t_t*
                           (_Hcoeff[_Hkdvec[j]*N_time_points*4+24+2]+t_t*
                                _Hcoeff[_Hkdvec[j]*N_time_points*4+24+3])))*_Hkdval[j];
                den = den*den*den*den ;
                den += (num + 1.0);
                ydot[0] = (*_Hr0)*(num/den)-(*_Hd)*y[0];
             }
             else if ( (*t-*_Htau > 6.0) && (*t-*_Htau <= 7.0) )
             {
                t_t = *t - *_Htau - 6.0;
                num = 0.0;
                den = 0.0;
                for ( int j = 0; j < *_Hnka; j++)
                    num += (_Hcoeff[_Hkavec[j]*N_time_points*4+28]+t_t*
                           (_Hcoeff[_Hkavec[j]*N_time_points*4+28+1]+t_t*
                           (_Hcoeff[_Hkavec[j]*N_time_points*4+28+2]+t_t*
                                _Hcoeff[_Hkavec[j]*N_time_points*4+28+3])))*_Hkaval[j];
                num *= num*num*num;
                num += (*_Hea)*(*_Hea)*(*_Hea)*(*_Hea);
                for ( int j = 0; j < *_Hnkd; j++)
                    den += (_Hcoeff[_Hkdvec[j]*N_time_points*4+28]+t_t*
                           (_Hcoeff[_Hkdvec[j]*N_time_points*4+28+1]+t_t*
                           (_Hcoeff[_Hkdvec[j]*N_time_points*4+28+2]+t_t*
                                _Hcoeff[_Hkdvec[j]*N_time_points*4+28+3])))*_Hkdval[j];
                den = den*den*den*den ;
                den += (num + 1.0);
                ydot[0] = (*_Hr0)*(num/den)-(*_Hd)*y[0];
             }
             else if ( (*t-*_Htau > 7.0) && (*t-*_Htau <= 8.0) )
             {
                t_t = *t - *_Htau - 7.0;
                num = 0.0;
                den = 0.0;
                for ( int j = 0; j < *_Hnka; j++)
                    num += (_Hcoeff[_Hkavec[j]*N_time_points*4+32]+t_t*
                           (_Hcoeff[_Hkavec[j]*N_time_points*4+32+1]+t_t*
                           (_Hcoeff[_Hkavec[j]*N_time_points*4+32+2]+t_t*
                                _Hcoeff[_Hkavec[j]*N_time_points*4+32+3])))*_Hkaval[j];
                num *= num*num*num;
                num += (*_Hea)*(*_Hea)*(*_Hea)*(*_Hea);
                for ( int j = 0; j < *_Hnkd; j++)
                    den += (_Hcoeff[_Hkdvec[j]*N_time_points*4+32]+t_t*
                           (_Hcoeff[_Hkdvec[j]*N_time_points*4+32+1]+t_t*
                           (_Hcoeff[_Hkdvec[j]*N_time_points*4+32+2]+t_t*
                                _Hcoeff[_Hkdvec[j]*N_time_points*4+32+3])))*_Hkdval[j];
                den = den*den*den*den ;
                den += (num + 1.0);
                ydot[0] = (*_Hr0)*(num/den)-(*_Hd)*y[0];
             }
             else if ( (*t-*_Htau > 8.0) && (*t-*_Htau <=9.0) )
             {
                t_t = *t - *_Htau - 8.0;
                num = 0.0;
                den = 0.0;
                for ( int j = 0; j < *_Hnka; j++)
                    num += (_Hcoeff[_Hkavec[j]*N_time_points*4+36]+t_t*
                           (_Hcoeff[_Hkavec[j]*N_time_points*4+36+1]+t_t*
                           (_Hcoeff[_Hkavec[j]*N_time_points*4+36+2]+t_t*
                                _Hcoeff[_Hkavec[j]*N_time_points*4+36+3])))*_Hkaval[j];
                num *= num*num*num;
                num += (*_Hea)*(*_Hea)*(*_Hea)*(*_Hea);
                for ( int j = 0; j < *_Hnkd; j++)
                    den += (_Hcoeff[_Hkdvec[j]*N_time_points*4+36]+t_t*
                           (_Hcoeff[_Hkdvec[j]*N_time_points*4+36+1]+t_t*
                           (_Hcoeff[_Hkdvec[j]*N_time_points*4+36+2]+t_t*
                                _Hcoeff[_Hkdvec[j]*N_time_points*4+36+3])))*_Hkdval[j];
                den = den*den*den*den ;
                den += (num + 1.0);
                ydot[0] = (*_Hr0)*(num/den)-(*_Hd)*y[0];
             }
             else if ( (*t-*_Htau >9.0) && (*t-*_Htau <=10.0) )
             {
                t_t = *t - *_Htau - 9.0;
                num = 0.0;
                den = 0.0;
                for ( int j = 0; j < *_Hnka; j++)
                    num += (_Hcoeff[_Hkavec[j]*N_time_points*4+40]+t_t*
                           (_Hcoeff[_Hkavec[j]*N_time_points*4+40+1]+t_t*
                           (_Hcoeff[_Hkavec[j]*N_time_points*4+40+2]+t_t*
                                _Hcoeff[_Hkavec[j]*N_time_points*4+40+3])))*_Hkaval[j];
                num *= num*num*num;
                num += (*_Hea)*(*_Hea)*(*_Hea)*(*_Hea);
                for ( int j = 0; j < *_Hnkd; j++)
                    den += (_Hcoeff[_Hkdvec[j]*N_time_points*4+40]+t_t*
                           (_Hcoeff[_Hkdvec[j]*N_time_points*4+40+1]+t_t*
                           (_Hcoeff[_Hkdvec[j]*N_time_points*4+40+2]+t_t*
                                _Hcoeff[_Hkdvec[j]*N_time_points*4+40+3])))*_Hkdval[j];
                den = den*den*den*den ;
                den += (num + 1.0);
                ydot[0] = (*_Hr0)*(num/den)-(*_Hd)*y[0];
             }
             else if ( (*t-*_Htau >10.0) && (*t-*_Htau <= 11.0  ) )
             {
                t_t = *t - *_Htau -10.0;
                num = 0.0;
                den = 0.0;
                for ( int j = 0; j < *_Hnka; j++)
                    num += (_Hcoeff[_Hkavec[j]*N_time_points*4+44]+t_t*
                           (_Hcoeff[_Hkavec[j]*N_time_points*4+44+1]+t_t*
                           (_Hcoeff[_Hkavec[j]*N_time_points*4+44+2]+t_t*
                                _Hcoeff[_Hkavec[j]*N_time_points*4+44+3])))*_Hkaval[j];
                num *= num*num*num;
                num += (*_Hea)*(*_Hea)*(*_Hea)*(*_Hea);
                for ( int j = 0; j < *_Hnkd; j++)
                    den += (_Hcoeff[_Hkdvec[j]*N_time_points*4+44]+t_t*
                           (_Hcoeff[_Hkdvec[j]*N_time_points*4+44+1]+t_t*
                           (_Hcoeff[_Hkdvec[j]*N_time_points*4+44+2]+t_t*
                                _Hcoeff[_Hkdvec[j]*N_time_points*4+44+3])))*_Hkdval[j];
                den = den*den*den*den ;
                den += (num + 1.0);
                ydot[0] = (*_Hr0)*(num/den)-(*_Hd)*y[0];
             }
             else if (*t-*_Htau > 11.0  )
             {
                t_t = *t - *_Htau - 11.0;
                num = 0.0;
                den = 0.0;
                for ( int j = 0; j < *_Hnka; j++)
                    num += (_Hcoeff[_Hkavec[j]*N_time_points*4+48]+t_t*
                           (_Hcoeff[_Hkavec[j]*N_time_points*4+48+1]+t_t*
                           (_Hcoeff[_Hkavec[j]*N_time_points*4+48+2]+t_t*
                                _Hcoeff[_Hkavec[j]*N_time_points*4+48+3])))*_Hkaval[j];
                num *= num*num*num;
                num += (*_Hea)*(*_Hea)*(*_Hea)*(*_Hea);
                for ( int j = 0; j < *_Hnkd; j++)
                    den += (_Hcoeff[_Hkdvec[j]*N_time_points*4+48]+t_t*
                           (_Hcoeff[_Hkdvec[j]*N_time_points*4+48+1]+t_t*
                           (_Hcoeff[_Hkdvec[j]*N_time_points*4+48+2]+t_t*
                                _Hcoeff[_Hkdvec[j]*N_time_points*4+48+3])))*_Hkdval[j];
                den = den*den*den*den ;
                den += (num + 1.0);
                ydot[0] = (*_Hr0)*(num/den)-(*_Hd)*y[0];
             }
	}
};

struct myJex_single
{
	void operator()(int *neq, value_type *t, value_type *y, int ml, int mu, value_type *pd, int nrowpd)
	{
		return;
	}
};
#endif
