/*
 *  cuLsoda.h
 *
 *	File Notes:	This file is a conversion of the double precision Livermore Solver for
 *	Ordinary Differential Equations with automatic switching for stiff and non-stiff
 *	problems (DLSODA)
 */

/**  barf  [ba:rf]  2.  "He suggested using FORTRAN, and everybody barfed."
 
 - From The Shogakukan DICTIONARY OF NEW ENGLISH (Second edition) */
# ifdef __CUDACC__
# define HOST __host__
# define DEVICE __device__
# else
# define HOST 
# define DEVICE
# endif 
#ifndef CULSODA_FEX_H
#define CULSODA_FEX_H
#include "../common.h"

#include <thrust/device_vector.h>

using namespace std;

#define Fex_and_Jex_definition
// The model class is basically Fex
struct myFex
{
        private:

        double *_Dcoeff,*r0_ptr, *d_ptr, *ea_ptr, *_Dkaval, *_Dkdval;
        int *_Dnka, *_Dkavec, *_Dkastart, *_Dnkd, *_Dkdvec, *_Dkdstart;
        double t_t, num, den;
//    thrust::device_vector<double> t_d;
    double* t_d_array;
    int t_d_size;
        public:
    void set_t_d(const double * t_d_in, int t_d_size_in) {
        size_t size = sizeof(double) * t_d_size_in;
        cudaMalloc((void**)&t_d_array, size);
        cudaMemcpy(t_d_array, t_d_in, size, cudaMemcpyHostToDevice);
        t_d_size = t_d_size_in;
    }
    void set_t_d_free() {
        cudaFree(t_d_array);
    }
        void set_r0 ( const double *r0 ) 
        {
	     cudaMalloc((void**)&r0_ptr,sizeof(double)*probSize);
	     cudaMemcpy(r0_ptr,r0,sizeof(double)*probSize,cudaMemcpyHostToDevice);
        };
        void get_r0_ptr ( thrust::device_ptr<double> &r0ptr ) 
        {
             r0ptr = thrust::device_pointer_cast(r0_ptr);
        };
        /*void set_r0_ptr ( thrust::device_ptr<double> &r0_opt ) 
        {
             thrust::device_ptr<double> r0ptr = thrust::device_pointer_cast(r0_ptr);
             thrust::copy(r0_opt,r0_opt+probSize,r0ptr);
        };*/
        void set_r0_free () { cudaFree(r0_ptr); };
        void set_d ( const double *d ) 
        {
	     cudaMalloc((void**)&d_ptr,sizeof(double)*probSize);
	     cudaMemcpy(d_ptr,d,sizeof(double)*probSize,cudaMemcpyHostToDevice);
        };
        void get_d_ptr ( thrust::device_ptr<double> &dptr ) 
        {
             dptr = thrust::device_pointer_cast(d_ptr);
        };
        /*void set_d_ptr ( thrust::device_ptr<double> &d_opt ) 
        {
             thrust::device_ptr<double> dptr = thrust::device_pointer_cast(d_ptr);
             thrust::copy(d_opt,d_opt+probSize,dptr);
        };*/
        void set_d_free () { cudaFree(d_ptr); };
        /*void set_d_array ( const double d[] ) 
        {
	     cudaMalloc((void**)&d_ptr,sizeof(double)*probSize);
	     cudaMemcpy(d_ptr,d,sizeof(double)*probSize,cudaMemcpyHostToDevice);
        };
        void set_d_free () { cudaFree(d_ptr); };*/
        void set_ea ( const double *ea ) 
        {
	     cudaMalloc((void**)&ea_ptr,sizeof(double)*probSize);
	     cudaMemcpy(ea_ptr,ea,sizeof(double)*probSize,cudaMemcpyHostToDevice);
        };
        void get_ea_ptr ( thrust::device_ptr<double> &eaptr ) 
        {
             eaptr = thrust::device_pointer_cast(ea_ptr);
        };
        /*void set_ea_ptr ( thrust::device_ptr<double> &ea_opt ) 
        {
             thrust::device_ptr<double> eaptr = thrust::device_pointer_cast(ea_ptr);
             thrust::copy(ea_opt,ea_opt+probSize,eaptr);
        };*/
        void set_ea_free () { cudaFree(ea_ptr); };
        void set_n_ka ( const int *nka ) 
        {
	     cudaMalloc((void**)&_Dnka,sizeof(int)*probSize);
	     cudaMemcpy(_Dnka,nka,sizeof(int)*probSize,cudaMemcpyHostToDevice);
        };
        void set_n_ka_free () { cudaFree(_Dnka); }; 
        /*void set_n_ka ( const int nka[] ) 
        {
	     cudaMalloc((void**)&_Dnka,sizeof(int)*probSize);
	     cudaMemcpy(_Dnka,nka,sizeof(int)*probSize,cudaMemcpyHostToDevice);
        };
        void set_n_ka_free () { cudaFree(_Dnka); }; */
        void set_ka_start ( const int *kastart ) 
        {
	     cudaMalloc((void**)&_Dkastart,sizeof(int)*probSize);
	     cudaMemcpy(_Dkastart,kastart,sizeof(int)*probSize,cudaMemcpyHostToDevice);
        };
        void set_ka_start_free () { cudaFree(_Dkastart); }; 
        /*void set_ka_start ( const int kastart[] ) 
        {
	     cudaMalloc((void**)&_Dkastart,sizeof(int)*probSize);
	     cudaMemcpy(_Dkastart,kastart,sizeof(int)*probSize,cudaMemcpyHostToDevice);
        };
        void set_ka_start_free () { cudaFree(_Dkastart); };*/ 
        void set_ka_vec ( const int *kavec, const int size ) 
        {
	     cudaMalloc((void**)&_Dkavec,sizeof(int)*size);
	     cudaMemcpy(_Dkavec,kavec,sizeof(int)*size,cudaMemcpyHostToDevice);
        };
        void get_kavec_vec ( thrust::host_vector<int> &kavec, const int size ) 
        {
             thrust::device_ptr<int> kaptr = thrust::device_pointer_cast(_Dkavec);
             thrust::copy(kaptr,kaptr+size,kavec.begin());
        };  
        void set_ka_vec_free () { cudaFree(_Dkavec); }; 
        /*void set_ka_vec ( const int kavec[], const int size ) 
        {
	     cudaMalloc((void**)&_Dkavec,sizeof(int)*size);
	     cudaMemcpy(_Dkavec,kavec,sizeof(int)*size,cudaMemcpyHostToDevice);
        };
        void set_ka_vec_free () { cudaFree(_Dkavec); };*/ 
        void set_ka_val ( const double *ka, const int size ) 
        {
	     cudaMalloc((void**)&_Dkaval,sizeof(double)*size);
	     cudaMemcpy(_Dkaval,ka,sizeof(double)*size,cudaMemcpyHostToDevice);
        };
        void get_kaval_ptr ( thrust::device_ptr<double> &kaptr ) 
        {
             kaptr = thrust::device_pointer_cast(_Dkaval);
        };
        void get_kaval_vec ( host_type &kavec, const int size ) 
        {
             thrust::device_ptr<double> kaptr = thrust::device_pointer_cast(_Dkaval);
             thrust::copy(kaptr,kaptr+size,kavec.begin());
        };  
        /*void set_ka_ptr ( thrust::device_ptr<double> &ka_opt, const int ka_size ) 
        {
             thrust::device_ptr<double> kaptr = thrust::device_pointer_cast(_Dkaval);
             thrust::copy(ka_opt,ka_opt+ka_size,kaptr);
        };*/
        void set_ka_val_free () { cudaFree(_Dkaval); }; 
        /*void set_ka_val ( const double ka[], const int size ) 
        {
	     cudaMalloc((void**)&_Dkaval,sizeof(double)*size);
	     cudaMemcpy(_Dkaval,ka,sizeof(double)*size,cudaMemcpyHostToDevice);
        };
        void set_ka_val_free () { cudaFree(_Dkaval); };*/
        void set_n_kd ( const int *nkd ) 
        {
	     cudaMalloc((void**)&_Dnkd,sizeof(int)*probSize);
	     cudaMemcpy(_Dnkd,nkd,sizeof(int)*probSize,cudaMemcpyHostToDevice);
        };
        void set_n_kd_free () { cudaFree(_Dnkd); }; 
        /*void set_n_kd ( const int nkd[] ) 
        {
	     cudaMalloc((void**)&_Dnkd,sizeof(int)*probSize);
	     cudaMemcpy(_Dnkd,nkd,sizeof(int)*probSize,cudaMemcpyHostToDevice);
        };
        void set_n_kd_free () { cudaFree(_Dnkd); };*/ 
        void set_kd_start ( const int *kdstart ) 
        {
	     cudaMalloc((void**)&_Dkdstart,sizeof(int)*probSize);
	     cudaMemcpy(_Dkdstart,kdstart,sizeof(int)*probSize,cudaMemcpyHostToDevice);
        };
        void set_kd_start_free () { cudaFree(_Dkdstart); }; 
        /*void set_kd_start ( const int kdstart[] ) 
        {
	     cudaMalloc((void**)&_Dkdstart,sizeof(int)*probSize);
	     cudaMemcpy(_Dkdstart,kdstart,sizeof(int)*probSize,cudaMemcpyHostToDevice);
        };
        void set_kd_start_free () { cudaFree(_Dkdstart); };*/ 
        void set_kd_vec ( const int *kdvec, const int size ) 
        {
	     cudaMalloc((void**)&_Dkdvec,sizeof(int)*size);
	     cudaMemcpy(_Dkdvec,kdvec,sizeof(int)*size,cudaMemcpyHostToDevice);
        };
        void set_kd_vec_free () { cudaFree(_Dkdvec); }; 
        void get_kdvec_vec ( thrust::host_vector<int> &kdvec, const int size ) 
        {
             thrust::device_ptr<int> kdptr = thrust::device_pointer_cast(_Dkdvec);
             thrust::copy(kdptr,kdptr+size,kdvec.begin());
        };
        /*void set_kd_vec ( const int kdvec[], const int size ) 
        {
	     cudaMalloc((void**)&_Dkdvec,sizeof(int)*size);
	     cudaMemcpy(_Dkdvec,kdvec,sizeof(int)*size,cudaMemcpyHostToDevice);
        };
        void set_kd_vec_free () { cudaFree(_Dkdvec); };*/ 
        void set_kd_val ( const double *kd, const int size ) 
        {
	     cudaMalloc((void**)&_Dkdval,sizeof(double)*size);
	     cudaMemcpy(_Dkdval,kd,sizeof(double)*size,cudaMemcpyHostToDevice);
        };
        void get_kdval_ptr ( thrust::device_ptr<double> &kdptr ) 
        {
             kdptr = thrust::device_pointer_cast(_Dkdval);
        };
        void get_kdval_vec ( host_type &kdvec, const int size ) 
        {
             thrust::device_ptr<double> kdptr = thrust::device_pointer_cast(_Dkdval);
             thrust::copy(kdptr,kdptr+size,kdvec.begin());
        };
        /*void set_kd_ptr ( thrust::device_ptr<double> &kd_opt, const int kd_size ) 
        {
             thrust::device_ptr<double> kdptr = thrust::device_pointer_cast(_Dkdval);
             thrust::copy(kd_opt,kd_opt+kd_size,kdptr);
        };*/
        void set_kd_val_free () { cudaFree(_Dkdval); }; 
        /*void set_kd_val ( const double kd[], const int size ) 
        {
	     cudaMalloc((void**)&_Dkdval,sizeof(double)*size);
	     cudaMemcpy(_Dkdval,kd,sizeof(double)*size,cudaMemcpyHostToDevice);
        };
        void set_kd_val_free () { cudaFree(_Dkdval); };*/ 
        void set_coeff ( const double c0_v[], const int size) 
        {
	     cudaMalloc((void**)&_Dcoeff,sizeof(double)*size);
	     cudaMemcpy(_Dcoeff,c0_v,sizeof(double)*size,cudaMemcpyHostToDevice);
        };
        void set_coeff_free () { cudaFree(_Dcoeff); }; 
        HOST DEVICE
        double get_co0 (int i, int j) const {return _Dcoeff[i*N_time_points*4+j*4];}; 
        HOST DEVICE
        double get_co1 (int i, int j) const {return _Dcoeff[i*N_time_points*4+j*4+1];};
        HOST DEVICE
        double get_co2 (int i, int j) const {return _Dcoeff[i*N_time_points*4+j*4+2];};
        HOST DEVICE
        double get_co3 (int i, int j) const {return _Dcoeff[i*N_time_points*4+j*4+3];};

    __device__ int t_d_lower_bound(double t_in) {
        int low, high, mid;
        low = 0;
        high = t_d_size - 1;
        while(low <= high) {
            mid = (low + high) / 2;
            if (t_in < t_d_array[mid]) {
                high = mid - 1;
            } else if (t_in > t_d_array[mid]) {
                low = mid + 1;
            } else {
                return mid;
            }
        }
        return high;
    }
    
	__device__ void operator()(int *neq, double *t, double *y, double *ydot)
	{



             t_t = 0.0;
	     int thread_id = threadIdx.x + blockIdx.x * blockDim.x;
             int ka_index = 0;
             int kd_index = 0;
             num = 0.0;
             den = 0.0; 
             
             if ( *t - 1.0 <= 0.0 )
             {
                while ( ka_index < _Dnka[thread_id] )
                { 
                       num += (get_co0(_Dkavec[_Dkastart[thread_id]+ka_index],0) * _Dkaval[_Dkastart[thread_id]+ka_index]);
                       ka_index++;
                }
                num *= num*num*num;
                num += ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id];
                while ( kd_index < _Dnkd[thread_id] )
                { 
                       den += (get_co0(_Dkdvec[_Dkdstart[thread_id]+kd_index],0) * _Dkdval[_Dkdstart[thread_id]+kd_index]);
                       kd_index++;
                }
                den *= den*den*den; 
                den += (num + 1.0);
                ydot[0] = r0_ptr[thread_id]*(num/den)-d_ptr[thread_id]*y[0];
	     }
             else {
                 int lower_index = t_d_lower_bound(*t);
                 t_t = *t - 1.0;
                while ( ka_index < _Dnka[thread_id] )
                { 
                       num += ((get_co0(_Dkavec[_Dkastart[thread_id]+ka_index],lower_index) + 
                              (t_t)*(get_co1(_Dkavec[_Dkastart[thread_id]+ka_index],lower_index)+
                              ((t_t)*(get_co2(_Dkavec[_Dkastart[thread_id]+ka_index],lower_index)+get_co3(_Dkavec[_Dkastart[thread_id]+ka_index],lower_index)*(t_t)))))
                              *_Dkaval[_Dkastart[thread_id]+ka_index]); 
                       ka_index++;
                }
                num *= num*num*num;
                num += ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id];
                while ( kd_index < _Dnkd[thread_id] )
                { 
                       den += ((get_co0(_Dkdvec[_Dkdstart[thread_id]+kd_index],lower_index) + 
                              (t_t)*(get_co1(_Dkdvec[_Dkdstart[thread_id]+kd_index],lower_index)+
                              ((t_t)*(get_co2(_Dkdvec[_Dkdstart[thread_id]+kd_index],lower_index)+get_co3(_Dkdvec[_Dkdstart[thread_id]+kd_index],lower_index)*(t_t)))))
                              *_Dkdval[_Dkdstart[thread_id]+kd_index]); 
                       kd_index++;
                }
                den *= den*den*den; 
                den += (num + 1.0);
                ydot[0] = r0_ptr[thread_id]*(num/den)-d_ptr[thread_id]*y[0];
             }
             
             /*
             if ( *t - 1.0 <= 0.0 )
             {
                while ( ka_index < _Dnka[thread_id] )
                { 
                       num += (get_co0(_Dkavec[_Dkastart[thread_id]+ka_index],0) * _Dkaval[_Dkastart[thread_id]+ka_index]);
                       ka_index++;
                }
                num *= num*num*num;
                num += ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id];
                while ( kd_index < _Dnkd[thread_id] )
                { 
                       den += (get_co0(_Dkdvec[_Dkdstart[thread_id]+kd_index],0) * _Dkdval[_Dkdstart[thread_id]+kd_index]);
                       kd_index++;
                }
                den *= den*den*den; 
                den += (num + 1.0);
                ydot[0] = r0_ptr[thread_id]*(num/den)-d_ptr[thread_id]*y[0];
	     }
             else if ( (*t-1.0 > 0.0) && (*t-1.0 <= 1.0) )
             {
                t_t = *t - 1.0;
                while ( ka_index < _Dnka[thread_id] )
                { 
                       num += ((get_co0(_Dkavec[_Dkastart[thread_id]+ka_index],1) + 
                              (t_t)*(get_co1(_Dkavec[_Dkastart[thread_id]+ka_index],1)+
                              ((t_t)*(get_co2(_Dkavec[_Dkastart[thread_id]+ka_index],1)+get_co3(_Dkavec[_Dkastart[thread_id]+ka_index],1)*(t_t)))))
                              *_Dkaval[_Dkastart[thread_id]+ka_index]); 
                       ka_index++;
                }
                num *= num*num*num;
                num += ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id];
                while ( kd_index < _Dnkd[thread_id] )
                { 
                       den += ((get_co0(_Dkdvec[_Dkdstart[thread_id]+kd_index],1) + 
                              (t_t)*(get_co1(_Dkdvec[_Dkdstart[thread_id]+kd_index],1)+
                              ((t_t)*(get_co2(_Dkdvec[_Dkdstart[thread_id]+kd_index],1)+get_co3(_Dkdvec[_Dkdstart[thread_id]+kd_index],1)*(t_t)))))
                              *_Dkdval[_Dkdstart[thread_id]+kd_index]); 
                       kd_index++;
                }
                den *= den*den*den; 
                den += (num + 1.0);
                ydot[0] = r0_ptr[thread_id]*(num/den)-d_ptr[thread_id]*y[0];
             }
             else if ( (*t-1.0 > 1.0) && (*t-1.0 <= 2.0) )
             {
                t_t = *t - 1.0 - 1.0;
                while ( ka_index < _Dnka[thread_id] )
                { 
                       num += ((get_co0(_Dkavec[_Dkastart[thread_id]+ka_index],2) + 
                              (t_t)*(get_co1(_Dkavec[_Dkastart[thread_id]+ka_index],2)+
                              ((t_t)*(get_co2(_Dkavec[_Dkastart[thread_id]+ka_index],2)+get_co3(_Dkavec[_Dkastart[thread_id]+ka_index],2)*(t_t)))))
                              *_Dkaval[_Dkastart[thread_id]+ka_index]); 
                       ka_index++;
                }
                num *= num*num*num;
                num += ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id];
                while ( kd_index < _Dnkd[thread_id] )
                { 
                       den += ((get_co0(_Dkdvec[_Dkdstart[thread_id]+kd_index],2) + 
                              (t_t)*(get_co1(_Dkdvec[_Dkdstart[thread_id]+kd_index],2)+
                              ((t_t)*(get_co2(_Dkdvec[_Dkdstart[thread_id]+kd_index],2)+get_co3(_Dkdvec[_Dkdstart[thread_id]+kd_index],2)*(t_t)))))
                              *_Dkdval[_Dkdstart[thread_id]+kd_index]); 
                       kd_index++;
                }
                den *= den*den*den; 
                den += (num + 1.0);
                ydot[0] = r0_ptr[thread_id]*(num/den)-d_ptr[thread_id]*y[0];
             }
             else if ( (*t-1.0 > 2.0) && (*t-1.0 <= 3.0) )
             {
                t_t = *t - 1.0 - 2.0;
                while ( ka_index < _Dnka[thread_id] )
                { 
                       num += ((get_co0(_Dkavec[_Dkastart[thread_id]+ka_index],3) + 
                              (t_t)*(get_co1(_Dkavec[_Dkastart[thread_id]+ka_index],3)+
                              ((t_t)*(get_co2(_Dkavec[_Dkastart[thread_id]+ka_index],3)+get_co3(_Dkavec[_Dkastart[thread_id]+ka_index],3)*(t_t)))))
                              *_Dkaval[_Dkastart[thread_id]+ka_index]); 
                       ka_index++;
                }
                num *= num*num*num;
                num += ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id];
                while ( kd_index < _Dnkd[thread_id] )
                { 
                       den += ((get_co0(_Dkdvec[_Dkdstart[thread_id]+kd_index],3) + 
                              (t_t)*(get_co1(_Dkdvec[_Dkdstart[thread_id]+kd_index],3)+
                              ((t_t)*(get_co2(_Dkdvec[_Dkdstart[thread_id]+kd_index],3)+get_co3(_Dkdvec[_Dkdstart[thread_id]+kd_index],3)*(t_t)))))
                              *_Dkdval[_Dkdstart[thread_id]+kd_index]); 
                       kd_index++;
                }
                den *= den*den*den; 
                den += (num + 1.0);
                ydot[0] = r0_ptr[thread_id]*(num/den)-d_ptr[thread_id]*y[0];
             }
             else if ( (*t-1.0 > 3.0) && (*t-1.0 <= 4.0) )
             {
                t_t = *t - 1.0 - 3.0;
                while ( ka_index < _Dnka[thread_id] )
                { 
                       num += ((get_co0(_Dkavec[_Dkastart[thread_id]+ka_index],4) + 
                              (t_t)*(get_co1(_Dkavec[_Dkastart[thread_id]+ka_index],4)+
                              ((t_t)*(get_co2(_Dkavec[_Dkastart[thread_id]+ka_index],4)+get_co3(_Dkavec[_Dkastart[thread_id]+ka_index],4)*(t_t)))))
                              *_Dkaval[_Dkastart[thread_id]+ka_index]); 
                       ka_index++;
                }
                num *= num*num*num;
                num += ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id];
                while ( kd_index < _Dnkd[thread_id] )
                { 
                       den += ((get_co0(_Dkdvec[_Dkdstart[thread_id]+kd_index],4) + 
                              (t_t)*(get_co1(_Dkdvec[_Dkdstart[thread_id]+kd_index],4)+
                              ((t_t)*(get_co2(_Dkdvec[_Dkdstart[thread_id]+kd_index],4)+get_co3(_Dkdvec[_Dkdstart[thread_id]+kd_index],4)*(t_t)))))
                              *_Dkdval[_Dkdstart[thread_id]+kd_index]); 
                       kd_index++;
                }
                den *= den*den*den; 
                den += (num + 1.0);
                ydot[0] = r0_ptr[thread_id]*(num/den)-d_ptr[thread_id]*y[0];
             }
             else if ( (*t-1.0 > 4.0) && (*t-1.0 <= 5.0) )
             {
                t_t = *t - 1.0 - 4.0;
                while ( ka_index < _Dnka[thread_id] )
                { 
                       num += ((get_co0(_Dkavec[_Dkastart[thread_id]+ka_index],5) + 
                              (t_t)*(get_co1(_Dkavec[_Dkastart[thread_id]+ka_index],5)+
                              ((t_t)*(get_co2(_Dkavec[_Dkastart[thread_id]+ka_index],5)+get_co3(_Dkavec[_Dkastart[thread_id]+ka_index],5)*(t_t)))))
                              *_Dkaval[_Dkastart[thread_id]+ka_index]); 
                       ka_index++;
                }
                num *= num*num*num;
                num += ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id];
                while ( kd_index < _Dnkd[thread_id] )
                { 
                       den += ((get_co0(_Dkdvec[_Dkdstart[thread_id]+kd_index],5) + 
                              (t_t)*(get_co1(_Dkdvec[_Dkdstart[thread_id]+kd_index],5)+
                              ((t_t)*(get_co2(_Dkdvec[_Dkdstart[thread_id]+kd_index],5)+get_co3(_Dkdvec[_Dkdstart[thread_id]+kd_index],5)*(t_t)))))
                              *_Dkdval[_Dkdstart[thread_id]+kd_index]); 
                       kd_index++;
                }
                den *= den*den*den; 
                den += (num + 1.0);
                ydot[0] = r0_ptr[thread_id]*(num/den)-d_ptr[thread_id]*y[0];
             }
             else if ( (*t-1.0 > 5.0) && (*t-1.0 <= 6.0) )
             {
                t_t = *t - 1.0 - 5.0;
                while ( ka_index < _Dnka[thread_id] )
                { 
                       num += ((get_co0(_Dkavec[_Dkastart[thread_id]+ka_index],6) + 
                              (t_t)*(get_co1(_Dkavec[_Dkastart[thread_id]+ka_index],6)+
                              ((t_t)*(get_co2(_Dkavec[_Dkastart[thread_id]+ka_index],6)+get_co3(_Dkavec[_Dkastart[thread_id]+ka_index],6)*(t_t)))))
                              *_Dkaval[_Dkastart[thread_id]+ka_index]); 
                       ka_index++;
                }
                num *= num*num*num;
                num += ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id];
                while ( kd_index < _Dnkd[thread_id] )
                { 
                       den += ((get_co0(_Dkdvec[_Dkdstart[thread_id]+kd_index],6) + 
                              (t_t)*(get_co1(_Dkdvec[_Dkdstart[thread_id]+kd_index],6)+
                              ((t_t)*(get_co2(_Dkdvec[_Dkdstart[thread_id]+kd_index],6)+get_co3(_Dkdvec[_Dkdstart[thread_id]+kd_index],6)*(t_t)))))
                              *_Dkdval[_Dkdstart[thread_id]+kd_index]); 
                       kd_index++;
                }
                den *= den*den*den; 
                den += (num + 1.0);
                ydot[0] = r0_ptr[thread_id]*(num/den)-d_ptr[thread_id]*y[0];
             }
             else if ( (*t-1.0 > 6.0) && (*t-1.0 <= 7.0) )
             {
                t_t = *t - 1.0 - 6.0;
                while ( ka_index < _Dnka[thread_id] )
                { 
                       num += ((get_co0(_Dkavec[_Dkastart[thread_id]+ka_index],7) + 
                              (t_t)*(get_co1(_Dkavec[_Dkastart[thread_id]+ka_index],7)+
                              ((t_t)*(get_co2(_Dkavec[_Dkastart[thread_id]+ka_index],7)+get_co3(_Dkavec[_Dkastart[thread_id]+ka_index],7)*(t_t)))))
                              *_Dkaval[_Dkastart[thread_id]+ka_index]); 
                       ka_index++;
                }
                num *= num*num*num;
                num += ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id];
                while ( kd_index < _Dnkd[thread_id] )
                { 
                       den += ((get_co0(_Dkdvec[_Dkdstart[thread_id]+kd_index],7) + 
                              (t_t)*(get_co1(_Dkdvec[_Dkdstart[thread_id]+kd_index],7)+
                              ((t_t)*(get_co2(_Dkdvec[_Dkdstart[thread_id]+kd_index],7)+get_co3(_Dkdvec[_Dkdstart[thread_id]+kd_index],7)*(t_t)))))
                              *_Dkdval[_Dkdstart[thread_id]+kd_index]); 
                       kd_index++;
                }
                den *= den*den*den; 
                den += (num + 1.0);
                ydot[0] = r0_ptr[thread_id]*(num/den)-d_ptr[thread_id]*y[0];
             }
             else if ( (*t-1.0 > 7.0) && (*t-1.0 <= 8.0) )
             {
                t_t = *t - 1.0 - 7.0;
                while ( ka_index < _Dnka[thread_id] )
                { 
                       num += ((get_co0(_Dkavec[_Dkastart[thread_id]+ka_index],8) + 
                              (t_t)*(get_co1(_Dkavec[_Dkastart[thread_id]+ka_index],8)+
                              ((t_t)*(get_co2(_Dkavec[_Dkastart[thread_id]+ka_index],8)+get_co3(_Dkavec[_Dkastart[thread_id]+ka_index],8)*(t_t)))))
                              *_Dkaval[_Dkastart[thread_id]+ka_index]); 
                       ka_index++;
                }
                num *= num*num*num;
                num += ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id];
                while ( kd_index < _Dnkd[thread_id] )
                { 
                       den += ((get_co0(_Dkdvec[_Dkdstart[thread_id]+kd_index],8) + 
                              (t_t)*(get_co1(_Dkdvec[_Dkdstart[thread_id]+kd_index],8)+
                              ((t_t)*(get_co2(_Dkdvec[_Dkdstart[thread_id]+kd_index],8)+get_co3(_Dkdvec[_Dkdstart[thread_id]+kd_index],8)*(t_t)))))
                              *_Dkdval[_Dkdstart[thread_id]+kd_index]); 
                       kd_index++;
                }
                den *= den*den*den; 
                den += (num + 1.0);
                ydot[0] = r0_ptr[thread_id]*(num/den)-d_ptr[thread_id]*y[0];
             }
             else if ( (*t-1.0 > 8.0) && (*t-1.0 <= 9.0) )
             {
                t_t = *t - 1.0 - 8.0;
                while ( ka_index < _Dnka[thread_id] )
                { 
                       num += ((get_co0(_Dkavec[_Dkastart[thread_id]+ka_index],9) + 
                              (t_t)*(get_co1(_Dkavec[_Dkastart[thread_id]+ka_index],9)+
                              ((t_t)*(get_co2(_Dkavec[_Dkastart[thread_id]+ka_index],9)+get_co3(_Dkavec[_Dkastart[thread_id]+ka_index],9)*(t_t)))))
                              *_Dkaval[_Dkastart[thread_id]+ka_index]); 
                       ka_index++;
                }
                num *= num*num*num;
                num += ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id];
                while ( kd_index < _Dnkd[thread_id] )
                { 
                       den += ((get_co0(_Dkdvec[_Dkdstart[thread_id]+kd_index],9) + 
                              (t_t)*(get_co1(_Dkdvec[_Dkdstart[thread_id]+kd_index],9)+
                              ((t_t)*(get_co2(_Dkdvec[_Dkdstart[thread_id]+kd_index],9)+get_co3(_Dkdvec[_Dkdstart[thread_id]+kd_index],9)*(t_t)))))
                              *_Dkdval[_Dkdstart[thread_id]+kd_index]); 
                       kd_index++;
                }
                den *= den*den*den; 
                den += (num + 1.0);
                ydot[0] = r0_ptr[thread_id]*(num/den)-d_ptr[thread_id]*y[0];
             }
             else if ( (*t-1.0 > 9.0) && (*t-1.0 <= 10.0) )
             {
                t_t = *t - 1.0 - 9.0;
                while ( ka_index < _Dnka[thread_id] )
                { 
                       num += ((get_co0(_Dkavec[_Dkastart[thread_id]+ka_index],10) + 
                              (t_t)*(get_co1(_Dkavec[_Dkastart[thread_id]+ka_index],10)+
                              ((t_t)*(get_co2(_Dkavec[_Dkastart[thread_id]+ka_index],10)+get_co3(_Dkavec[_Dkastart[thread_id]+ka_index],10)*(t_t)))))
                              *_Dkaval[_Dkastart[thread_id]+ka_index]); 
                       ka_index++;
                }
                num *= num*num*num;
                num += ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id];
                while ( kd_index < _Dnkd[thread_id] )
                { 
                       den += ((get_co0(_Dkdvec[_Dkdstart[thread_id]+kd_index],10) + 
                              (t_t)*(get_co1(_Dkdvec[_Dkdstart[thread_id]+kd_index],10)+
                              ((t_t)*(get_co2(_Dkdvec[_Dkdstart[thread_id]+kd_index],10)+get_co3(_Dkdvec[_Dkdstart[thread_id]+kd_index],10)*(t_t)))))
                              *_Dkdval[_Dkdstart[thread_id]+kd_index]); 
                       kd_index++;
                }
                den *= den*den*den; 
                den += (num + 1.0);
                ydot[0] = r0_ptr[thread_id]*(num/den)-d_ptr[thread_id]*y[0];
             }
             else if ( (*t-1.0 > 10.0) && (*t-1.0 <= 11.0) )
             {
                t_t = *t - 1.0 - 10.0;
                while ( ka_index < _Dnka[thread_id] )
                { 
                       num += ((get_co0(_Dkavec[_Dkastart[thread_id]+ka_index],11) + 
                              (t_t)*(get_co1(_Dkavec[_Dkastart[thread_id]+ka_index],11)+
                              ((t_t)*(get_co2(_Dkavec[_Dkastart[thread_id]+ka_index],11)+get_co3(_Dkavec[_Dkastart[thread_id]+ka_index],11)*(t_t)))))
                              *_Dkaval[_Dkastart[thread_id]+ka_index]); 
                       ka_index++;
                }
                num *= num*num*num;
                num += ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id];
                while ( kd_index < _Dnkd[thread_id] )
                { 
                       den += ((get_co0(_Dkdvec[_Dkdstart[thread_id]+kd_index],11) + 
                              (t_t)*(get_co1(_Dkdvec[_Dkdstart[thread_id]+kd_index],11)+
                              ((t_t)*(get_co2(_Dkdvec[_Dkdstart[thread_id]+kd_index],11)+get_co3(_Dkdvec[_Dkdstart[thread_id]+kd_index],11)*(t_t)))))
                              *_Dkdval[_Dkdstart[thread_id]+kd_index]); 
                       kd_index++;
                }
                den *= den*den*den; 
                den += (num + 1.0);
                ydot[0] = r0_ptr[thread_id]*(num/den)-d_ptr[thread_id]*y[0];
             }
             else if (*t-1.0 > 11.0)
             {
                t_t = *t - 1.0 - 11.0;
                while ( ka_index < _Dnka[thread_id] )
                { 
                       num += ((get_co0(_Dkavec[_Dkastart[thread_id]+ka_index],12) + 
                              (t_t)*(get_co1(_Dkavec[_Dkastart[thread_id]+ka_index],12)+
                              ((t_t)*(get_co2(_Dkavec[_Dkastart[thread_id]+ka_index],12)+get_co3(_Dkavec[_Dkastart[thread_id]+ka_index],12)*(t_t)))))
                              *_Dkaval[_Dkastart[thread_id]+ka_index]); 
                       ka_index++;
                }
                num *= num*num*num;
                num += ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id];
                while ( kd_index < _Dnkd[thread_id] )
                { 
                       den += ((get_co0(_Dkdvec[_Dkdstart[thread_id]+kd_index],12) + 
                              (t_t)*(get_co1(_Dkdvec[_Dkdstart[thread_id]+kd_index],12)+
                              ((t_t)*(get_co2(_Dkdvec[_Dkdstart[thread_id]+kd_index],12)+get_co3(_Dkdvec[_Dkdstart[thread_id]+kd_index],12)*(t_t)))))
                              *_Dkdval[_Dkdstart[thread_id]+kd_index]); 
                       kd_index++;
                }
                den *= den*den*den; 
                den += (num + 1.0);
                ydot[0] = r0_ptr[thread_id]*(num/den)-d_ptr[thread_id]*y[0];
                }*/
             /*else if (*t-1.0 > double(60/3))
             {
                t_t = *t - 1.0 - double(60/3);
                while ( ka_index < _Dnka[thread_id] )
                { 
                       num += ((get_co0(_Dkavec[_Dkastart[thread_id]+ka_index],13) + 
                              (t_t)*(get_co1(_Dkavec[_Dkastart[thread_id]+ka_index],13)+
                              ((t_t)*(get_co2(_Dkavec[_Dkastart[thread_id]+ka_index],13)+get_co3(_Dkavec[_Dkastart[thread_id]+ka_index],13)*(t_t)))))
                              *_Dkaval[_Dkastart[thread_id]+ka_index]); 
                       ka_index++;
                }
                num *= num*num*num;
                num += ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id]*ea_ptr[thread_id];
                while ( kd_index < _Dnkd[thread_id] )
                { 
                       den += ((get_co0(_Dkdvec[_Dkdstart[thread_id]+kd_index],13) + 
                              (t_t)*(get_co1(_Dkdvec[_Dkdstart[thread_id]+kd_index],13)+
                              ((t_t)*(get_co2(_Dkdvec[_Dkdstart[thread_id]+kd_index],13)+get_co3(_Dkdvec[_Dkdstart[thread_id]+kd_index],13)*(t_t)))))
                              *_Dkdval[_Dkdstart[thread_id]+kd_index]); 
                       kd_index++;
                }
                den *= den*den*den; 
                den += (num + 1.0);
                ydot[0] = r0_ptr[thread_id]*(num/den)-d_ptr[thread_id]*y[0];
             }*/
	}
};

struct myJex
{
	__device__ void operator()(int *neq, double *t, double *y, int ml, int mu, double *pd, int nrowpd)
	{
		return;
	}
};

#endif
