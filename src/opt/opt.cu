# include "../common.h"
using namespace std;
# include "../lsoda/cuLsoda_all.cu"
# include "../lsoda/cuLsoda.hpp"
# include "./opt.hpp"
// A guassian distribution with mean 0 standard deviation of 1
double guassrand()
{
      static double V1, V2, S;
      static int phase = 0;
      double X;

      srand(time(NULL));
      if (phase == 0) 
      {
      	 do {
	    double U1 = (double)rand() / RAND_MAX;
	    double U2 = (double)rand() / RAND_MAX;
            
            V1 = 2 * U1 - 1;
	    V2 = 2 * U2 - 1;
	    S = V1 * V1 + V2 * V2;
	    } while(S >= 1 || S == 0);

	    X = V1 * sqrt(-2 * log(S) / S);
      } else
      X = V2 * sqrt(-2 * log(S) / S);
      
      phase = 1 - phase;

      return X;
}
// LSODE integrator function
//void integrate_lsoda_ode ( const vector <double> &x_d, const vector <double> &t_d, const vector <double> &sd_d, const myFex &fex, const myJex &jex, state_type &error_sim_d)
void integrate_lsoda_ode ( const int &gene_ind, const vector <double> &x_d, const vector <double> &t_d, const double &mean_xd, const myFex &fex, const myJex &jex, state_type &error_sim_d)
{
   bool check_int = false; 
   host_type error_sim_h(probSize);
   //state_type error_sim_d(probSize);
   // Studying gene 20 for the N= 25 set
   /* Local variables: input arguments for Lsoda.
      For a more detailed description see cuLsoda.cu L 171*/
   // initial value of independent variable t
   double *t = (double*)malloc(sizeof(double)*probSize);
   // initial value of dependent variable Y.SIZE() = NEQ 
   double *y = (double*)malloc(sizeof(double)*probSize);
   // JT is the jacobian type indicator
   int *jt = (int*)malloc(sizeof(int)*probSize);
   // NEQ is the number of equations
   int *neq = (int*)malloc(sizeof(int)*probSize);
   // lengths of IWORK and RWORK
   int *liw = (int*)malloc(sizeof(int)*probSize);
   int *lrw = (int*)malloc(sizeof(int)*probSize);
   // ATOL is the absolute tolerance parameter
   double *atol = (double*)malloc(sizeof(double)*probSize);
   // ITOL size of ATOL, can be same or different for each EQ in NEQ
   int *itol =(int*) malloc(sizeof(int)*probSize);
   // IOPT optional inputs argument
   int *iopt =(int*) malloc(sizeof(int)*probSize);
   // RTOL relative tolerance parameter
   double *rtol = (double*)malloc(sizeof(double)*probSize);
   // IOUT forward step iterations   
   //int *iout =(int*) malloc(sizeof(int)*probSize);
   // TOUT time point where output is desired    
   double *tout =(double*) malloc(sizeof(double)*probSize);
   // ITASK determines normal computation of Y at TOUT
   int *itask = (int*)malloc(sizeof(int)*probSize);
   // IWORK int array of length of at least 20 + NEQ
   int *iwork =(int*) malloc(sizeof(int)*21*probSize);
   // RWORK array of length of at least 22 + NEQ * MAX(16, NEQ + 9)
   double *rwork = (double*)malloc(sizeof(double)*38*probSize);
   // ISTATE input flag  
   int *istate = (int*)malloc(sizeof(int)*probSize);
   // COMMON BLOCK DECLARATIONS
   struct cuLsodaCommonBlock common[probSize];
   struct cuLsodaCommonBlock *Hcommon = common;
   int *err = (int*)malloc(sizeof(int)*probSize);
   //  End Local Block 

   // Pointers to Device versions of Local variables 
   double *_Dt;
   double *_Dy;	// [3]
   int *_Djt;
   int *_Dneq;
   int *_Dliw;
   int *_Dlrw;
   double *_Datol;	//[3]
   int *_Ditol;
   int *_Diopt;
   double *_Drtol;
   double *_Dtout;
   int *_Ditask;
   int *_Diwork;	// [23]
   double *_Drwork;	// [70]
   int *_Distate;
   struct cuLsodaCommonBlock *_Dcommon;
   int *_Derr;
   // End Pointer Block 
   // Transfer other data to device
   state_type x_d_d = x_d;
   state_type t_d_d = t_d;
   //state_type sd_d_d = sd_d;
   host_type xd_h, td_h;
   //host_type y_val_h(N_time_points*probSize);
	
   // Method instantiations for Derivative and Jacobian functions to send to template 
   thrust::device_ptr<double> t0_ptr;
   thrust::device_ptr<int> index_ptr;
   thrust::fill(error_sim_d.begin(), error_sim_d.end(), 0.0);
   // Assignment of initial values to locals 
   for (int i = 0; i < probSize; i++)
   {
      *(neq+i) = 1;
      *(y+0+i) = x_d[gene_ind*N_time_points];
      //*(y+0+i) = 1.0;
      *(t+i) = (double)0.;
      *(tout+i) = 1.0;
      *(itol+i) = 1;
      *(rtol+i) = (double)1e-3;
      *(atol+i) = (double)1e-6;
      *(itask+i) = 1;
      *(istate+i) = 1;
      *(iopt+i) = 0;
      *(lrw+i) = 38;
      *(liw+i) = 21;
      *(jt+i) = 2;
      cuLsodaCommonBlockInit(&Hcommon[i]);
      *(err +i) = -1;
   } 
   // Allocate device memory for each of the pointers, and copy the values from local to device
   cudaMalloc((void**)&_Dt,sizeof(double)*probSize);
   cudaMemcpy(_Dt,t,sizeof(double)*probSize,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&_Dy,sizeof(double)*probSize);							
   cudaMemcpy(_Dy,y,sizeof(double)*probSize,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&_Djt,sizeof(int)*probSize);
   cudaMemcpy(_Djt,jt,sizeof(int)*probSize,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&_Dneq,sizeof(int)*probSize);
   cudaMemcpy(_Dneq,neq,sizeof(int)*probSize,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&_Dliw,sizeof(int)*probSize);
   cudaMemcpy(_Dliw,liw,sizeof(int)*probSize,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&_Dlrw,sizeof(int)*probSize);
   cudaMemcpy(_Dlrw,lrw,sizeof(int)*probSize,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&_Datol,sizeof(double)*probSize);
   cudaMemcpy(_Datol,atol,sizeof(double)*probSize,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&_Ditol,sizeof(int)*probSize);							
   cudaMemcpy(_Ditol,itol,sizeof(int)*probSize,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&_Diopt,sizeof(int)*probSize);							
   cudaMemcpy(_Diopt,iopt,sizeof(int)*probSize,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&_Drtol,sizeof(double)*probSize);							
   cudaMemcpy(_Drtol,rtol,sizeof(double)*probSize,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&_Dtout,sizeof(double)*probSize);
   cudaMemcpy(_Dtout,tout,sizeof(double)*probSize,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&_Ditask,sizeof(int)*probSize);
   cudaMemcpy(_Ditask,itask,sizeof(int)*probSize,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&_Diwork,sizeof(int)*21*probSize);
   cudaMemcpy(_Diwork,iwork,sizeof(int)*21*probSize,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&_Drwork,sizeof(double)*38*probSize);
   cudaMemcpy(_Drwork,rwork,sizeof(double)*38*probSize,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&_Distate,sizeof(int)*probSize);							
   cudaMemcpy(_Distate,istate,sizeof(int)*probSize,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&_Dcommon,sizeof(struct cuLsodaCommonBlock)*probSize);	
   cudaMemcpy(_Dcommon,Hcommon,sizeof(struct cuLsodaCommonBlock)*probSize, cudaMemcpyHostToDevice);
   cudaMalloc((void**)&_Derr,sizeof(double)*probSize);
   cudaMemcpy(_Derr,istate,sizeof(double)*probSize,cudaMemcpyHostToDevice);
   // End Allocation and Copy Block 
   thrust::device_ptr<double> y_dev_p = thrust::device_pointer_cast(_Dy);
   thrust::device_ptr<double> t_dev_p = thrust::device_pointer_cast(_Dt);
   thrust::device_ptr<double> tout_dev_p = thrust::device_pointer_cast(_Dtout);
   cuLsoda<<<blocksPerGrid,threadsPerBlock>>>(fex, _Dneq, _Dy, _Dt, _Dtout, _Ditol, _Drtol, _Datol, _Ditask, _Distate, _Diopt, _Drwork, _Dlrw, _Diwork, _Dliw, jex, _Djt, _Dcommon, _Derr, probSize);
   //thrust::copy(y_dev_p,y_dev_p+probSize,y_val_h.begin());
   //xd_h.push_back(y_dev_p[0]);
   //td_h.push_back(tout_dev_p[0]);
   thrust::fill(tout_dev_p, tout_dev_p+probSize, t_d_d[1]);
   thrust::fill(y_dev_p, y_dev_p+probSize, x_d_d[gene_ind*N_time_points]);
   //thrust::fill(y_dev_p, y_dev_p+probSize, 1.0);
   thrust::fill(t_dev_p, t_dev_p+probSize, 0.0);
   //for ( int t_ind = 1; t_ind <= N_time_points; t_ind++)
   int t_ind = 1; 
   while ( t_ind < N_time_points )
   {
       cuLsoda<<<blocksPerGrid,threadsPerBlock>>>(fex, _Dneq, _Dy, _Dt, _Dtout, _Ditol, _Drtol, _Datol, _Ditask, _Distate, _Diopt, _Drwork, _Dlrw, _Diwork, _Dliw, jex, _Djt, _Dcommon, _Derr, probSize);
       if (check_int)
       {
          xd_h.push_back(y_dev_p[0]);
          td_h.push_back(tout_dev_p[0]);
       }
       thrust::fill(tout_dev_p, tout_dev_p+probSize, t_d_d[t_ind+1]);
       thrust::transform(error_sim_d.begin(), error_sim_d.end(), y_dev_p, 
                         error_sim_d.begin(), error_functor(x_d_d[gene_ind*N_time_points+t_ind],mean_xd));
       //thrust::fill(y_dev_p, y_dev_p+probSize, 1.0);
       //thrust::fill(t_dev_p, t_dev_p+probSize, 0.0);
       t_ind++;
   }
   error_sim_h = error_sim_d;
   // Copy memory back from Device to Host 
   cudaMemcpy(t,_Dt,sizeof(double)*probSize,cudaMemcpyDeviceToHost);
   cudaMemcpy(y,_Dy,sizeof(double)*probSize,cudaMemcpyDeviceToHost);
   cudaMemcpy(jt,_Djt,sizeof(int)*probSize,cudaMemcpyDeviceToHost);
   cudaMemcpy(neq,_Dneq,sizeof(int)*probSize,cudaMemcpyDeviceToHost);
   cudaMemcpy(liw,_Dliw,sizeof(int)*probSize,cudaMemcpyDeviceToHost);
   cudaMemcpy(lrw,_Dlrw,sizeof(int)*probSize,cudaMemcpyDeviceToHost);
   cudaMemcpy(atol,_Datol,sizeof(double)*probSize,cudaMemcpyDeviceToHost);
   cudaMemcpy(itol,_Ditol,sizeof(int)*probSize,cudaMemcpyDeviceToHost);
   cudaMemcpy(iopt,_Diopt,sizeof(int)*probSize,cudaMemcpyDeviceToHost);
   cudaMemcpy(rtol,_Drtol,sizeof(double)*probSize,cudaMemcpyDeviceToHost);
   cudaMemcpy(tout,_Dtout,sizeof(double)*probSize,cudaMemcpyDeviceToHost);
   cudaMemcpy(itask,_Ditask,sizeof(int)*probSize,cudaMemcpyDeviceToHost);
   cudaMemcpy(iwork,_Diwork,sizeof(int)*21*probSize,cudaMemcpyDeviceToHost);
   cudaMemcpy(rwork,_Drwork,sizeof(double)*38*probSize,cudaMemcpyDeviceToHost);
   cudaMemcpy(istate,_Distate,sizeof(int)*probSize,cudaMemcpyDeviceToHost);
   cudaMemcpy(Hcommon,_Dcommon,sizeof(struct cuLsodaCommonBlock)*probSize, cudaMemcpyDeviceToHost);
   cudaMemcpy(err,_Derr,sizeof(int)*probSize,cudaMemcpyDeviceToHost);
   // Free memory on Device 
   cudaFree(_Dt);
   cudaFree(_Dy);
   cudaFree(_Djt);
   cudaFree(_Dneq);
   cudaFree(_Dliw);
   cudaFree(_Dlrw);
   cudaFree(_Datol);
   cudaFree(_Ditol);
   cudaFree(_Diopt);
   cudaFree(_Drtol);
   cudaFree(_Dtout);
   cudaFree(_Ditask);
   cudaFree(_Diwork);
   cudaFree(_Drwork);
   cudaFree(_Distate);
   cudaFree(_Dcommon);
   cudaFree(_Derr);
   free(t);
   free(y);
   free(jt);
   free(neq);
   free(liw);
   free(lrw);
   free(atol);
   free(itol);
   free(iopt);
   free(rtol);
   free(tout);
   free(itask);
   free(iwork);
   free(rwork);
   free(istate);
   //free(common);
   free(err);
   //for ( int i = 0; i < td_h.size(); i++) cout << " " << td_h[i] << "   " << xd_h[i] << endl;
   if ( check_int )
   {
        //if ( (error_sim_h[j] < 0.01) && (check == false) )
        if (error_sim_h[0] < 0.003) 
        {  
           for ( int i = 0; i < td_h.size(); i++) cout << " " << td_h[i] << "   " << xd_h[i] << endl;
        }
   }
}
// Functor used to update parameters
template < class Tuple > 
void update_functor::operator()( Tuple step )
{
   temp_p = thrust::get<0>(step);
   temp_error = thrust::get<1>(step);
   temp_err = thrust::get<2>(step);
   //if (m_t_rand < temp_p)
   //if ( (m_t_rand < thrust::get<0>(step)) && (thrust::get<1>(step) > thrust::get<2>(step)) )    
   //if ( (m_t_rand < temp_p) && (temp_error > temp_err) )    
   if (temp_error > temp_err)    
   {
       //update error, r0, d, ea
       thrust::get<1>(step) = temp_err;
       thrust::get<3>(step) = thrust::get<4>(step);
       thrust::get<5>(step) = thrust::get<6>(step);
       thrust::get<7>(step) = thrust::get<8>(step);
       thrust::get<9>(step) = 1;
   }
} 
//MC Simulation
//void MC_sim ( const vector <double> &x_d, const vector <double> &t_d, const vector <double> &sd_d, const int *n_ka,  myFex &fex, const myJex &jex, state_type &error_opt )
void MC_sim ( const int &gene_ind, const vector <double> &x_d, const vector <double> &t_d, const double &mean_xd, const int *n_ka, const int *n_kd, const int &size_ka, const int &size_kd,  myFex &fex, const myJex &jex, state_type &error_opt )
{
   const int N_steps = 100;
   state_type r0_opt(probSize), d_opt(probSize), ea_opt(probSize), kaval_opt(size_ka), kdval_opt(size_kd);
   thrust::device_ptr<double> _Dr0_ptr, _Dd_ptr, _Dea_ptr, _Dkaval_ptr, _Dkdval_ptr, d_pt;
   thrust::counting_iterator<int> sequence_begin(0);
   fex.get_r0_ptr(_Dr0_ptr);
   thrust::copy(_Dr0_ptr, _Dr0_ptr+probSize, r0_opt.begin()); 
   fex.get_d_ptr(_Dd_ptr);
   thrust::copy(_Dd_ptr, _Dd_ptr+probSize, d_opt.begin());
   fex.get_ea_ptr(_Dea_ptr);
   thrust::copy(_Dea_ptr, _Dea_ptr+probSize, ea_opt.begin()); 
   fex.get_kaval_ptr(_Dkaval_ptr);
   thrust::copy(_Dkaval_ptr, _Dkaval_ptr+size_ka, kaval_opt.begin()); 
   fex.get_kdval_ptr(_Dkdval_ptr);
   thrust::copy(_Dkdval_ptr, _Dkdval_ptr+size_kd, kdval_opt.begin()); 
   state_type error_sim(probSize);
   srand(time(NULL));
   for ( int k =0; k < N_steps; k++)
   {
      // update param values
      thrust::transform(sequence_begin,sequence_begin + probSize,r0_opt.begin(),_Dr0_ptr,move_functor());
      thrust::transform(sequence_begin,sequence_begin + probSize,d_opt.begin(),_Dd_ptr,move_functor());
      thrust::transform(sequence_begin,sequence_begin + probSize,ea_opt.begin(),_Dea_ptr,move_functor());
      thrust::transform(sequence_begin,sequence_begin + size_ka,kaval_opt.begin(),_Dkaval_ptr,move_functor());
      thrust::transform(sequence_begin,sequence_begin + size_kd,kdval_opt.begin(),_Dkdval_ptr,move_functor());
      //fex.set_r0_ptr (_Dr0_ptr); 
      //fex.set_d_ptr (_Dd_ptr); 
      //fex.set_ea_ptr (_Dea_ptr); 
      //fex.set_ka_ptr (_Dkaval_ptr,size_ka); 
      //fex.set_kd_ptr (_Dkdval_ptr,size_kd); 
      integrate_lsoda_ode (gene_ind, x_d, t_d, mean_xd, fex, jex, error_sim);
      // Simulated Annealing calculation
      state_type p_MC(probSize); 
      double S_MC = 0.1 * exp( -2.0 * double(k+1)/double(N_steps));
      thrust::transform(error_opt.begin(),error_opt.end(),error_sim.begin(),p_MC.begin(),p_functor(S_MC));
      double temp_rand = abs(double(rand()%100)/100.0);
      // update_track 0 = no update, 1 = update
      thrust::device_vector<int> update_track(probSize);
      thrust::fill( update_track.begin(), update_track.end(), 0 );
      thrust::for_each( thrust::make_zip_iterator(
      thrust::make_tuple(
      p_MC.begin(),error_opt.begin(),error_sim.begin(),
      r0_opt.begin(),_Dr0_ptr,ea_opt.begin(),_Dea_ptr,d_opt.begin(),_Dd_ptr,update_track.begin())), 
      thrust::make_zip_iterator(
      thrust::make_tuple(
      p_MC.end(),error_opt.end(),error_sim.end(),
      r0_opt.end(),_Dr0_ptr+probSize,ea_opt.end(),_Dea_ptr+probSize,d_opt.end(),_Dd_ptr+probSize,update_track.end())),
              update_functor(temp_rand) );
      // update kavals now
      int sum_check = thrust::reduce(update_track.begin(),update_track.end());
      int ka_index = 0;
      int kd_index = 0;
      //host_type error_sim_h;
      //thrust::host_vector <int> update_track_h;
      //update_track_h =  update_track;
      //error_sim_h = error_sim;
      if ( sum_check > 0)
      {
         for (int k = 0; k < probSize; k++)
         {
             if ( update_track[k] == 1 )
             {
                //if (error_sim_h[k] < 0.1)cout << "Update for thread " << k << " with error of: " << error_sim_h[k] << endl;
                for (int j = 0; j < n_ka[k]; j++ )
                {
                    kaval_opt[ka_index] = _Dkaval_ptr[ka_index];
                    ka_index++;
                }
                for (int j = 0; j < n_kd[k]; j++ )
                {
                    kdval_opt[kd_index] = _Dkdval_ptr[kd_index];
                    kd_index++;
                }
             }  
             else 
             { 
                ka_index += n_ka[k];
                kd_index += n_kd[k];
             }
         }
      }
   }
   // update param values to optimized parameters in class fex
   thrust::copy(r0_opt.begin(),r0_opt.end(),_Dr0_ptr); 
   thrust::copy(d_opt.begin(),d_opt.end(),_Dd_ptr); 
   thrust::copy(ea_opt.begin(),ea_opt.end(),_Dea_ptr); 
   thrust::copy(kaval_opt.begin(),kaval_opt.end(),_Dkaval_ptr); 
   thrust::copy(kdval_opt.begin(),kdval_opt.end(),_Dkdval_ptr); 
}
