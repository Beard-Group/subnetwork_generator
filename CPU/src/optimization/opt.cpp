# include "../common.h"
using namespace std;
# include "../lsoda/cuLsoda_all.cpp"
# include "../lsoda/cuLsoda.hpp"
# include "./opt.hpp"
// A guassian distribution with mean 0 standard deviation of 1
double guassrand()
{
      static double V1, V2, S;
      static int phase = 0;
      double X;

      mt19937 mt;
      mt.seed(time(NULL));
      std::uniform_int_distribution<> dist_int(0,RAND_MAX);
      if (phase == 0) 
      {
      	 do {
	    double U1 = (double)dist_int(mt) / RAND_MAX;
	    double U2 = (double)dist_int(mt) / RAND_MAX;
            
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
void integrate_lsoda_ode ( const state_type &x_d, const state_type &t_d, const state_type &mean_xd, const myFex_single &fex, const myJex_single &jex, const int &gene_index, const bool out_val, double &error)
{
    error = 0.0;
    // Local variables 
    value_type *t  = (value_type*)malloc(sizeof(value_type));
    value_type *y  = (value_type*)malloc(sizeof(value_type));
    int *jt  = (int*)malloc(sizeof(int));
    int *neq  = (int*)malloc(sizeof(int));
    int *liw  = (int*)malloc(sizeof(int));
    int *lrw  = (int*)malloc(sizeof(int));
    value_type *atol  = (value_type*)malloc(sizeof(value_type));
    int *itol  =(int*) malloc(sizeof(int));
    int *iopt  =(int*) malloc(sizeof(int));
    value_type *rtol  = (value_type*)malloc(sizeof(value_type));
    int *iout  =(int*) malloc(sizeof(int));
    value_type *tout  =(value_type*) malloc(sizeof(value_type));
    int *itask  = (int*)malloc(sizeof(int));
    int *iwork  =(int*) malloc(sizeof(int)*21);
    value_type *rwork  = (value_type*)malloc(sizeof(value_type)*(38));//rwork is 22+NEQ*MAX(16,NEQ+9)
    int *istate  = (int*)malloc(sizeof(int));
    struct cuLsodaCommonBlock common ;
    struct cuLsodaCommonBlock *Hcommon  = &common ;
    int *err  = (int*)malloc(sizeof(int));
    // End Local Block 
    // Assignment of initial values to locals 
    *neq  = 1;
    *y  = x_d[gene_index*N_time_points];
    *t  = (value_type)0.;
    *tout  = 1.0;
    *itol  = 1;
    *rtol  = (value_type)1E-3;
    *atol  = (value_type)1E-6;
    *itask  = 1;
    *istate  = 1;
    *iopt  = 0;
    *lrw  = (38);
    *liw  = (21);
    *jt  = 2;
    cuLsodaCommonBlockInit(Hcommon );
    *err  = -1;
    // ODE integration
    Sanders(fex, neq, y, t, tout, itol, rtol, atol, itask, istate, iopt, rwork, lrw, iwork, liw, jex, jt, Hcommon, err);
    if (out_val)
    {
       cout << setw(16) << setprecision(8) << *tout << "   "; 
       cout << setw(16) << setprecision(8) << y[0] << endl; 
    }
    for ( int t_ind = 1; t_ind < N_time_points; t_ind++) 
    {
        *y  = x_d[gene_index*N_time_points];
        *t  = (value_type)0.;
        *tout  = t_d[t_ind]; 
        *itask  = 1;
        *istate  = 1;
        *iopt  = 0;
        Sanders(fex , neq , y , t , tout , itol , rtol , atol , itask , istate , iopt , rwork , lrw , iwork , liw , jex , jt , Hcommon , err );
        error += ((((y[0]-x_d[gene_index*N_time_points+t_ind])/mean_xd[gene_index])
                   *((y[0]-x_d[gene_index*N_time_points+t_ind])/mean_xd[gene_index]))/N_time_points); 
        if (out_val)
        {
           cout << setw(16) << setprecision(8) << *tout << "   "; 
           cout << setw(16) << setprecision(8) << y[0] << endl; 
        }
    }
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
    free(err);
}
// LSODE integrator function
void integrate_lsoda_nm ( const state_type &x_d, const state_type &t_d, const state_type &mean_xd, myFex_single &fex, myJex_single &jex, const int &gene_index, const int nka, const int nkd, double &error, const state_type &param)
{
    error = 0.0;
    fex.set_r0(pow(10,param[0]));
    fex.set_d(pow(10,param[1]));
    fex.set_ea(pow(10,param[2]));
    state_type kaval(nka), kdval(nkd);
    copy( param.begin()+3, param.begin()+3+nka, kaval.begin());
    copy( param.begin()+3+nka, param.begin()+3+nka+nkd, kdval.begin());
    transform(kaval.begin(),kaval.end(),kaval.begin(),powof10());
    transform(kdval.begin(),kdval.end(),kdval.begin(),powof10());
    fex.set_ka_val(kaval);
    fex.set_kd_val(kdval);
    // Local variables 
    value_type *t  = (value_type*)malloc(sizeof(value_type));
    value_type *y  = (value_type*)malloc(sizeof(value_type));
    int *jt  = (int*)malloc(sizeof(int));
    int *neq  = (int*)malloc(sizeof(int));
    int *liw  = (int*)malloc(sizeof(int));
    int *lrw  = (int*)malloc(sizeof(int));
    value_type *atol  = (value_type*)malloc(sizeof(value_type));
    int *itol  =(int*) malloc(sizeof(int));
    int *iopt  =(int*) malloc(sizeof(int));
    value_type *rtol  = (value_type*)malloc(sizeof(value_type));
    int *iout  =(int*) malloc(sizeof(int));
    value_type *tout  =(value_type*) malloc(sizeof(value_type));
    int *itask  = (int*)malloc(sizeof(int));
    int *iwork  =(int*) malloc(sizeof(int)*21);
    value_type *rwork  = (value_type*)malloc(sizeof(value_type)*(38));//rwork is 22+NEQ*MAX(16,NEQ+9)
    int *istate  = (int*)malloc(sizeof(int));
    struct cuLsodaCommonBlock common ;
    struct cuLsodaCommonBlock *Hcommon  = &common ;
    int *err  = (int*)malloc(sizeof(int));
    // End Local Block 
    // Assignment of initial values to locals 
    *neq  = 1;
    *y  = x_d[gene_index*N_time_points];
    *t  = (value_type)0.;
    *tout  = 1.0;
    *itol  = 1;
    *rtol  = (value_type)1E-3;
    *atol  = (value_type)1E-6;
    *itask  = 1;
    *istate  = 1;
    *iopt  = 0;
    *lrw  = (38);
    *liw  = (21);
    *jt  = 2;
    cuLsodaCommonBlockInit(Hcommon );
    *err  = -1;
    // ODE integration
    Sanders(fex, neq, y, t, tout, itol, rtol, atol, itask, istate, iopt, rwork, lrw, iwork, liw, jex, jt, Hcommon, err);
    for ( int t_ind = 1; t_ind < N_time_points; t_ind++) 
    {
        *y  = x_d[gene_index*N_time_points];
        *t  = (value_type)0.;
        *tout  = t_d[t_ind]; 
        *itask  = 1;
        *istate  = 1;
        *iopt  = 0;
        Sanders(fex , neq , y , t , tout , itol , rtol , atol , itask , istate , iopt , rwork , lrw , iwork , liw , jex , jt , Hcommon , err );
        error += ((((y[0]-x_d[gene_index*N_time_points+t_ind])/mean_xd[gene_index])
                   *((y[0]-x_d[gene_index*N_time_points+t_ind])/mean_xd[gene_index]))/N_time_points); 
    }
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
    free(err);
}
//MC Simulation New
void MC_sim_new ( const state_type &x_d, const state_type &t_d, const state_type &mean_xd, myFex_single &fex, myJex_single &jex, const int &gene_ind, const int nka, const int nkd, double &error_ode, state_type &param)
{
    const int N_steps = 100;
    double S0 = 0.1;
    double S = 0.0;
    double error_mc, p;
    state_type param_opt;
    param_opt = param;
    double g = 0.1;
    for ( int k =0; k < N_steps; k++)
    {
       // update param values
       srand(time(NULL)*(k+1));
       transform(param_opt.begin(),param_opt.end(),param.begin(),move_functor_new());
       integrate_lsoda_nm ( x_d, t_d, mean_xd, fex, jex, gene_ind, nka, nkd, error_mc, param );
       // Simulated Annealing calculation
       S = S0 * exp( -2.0 * double(k+1)/double(N_steps));
       p = exp( (error_ode - error_mc) / S );
       if ( ( error_ode > error_mc ) )
       {
          error_ode = error_mc;
          param_opt = param;
       }
    }
    param = param_opt;
}
