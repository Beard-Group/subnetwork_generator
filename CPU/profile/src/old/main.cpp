/*
 *  main.cpp
 */

#include "./common.h"
#include "./io/io.hpp"
#include "./spline/spline.hpp"
#include "./lsoda/cuLsoda.hpp"
#include "./lsoda/cuLsoda_all.cpp"
//#include "./lsoda/cuLsoda_kernel.cpp"

using namespace std;
// LSODE integrator function
void integrate_lsoda_ode ( const state_type &x_d, const state_type &t_d, const state_type &mean_xd, myFex_single &fex, myJex_single &jex, const int &gene_index, double &error)
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
    for ( int t_ind = 1; t_ind < N_time_points; t_ind++) 
    {
        *y  = x_d[gene_index*N_time_points];
        *t  = (value_type)0.;
        *tout  = t_d[t_ind]; 
        *itask  = 1;
        *istate  = 1;
        *iopt  = 0;
        Sanders(fex , neq , y , t , tout , itol , rtol , atol , itask , istate , iopt , rwork , lrw , iwork , liw , jex , jt , Hcommon , err );
        error += (((y[0]-x_d[gene_index*N_time_points+t_ind])/mean_xd[gene_index])
                   *((y[0]-x_d[gene_index*N_time_points+t_ind])/mean_xd[gene_index]))/N_time_points; 
    }
}
// powof10_functor
struct powof10 : public unary_function<value_type,value_type>
{
       value_type operator()(value_type &x) {return pow(10,x);} 
};
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
        error += (((y[0]-x_d[gene_index*N_time_points+t_ind])/mean_xd[gene_index])
                   *((y[0]-x_d[gene_index*N_time_points+t_ind])/mean_xd[gene_index]))/N_time_points; 
    }
}
// A guassian distribution with mean 0 standard deviation of 1
double guassrand()
{
      static double V1, V2, S;
      static int phase = 0;
      double X;

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
//
struct move_functor
{
   value_type operator()(const value_type &param)
   {
       double g = 0.1;
       return exp10(fmod(log10(param)+g*guassrand(),2.0));
   }
};
// log10_functor
struct logb10 : public unary_function<value_type,value_type>
{
       value_type operator()(value_type &x) {return log10(x);} 
};
//MC Simulation
void MC_sim ( const state_type &x_d, const state_type &t_d, const state_type &mean_xd, myFex_single &fex, myJex_single &jex, const int &gene_ind, const int nka, const int nkd, double &error_ode)
{
    //error_ode = 0.0;
    const int N_steps = 100;
    double S0 = 0.1;
    double S = 0.0;
    double error_mc, p;
    double r0_mc, d_mc, ea_mc;
    state_type kaval_mc(nka), kdval_mc(nkd);
    double r0_mc_opt, d_mc_opt, ea_mc_opt;
    state_type kaval_mc_opt(nka), kdval_mc_opt(nkd);
    fex.get_r0(r0_mc_opt);
    fex.get_d(d_mc_opt); 
    fex.get_ea(ea_mc_opt); 
    fex.get_ka_val(kaval_mc_opt); 
    fex.get_kd_val(kdval_mc_opt); 
    double g = 0.1;
    cout << "-----------MC_Sim_begin------------- " << endl;
    for ( int k =0; k < N_steps; k++)
    {
       // update r0, ea, d, ka values
       r0_mc =  exp10(fmod(log10(r0_mc_opt)+g*guassrand(),2.0));
       d_mc =  exp10(fmod(log10(d_mc_opt)+g*guassrand(),2.0));
       ea_mc =  exp10(fmod(log10(ea_mc_opt)+g*guassrand(),2.0));
       transform(kaval_mc_opt.begin(),kaval_mc_opt.end(),kaval_mc.begin(),move_functor());
       transform(kdval_mc_opt.begin(),kdval_mc_opt.end(),kdval_mc.begin(),move_functor());
       fex.set_r0(r0_mc);
       fex.set_d(d_mc);
       fex.set_ea(ea_mc);
       fex.set_ka_val(kaval_mc);
       fex.set_kd_val(kdval_mc);
       integrate_lsoda_ode ( x_d, t_d, mean_xd, fex, jex, gene_ind, error_mc );
       // Simulated Annealing calculation
       S = S0 * exp( -2.0 * double(k+1)/double(N_steps));
       p = exp( (error_ode - error_mc) / S );
       if ( (guassrand() < p ) && ( error_ode > error_mc ) )
       //if ( ( error_ode > error_mc ) )
       {
          //cout << "----------------MC UPDATE----------------" << endl;
          error_ode = error_mc;
          r0_mc_opt = r0_mc;
          d_mc_opt = d_mc;
          ea_mc_opt = ea_mc;
          for (int ind_a = 0; ind_a < nka; ind_a++) 
              kaval_mc_opt[ind_a] = kaval_mc[ind_a];
          for (int ind_d = 0; ind_d < nkd; ind_d++)
              kdval_mc_opt[ind_d] = kdval_mc[ind_d];
       }
    }
    fex.set_r0(r0_mc_opt);
    fex.set_d(d_mc_opt); 
    fex.set_ea(ea_mc_opt); 
    fex.set_ka_val(kaval_mc_opt); 
    fex.set_kd_val(kdval_mc_opt); 
};
//
int main(void)    
{
    cout << "----Reading in data----" << endl; 	
    vector <double> t_d(N_time_points);
    vector <double> x_d(N_gene*N_time_points);
    read_data ( N_time_points, N_gene, t_d, x_d);
    cout << "----Setting up spline co-efficients----" << endl; 	
    vector <double> deriv_time(N_time_points*N_gene);
    cout << "  Setting up co-fficients of cubic spline" << endl;
    for ( int i = 0 ; i < N_gene ; i++ )
    {
        spline_deriv ( i, N_time_points, t_d, x_d, deriv_time);
    }
    // Set up the co-effients matrix 
    vector <double> cub_coeff_spline (N_gene*N_time_points*4);
    spline_coeff ( N_gene, N_time_points, t_d, x_d, deriv_time, cub_coeff_spline);
    // Parameters for the governing ODE
    vector <double> r0(N_gene); 
    vector <double> d(N_gene); 
    vector <double> ea(N_gene);
    vector<int> n_ka(N_gene); 
    vector<int> n_kd(N_gene);
    vector <double> mean_xd(N_gene);
    /*r0[0] = 6.5525;
    r0[1] = 2.8578;
    r0[2] = 6.3162;
    r0[3] = 7.2104;*/
    
    r0[0] = fabs(guassrand());
    r0[1] = fabs(guassrand());
    r0[2] = fabs(guassrand());
    r0[3] = fabs(guassrand());
    
    /*d[0] = 0.3609;
    d[1] = 0.4201;
    d[2] = 0.4323;
    d[3] = 0.4366;*/
 
    d[0] = fabs(guassrand());
    d[1] = fabs(guassrand());
    d[2] = fabs(guassrand());
    d[3] = fabs(guassrand());
 
    ea[0] = 0.3351;
    ea[1] = 1.1391;
    ea[2] = 0.3420;
    ea[3] = 1.7558;
 
    n_ka[0] = 2; 
    n_ka[1] = 1; 
    n_ka[2] = 0; 
    n_ka[3] = 1; 
    
    n_kd[0] = 1; 
    n_kd[1] = 0; 
    n_kd[2] = 3; 
    n_kd[3] = 2; 
    
    vector <int> ka_vec(4);
    ka_vec[0]=1; 
    ka_vec[1]=3; 
    ka_vec[2]=2; 
    ka_vec[3]=0; 
    vector <int> kd_vec(6);
    kd_vec[0]=2; 
    kd_vec[1]=0; 
    kd_vec[2]=1; 
    kd_vec[3]=3; 
    kd_vec[4]=1; 
    kd_vec[5]=2; 
    vector <double> ka_val(4);
    ka_val[0]=1.6975;
    ka_val[1]=1.2929;
    ka_val[2]=18.0991;
    ka_val[3]=0.4726;
    vector <double> kd_val(6);
    kd_val[0]=5.0678;
    kd_val[1]=3.5361;
    kd_val[2]=0.1899;
    kd_val[3]=7.7584;
    kd_val[4]=0.9993;
    kd_val[5]=0.9606;
    // Calculate Mean
    for ( int i = 0; i < N_gene; i++)
    {
        for ( int j = 0; j < N_time_points; j++ ) mean_xd[i] += x_d[i*N_time_points+j]; 
    }
 
    for ( int i = 0; i < N_gene; i++)
    {
        mean_xd[i] /= double(N_time_points);
    }
    state_type r0_opt(N_gene), d_opt(N_gene), ea_opt(N_gene), kaval_opt(4), kdval_opt(6); 
    int ka_index = 0;
    int kd_index = 0;
    // Engine for optimization of parameters  
    for (int gene_ind = 0; gene_ind < N_gene; gene_ind++)
    {
        //Preprocess to set up data for per gene
        const int param_size = 3 + n_ka[gene_ind] + n_kd[gene_ind];
        vector <double> param;
        vector <int> kavec_temp(n_ka[gene_ind]), kdvec_temp(n_kd[gene_ind]);
        vector <double> kaval_temp(n_ka[gene_ind]), kdval_temp(n_kd[gene_ind]);
        myFex_single fex_nm;
        myJex_single jex_nm;
        fex_nm.set_tau(1.0);
        fex_nm.set_r0(r0[gene_ind]);
        fex_nm.set_d(d[gene_ind]);
        fex_nm.set_ea(ea[gene_ind]);
        fex_nm.set_n_ka(n_ka[gene_ind]);
        fex_nm.set_n_kd(n_kd[gene_ind]);
        copy( ka_vec.begin()+ka_index, ka_vec.begin()+ka_index+n_ka[gene_ind], kavec_temp.begin());
        copy( kd_vec.begin()+kd_index, kd_vec.begin()+kd_index+n_kd[gene_ind], kdvec_temp.begin());
        fex_nm.set_ka_vec(kavec_temp);
        fex_nm.set_kd_vec(kdvec_temp);
        copy( ka_val.begin()+ka_index, ka_val.begin()+ka_index+n_ka[gene_ind], kaval_temp.begin());
        copy( kd_val.begin()+kd_index, kd_val.begin()+kd_index+n_kd[gene_ind], kdval_temp.begin());
        fex_nm.set_ka_val(kaval_temp);
        fex_nm.set_kd_val(kdval_temp);
        fex_nm.set_coeff(cub_coeff_spline);
        // ODE calculation
        double error_ode; 
        integrate_lsoda_ode ( x_d, t_d, mean_xd, fex_nm, jex_nm, gene_ind, error_ode );
        cout << "Value of error for the first run of ODE is: " << error_ode << endl;
        /*// MC Simulation 
        double r0_mc_opt, d_mc_opt, ea_mc_opt;
        state_type kaval_mc_opt(nka), kdval_mc_opt(nkd);
        const int N_steps = 100;
        double S0 = 0.1;
        double S = 0.0;
        double error_mc, p;
        double r0_mc, d_mc, ea_mc;
        state_type kaval_mc(n_ka[gene_ind]), kdval_mc(n_kd[gene_ind]);
        double g = 0.1;
        cout << "-----------MC_Sim_begin------------- " << endl;
        for ( int k =0; k < N_steps; k++)
        {
           // update r0, ea, d, ka values
           r0_mc =  exp10(fmod(log10(r0_mc_opt)+g*guassrand(),2.0));
           d_mc =  exp10(fmod(log10(d_mc_opt)+g*guassrand(),2.0));
           ea_mc =  exp10(fmod(log10(ea_mc_opt)+g*guassrand(),2.0));
           transform(kaval_mc_opt.begin(),kaval_mc_opt.end(),kaval_mc.begin(),move_functor());
           transform(kdval_mc_opt.begin(),kdval_mc_opt.end(),kdval_mc.begin(),move_functor());
           fex_nm.set_r0(r0_mc);
           fex_nm.set_d(d_mc);
           fex_nm.set_ea(ea_mc);
           fex_nm.set_ka_val(kaval_mc);
           fex_nm.set_kd_val(kdval_mc);
           integrate_lsoda_ode ( x_d, t_d, mean_xd, fex_nm, jex_nm, gene_ind, error_mc );
           // Simulated Annealing calculation
           S = S0 * exp( -2.0 * double(k+1)/double(N_steps));
           p = exp( (error_ode - error_mc) / S );
           if ( (guassrand() < p ) && ( error_ode > error_mc ) )
           //if ( ( error_ode > error_mc ) )
           {
              error_ode = error_mc;
              r0_mc_opt = r0_mc;
              d_mc_opt = d_mc;
              ea_mc_opt = ea_mc;
              for (int ind_a = 0; ind_a < n_ka[gene_ind]; ind_a++) 
                  kaval_mc_opt[ind_a] = kaval_mc[ind_a];
              for (int ind_d = 0; ind_d < n_kd[gene_ind]; ind_d++)
                  kdval_mc_opt[ind_d] = kdval_mc[ind_d];
           }
        }*/
        MC_sim ( x_d, t_d, mean_xd, fex_nm, jex_nm, gene_ind, n_ka[gene_ind], n_kd[gene_ind], error_ode);
        cout << "The optimized value of error from MC calculation is: " << error_ode << endl;
        cout << "-----------MC_Sim_end------------- " << endl;
        /*param.push_back(r0_mc_opt);
        param.push_back(d_mc_opt);
        param.push_back(ea_mc_opt);
        param.insert(param.end(),kaval_mc_opt.begin(),kaval_mc_opt.end());
        param.insert(param.end(),kdval_mc_opt.begin(),kdval_mc_opt.end());
        double error_nm;
        // Local optimization search in log10 space
        transform(param.begin(),param.end(),param.begin(),logb10());
        // Nelder Mead optimization
        double reqmin = 1.0E-08;
        vector <double> step(param_size,1.0);
        int konvge = 10;
        int kcount = 100;
        double ccoeff = 0.5;
        double del;
        double dn;
        double dnn;
        double ecoeff = 2.0;
        double eps = 0.001;
        int ihi;
        int ilo;
        int jcount;
        int l_min;
        int nn;
        int ifault_nm, icount_nm, numres_nm;
        double ynewlo;
        vector <double> xmin(param_size);
        double rcoeff = 1.0;
        double rq;
        double x_min;
        vector <double> y_min(param_size+1);
        double y2star;
        double ylo;
        double ystar;
        double z_min;
        double err = 0.0;
///  /  Check the input parameters.
        if ( reqmin <= 0.0 )
        {
          ifault_nm = 1;
         // return;
        }
        if ( param_size < 1 )
        {
          ifault_nm = 1;
          //return;
        }
        if ( konvge < 1 )
        {
          ifault_nm = 1;
          //return;
        }
        //
        vector <value_type> p_min(param_size*(param_size+1));
        vector <value_type> p2star(param_size);
        vector <value_type> pbar(param_size);
        vector <value_type> pstar(param_size);
        icount_nm = 0;
        numres_nm = 0;
        jcount = konvge; 
        dn = ( double ) ( param_size );
        nn = param_size + 1;
        dnn = ( double ) ( nn );
        del = 0.01;
        rq = reqmin * dn;
///  /  Initial or restarted loop.
        for ( ; ; )
        {
          for ( int i = 0; i < param_size; i++ )
          { 
            p_min[i+param_size*param_size] = param[i];
          }
          y_min[param_size] = error_ode;
  
          for ( int j = 0; j < param_size; j++ )
          {
            x_min = param[j];
            param[j] += step[j] * del;
            for ( int i = 0; i < param_size; i++ )
            {
              p_min[i+j*param_size] = param[i];
            }
            integrate_lsoda_nm ( x_d, t_d, mean_xd, fex_nm, jex_nm, gene_ind, n_ka[gene_ind], n_kd[gene_ind], err, param );
            y_min[j] = err;
            icount_nm += 1;
            param[j] = x_min;
          }
////////  The simplex construction is complete.
////////  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
////////  the vertex of the simplex to be replaced.
          ylo = y_min[0];
          ilo = 0;
  
          for ( int i = 1; i < nn; i++ )
          {
            if ( y_min[i] < ylo )
            {
              ylo = y_min[i];
              ilo = i;
            }
          }
////////  Inner loop.
          for ( ; ; )
          {
            if ( kcount <= icount_nm )
            {
              break;
            }
            ynewlo = y_min[0];
            ihi = 0;
  
            for ( int i = 1; i < nn; i++ )
            {
              if ( ynewlo < y_min[i] )
              {
                ynewlo = y_min[i];
                ihi = i;
              }
            }
////////  Calculate PBAR, the centroid of the simplex vertices
////////  excepting the vertex with Y value YNEWLO.
            for ( int i = 0; i < param_size; i++ )
            {
              z_min = 0.0;
              for ( int j = 0; j < nn; j++ )
              { 
                z_min = z_min + p_min[i+j*param_size];
              }
              z_min = z_min - p_min[i+ihi*param_size];  
              pbar[i] = z_min / dn;
            }
////////  Reflection through the centroid.
            for ( int i = 0; i < param_size; i++ )
            {
              pstar[i] = pbar[i] + rcoeff * ( pbar[i] - p_min[i+ihi*param_size] );
            }
            // Calculate error with new parameters
            integrate_lsoda_nm ( x_d, t_d, mean_xd, fex_nm, jex_nm, gene_ind, n_ka[gene_ind], n_kd[gene_ind], err, pstar );
            ystar = err;
            icount_nm += 1;
////////  Successful reflection, so extension.
            if ( ystar < ylo )
            {
              for ( int i = 0; i < param_size; i++ )
              {
                p2star[i] = pbar[i] + ecoeff * ( pstar[i] - pbar[i] );
              }
              // Calculate error with new parameters
              integrate_lsoda_nm ( x_d, t_d, mean_xd, fex_nm, jex_nm, gene_ind, n_ka[gene_ind], n_kd[gene_ind], err, p2star );
              y2star = err;
              icount_nm += 1;
////////  Check extension.
              if ( ystar < y2star )
              {
                for ( int i = 0; i < param_size; i++ )
                {
                  p_min[i+ihi*param_size] = pstar[i];
                }
                y_min[ihi] = ystar;
              }
////////  Retain extension or contraction.
              else
              {
                for ( int i = 0; i < param_size; i++ )
                {
                  p_min[i+ihi*param_size] = p2star[i];
                }
                y_min[ihi] = y2star;
              }
            }
////////  No extension.
            else
            {
              l_min = 0;
              for ( int i = 0; i < nn; i++ )
              {
                if ( ystar < y_min[i] )
                {
                  l_min += 1;
                }
              }
  
              if ( 1 < l_min )
              {
                for ( int i = 0; i < param_size; i++ )
                {
                  p_min[i+ihi*param_size] = pstar[i];
                }
                y_min[ihi] = ystar;
              }
////////  Contraction on the Y(IHI) side of the centroid.
              else if ( l_min == 0 )
              {
                for ( int i = 0; i < param_size; i++ )
                {
                  p2star[i] = pbar[i] + ccoeff * ( p_min[i+ihi*param_size] - pbar[i] );
                }
                integrate_lsoda_nm ( x_d, t_d, mean_xd, fex_nm, jex_nm, gene_ind, n_ka[gene_ind], n_kd[gene_ind], err, p2star );
                y2star = err;
                icount_nm += 1;
////////  Contract the whole simplex.
                if ( y_min[ihi] < y2star )
                {
                  for ( int j = 0; j < nn; j++ )
                  {
                    for ( int i = 0; i < param_size; i++ )
                    {
                      p_min[i+j*param_size] = ( p_min[i+j*param_size] + p_min[i+ilo*param_size] ) * 0.5;
                      xmin[i] = p_min[i+j*param_size];
                    }
                     integrate_lsoda_nm ( x_d, t_d, mean_xd, fex_nm, jex_nm, gene_ind, n_ka[gene_ind], n_kd[gene_ind], err, xmin );
                     y2star = err;
                     icount_nm += 1;
                  }
                  ylo = y_min[0];
                  ilo = 0;
  
                  for ( int i = 1; i < nn; i++ )
                  {
                    if ( y_min[i] < ylo )
                    {
                      ylo = y_min[i];
                      ilo = i;
                    }
                  }
                  continue;
                }
////////  Retain contraction.
                else
                {
                  for ( int i = 0; i < param_size; i++ )
                  {
                    p_min[i+ihi*param_size] = p2star[i];
                  }
                  y_min[ihi] = y2star;
                }
              }
////////  Contraction on the reflection side of the centroid.
              else if ( l_min == 1 )
              {
                for ( int i = 0; i < param_size; i++ )
                {
                  p2star[i] = pbar[i] + ccoeff * ( pstar[i] - pbar[i] );
                }
                integrate_lsoda_nm ( x_d, t_d, mean_xd, fex_nm, jex_nm, gene_ind, n_ka[gene_ind], n_kd[gene_ind], err, p2star );
                y2star = err;
                icount_nm += 1;
////////  Retain reflection?
                if ( y2star <= ystar )
                {
                  for ( int i = 0; i < param_size ; i++ )
                  {
                    p_min[i+ihi*param_size] = p2star[i];
                  }
                  y_min[ihi] = y2star;
                }
                else
                {
                  for ( int i = 0; i < param_size; i++ )
                  {
                    p_min[i+ihi*param_size] = pstar[i];
                  }
                  y_min[ihi] = ystar;
                }
              }
            }
////////  Check if YLO improved.
            if ( y_min[ihi] < ylo )
            {
              ylo = y_min[ihi];
              ilo = ihi;
            }
            jcount -= 1;
  
            if ( 0 < jcount )
            {
              continue;
            }
////////  Check to see if minimum reached.
            if ( icount_nm <= kcount )
            {
              jcount = konvge;
  
              z_min = 0.0;
              for ( int i = 0; i < nn; i++ )
              {
                z_min = z_min + y_min[i];
              }
              x_min = z_min / dnn;
  
              z_min = 0.0;
              for ( int i = 0; i < nn; i++ )
              {
                z_min = z_min + pow ( y_min[i] - x_min, 2 );
              }
  
              if ( z_min <= rq )
              {
                break;
              }
            }
          }
////////  Factorial tests to check that YNEWLO is a local minimum.
          for ( int i = 0; i < param_size; i++ )
          {
            xmin[i] = p_min[i+ilo*param_size];
          }
          ynewlo = y_min[ilo];
  
          if ( kcount < icount_nm )
          {
            ifault_nm = 2;
            break;
          }
  
          ifault_nm = 0;
  
          for ( int i = 0; i < param_size; i++ )
          {
            del = step[i] * eps;
            xmin[i] = xmin[i] + del;
            integrate_lsoda_nm ( x_d, t_d, mean_xd, fex_nm, jex_nm, gene_ind, n_ka[gene_ind], n_kd[gene_ind], err, xmin );
            z_min = err;
            icount_nm += 1;
            if ( z_min < ynewlo )
            {
              ifault_nm = 2;
              break;
            }
            xmin[i] = xmin[i] - del - del;
            integrate_lsoda_nm ( x_d, t_d, mean_xd, fex_nm, jex_nm, gene_ind, n_ka[gene_ind], n_kd[gene_ind], err, xmin );
            z_min = err;
            icount_nm += 1;
            if ( z_min < ynewlo )
            {
              ifault_nm = 2;
              break;
            }
            xmin[i] = xmin[i] + del;
          }
  
          if ( ifault_nm == 0 )
          {
            break;
          }
////////  Restart the procedure.
          for ( int i = 0; i < param_size; i++ )
          {
            param[i] = xmin[i];
          }
          del = eps;
          numres_nm += 1;
        }
        integrate_lsoda_nm ( x_d, t_d, mean_xd, fex_nm, jex_nm, gene_ind, n_ka[gene_ind], n_kd[gene_ind], err, xmin );
        cout << "The optimized value of error after NM is: " << err << endl;
        cout << "-----------NM_Sim_end------------- " << endl;
        r0_opt[gene_ind] = xmin[0];
        d_opt[gene_ind] = xmin[1];
        ea_opt[gene_ind] = xmin[2];
        copy( xmin.begin()+3, xmin.begin()+3+n_ka[gene_ind], kaval_opt.begin()+ka_index);
        copy( xmin.begin()+3+n_ka[gene_ind], xmin.begin()+3+n_ka[gene_ind]+n_kd[gene_ind], kdval_opt.begin()+kd_index);*/
        ka_index += n_ka[gene_ind];
        kd_index += n_kd[gene_ind];
    }
    cout << "-------------Nelder-Mead-End--------------" << endl;

    transform(r0_opt.begin(),r0_opt.end(),r0_opt.begin(),powof10());
    transform(d_opt.begin(),d_opt.end(),d_opt.begin(),powof10());
    transform(ea_opt.begin(),ea_opt.end(),ea_opt.begin(),powof10());
    transform(kaval_opt.begin(),kaval_opt.end(),kaval_opt.begin(),powof10());
    transform(kdval_opt.begin(),kdval_opt.end(),kdval_opt.begin(),powof10());
    
    return 0;
} 
