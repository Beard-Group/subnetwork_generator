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

      //srand(time(NULL));
      mt19937 mt;
      mt.seed(time(NULL));
      std::uniform_int_distribution<> dist_int(0,RAND_MAX);
      if (phase == 0) 
      {
      	 do {
	    //double U1 = (double)rand() / RAND_MAX;
	    //double U2 = (double)rand() / RAND_MAX;
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
void integrate_lsoda_ode ( const state_type &x_d, const state_type &t_d, const state_type &mean_xd, const myFex_single &fex, const myJex_single &jex, const int &gene_index, const bool out_val, double &error, ofstream& outfile)
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
       cout << setw(16) << setprecision(8) << double(0.0) << "   "; 
       cout << setw(16) << setprecision(8) << double(1.0) << endl; 
       outfile << setw(16) << setprecision(8) << double(0.0) << "   "; 
       outfile << setw(16) << setprecision(8) << double(1.0) << endl; 
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
           outfile << setw(16) << setprecision(8) << *tout << "   "; 
           outfile << setw(16) << setprecision(8) << y[0] << endl; 
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
    //free(common);
    free(err);
    //cout << "11111111111111111111111111111111111111" << endl;
}
//MC Simulation
/*void MC_sim ( const state_type &x_d, const state_type &t_d, const state_type &mean_xd, myFex_single &fex, myJex_single &jex, const int &gene_ind, const int nka, const int nkd, double &error_ode, state_type &param)
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
       //cout << "-----------Simulated Annealing begin------------- " << endl;
       transform(kaval_mc_opt.begin(),kaval_mc_opt.end(),kaval_mc.begin(),move_functor());
       transform(kdval_mc_opt.begin(),kdval_mc_opt.end(),kdval_mc.begin(),move_functor());
       //cout << "-----------Simulated Annealing begin------------- " << endl;
       fex.set_r0(r0_mc);
       fex.set_d(d_mc);
       fex.set_ea(ea_mc);
       fex.set_ka_val(kaval_mc);
       fex.set_kd_val(kdval_mc);
       //cout << "-----------Simulated Annealing begin------------- " << endl;
       integrate_lsoda_ode ( x_d, t_d, mean_xd, fex, jex, gene_ind, error_mc );
       //cout << "-----------Simulated Annealing begin------------- " << endl;
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
    cout << "-----------MC_END------------- " << endl;
    fex.set_r0(r0_mc_opt);
    fex.set_d(d_mc_opt); 
    fex.set_ea(ea_mc_opt); 
    fex.set_ka_val(kaval_mc_opt); 
    fex.set_kd_val(kdval_mc_opt); 
    param.push_back(r0_mc_opt);
    param.push_back(d_mc_opt);
    param.push_back(ea_mc_opt);
    param.insert(param.end(),kaval_mc_opt.begin(),kaval_mc_opt.end());
    param.insert(param.end(),kdval_mc_opt.begin(),kdval_mc_opt.end());
}*/
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
    //free(common);
    free(err);
    //cout << "1111111111111111111111111111111111111111111" << endl;
}
//MC Simulation New
void MC_sim_new ( const state_type &x_d, const state_type &t_d, const state_type &mean_xd, myFex_single &fex, myJex_single &jex, const int &gene_ind, const int nka, const int nkd, double &error_ode, state_type &param)
{
    // Local optimization search in log10 space
    //transform(param.begin(),param.end(),param.begin(),logb10());
    const int N_steps = 100;
    double S0 = 0.1;
    double S = 0.0;
    double error_mc, p;
    //double r0_mc, d_mc, ea_mc;
    //state_type kaval_mc(nka), kdval_mc(nkd);
    state_type param_opt;
    param_opt = param;
    double g = 0.1;
    //cout << "-----------MC_Sim_begin------------- " << endl;
    for ( int k =0; k < N_steps; k++)
    {
       // update param values
       srand(time(NULL)*(k+1));
       transform(param_opt.begin(),param_opt.begin()+3,param.begin(),move_functor_new());
       //transform(param_opt.begin(),param_opt.end(),param.begin(),move_functor_new());
       integrate_lsoda_nm ( x_d, t_d, mean_xd, fex, jex, gene_ind, nka, nkd, error_mc, param );
       // Simulated Annealing calculation
       S = S0 * exp( -2.0 * double(k+1)/double(N_steps));
       p = exp( (error_ode - error_mc) / S );
       //double temp_rand = abs(double(rand()%100)/100.0);
       //if ( (temp_rand < p ) && ( error_ode > error_mc ) )
       if ( ( error_ode > error_mc ) )
       {
          error_ode = error_mc;
          param_opt = param;
       }
    }
    param = param_opt;
    //cout << "-----------MC_END------------- " << endl;
}
//Nelder-Mead Optimization
/*void NM_sim( const state_type &x_d, const state_type &t_d, const state_type &mean_xd, myFex_single &fex_nm, myJex_single &jex_nm, const int &gene_ind, const int nka, const int nkd, double &error_ode, state_type &param)
{
    const int param_size = 3 + nka + nkd;
    // Local optimization search in log10 space
    //transform(param.begin(),param.end(),param.begin(),logb10());
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
//  Check the input parameters.
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
//  Initial or restarted loop.
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
        integrate_lsoda_nm ( x_d, t_d, mean_xd, fex_nm, jex_nm, gene_ind, nka, nkd, err, param );
        y_min[j] = err;
        icount_nm += 1;
        param[j] = x_min;
      }
////  The simplex construction is complete.
////  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
////  the vertex of the simplex to be replaced.
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
////  Inner loop.
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
////  Calculate PBAR, the centroid of the simplex vertices
////  excepting the vertex with Y value YNEWLO.
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
////  Reflection through the centroid.
        for ( int i = 0; i < param_size; i++ )
        {
          pstar[i] = pbar[i] + rcoeff * ( pbar[i] - p_min[i+ihi*param_size] );
        }
        // Calculate error with new parameters
        integrate_lsoda_nm ( x_d, t_d, mean_xd, fex_nm, jex_nm, gene_ind, nka, nkd, err, pstar );
        ystar = err;
        icount_nm += 1;
////  Successful reflection, so extension.
        if ( ystar < ylo )
        {
          for ( int i = 0; i < param_size; i++ )
          {
            p2star[i] = pbar[i] + ecoeff * ( pstar[i] - pbar[i] );
          }
          // Calculate error with new parameters
          integrate_lsoda_nm ( x_d, t_d, mean_xd, fex_nm, jex_nm, gene_ind, nka, nkd, err, p2star );
          y2star = err;
          icount_nm += 1;
////  Check extension.
          if ( ystar < y2star )
          {
            for ( int i = 0; i < param_size; i++ )
            {
              p_min[i+ihi*param_size] = pstar[i];
            }
            y_min[ihi] = ystar;
          }
////  Retain extension or contraction.
          else
          {
            for ( int i = 0; i < param_size; i++ )
            {
              p_min[i+ihi*param_size] = p2star[i];
            }
            y_min[ihi] = y2star;
          }
        }
////  No extension.
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
////  Contraction on the Y(IHI) side of the centroid.
          else if ( l_min == 0 )
          {
            for ( int i = 0; i < param_size; i++ )
            {
              p2star[i] = pbar[i] + ccoeff * ( p_min[i+ihi*param_size] - pbar[i] );
            }
            integrate_lsoda_nm ( x_d, t_d, mean_xd, fex_nm, jex_nm, gene_ind, nka, nkd, err, p2star );
            y2star = err;
            icount_nm += 1;
////  Contract the whole simplex.
            if ( y_min[ihi] < y2star )
            {
              for ( int j = 0; j < nn; j++ )
              {
                for ( int i = 0; i < param_size; i++ )
                {
                  p_min[i+j*param_size] = ( p_min[i+j*param_size] + p_min[i+ilo*param_size] ) * 0.5;
                  xmin[i] = p_min[i+j*param_size];
                }
                 integrate_lsoda_nm ( x_d, t_d, mean_xd, fex_nm, jex_nm, gene_ind, nka, nkd, err, xmin );
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
////  Retain contraction.
            else
            {
              for ( int i = 0; i < param_size; i++ )
              {
                p_min[i+ihi*param_size] = p2star[i];
              }
              y_min[ihi] = y2star;
            }
          }
////  Contraction on the reflection side of the centroid.
          else if ( l_min == 1 )
          {
            for ( int i = 0; i < param_size; i++ )
            {
              p2star[i] = pbar[i] + ccoeff * ( pstar[i] - pbar[i] );
            }
            integrate_lsoda_nm ( x_d, t_d, mean_xd, fex_nm, jex_nm, gene_ind, nka, nkd, err, p2star );
            y2star = err;
            icount_nm += 1;
////  Retain reflection?
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
////  Check if YLO improved.
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
////  Check to see if minimum reached.
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
////  Factorial tests to check that YNEWLO is a local minimum.
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
        integrate_lsoda_nm ( x_d, t_d, mean_xd, fex_nm, jex_nm, gene_ind, nka, nkd, err, xmin );
        z_min = err;
        icount_nm += 1;
        if ( z_min < ynewlo )
        {
          ifault_nm = 2;
          break;
        }
        xmin[i] = xmin[i] - del - del;
        integrate_lsoda_nm ( x_d, t_d, mean_xd, fex_nm, jex_nm, gene_ind, nka, nkd, err, xmin );
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
////  Restart the procedure.
      for ( int i = 0; i < param_size; i++ )
      {
        param[i] = xmin[i];
      }
      del = eps;
      numres_nm += 1;
    }
    //integrate_lsoda_nm ( x_d, t_d, mean_xd, fex_nm, jex_nm, gene_ind, nka, nkd, err, xmin );
    integrate_lsoda_nm ( x_d, t_d, mean_xd, fex_nm, jex_nm, gene_ind, nka, nkd, err, param );
    //cout << "The optimized value of error after NM is: " << err << endl;
    error_ode = err;
    //cout << "-----------NM_Sim_end------------- " << endl;
    param = xmin;
};*/
