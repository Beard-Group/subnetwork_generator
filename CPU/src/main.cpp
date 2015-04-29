#include "./common.h"
#include "./io/io.hpp"
#include "./spline/spline.hpp"
#include "./lsoda/cuLsoda.hpp"
#include "./optimization/opt.hpp"
#include "./subnetsearch/search.hpp"

using namespace std;
//
int main(void)    
{
    cout << "----Reading in data----" << endl; 	
    vector <double> t_d(N_time_points);
    vector <double> x_d(N_gene*N_time_points);
    vector <double> sd_d(N_gene*N_time_points);
    read_data ( N_time_points, N_gene, t_d, x_d);
    cout << "----Setting up spline co-efficients----" << endl;
    vector <double> deriv_time(N_time_points*N_gene);
    for ( int i = 0 ; i < N_gene ; i++ )
    {
        spline_deriv ( i, N_time_points, t_d, x_d, deriv_time);
    }
    // Set up the co-effients matrix 
    vector <double> cub_coeff_spline (N_gene*N_time_points*4);
    spline_coeff ( N_gene, N_time_points, t_d, x_d, deriv_time, cub_coeff_spline);
    // Determine mean, standard deviation, and correlation statistics
    vector <double> mean_xd, sd_xd, x_spline ;
    for ( int ind = 0; ind < N_gene; ind++)
    { 
        vector <double> x_spline_temp(N_time_points);
        copy( x_d.begin()+ind*N_time_points, x_d.begin()+(ind+1)*N_time_points, x_spline_temp.begin());
        calc_stats( x_spline_temp, mean_xd, sd_xd);
        vector<double>().swap(x_spline_temp);
    }
    double E_0 = 5.0E-03;
    boost::timer timer;
    double proc_time = 0.0;
    timer.restart();
    int acc_subnet = 0;
    for (int gene_ind = 644; gene_ind < 645; gene_ind++)
    {
        vector <double> act_vec, inh_vec;
        //cout << "--------START FOR GENE "<< gene_ind+1 << " -----------" << endl;
        //cout << " Values arranged as error, n_ka, gene ids for ka, n_kd, gene ids for kd," << endl;
        //cout<< " r0, d, ea, ka_values, kd_values for each subnet. " << endl;
        int subnet_attempts = 0;
        mt19937 mt_s;
        mt_s.seed(time(NULL));
        std::uniform_int_distribution<> dist_in(50,200);
        while (acc_subnet < probSize) 
        {
           int n_ka_temp = 0;
           int n_kd_temp = 0;
           vector <int> kavec_temp, kdvec_temp;
           vector <double> kaval_temp, kdval_temp;
           for (unsigned int k = 0; k < dist_in(mt_s); k++) 
               random_select(gene_ind, k+subnet_attempts, mean_xd, n_ka_temp, n_kd_temp, kavec_temp, kdvec_temp, kaval_temp, kdval_temp);
           subnet_attempts++;
           myFex_single fex_nm;
           myJex_single jex_nm;
           fex_nm.set_tau(1.0);
           double r0_temp = mean_xd[gene_ind];
           fex_nm.set_r0(r0_temp);
           double d_temp = 1.0E-01;
           fex_nm.set_d(d_temp);
           double ea_temp = 1.0;
           fex_nm.set_ea(ea_temp);
           fex_nm.set_n_ka(n_ka_temp);
           fex_nm.set_n_kd(n_kd_temp);
           fex_nm.set_ka_vec(kavec_temp);
           fex_nm.set_kd_vec(kdvec_temp);
           fex_nm.set_ka_val(kaval_temp);
           fex_nm.set_kd_val(kdval_temp);
           fex_nm.set_coeff(cub_coeff_spline);
           // ODE calculation
           timer.restart();
           double error_ode = 0.0;
           bool out_val = false; 
           integrate_lsoda_ode ( x_d, t_d, mean_xd, fex_nm, jex_nm, gene_ind, out_val, error_ode );
           state_type param;
           param.push_back(r0_temp);
           param.push_back(d_temp);
           param.push_back(ea_temp);
           param.insert(param.end(),kaval_temp.begin(),kaval_temp.end());
           param.insert(param.end(),kdval_temp.begin(),kdval_temp.end());
           // MC and NM search conducted in log space so transform param to log10 space
           transform(param.begin(),param.end(),param.begin(),logb10());
           MC_sim_new ( x_d, t_d, mean_xd, fex_nm, jex_nm, gene_ind, n_ka_temp, n_kd_temp, error_ode, param);
           proc_time += timer.elapsed();
	   fex_nm.clear_mem();
           if (error_ode < E_0)
           {
                  bool acc_sub = false;
                  out_val = true;
                  int n_ka_val = 0;
                  int n_kd_val = 0;
                  vector <int> kavec_val, kdvec_val;
                  vector <double> kaval_val, kdval_val;
                  bool acc_sub_i = false;
                  for ( int k = 0; k < n_ka_temp; k++) 
                  {
                     n_ka_val++;
                     kavec_val.push_back(kavec_temp[k]);
                     kaval_val.push_back(pow(10,param[3+k]));
                     if ( pow(10,param[3+k]) > 0.05 ) 
                     {
                        act_vec.push_back(kavec_temp[k]);
                        acc_sub = true; 
                     } 
                  }
                  for ( int k = 0; k < n_kd_temp; k++) 
                  {
                     n_kd_val++;
                     kdvec_val.push_back(kdvec_temp[k]);
                     kdval_val.push_back(pow(10,param[3+n_ka_temp+k]));
                     if ( pow(10,param[3+n_ka_temp+k]) > 0.05) 
                     {  
                        inh_vec.push_back(kdvec_temp[k]);   
                        acc_sub = true; 
                     }
                  }
                  if (acc_sub) 
                  {  
                     acc_subnet++;
                  }
                  myFex_single fex_val;
                  myJex_single jex_val;
                  fex_val.set_tau(1.0);
                  double r0_val = pow(10,param[0]);
                  fex_val.set_r0(r0_val);
                  double d_val = pow(10,param[1]);
                  fex_val.set_d(d_val);
                  double ea_val = pow(10,param[2]);
                  fex_val.set_ea(ea_val);
                  fex_val.set_n_ka(n_ka_val);
                  fex_val.set_n_kd(n_kd_val);
                  fex_val.set_ka_vec(kavec_val);
                  fex_val.set_kd_vec(kdvec_val);
                  fex_val.set_ka_val(kaval_val);
                  fex_val.set_kd_val(kdval_val);
                  fex_val.set_coeff(cub_coeff_spline);
                  double error_val = 0.0;
                  integrate_lsoda_ode ( x_d, t_d, mean_xd, fex_val, jex_val, gene_ind, out_val, error_val );
                  cout << "Value of error for an acceptanle subnet is: " << error_val << endl;
           } 
           vector<double>().swap(param);
           vector<int>().swap(kavec_temp);
           vector<int>().swap(kdvec_temp);
           vector<double>().swap(kaval_temp);
           vector<double>().swap(kdval_temp);
        }
        cout << "------------ START FOR ACTIVATORS --------------" << endl;
        for ( int ind_con = 0; ind_con < N_gene; ind_con++ )
        {
           int act_count = count(act_vec.begin(),act_vec.end(),ind_con);
           cout << setw(3) << ind_con+1 << "   " << setw(3) << act_count << endl; 
        }
        cout << "-------------- END FOR ACTIVATORS --------------" << endl;
        cout << "------------ START FOR INHIBITORS --------------" << endl;
        for ( int ind_con = 0; ind_con < N_gene; ind_con++ )
        {
           int inh_count = count(inh_vec.begin(),inh_vec.end(),ind_con); 
           cout << setw(3) << ind_con+1 << "   " << setw(3) << inh_count << endl; 
           //cout << "Number of times gene " << ind_con+1 << " appears as an inhibitor is: " << inh_count << endl; 
        }
        cout << "------------ END FOR INHIBITORS --------------" << endl;
    }
    cout << "----------Time taken for the process is: "<< setw(18) << setprecision(8)<< proc_time << " -----------" << endl;
    cout << "----------END FOR GENE 4 -----------" << endl;
    return 0;
} 
