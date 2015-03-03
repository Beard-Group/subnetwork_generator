#include "./common.h"
#include "./io/io.hpp"
#include "./spline/spline.hpp"
#include "./lsoda/cuLsoda.hpp"
#include "./optimization/opt.hpp"
#include "./subnetsearch/search.hpp"

using namespace std;
//
const char * const in_file_path_c = "../../../out/N1000/";

int main(void)    
{
    //cout << "----Reading in data----" << endl; 	
    vector <double> t_d(N_time_points);
    vector <double> x_d(N_gene*N_time_points);
    read_data ( N_time_points, N_gene, t_d, x_d);
    //cout << "----Setting up spline co-efficients----" << endl;
    vector <double> deriv_time(N_time_points*N_gene);
    //cout << "  Setting up co-fficients of cubic spline" << endl;
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
    //const int alpha = 0.035;
    //boost::timer timer;
    double proc_time = 0.0;
    //double init_time = 0.0;
    //timer.restart();
    vector <double> act_vec, inh_vec;
    for (int gene_ind = 644; gene_ind < 645; gene_ind++)
    {
        cout << "--------START FOR GENE "<< gene_ind+1 << " -----------" << endl;
        //cout << " Values arranged as error, n_ka, gene ids for ka, n_kd, gene ids for kd," << endl;
        //cout<< " r0, d, ea, ka_values, kd_values for each subnet. " << endl;
        //vector <int> act_vec, inh_vec;
        int subnet_attempts = 0;
        //while (subnet_attempts < probSize) 
        //mt19937 mt_s;
        //mt_s.seed(time(NULL));
        //std::uniform_int_distribution<> dist_in(50,200);
        stringstream ss;
        ss << gene_ind + 1;
        string in_r0_file = in_file_path_c + ss.str() + ".r0";
        string in_ea_file = in_file_path_c + ss.str() + ".ea";
        string in_d_file = in_file_path_c + ss.str() + ".d";
        string in_ka_file = in_file_path_c + ss.str() + ".nka";
        string in_kd_file = in_file_path_c + ss.str() + ".nkd";
        string in_kavec_file = in_file_path_c + ss.str() + ".kavec";
        string in_kdvec_file = in_file_path_c + ss.str() + ".kdvec";
        string in_kaval_file = in_file_path_c + ss.str() + ".kaval";
        string in_kdval_file = in_file_path_c + ss.str() + ".kdval";
        ss.clear();

        vector <int> nka_vec, nkd_vec, kavec, kdvec;
        vector <double> r0_, ea_, d_, kaval, kdval;
        ifstream file_in_r0;
        file_in_r0.open( in_r0_file.c_str(), ios_base::in );
        if (file_in_r0.is_open())
        {
           double a;
           string line;
           while (!file_in_r0.eof())
           {
              file_in_r0 >> a;
              r0_.push_back(a);
           }
           r0_.pop_back();
           file_in_r0.close();
        }
        else cout << " open() failed" << endl;
        ifstream file_in_ea;
        file_in_ea.open( in_ea_file.c_str(), ios_base::in );
        if (file_in_ea.is_open())
        {
           double a;
           string line;
           while (!file_in_ea.eof())
           {
              file_in_ea >> a;
              ea_.push_back(a);
           }
           ea_.pop_back();
           file_in_ea.close();
        }
        else cout << " open() failed" << endl;
        ifstream file_in_d;
        file_in_d.open( in_d_file.c_str(), ios_base::in );
        if (file_in_d.is_open())
        {
           double a;
           string line;
           while (!file_in_d.eof())
           {
              file_in_d >> a;
              d_.push_back(a);
           }
           d_.pop_back();
           file_in_d.close();
        }
        else cout << " open() failed" << endl;
        ifstream file_in_ka;
        file_in_ka.open( in_ka_file.c_str(), ios_base::in );
        if (file_in_ka.is_open())
        {
           int a;
           string line;
           while (!file_in_ka.eof())
           {
              file_in_ka >> a;
              nka_vec.push_back(a);
           }
           nka_vec.pop_back();
           file_in_ka.close();
        }
        else cout << " open() failed" << endl;
        ifstream file_in_kd;
        file_in_kd.open( in_kd_file.c_str(), ios_base::in );
        if (file_in_kd.is_open())
        {
           int a;
           //string line;
           //while ((getline(file_in_kd,line)))
           while (!file_in_kd.eof())
           {
              file_in_kd >> a;
              nkd_vec.push_back(a);
           }
           nkd_vec.pop_back();
           file_in_kd.close();
        }
        else cout << " open() failed" << endl;
        if (nka_vec.size() != nkd_vec.size())  cout << "Error in vector sizes of nka and nkd" << endl;
        else  cout << "Copied nka and nkd files perfectly" << endl;
        const int N_subnets = nka_vec.size();
        cout << "Number of subnets found is: " << N_subnets << endl;
        ifstream file_in_kavec;
        file_in_kavec.open( in_kavec_file.c_str(), ios_base::in );
        if (file_in_kavec.is_open())
        {
           int a;
           //string line;
           //while ((getline(file_in_kavec,line)))
           while(!file_in_kavec.eof())
           {
              file_in_kavec >> a;
              kavec.push_back(a);
           }
           kavec.pop_back(); 
           file_in_kavec.close();
        }
        else cout << " open() failed" << endl;
        //cout << "Size of kavec is: " << kavec.size() << endl;
        ifstream file_in_kdvec;
        file_in_kdvec.open( in_kdvec_file.c_str(), ios_base::in );
        if (file_in_kdvec.is_open())
        {
           int a;
           //string line;
           //while ((getline(file_in_kdvec,line)))
	   while(!file_in_kdvec.eof())
           {
              file_in_kdvec >> a;
              kdvec.push_back(a);
           }
           kdvec.pop_back(); 
           file_in_kdvec.close();
        }
        else cout << " open() failed" << endl;
        //cout << "Size of kdvec is: " << kdvec.size() << endl;
        ifstream file_in_kaval;
        file_in_kaval.open( in_kaval_file.c_str(), ios_base::in );
        if (file_in_kaval.is_open())
        {
           double a;
           //string line;
           //while ((getline(file_in_kaval,line)))
	   while(!file_in_kaval.eof())
           {
              file_in_kaval >> a;
              kaval.push_back(a);
           }
           kaval.pop_back();
           file_in_kaval.close();
        }
        else cout << " open() failed" << endl;
        //cout << "Size of kaval is: " << kaval.size() << endl;
        ifstream file_in_kdval;
        file_in_kdval.open( in_kdval_file.c_str(), ios_base::in );
        if (file_in_kdval.is_open())
        {
           double a;
           //string line;
           //while ((getline(file_in_kdval,line)))
	   while(!file_in_kdval.eof())
           {
              file_in_kdval >> a;
              kdval.push_back(a);
           }
           kdval.pop_back();
           file_in_kdval.close();
        }
        else cout << " open() failed" << endl;
        int sum_act = accumulate(nka_vec.begin(),nka_vec.end(),0);
        cout << "Value of sum_act is: " << sum_act << endl;
        cout << "Value of kaval size is: " << kaval.size() << endl;
        if (sum_act != kaval.size())  cout << "Error in vector sizes of nka and kaval" << endl;
        else  cout << "Copied nka and kaval files perfectly" << endl;
        int sum_inh = accumulate(nkd_vec.begin(),nkd_vec.end(),0);
        if (sum_inh != kdval.size())  cout << "Error in vector sizes of nkd and kdval" << endl;
        else  cout << "Copied nkd and kdval files perfectly" << endl;
        //
        
        int acc_subnet = 0;
        //for (unsigned int k = 0; k < dist_in(mt_s); k++) 
        //    random_select(gene_ind, k+subnet_attempts, mean_xd, n_ka_temp, n_kd_temp, kavec_temp, kdvec_temp, kaval_temp, kdval_temp);
        //subnet_attempts++;
        while (acc_subnet < probSize) 
        //while (subnet_attempts < probSize) 
        {
           int counter_ka = 0;
           int counter_kd = 0;
           for ( int ind = 0; ind < probSize; ind++)
           {
              int n_ka_temp = 0;
              int n_kd_temp = 0;
              vector <int> kavec_temp, kdvec_temp;
              vector <double> kaval_temp, kdval_temp;
              myFex_single fex_nm;
              myJex_single jex_nm;
              fex_nm.set_tau(1.0);
              double r0_temp = r0_[ind];
              //double r0_temp = pow(10,guassrand());
              fex_nm.set_r0(r0_temp);
              double d_temp = d_[ind];
              //double d_temp = pow(10,guassrand());
              fex_nm.set_d(d_temp);
              double ea_temp = ea_[ind];
              //double ea_temp = pow(10,guassrand());
              fex_nm.set_ea(ea_temp);
              n_ka_temp = nka_vec[ind];
              n_kd_temp = nkd_vec[ind];
              fex_nm.set_n_ka(n_ka_temp);
              fex_nm.set_n_kd(n_kd_temp);
              for ( int k = 0; k < nka_vec[ind]; k++) 
              {
                  int ind_ka = counter_ka + k;
                  kavec_temp.push_back(kavec[ind_ka]);
                  kaval_temp.push_back(kaval[ind_ka]);
              } 
              counter_ka += nka_vec[ind];
              for ( int k = 0; k < nkd_vec[ind]; k++) 
              {
                  int ind_kd = counter_kd + k;
                  kdvec_temp.push_back(kdvec[ind_kd]);
                  kdval_temp.push_back(kdval[ind_kd]);
              } 
              counter_kd += nkd_vec[ind];
              fex_nm.set_ka_vec(kavec_temp);
              fex_nm.set_kd_vec(kdvec_temp);
              fex_nm.set_ka_val(kaval_temp);
              fex_nm.set_kd_val(kdval_temp);
              fex_nm.set_coeff(cub_coeff_spline);
              cout << "------------ Starting ODE solver -------------"  << endl;
              // ODE calculation
              //timer.restart();
              double error_ode = 0.0;
              bool out_val = true; 
              integrate_lsoda_ode ( x_d, t_d, mean_xd, fex_nm, jex_nm, gene_ind, out_val, error_ode );
              cout << "Value of error for the first run of ODE is: " << error_ode << endl;
              acc_subnet++; 
	      fex_nm.clear_mem();
              vector<int>().swap(kavec_temp);
              vector<int>().swap(kdvec_temp);
              vector<double>().swap(kaval_temp);
              vector<double>().swap(kdval_temp);
           }
        }
        cout << "--------END FOR GENE "<< gene_ind+1 << " -----------" << endl;
    }
    cout << "----------Time taken for the process is: "<< setw(18) << setprecision(8)<< proc_time << " -----------" << endl;
    
    return 0;
} 
