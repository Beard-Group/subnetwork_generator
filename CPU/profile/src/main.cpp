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
    boost::timer timer;
    double proc_time = 0.0;
    //double init_time = 0.0;
    timer.restart();
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
        char in_r0_file[32] = ".r0";
        char in_ea_file[32] = ".ea";
        char in_d_file[32] = ".d";
        char in_ka_file[32] = ".nka";
        char in_kd_file[32] = ".nkd";
        char in_kavec_file[32] = ".kavec";
        char in_kdvec_file[32] = ".kdvec";
        char in_kaval_file[32] = ".kaval";
        char in_kdval_file[32] = ".kdval";
        char index_string_r0[32];
        char index_string_ea[32];
        char index_string_d[32];
        char index_string_ka[32];
        char index_string_kd[32];
        char index_string_kavec[32];
        char index_string_kdvec[32];
        char index_string_kaval[32];
        char index_string_kdval[32];
        char in_r0_path[64];
        char in_ea_path[64];
        char in_d_path[64];
        char in_ka_path[64];
        char in_kd_path[64];
        char in_kavec_path[64];
        char in_kdvec_path[64];
        char in_kaval_path[64];
        char in_kdval_path[64];
        sprintf(index_string_r0, "%d", gene_ind+1);
        sprintf(index_string_ea, "%d", gene_ind+1);
        sprintf(index_string_d, "%d", gene_ind+1);
        sprintf(index_string_ka, "%d", gene_ind+1);
        sprintf(index_string_kd, "%d", gene_ind+1);
        sprintf(index_string_kavec, "%d", gene_ind+1);
        sprintf(index_string_kdvec, "%d", gene_ind+1);
        sprintf(index_string_kaval, "%d", gene_ind+1);
        sprintf(index_string_kdval, "%d", gene_ind+1);
        strcat(index_string_r0,in_r0_file);
        strcat(index_string_ea,in_ea_file);
        strcat(index_string_d,in_d_file);
        strcat(index_string_ka,in_ka_file);
        strcat(index_string_kd,in_kd_file);
        strcat(index_string_kavec,in_kavec_file);
        strcat(index_string_kdvec,in_kdvec_file);
        strcat(index_string_kaval,in_kaval_file);
        strcat(index_string_kdval,in_kdval_file);
        char file_path_r0[32] = "../in/";
        char file_path_ea[32] = "../in/";
        char file_path_d[32] = "../in/";
        char file_path_ka[32] = "../in/";
        char file_path_kd[32] = "../in/";
        char file_path_kavec[32] = "../in/";
        char file_path_kdvec[32] = "../in/";
        char file_path_kaval[32] = "../in/";
        char file_path_kdval[32] = "../in/";
        strcat(file_path_r0,index_string_r0);
        strcat(file_path_ea,index_string_ea);
        strcat(file_path_d,index_string_d);
        strcat(file_path_ka,index_string_ka);
        strcat(file_path_kd,index_string_kd);
        strcat(file_path_kavec,index_string_kavec);
        strcat(file_path_kdvec,index_string_kdvec);
        strcat(file_path_kaval,index_string_kaval);
        strcat(file_path_kdval,index_string_kdval);
        strncpy(in_r0_path,file_path_r0,sizeof(file_path_r0));
        strncpy(in_ea_path,file_path_ea,sizeof(file_path_ea));
        strncpy(in_d_path,file_path_d,sizeof(file_path_d));
        strncpy(in_ka_path,file_path_ka,sizeof(file_path_ka));
        strncpy(in_kd_path,file_path_kd,sizeof(file_path_kd));
        strncpy(in_kavec_path,file_path_kavec,sizeof(file_path_kavec));
        strncpy(in_kdvec_path,file_path_kdvec,sizeof(file_path_kdvec));
        strncpy(in_kaval_path,file_path_kaval,sizeof(file_path_kaval));
        strncpy(in_kdval_path,file_path_kdval,sizeof(file_path_kdval));
        vector <int> nka_vec, nkd_vec, kavec, kdvec;
        vector <double> r0_, ea_, d_, kaval, kdval;
        ifstream file_in_r0;
        file_in_r0.open( in_r0_path, ios_base::in );
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
        file_in_ea.open( in_ea_path, ios_base::in );
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
        file_in_d.open( in_d_path, ios_base::in );
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
        file_in_ka.open( in_ka_path, ios_base::in );
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
        file_in_kd.open( in_kd_path, ios_base::in );
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
        file_in_kavec.open( in_kavec_path, ios_base::in );
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
        file_in_kdvec.open( in_kdvec_path, ios_base::in );
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
        file_in_kaval.open( in_kaval_path, ios_base::in );
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
        file_in_kdval.open( in_kdval_path, ios_base::in );
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
              timer.restart();
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
        /*while (acc_subnet < probSize) 
        //while (subnet_attempts < probSize) 
        {
           int counter_ka = 0;
           int counter_kd = 0;
           for ( int ind = 0; ind < N_subnets; ind++)
           {
              int n_ka_temp = 0;
              int n_kd_temp = 0;
              vector <int> kavec_temp, kdvec_temp;
              vector <double> kaval_temp, kdval_temp;
              myFex_single fex_nm;
              myJex_single jex_nm;
              fex_nm.set_tau(1.0);
              //double r0_temp = mean_xd[gene_ind];
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
              for ( int k = 0; k < nkd_vec[ind]; k++) 
              {
                  int ind_kd = counter_kd + k;
                  kdvec_temp.push_back(kdvec[ind_kd]);
                  kdval_temp.push_back(kdval[ind_kd]);
              } 
              fex_nm.set_ka_vec(kavec_temp);
              fex_nm.set_kd_vec(kdvec_temp);
              fex_nm.set_ka_val(kaval_temp);
              fex_nm.set_kd_val(kdval_temp);
              fex_nm.set_coeff(cub_coeff_spline);
              cout << "------------ Starting ODE solver -------------"  << endl;
              // ODE calculation
              timer.restart();
              double error_ode = 0.0;
              bool out_val = false; 
              integrate_lsoda_ode ( x_d, t_d, mean_xd, fex_nm, jex_nm, gene_ind, out_val, error_ode );
              cout << "Value of error for the first run of ODE is: " << error_ode << endl;
              state_type param;
              param.push_back(r0_temp);
              param.push_back(d_temp);
              param.push_back(ea_temp);
              param.insert(param.end(),kaval_temp.begin(),kaval_temp.end());
              param.insert(param.end(),kdval_temp.begin(),kdval_temp.end());
              // MC and NM search conducted in log space so transform param to log10 space
	      cout << "------------ Starting MC solver -------------"  << endl;
              transform(param.begin(),param.end(),param.begin(),logb10());
              MC_sim_new ( x_d, t_d, mean_xd, fex_nm, jex_nm, gene_ind, n_ka_temp, n_kd_temp, error_ode, param);
              cout << "Value of error for the MC run is: " << error_ode << endl;
              proc_time += timer.elapsed();
	      cout << "------------ Ending MC solver -------------"  << endl;
              //error_vec.push_back(error_ode); 
              //if (N_itr==0) time_vec.push_back(timer.elapsed());
              //if ( error_vec[ind] < E_0 ) 
              //if (( error_ode < E_0 ) && ((param.size()<8) && (param.size()>3)))
	      fex_nm.clear_mem();
              //if ( (error_ode < E_0) && ((param.size()>3) && (param.size()<8)) )
              //if ((error_ode < E_0) && (param.size()<10))
              if (error_ode < E_0)
              {
                     //act_vec.insert(act_vec.end(),temp_ka_vec.begin()+*(ka_start+ind),temp_ka_vec.begin()+*(ka_start+ind)+*(n_ka+ind)); 
                     //inh_vec.insert(inh_vec.end(),temp_kd_vec.begin()+*(kd_start+ind),temp_kd_vec.begin()+*(kd_start+ind)+*(n_kd+ind)); 
                     //act_vec.insert(act_vec.end(),kavec_temp.begin(),kavec_temp.end()); 
                     //inh_vec.insert(inh_vec.end(),kdvec_temp.begin(),kdvec_temp.end());
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
                        cout << "-------- Acceptable subnet number: " << acc_subnet << "  ------------" << endl; 
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
                     cout << "Value of error for the first run of ODE is: " << error_val << endl;
              } 
              vector<double>().swap(param);
              vector<int>().swap(kavec_temp);
              vector<int>().swap(kdvec_temp);
              vector<double>().swap(kaval_temp);
              vector<double>().swap(kdval_temp);
           }
        }*/
	//cout << "------------ ADD MC solver -------------"  << endl;
    }
    //}
    /*cout << "------------ START FOR ACTIVATORS --------------" << endl;
    for ( int ind_con = 0; ind_con < N_gene; ind_con++ )
    {
       int act_count = count(act_vec.begin(),act_vec.end(),ind_con);
       cout << setw(3) << ind_con+1 << "   " << setw(3) << act_count << endl; 
    }
    cout << "-------------- END FOR ACTIVATORS --------------" << endl;
    for ( int ind_con = 0; ind_con < N_gene; ind_con++ )
    {
       int inh_count = count(inh_vec.begin(),inh_vec.end(),ind_con); 
       cout << setw(3) << ind_con+1 << "   " << setw(3) << inh_count << endl; 
       //cout << "Number of times gene " << ind_con+1 << " appears as an inhibitor is: " << inh_count << endl; 
    }*/
    //cout << "----------Size of error vector is: "<< error_vec.size() << " --------------------" << endl;
    //cout << "----------Number of Acceptable Subnets found: "<< acc_subnet << " -----------" << endl;
    //cout << "----------Time taken for the process is: "<< setw(18) << setprecision(8)<< double(proc_time/double(N_tot)) << " -----------" << endl;
    cout << "----------Time taken for the process is: "<< setw(18) << setprecision(8)<< proc_time << " -----------" << endl;
    //for (int i = 0; i < time_vec.size(); i++) cout << setw(4) << i << "  " << setw(18) << setprecision(8) << time_vec[i] << endl;
    return 0;
} 
