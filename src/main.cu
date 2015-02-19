#include "./common.h"
#include "./io/io.hpp"
#include "./spline/spline.hpp"
#include "./subnetsearch/search.hpp"
#include "./lsoda/cuLsoda.hpp"
#include "./opt/opt.hpp"

using namespace std;
//
int main(void)    
{  
   cudaSetDevice(2);
   cout << "----Starting simulation----" << endl; 	
   cout << "----Reading in data----" << endl; 	
   vector <double> t_d(N_time_points);
   vector <double> x_d(N_gene*N_time_points);
   read_data (N_time_points, N_gene, t_d, x_d);
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
   vector <double> mean_xd, sd_xd;
   for ( int ind = 0; ind < N_gene; ind++)
   { 
       //x_spline represented in terms of dt = 1. Used to determine mean, sd, and correlation
       // N=1000 case dt=1 for the data anyway so not necessary
       vector <double> x_spline_temp(N_time_points);
       copy( x_d.begin()+ind*N_time_points, x_d.begin()+(ind+1)*N_time_points, x_spline_temp.begin());
       //calc_spline( x_d, t_d, cub_coeff_spline, ind, x_spline_temp);
       calc_stats( x_spline_temp, mean_xd, sd_xd);
       vector<double>().swap(x_spline_temp);
   }
   cout << "----END spline setup----" << endl; 	
   //Error
   double E0 = 3.0E-03;
   for (int id_gene = 644; id_gene < 645; id_gene++)
   {
      const int gene_ind = id_gene;
      vector <double> act_vec, inh_vec;
      vector <int> subnet_size;
      bool r0_out = false;
      bool ea_out = false;
      bool d_out = false;
      bool nka_out = false;
      bool nkd_out = false;
      bool kavec_out = false;
      bool kdvec_out = false;
      bool kaval_out = false;
      bool kdval_out = false;
      char r0_out_file[32] = ".r0";
      char ea_out_file[32] = ".ea";
      char d_out_file[32] = ".d";
      char nka_out_file[32] = ".nka";
      char nkd_out_file[32] = ".nkd";
      char kavec_out_file[32] = ".kavec";
      char kdvec_out_file[32] = ".kdvec";
      char kaval_out_file[32] = ".kaval";
      char kdval_out_file[32] = ".kdval";
      char r0_index_string[64];
      char r0_out_path[96] = "../out/N1000/";
      char ea_index_string[64];
      char ea_out_path[96] = "../out/N1000/";
      char d_index_string[64];
      char d_out_path[96] = "../out/N1000/";
      char nka_index_string[64];
      char nka_out_path[96] = "../out/N1000/";
      char nkd_index_string[64];
      char nkd_out_path[96] = "../out/N1000/";
      char kavec_index_string[64];
      char kavec_out_path[96] = "../out/N1000/";
      char kdvec_index_string[64];
      char kdvec_out_path[96] = "../out/N1000/";
      char kaval_index_string[64];
      char kaval_out_path[96] = "../out/N1000/";
      char kdval_index_string[64];
      char kdval_out_path[96] = "../out/N1000/";
      if (!r0_out)
      {
         r0_out = true;
         sprintf(r0_index_string, "%d", gene_ind+1);
         strcat(r0_index_string,r0_out_file);
         strcat(r0_out_path,r0_index_string);
      }
      string out_r0(r0_out_path);
      if (!ea_out)
      {
         ea_out = true;
         sprintf(ea_index_string, "%d", gene_ind+1);
         strcat(ea_index_string,ea_out_file);
         strcat(ea_out_path,ea_index_string);
      }
      string out_ea(ea_out_path);
      if (!d_out)
      {
         d_out = true;
         sprintf(d_index_string, "%d", gene_ind+1);
         strcat(d_index_string,d_out_file);
         strcat(d_out_path,d_index_string);
      }
      string out_d(d_out_path);
      if (!nka_out)
      {
         nka_out = true;
         sprintf(nka_index_string, "%d", gene_ind+1);
         strcat(nka_index_string,nka_out_file);
         strcat(nka_out_path,nka_index_string);
      }
      string out_nka(nka_out_path);
      if (!nkd_out)
      {
         nkd_out = true;
         sprintf(nkd_index_string, "%d", gene_ind+1);
         strcat(nkd_index_string,nkd_out_file);
         strcat(nkd_out_path,nkd_index_string);
      }
      string out_nkd(nkd_out_path);
      if (!kavec_out)
      {
         kavec_out = true;
         sprintf(kavec_index_string, "%d", gene_ind+1);
         strcat(kavec_index_string,kavec_out_file);
         strcat(kavec_out_path,kavec_index_string);
      }
      string out_kavec(kavec_out_path);
      if (!kdvec_out)
      {
         kdvec_out = true;
         sprintf(kdvec_index_string, "%d", gene_ind+1);
         strcat(kdvec_index_string,kdvec_out_file);
         strcat(kdvec_out_path,kdvec_index_string);
      }
      string out_kdvec(kdvec_out_path);
      if (!kaval_out)
      {
         kaval_out = true;
         sprintf(kaval_index_string, "%d", gene_ind+1);
         strcat(kaval_index_string,kaval_out_file);
         strcat(kaval_out_path,kaval_index_string);
      }
      string out_kaval(kaval_out_path);
      if (!kdval_out)
      {
         kdval_out = true;
         sprintf(kdval_index_string, "%d", gene_ind+1);
         strcat(kdval_index_string,kdval_out_file);
         strcat(kdval_out_path,kdval_index_string);
      }
      string out_kdval(kdval_out_path);
      // END DECLARATIONS 
      int acc_subnet = 0;
      int subnet_att = 0;
      cout << "------------ Start for gene " << gene_ind+1 << "  ------" << endl;
      while ( (acc_subnet < N_att) && (subnet_att<50) )
      {
         if (subnet_att <= 500) E0=5.0E-03;
         else E0=7.0E-03;
         cout << "------------ Attempt number " << subnet_att << "  ------" << endl;
         subnet_att++;
         int *n_ka = (int*)malloc(sizeof(int)*probSize);
         int *n_kd = (int*)malloc(sizeof(int)*probSize);
         int *ka_start = (int*)malloc(sizeof(int)*probSize);
         int *kd_start = (int*)malloc(sizeof(int)*probSize);
         int ka_ind = 0;
         int kd_ind = 0;
         vector <int> temp_ka_vec, temp_kd_vec;
         vector <double> temp_ka_val, temp_kd_val;
         srand(time(NULL));
         // generate integers between 50 to 200
         thrust::default_random_engine rng_gen(time(NULL)*rand());
         thrust::uniform_int_distribution<int> u_l(50,200); 
         for ( unsigned int ind = 0; ind < probSize; ind++ )
         {
            int n_ka_temp = 0;
            int n_kd_temp = 0;
            // How to declare kavec_temp and kd_vec_temp
            vector <int> kavec_temp, kdvec_temp;
            vector <double> kaval_temp, kdval_temp;
            for (unsigned int k = 0; k < u_l(rng_gen); k++) 
                random_conn(gene_ind, k+ind, mean_xd, n_ka_temp, n_kd_temp, kavec_temp, kdvec_temp, kaval_temp, kdval_temp);
            *(n_ka+ind) = n_ka_temp;
            *(n_kd+ind) = n_kd_temp;
            *(ka_start+ind) = ka_ind;
            *(kd_start+ind) = kd_ind;
            ka_ind += n_ka_temp;
            kd_ind += n_kd_temp;
            temp_ka_vec.insert(temp_ka_vec.end(),kavec_temp.begin(),kavec_temp.end()); 
            temp_kd_vec.insert(temp_kd_vec.end(),kdvec_temp.begin(),kdvec_temp.end()); 
            temp_ka_val.insert(temp_ka_val.end(),kaval_temp.begin(),kaval_temp.end()); 
            temp_kd_val.insert(temp_kd_val.end(),kdval_temp.begin(),kdval_temp.end());
            vector <int>().swap(kavec_temp);
            vector <int>().swap(kdvec_temp);
            vector <double>().swap(kaval_temp);
            vector <double>().swap(kdval_temp);
         }
         //
         const int size_ka = temp_ka_vec.size();
         const int size_kd = temp_kd_vec.size();
         host_type error_sim_h;
         int *ka_vec = (int*)malloc(sizeof(int)*size_ka);
         double *ka_val = (double*)malloc(sizeof(double)*size_ka);
         int l = 0;
         while ( l < size_ka )
         {
             *(ka_vec+l) = temp_ka_vec[l];
             *(ka_val+l) = temp_ka_val[l];
             l++;
         } 
         l = 0;
         int *kd_vec = (int*)malloc(sizeof(int)*size_kd);
         double *kd_val = (double*)malloc(sizeof(double)*size_kd);
         while ( l < size_kd )
         {
             *(kd_vec+l) = temp_kd_vec[l];
             *(kd_val+l) = temp_kd_val[l];
             l++;
         } 
         double *r0 = (double*)malloc(sizeof(double)*probSize);
         double *d = (double*)malloc(sizeof(double)*probSize);
         double *ea = (double*)malloc(sizeof(double)*probSize);
         srand(time(NULL));
         for ( int i = 0; i < probSize; i++)
         {
               *(r0+i) = mean_xd[gene_ind];
               *(d+i) = 1.0;
               *(ea+i) = 1.0;
         }
         double cub_coeff[N_gene*N_time_points*4];
         for ( int i = 0; i < cub_coeff_spline.size(); i++) cub_coeff[i] = cub_coeff_spline[i];
         myFex fex;
         fex.set_r0(r0);
         fex.set_d(d);
         fex.set_ea(ea);
         fex.set_n_ka(n_ka);
         fex.set_ka_vec(ka_vec,size_ka);
         fex.set_ka_start(ka_start);
         fex.set_ka_val(ka_val,size_ka);
         fex.set_n_kd(n_kd);
         fex.set_kd_start(kd_start);
         fex.set_kd_vec(kd_vec,size_kd);
         fex.set_kd_val(kd_val,size_kd);
         int size_coeff = N_gene*N_time_points*4;
         fex.set_coeff(cub_coeff, size_coeff); 
         myJex jex;
         state_type error_sim_d(probSize);
         // Integrate ODEs
         integrate_lsoda_ode (gene_ind, x_d, t_d, mean_xd[gene_ind], fex, jex, error_sim_d);
         // MC Simulation
         MC_sim(gene_ind, x_d, t_d, mean_xd[gene_ind], n_ka, n_kd, size_ka, size_kd, fex, jex, error_sim_d);
         host_type ka_val_h(size_ka), kd_val_h(size_kd);
         host_type r0_mc(probSize), d_mc(probSize), ea_mc(probSize);
         thrust::device_ptr<double> r0_mc_d, d_mc_d, ea_mc_d;
         fex.get_r0_ptr(r0_mc_d);
         fex.get_ea_ptr(ea_mc_d);
         fex.get_d_ptr(d_mc_d);
         thrust::copy(r0_mc_d, r0_mc_d+probSize, r0_mc.begin()); 
         thrust::copy(ea_mc_d, ea_mc_d+probSize, ea_mc.begin()); 
         thrust::copy(d_mc_d, d_mc_d+probSize, d_mc.begin()); 
         fex.get_kaval_vec(ka_val_h,size_ka);
         fex.get_kdval_vec(kd_val_h,size_kd);
         error_sim_h = error_sim_d;
         fex.set_r0_free();
         fex.set_d_free();
         fex.set_ea_free();
         fex.set_n_ka_free();
         fex.set_ka_vec_free();
         fex.set_ka_start_free();
         fex.set_ka_val_free();
         fex.set_n_kd_free();
         fex.set_kd_start_free();
         fex.set_kd_vec_free();
         fex.set_kd_val_free();
         fex.set_coeff_free(); 
         for ( int ind = 0; ind < probSize; ind++)
         {
             if ( error_sim_h[ind] < E0 ) 
             {
                    bool acc_sub = false; 
                    int act_count = 0;
                    int inh_count = 0;
                    int start_ka = *(ka_start+ind);
                    int end_ka = *(ka_start+ind)+*(n_ka+ind);
                    int start_kd = *(kd_start+ind);
                    int end_kd = *(kd_start+ind)+*(n_kd+ind);
                    double r0_ = r0_mc[ind];
                    double ea_ = ea_mc[ind];
                    double d_ = d_mc[ind];
                    output_data ( gene_ind, r0_, ea_, d_, out_r0, out_ea, out_d, start_ka, end_ka, start_kd, end_kd, temp_ka_vec, temp_kd_vec, ka_val_h, kd_val_h, out_nka, out_nkd, out_kavec, out_kdvec, out_kaval, out_kdval);
                    for ( int k = *(ka_start+ind); k < (*(ka_start+ind)+*(n_ka+ind)); k++) 
                    {
                       if ( ka_val_h[k] > 0.05 ) 
                       {  
                          act_vec.push_back(temp_ka_vec[k]);   
                          acc_sub = true;
                          act_count++;
                       } 
                    }
                    for ( int k = *(kd_start+ind); k < (*(kd_start+ind)+*(n_kd+ind)); k++) 
                    {
                       if ( kd_val_h[k] > 0.05 ) 
                       {   
                          inh_vec.push_back(temp_kd_vec[k]);   
                          acc_sub = true;
                          inh_count++;
                       }
                    }
                    if (acc_sub) 
                    {
                       acc_subnet++;
                       subnet_size.push_back(act_count+inh_count);
                       cout << "-------- Acceptable subnet number: " << acc_subnet << "  ------------" << endl; 
                       cout << "-------- Size of subnetwork: " << act_count+inh_count << "  ------------" << endl; 
                       cout << "-------- Value of error is : " << error_sim_h[ind] << "  ------------" << endl; 
                       cout << "-------- Number of attempts made till now: " << subnet_att << "  ------------" << endl; 
                    }
             } 
         }
         error_sim_h.clear();
         error_sim_h.shrink_to_fit();
         vector <int>().swap(temp_ka_vec);
         vector <int>().swap(temp_kd_vec);
         vector <double>().swap(temp_ka_val);
         vector <double>().swap(temp_kd_val);
         free(r0);
         free(d);
         free(ea);
         free(n_ka);
         free(n_kd);
         free(ka_start);
         free(kd_start);
         free(ka_vec);
         free(kd_vec);
         free(ka_val);
         free(kd_val);
      }
      char out_file[32] = ".out";
      char index_string[64];
      sprintf(index_string, "%d", gene_ind+1);
      strcat(index_string,out_file);
      char out_path[96] = "../out/N1000/";
      strcat(out_path,index_string);
      ofstream fileout;
      fileout.open( out_path, ios_base::binary|ios_base::app|ios_base::out );
      if (fileout.is_open())
      {
         cout << "File open successful!" << endl;
         fileout << "------------ START FOR ACTIVATORS --------------" << endl;
         for ( int ind_con = 0; ind_con < N_gene; ind_con++ )
         {
            int act_count = count(act_vec.begin(),act_vec.end(),ind_con);
            fileout << setw(3) << ind_con+1 << "   " << setw(3) << act_count << endl; 
         }
         fileout << "-------------- END FOR ACTIVATORS --------------" << endl;
         for ( int ind_con = 0; ind_con < N_gene; ind_con++ )
         {
            int inh_count = count(inh_vec.begin(),inh_vec.end(),ind_con); 
            fileout << setw(3) << ind_con+1 << "   " << setw(3) << inh_count << endl; 
         }
         fileout << "------------ Size of subnetworks START----------- " << endl;
         for ( int ind_con = 0; ind_con < 20; ind_con++ )
         {
            int size_count = count(subnet_size.begin(),subnet_size.end(),ind_con); 
            fileout << setw(3) << ind_con << "   " << setw(3) << size_count << endl; 
         }
         fileout << "------------ Size of subnetworks END----------- " << endl;
         fileout << "------------ Number of acceptable subnetworks for gene " << gene_ind+1 <<
          " is: " << acc_subnet << " --------------" << endl;
         fileout << "------------ Number of attempts for gene " << gene_ind+1 <<
          " is: " << subnet_att*probSize << " --------------" << endl;
         fileout.close();
      }
      else cout << " open() failed" << endl;
      vector <double>().swap(act_vec);
      vector <double>().swap(inh_vec);
      vector <int>().swap(subnet_size);
      cout << "------------ End for gene " << gene_ind+1 << "  ------" << endl;
   }
   //
   return 0;
} /* MAIN */
