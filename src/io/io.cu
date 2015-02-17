#include "../common.h"
using namespace std;
# include "./io.hpp"
//
//void read_data ( const int N_time_points, const int N_gene, vector<double> &t_d, vector <double> &x_d, vector<double> &sd_d)
void read_data ( const int N_time_points, const int N_gene, vector<double> &t_d, vector <double> &x_d)
{
   // Read in time steps t_d
   ifstream FileInpt_d;
   char path1[100] = "../in/t_d_1000";
   //char path1[100] = "/tmp/raghut/GPU_g1g50/inp/t_d_1000";
   FileInpt_d.open( path1, ios_base::in);
   if (FileInpt_d.is_open())
   {
      cout << "File open successful. It contains: " << endl; 
      for ( int i = 0; i < N_time_points; i++ )
      {
            FileInpt_d >> t_d[i];
            //cout << t_d[i] << endl;  
      }
      cout << "Finished reading file, will close now." << endl;
      FileInpt_d.close();
   }
   else
        cout << "open() failed" << endl;
   // Read in time steps expression profile of all genes 
   ifstream FileInpx_dt;
   char path2[100] = "../in/x_d_1000_N8_200_e1";
   //char path2[100] = "/tmp/raghut/GPU_g1g50/inp/x_d_1000_N7_e1";
   FileInpx_dt.open( path2, ios_base::in);
   if (FileInpx_dt.is_open())
   {
      cout << "File x_dt open successful." << endl; 
      for ( int i = 0; i < N_time_points * N_gene ; i++ )
      {
           FileInpx_dt >> x_d[i];
           //cout << "Element "<< i << " is: "<< x_d[i] << endl;  
      }
      cout << "Finished reading file x_dt, will close now." << endl;
      FileInpx_dt.close();
   }
   else
        cout << "open() failed" << endl;
   // Read in standard deviation in expression profiles
   /*ifstream FileInpsd_d;
   FileInpsd_d.open("/home/rthiagarajan/workspace/lsoda/MCsim/update/in/sd_d", ios_base::in);
   if (FileInpsd_d.is_open())
   {
      cout << "File sd_d open successful." << endl; 
      for ( int i = 0; i < N_time_points * N_gene ; i++ )
      {
           FileInpsd_d >> sd_d[i];
           //cout << sd_d[i] << endl;  
      }
      cout << "Finished reading file sd_d, will close now." << endl;
      FileInpsd_d.close();
   }
   else
        cout << "open() failed" << endl;*/
}
void read_in_subnet ( const int gene_ind, vector <int> &nka_vec, vector <int> &nkd_vec, vector <int> &kavec, vector <int> &kdvec)
{
     char in_ka_file[32] = ".nkd";
     char in_kd_file[32] = ".nka";
     char in_kavec_file[32] = ".kavec";
     char in_kdvec_file[32] = ".kdvec";
     char index_string_ka[32];
     char index_string_kd[32];
     char index_string_kavec[32];
     char index_string_kdvec[32];
     char in_ka_path[64];
     char in_kd_path[64];
     char in_kavec_path[64];
     char in_kdvec_path[64];
     sprintf(index_string_ka, "%d", gene_ind+1);
     sprintf(index_string_kd, "%d", gene_ind+1);
     sprintf(index_string_kavec, "%d", gene_ind+1);
     sprintf(index_string_kdvec, "%d", gene_ind+1);
     strcat(index_string_ka,in_ka_file);
     strcat(index_string_kd,in_kd_file);
     strcat(index_string_kavec,in_kavec_file);
     strcat(index_string_kdvec,in_kdvec_file);
     char file_path_ka[32] = "../in/200eas/";
     char file_path_kd[32] = "../in/200eas/";
     char file_path_kavec[32] = "../in/200eas/";
     char file_path_kdvec[32] = "../in/200eas/";
     strcat(file_path_ka,index_string_ka);
     strcat(file_path_kd,index_string_kd);
     strcat(file_path_kavec,index_string_kavec);
     strcat(file_path_kdvec,index_string_kdvec);
     strncpy(in_ka_path,file_path_ka,sizeof(file_path_ka));
     strncpy(in_kd_path,file_path_kd,sizeof(file_path_kd));
     strncpy(in_kavec_path,file_path_kavec,sizeof(file_path_kavec));
     strncpy(in_kdvec_path,file_path_kdvec,sizeof(file_path_kdvec));
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
     //se  cout << "Copied nka and nkd files perfectly" << endl;
     //
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
   
} 
void output_data ( const int gene_ind, const double r0_, const double ea_, const double d_, const string out_r0, const string out_ea, const string out_d, const int start_ka, const int end_ka, const int start_kd, const int end_kd, const vector <int> kavec, const vector <int> kdvec, const host_type &kaval, const host_type &kdval, const string out_nka, const string out_nkd, const string kavec_out_path, const string kdvec_out_path, const string kaval_out_path, const string kdval_out_path)
{
   // Output r0 
   const char *r0_out = out_r0.c_str();
   ofstream r0_fileout;
   r0_fileout.open( r0_out, ios_base::binary|ios_base::app|ios_base::out );
   if (r0_fileout.is_open())
          r0_fileout << setw(20) << setprecision(12) << r0_ << endl; 
   else cout << " open() r0_fileout failed" << endl;
   // Output ea 
   const char *ea_out = out_ea.c_str();
   ofstream ea_fileout;
   ea_fileout.open( ea_out, ios_base::binary|ios_base::app|ios_base::out );
   if (ea_fileout.is_open())
          ea_fileout << setw(20) << setprecision(12) << ea_ << endl; 
   else cout << " open() ea_fileout failed" << endl;
   // Output d 
   const char *d_out = out_d.c_str();
   ofstream d_fileout;
   d_fileout.open( d_out, ios_base::binary|ios_base::app|ios_base::out );
   if (d_fileout.is_open())
          d_fileout << setw(20) << setprecision(12) << d_ << endl; 
   else cout << " open() d_fileout failed" << endl;
   // Output nka 
   const char *nka_out = out_nka.c_str();
   ofstream nka_fileout;
   nka_fileout.open( nka_out, ios_base::binary|ios_base::app|ios_base::out );
   if (nka_fileout.is_open())
          nka_fileout << setw(3) << (end_ka-start_ka) << endl; 
   else cout << " open() nka_fileout failed" << endl;
   // Output nkd 
   const char *nkd_out = out_nkd.c_str();
   ofstream nkd_fileout;
   nkd_fileout.open( nkd_out, ios_base::binary|ios_base::app|ios_base::out );
   if (nkd_fileout.is_open())
          nkd_fileout << setw(3) << (end_kd-start_kd) << endl; 
   else cout << " open() nkd_fileout failed" << endl;
   // Output ka vec
   const char *kavec_out = kavec_out_path.c_str();
   ofstream kavec_fileout;
   kavec_fileout.open( kavec_out, ios_base::binary|ios_base::app|ios_base::out );
   if (kavec_fileout.is_open())
   {
      for ( int ind = start_ka; ind < end_ka; ind++ )
          kavec_fileout << setw(3) << kavec[ind] << endl; 
      kavec_fileout.close();
   }
   else cout << " open() kavec_fileout failed" << endl;
   // Output kd vec
   const char *kdvec_out = kdvec_out_path.c_str();
   ofstream kdvec_fileout;
   kdvec_fileout.open( kdvec_out, ios_base::binary|ios_base::app|ios_base::out );
   if (kdvec_fileout.is_open())
   {
      for ( int ind = start_kd; ind < end_kd; ind++ )
          kdvec_fileout << setw(3) << kdvec[ind] << endl; 
      kdvec_fileout.close();
   }
   else cout << " open() kdvec_fileout failed" << endl;
   // Output ka val
   const char *kaval_out = kaval_out_path.c_str();
   ofstream kaval_fileout;
   kaval_fileout.open( kaval_out, ios_base::binary|ios_base::app|ios_base::out );
   if (kaval_fileout.is_open())
   {
      for ( int ind = start_ka; ind < end_ka; ind++ )
          kaval_fileout << setw(20) << setprecision(12) << kaval[ind] << endl; 
      kaval_fileout.close();
   }
   else cout << " open() kaval_fileout failed" << endl;
   // Output kd val
   const char *kdval_out = kdval_out_path.c_str();
   ofstream kdval_fileout;
   kdval_fileout.open( kdval_out, ios_base::binary|ios_base::app|ios_base::out );
   if (kdval_fileout.is_open())
   {
      for ( int ind = start_kd; ind < end_kd; ind++ )
          kdval_fileout << setw(20) << setprecision(12) << kdval[ind] << endl; 
      kdval_fileout.close();
   }
   else cout << " open() kdval_fileout failed" << endl;
}
