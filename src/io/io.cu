#include "../common.h"
using namespace std;
# include "./io.hpp"
//
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
}

void output_data ( const int gene_ind, const double r0_, const double ea_, const double d_, const int start_ka, const int end_ka, const int start_kd, const int end_kd, const vector <int> kavec, const vector <int> kdvec, const host_type &kaval, const host_type &kdval, const vector<string>& out_paths)
{
   // Output r0 
   ofstream r0_fileout( out_paths[0].c_str(), ios_base::binary|ios_base::app|ios_base::out );
   if (r0_fileout.is_open()) {
     r0_fileout << setw(20) << setprecision(12) << r0_ << endl; 
     r0_fileout.close();
   }
   else
     cout << " open() r0_fileout failed" << endl;


   // Output ea 
   ofstream ea_fileout( out_paths[1].c_str(), ios_base::binary|ios_base::app|ios_base::out );
   if (ea_fileout.is_open()) {
     ea_fileout << setw(20) << setprecision(12) << ea_ << endl; 
     ea_fileout.close();
   }
   else
     cout << " open() ea_fileout failed" << endl;


   // Output d 
   ofstream d_fileout( out_paths[2].c_str(), ios_base::binary|ios_base::app|ios_base::out );
   if (d_fileout.is_open()) {
     d_fileout << setw(20) << setprecision(12) << d_ << endl; 
     d_fileout.close();
   }
   else
     cout << " open() d_fileout failed" << endl;


   // Output nka 
   ofstream nka_fileout( out_paths[3].c_str(), ios_base::binary|ios_base::app|ios_base::out );
   if (nka_fileout.is_open()) {
     nka_fileout << setw(3) << (end_ka-start_ka) << endl; 
     nka_fileout.close();
   }
   else
     cout << " open() nka_fileout failed" << endl;


   // Output nkd 
   ofstream nkd_fileout( out_paths[4].c_str(), ios_base::binary|ios_base::app|ios_base::out );
   if (nkd_fileout.is_open()) {
     nkd_fileout << setw(3) << (end_kd-start_kd) << endl; 
     nkd_fileout.close();
   }
   else
     cout << " open() nkd_fileout failed" << endl;


   // Output ka vec
   ofstream kavec_fileout( out_paths[5].c_str(), ios_base::binary|ios_base::app|ios_base::out );
   if (kavec_fileout.is_open())
   {
      for ( int ind = start_ka; ind < end_ka; ind++ )
	kavec_fileout << setw(3) << kavec[ind] << endl; 
      kavec_fileout.close();
   }
   else
     cout << " open() kavec_fileout failed" << endl;

   // Output kd vec
   ofstream kdvec_fileout( out_paths[6].c_str(), ios_base::binary|ios_base::app|ios_base::out );
   if (kdvec_fileout.is_open())
   {
     for ( int ind = start_kd; ind < end_kd; ind++ )
       kdvec_fileout << setw(3) << kdvec[ind] << endl; 
     kdvec_fileout.close();
   }
   else
     cout << " open() kdvec_fileout failed" << endl;

   // Output ka val
   ofstream kaval_fileout( out_paths[7].c_str(), ios_base::binary|ios_base::app|ios_base::out );
   if (kaval_fileout.is_open())
   {
      for ( int ind = start_ka; ind < end_ka; ind++ )
          kaval_fileout << setw(20) << setprecision(12) << kaval[ind] << endl; 
      kaval_fileout.close();
   }
   else
     cout << " open() kaval_fileout failed" << endl;

   // Output kd val
   ofstream kdval_fileout( out_paths[8].c_str(), ios_base::binary|ios_base::app|ios_base::out );
   if (kdval_fileout.is_open())
   {
      for ( int ind = start_kd; ind < end_kd; ind++ )
          kdval_fileout << setw(20) << setprecision(12) << kdval[ind] << endl; 
      kdval_fileout.close();
   }
   else
     cout << " open() kdval_fileout failed" << endl;

}
