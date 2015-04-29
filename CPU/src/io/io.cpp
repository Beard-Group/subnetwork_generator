#include "../common.h"
using namespace std;
# include "./io.hpp"
//
void read_data ( const int N_time_points, const int N_gene, vector<double> &t_d, vector <double> &x_d )
//void read_data ( const int N_time_points, const int N_gene, vector<double> &t_d, vector <double> &x_d, vector<double> &sd_d)
{
   // Read in time steps t_d
   ifstream FileInpt_d;
   char path1[100] = "../in/t_d";
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
   char path2[100] = "../in/x_d_1000_N8_250_e1";
   FileInpx_dt.open( path2, ios_base::in);
   if (FileInpx_dt.is_open())
   {
      cout << "File x_dt open successful." << endl; 
      for ( int i = 0; i < N_time_points * N_gene ; i++ )
      {
           FileInpx_dt >> x_d[i];
          // cout << "Element "<< i << " is: "<< x_d[i] << endl;  
      }
      cout << "Finished reading file x_dt, will close now." << endl;
      FileInpx_dt.close();
   }
   else
        cout << "open() failed" << endl;
}
