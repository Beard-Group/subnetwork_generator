#include "./common.h"
#include "./io/io.hpp"
#include "./spline/spline.hpp"
#include "./subnetsearch/search.hpp"
#include "./lsoda/cuLsoda.hpp"
#include "./opt/opt.hpp"

#include <sstream>
using namespace std;

// Global constants, single point of maintenance
const char * const output_path_c = "../out/N1000_test/";

void setup_spline(vector<double>& t_d_in, vector<double> x_d_in,
                  vector<double>& cub_coeff_spline_out, vector<double>& mean_xd_out);

void setup_output(int gene_ind_n, vector<string>& out_files_out);

int main(void)
{
    cudaSetDevice(2);
    cout << "----Starting simulation----" << endl;
    cout << "----Reading in data----" << endl;
    vector <double> t_d(N_time_points);
    vector <double> x_d(N_gene*N_time_points);
    read_data (N_time_points, N_gene, t_d, x_d);
    vector <double> cub_coeff_spline (N_gene*N_time_points*4);
    vector<double> mean_xd;
    setup_spline(t_d, x_d, cub_coeff_spline, mean_xd);
    //Error
    double E0 = 3.0E-03;
    for (int id_gene = 644; id_gene < 645; id_gene++) {
        const int gene_ind = id_gene;
        vector <double> act_vec, inh_vec;
        vector <int> subnet_size;
        // Set up output file paths for this gene
        bool created_paths = false; // Why do we need this check? Paths only created once per gene index anyway.
        vector<string> out_files;
        if (!created_paths) {
            created_paths = true;
            setup_output(gene_ind, out_files);
        }
        // END DECLARATIONS
        int acc_subnet = 0; // count of accepted subnetworks
        int subnet_att = 0; // attempts
        cout << "------------ Start for gene " << gene_ind+1 << "  ------" << endl;
        while ( (acc_subnet < N_att) && (subnet_att<50) ) {
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
            // Problem size loop, 1600 subnetworks
            for ( unsigned int ind = 0; ind < probSize; ind++ ) {
                int n_ka_temp = 0;
                int n_kd_temp = 0;
                // How to declare kavec_temp and kd_vec_temp
                vector <int> kavec_temp, kdvec_temp;
                vector <double> kaval_temp, kdval_temp;
                for (unsigned int k = 0; k < u_l(rng_gen); k++)
                    random_conn(gene_ind, k+ind, mean_xd, n_ka_temp, n_kd_temp, kavec_temp,
                                kdvec_temp, kaval_temp, kdval_temp); // subnetsearch subfolder
                // Parameters
                *(n_ka+ind) = n_ka_temp; // n_ka number of activating genes
                *(n_kd+ind) = n_kd_temp; // n_kd number of deactivating genes
                *(ka_start+ind) = ka_ind;
                *(kd_start+ind) = kd_ind;
                ka_ind += n_ka_temp;
                kd_ind += n_kd_temp;
                // Here we select which are the activators and deactivators
                temp_ka_vec.insert(temp_ka_vec.end(),kavec_temp.begin(),kavec_temp.end());
                temp_kd_vec.insert(temp_kd_vec.end(),kdvec_temp.begin(),kdvec_temp.end());
                // The corresponding values go in these two
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
            while ( l < size_ka ) {
                *(ka_vec+l) = temp_ka_vec[l];
                *(ka_val+l) = temp_ka_val[l];
                l++;
            }
            l = 0;
            int *kd_vec = (int*)malloc(sizeof(int)*size_kd);
            double *kd_val = (double*)malloc(sizeof(double)*size_kd);
            while ( l < size_kd ) {
                *(kd_vec+l) = temp_kd_vec[l];
                *(kd_val+l) = temp_kd_val[l];
                l++;
            }
            double *r0 = (double*)malloc(sizeof(double)*probSize);
            double *d = (double*)malloc(sizeof(double)*probSize);
            double *ea = (double*)malloc(sizeof(double)*probSize);
            srand(time(NULL));
            for ( int i = 0; i < probSize; i++) {
                *(r0+i) = mean_xd[gene_ind];
                *(d+i) = 1.0;
                *(ea+i) = 1.0;
            }
            double cub_coeff[N_gene*N_time_points*4];
            for ( int i = 0; i < cub_coeff_spline.size();
                    i++) cub_coeff[i] = cub_coeff_spline[i];
            myFex fex; // ./lsoda/cuLSODA.hpp
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
            // Integrate ODEs, gives us an error value
            integrate_lsoda_ode (gene_ind, x_d, t_d, mean_xd[gene_ind], fex, jex,
                                 error_sim_d);
            // MC Simulation, changes parameters, generates new error value
            MC_sim(gene_ind, x_d, t_d, mean_xd[gene_ind], n_ka, n_kd, size_ka, size_kd, fex,
                   jex, error_sim_d);
            host_type ka_val_h(size_ka), kd_val_h(size_kd);
            host_type r0_mc(probSize), d_mc(probSize), ea_mc(probSize);
            thrust::device_ptr<double> r0_mc_d, d_mc_d, ea_mc_d;
            fex.get_r0_ptr(r0_mc_d);
            fex.get_ea_ptr(ea_mc_d);
            fex.get_d_ptr(d_mc_d);
            // Saving the error values and param values from device to host
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
            for ( int ind = 0; ind < probSize; ind++) {
                if ( error_sim_h[ind] < E0 ) {
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
                    output_data(gene_ind, r0_, ea_, d_, start_ka, end_ka, start_kd, end_kd,
                                temp_ka_vec, temp_kd_vec, ka_val_h, kd_val_h, out_files);
                    for ( int k = *(ka_start+ind); k < (*(ka_start+ind)+*(n_ka+ind)); k++) {
                        // Sensitivity for activators
                        if ( ka_val_h[k] > 0.05 ) {
                            act_vec.push_back(temp_ka_vec[k]);
                            acc_sub = true;
                            act_count++;
                        }
                    }
                    for ( int k = *(kd_start+ind); k < (*(kd_start+ind)+*(n_kd+ind)); k++) {
                        // Sensitivity for deactivators
                        if ( kd_val_h[k] > 0.05 ) {
                            inh_vec.push_back(temp_kd_vec[k]);
                            acc_sub = true;
                            inh_count++;
                        }
                    }
                    if (acc_sub) {
                        // Adding accepted candidates
                        acc_subnet++;
                        subnet_size.push_back(act_count+inh_count);
                        cout << "-------- Acceptable subnet number: " << acc_subnet << "  ------------"
                             << endl;
                        cout << "-------- Size of subnetwork: " << act_count+inh_count <<
                             "  ------------" << endl;
                        cout << "-------- Value of error is : " << error_sim_h[ind] << "  ------------"
                             << endl;
                        cout << "-------- Number of attempts made till now: " << subnet_att <<
                             "  ------------" << endl;
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
        stringstream ss;
        ss << (gene_ind + 1);
        string outfile(output_path_c);
        outfile += ss.str();
        outfile += ".out";
//        char out_file[32] = ".out";
//        char index_string[64];
//        sprintf(index_string, "%d", gene_ind+1);
//        strcat(index_string,out_file);
//        char out_path[96] = "../out/N1000/";
//        strcat(out_path,index_string);
        ofstream fileout;
        fileout.open( outfile.c_str(), ios_base::binary|ios_base::app|ios_base::out );
        if (fileout.is_open()) {
            cout << "File open successful!" << endl;
            fileout << "------------ START FOR ACTIVATORS --------------" << endl;
            for ( int ind_con = 0; ind_con < N_gene; ind_con++ ) {
                int act_count = count(act_vec.begin(),act_vec.end(),ind_con);
                fileout << setw(3) << ind_con+1 << "   " << setw(3) << act_count << endl;
            }
            fileout << "-------------- END FOR ACTIVATORS --------------" << endl;
            for ( int ind_con = 0; ind_con < N_gene; ind_con++ ) {
                int inh_count = count(inh_vec.begin(),inh_vec.end(),ind_con);
                fileout << setw(3) << ind_con+1 << "   " << setw(3) << inh_count << endl;
            }
            fileout << "------------ Size of subnetworks START----------- " << endl;
            for ( int ind_con = 0; ind_con < 20; ind_con++ ) {
                int size_count = count(subnet_size.begin(),subnet_size.end(),ind_con);
                fileout << setw(3) << ind_con << "   " << setw(3) << size_count << endl;
            }
            fileout << "------------ Size of subnetworks END----------- " << endl;
            fileout << "------------ Number of acceptable subnetworks for gene " << gene_ind
                    +1 <<
                    " is: " << acc_subnet << " --------------" << endl;
            fileout << "------------ Number of attempts for gene " << gene_ind+1 <<
                    " is: " << subnet_att*probSize << " --------------" << endl;
            fileout.close();
        } else cout << " open() failed" << endl;
        vector <double>().swap(act_vec);
        vector <double>().swap(inh_vec);
        vector <int>().swap(subnet_size);
        cout << "------------ End for gene " << gene_ind+1 << "  ------" << endl;
    }//for (int id_gene = 644; id_gene < 645; id_gene++)
    //
    return 0;
} /* MAIN */

void setup_spline(vector<double>& t_d_in, vector<double> x_d_in,
                  vector<double>& cub_coeff_spline_out, vector<double>& mean_xd_out)
{
    cout << "----Setting up spline co-efficients----" << endl;
    vector <double> deriv_time(N_time_points*N_gene);
    for ( int i = 0 ; i < N_gene ; i++ ) {
        spline_deriv ( i, N_time_points, t_d_in, x_d_in, deriv_time);
    }
    // Set up the co-effients matrix
    spline_coeff ( N_gene, N_time_points, t_d_in, x_d_in, deriv_time,
                   cub_coeff_spline_out);
    // Determine mean, standard deviation, and correlation statistics
    vector<double> sd_xd;
    for ( int ind = 0; ind < N_gene; ind++) {
        //x_spline represented in terms of dt = 1. Used to determine mean, sd, and correlation
        // N=1000 case dt=1 for the data anyway so not necessary
        vector <double> x_spline_temp(N_time_points);
        copy( x_d_in.begin()+ind*N_time_points, x_d_in.begin()+(ind+1)*N_time_points,
              x_spline_temp.begin());
        //calc_spline( x_d, t_d, cub_coeff_spline, ind, x_spline_temp);
        calc_stats( x_spline_temp, mean_xd_out, sd_xd);
        vector<double>().swap(x_spline_temp); // Why?
    }
    cout << "----END spline setup----" << endl;
}

void setup_output(int gene_ind_in, vector<string>& out_files_out) {
    stringstream ss;
    ss << (gene_ind_in + 1);
    string file_prefix(output_path_c);
    file_prefix += ss.str();
    out_files_out.push_back(file_prefix + ".r0");
    out_files_out.push_back(file_prefix + ".ea");
    out_files_out.push_back(file_prefix + ".d");
    out_files_out.push_back(file_prefix + ".nka");
    out_files_out.push_back(file_prefix + ".nkd");
    out_files_out.push_back(file_prefix + ".kavec");
    out_files_out.push_back(file_prefix + ".kdvec");
    out_files_out.push_back(file_prefix + ".kaval");
    out_files_out.push_back(file_prefix + ".kdval");
}

// Error bound
// Max number attempts per gene
// Acceptable subn for each genes

// First loop is which gene (specify start gene end gene)
// Second loop is acceptable subn for each gene and max num for each gene (AND condition)
// Final loop tries to go through generated error values, and sees which are accepted within bounds
// opt contains integrate and MC sim
// fex def is in lsoda/CUlsoda.hpp
