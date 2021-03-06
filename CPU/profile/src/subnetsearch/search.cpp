# include "../common.h"
using namespace std;
# include "../lsoda/cuLsoda.hpp"
# include "../optimization/opt.hpp"
# include "./search.hpp"
void random_select( const int &gene_ind, unsigned int seed_int, const state_type &mean, int &n_ka_temp, int &n_kd_temp, vector <int> &kavec_temp, vector <int> &kdvec_temp, state_type &kaval_temp, state_type &kdval_temp)
{
    mt19937 mt_r;
    mt_r.seed(time(NULL)*seed_int);
    std::uniform_int_distribution<> dist_r(0,RAND_MAX);
    srand(dist_r(mt_r));
    vector<bool> nodes_on(N_gene), free_nodes(N_gene);
    vector<int> free_nodes_vec;
    fill(nodes_on.begin(), nodes_on.end(), false);
    nodes_on[gene_ind] = true;	
    for ( int i = 0; i < n_ka_temp; i++ ) nodes_on[kavec_temp[i]] = true;	
    for ( int i = 0; i < n_kd_temp; i++ ) nodes_on[kdvec_temp[i]] = true;	
    free_nodes = nodes_on;
    free_nodes.flip();  
    /*for ( int j = 0; j < N_gene; j++) 
    {
        if (nodes_on[j]) cout << "###############  Node " << j << " is ON!  #############" << endl;
    }*/
    for ( int j = 0; j < N_gene; j++) 
    {
        if (free_nodes[j]) free_nodes_vec.push_back(j); 
    }
    // switch_vec represents 45% probability of adding(1)/deleting(2) a connection
    // and 10% probability of doing nothing
    int sum_elems, switch_call;
    sum_elems = accumulate(nodes_on.begin(),nodes_on.end(),0);
    
    if (sum_elems == N_gene)//Only remove connections or do nothing
    {
       // switch_vec represents 90% probability of deleting(2) a connection
       // and 10% probability of doing nothing
       vector<int> switch_vec;
       for (int i = 0; i < 9;i++) switch_vec.push_back(2);
       switch_vec.push_back(3);
       random_shuffle(switch_vec.begin(), switch_vec.end());
       switch_call = switch_vec[rand()/(RAND_MAX/(9+1) + 1)];  
    } 
    else if (sum_elems == 1)
    {
       // switch_vec represents 90% probability of adding(1) a connection
       // and 10% probability of doing nothing
       vector<int> switch_vec;
       for (int i = 0; i < 9;i++) switch_vec.push_back(1);
       switch_vec.push_back(3);
       random_shuffle(switch_vec.begin(), switch_vec.end());
       switch_call = switch_vec[rand()/(RAND_MAX/(9+1) + 1)];  
    }
    else
    {
       vector<int> switch_vec;
       for (int i = 0; i < 9;i++) switch_vec.push_back(1);
       for (int i = 0; i < 9;i++) switch_vec.push_back(2);
       switch_vec.push_back(3);
       switch_vec.push_back(3);
       random_shuffle(switch_vec.begin(), switch_vec.end());
       switch_call = switch_vec[rand()/(RAND_MAX/(19+1) + 1)];  
    }
    //cout << "---------- switch_call is: "<< switch_call << " -----------" << endl;
    srand(time(NULL));
    switch(switch_call)
    {
       case 1:
              {
                 random_shuffle(free_nodes_vec.begin(), free_nodes_vec.end());
                 int add_conn;
                 //add_conn = free_nodes_vec[rand()/(RAND_MAX/(free_nodes_vec.size()) + 1)];
                 add_conn = free_nodes_vec[dist_r(mt_r)/(RAND_MAX/(free_nodes_vec.size()) + 1)];
                 // ka and kd values in log space so no need of fabs for guassrand()   
                 //if ((rand()%2 == 0) && (n_ka_temp != 0))
                 //if (rand()%2 == 0) 
                 if (dist_r(mt_r)%2 == 0) 
                 {
                    kavec_temp.push_back(add_conn); 
                    //kaval_temp.push_back(pow(10,guassrand()));
                    kaval_temp.push_back(mean[add_conn]);
                    n_ka_temp = n_ka_temp + 1; 
                 }     
                 else 
                 {
                    kdvec_temp.push_back(add_conn); 
                    //kdval_temp.push_back(pow(10,guassrand())); 
                    kdval_temp.push_back(mean[add_conn]); 
                    n_kd_temp = n_kd_temp + 1; 
                 }                
                 break;
              }
       case 2:
              {
                 //if ((rand()%2 == 0) && (n_ka_temp != 0))
                 if ((dist_r(mt_r)%2 == 0) && (n_ka_temp != 0))
                 {
                    vector<int> indices(n_ka_temp);
                    for (int ind = 0; ind < n_ka_temp; ind++) indices[ind] = ind;
                    random_shuffle(indices.begin(), indices.end());
                    //int del_conn = indices[rand()/(RAND_MAX/(indices.size()) + 1)];
                    int del_conn = indices[dist_r(mt_r)/(RAND_MAX/(indices.size()) + 1)];
                    kavec_temp.erase(kavec_temp.begin()+del_conn);
                    kaval_temp.erase(kaval_temp.begin()+del_conn);
                    n_ka_temp = n_ka_temp - 1; 
                 }
                 else if (n_kd_temp !=0)
                 {
                    vector<int> indices(n_kd_temp);
                    for (int ind = 0; ind < n_kd_temp; ind++) indices[ind] = ind;
                    //random_shuffle(indices.begin(), indices.end());
                    //int del_conn = indices[rand()/(RAND_MAX/(indices.size()) + 1)];
                    int del_conn = indices[dist_r(mt_r)/(RAND_MAX/(indices.size()) + 1)];
                    kdvec_temp.erase(kdvec_temp.begin()+del_conn);
                    kdval_temp.erase(kdval_temp.begin()+del_conn);
                    n_kd_temp = n_kd_temp - 1; 
                 }         
                 break;
              }
       case 3:
              break;
    }
};
/*void biased_select( const int &gene_ind, const state_type &corr, int &n_ka_temp, int &n_kd_temp, vector <int> &kavec_temp, vector <int> &kdvec_temp, state_type &kaval_temp, state_type &kdval_temp)
{
    vector<bool> nodes_on(N_gene), free_nodes(N_gene);
    vector<int> free_nodes_vec;
    fill(nodes_on.begin(), nodes_on.end(), false);
    nodes_on[gene_ind] = true;	
    for ( int i = 0; i < n_ka_temp; i++ ) nodes_on[kavec_temp[i]] = true;	
    for ( int i = 0; i < n_kd_temp; i++ ) nodes_on[kdvec_temp[i]] = true;	
    free_nodes = nodes_on;
    free_nodes.flip();  
    for ( int j = 0; j < N_gene; j++) 
    {
        if (free_nodes[j]) free_nodes_vec.push_back(j); 
    }
    // switch_vec represents 45% probability of adding(1)/deleting(2) a connection
    // and 10% probability of doing nothing
    int sum_elems, switch_call;
    sum_elems = accumulate(nodes_on.begin(),nodes_on.end(),0);
    
    if (sum_elems == N_gene)//Only remove connections or do nothing
    {
       // switch_vec represents 90% probability of deleting(2) a connection
       // and 10% probability of doing nothing
       vector<int> switch_vec;
       for (int i = 0; i < 9;i++) switch_vec.push_back(2);
       switch_vec.push_back(3);
       random_shuffle(switch_vec.begin(), switch_vec.end());
       switch_call = switch_vec[rand()/(RAND_MAX/(9+1) + 1)];  
    } 
    else if (sum_elems == 1)
    {
       // switch_vec represents 90% probability of adding(1) a connection
       // and 10% probability of doing nothing
       vector<int> switch_vec;
       for (int i = 0; i < 9;i++) switch_vec.push_back(1);
       switch_vec.push_back(3);
       random_shuffle(switch_vec.begin(), switch_vec.end());
       switch_call = switch_vec[rand()/(RAND_MAX/(9+1) + 1)];  
    }
    else
    {
       vector<int> switch_vec;
       for (int i = 0; i < 9;i++) switch_vec.push_back(1);
       for (int i = 0; i < 9;i++) switch_vec.push_back(2);
       switch_vec.push_back(3);
       switch_vec.push_back(3);
       random_shuffle(switch_vec.begin(), switch_vec.end());
       switch_call = switch_vec[rand()/(RAND_MAX/(19+1) + 1)];  
    }
    //cout << "---------- switch_call is: "<< switch_call << " -----------" << endl;
    switch(switch_call)
    {
       case 1:
              {
                 //vector <double> free_nodes_weight(free_nodes_vec.size());
                 vector <double> free_nodes_weight;
                 // assign interval as indices to the free_node vector
                 vector <double> indices; 
                 // Determine the weights of each interval
                 indices.push_back(0.0);
                 for ( int ind = 0; ind < free_nodes_vec.size(); ind++)
                 {
                     //free_nodes_weight[ind] = fabs(corr[gene_ind*(N_gene-1)+free_nodes_vec[ind]]);  
                     free_nodes_weight.push_back(fabs(corr[gene_ind*(N_gene-1)+free_nodes_vec[ind]]));  
                     indices.push_back(double(ind+1));
                 }
                 piecewise_constant_distribution<>dist(indices.begin(), indices.end(), free_nodes_weight.begin());
                 mt19937 gen(time(0));
                 // Then determine the free node index using correlation
                 int add_conn_ind = int(dist(gen));
                 int add_conn = free_nodes_vec[add_conn_ind];
                 if ( corr[gene_ind*(N_gene-1)+free_nodes_vec[add_conn_ind]] > 0.0 ) 
                 {
                    kavec_temp.push_back(add_conn); 
                    kaval_temp.push_back(guassrand());
                    n_ka_temp = n_ka_temp + 1; 
                 }     
                 else
                 {
                    kdvec_temp.push_back(add_conn); 
                    kdval_temp.push_back(guassrand()); 
                    n_kd_temp = n_kd_temp + 1; 
                 }                
                 //add_conn = free_nodes_vec[int(dist(gen))];
                 // Use that to determine the free node vector
                 // Add a connection from the free nodes vector
                 break;
              }
       case 2:
              {
                 if ((rand()%2 == 0) && (n_ka_temp != 0))
                 {
                    vector<int> indices(n_ka_temp);
                    for (int ind = 0; ind < n_ka_temp; ind++) indices[ind] = ind;
                    random_shuffle(indices.begin(), indices.end());
                    int del_conn = indices[rand()/(RAND_MAX/(indices.size()) + 1)];
                    kavec_temp.erase(kavec_temp.begin()+del_conn);
                    kaval_temp.erase(kaval_temp.begin()+del_conn);
                    n_ka_temp = n_ka_temp - 1; 
                 }
                 else if (n_kd_temp !=0)
                 {
                    vector<int> indices(n_kd_temp);
                    for (int ind = 0; ind < n_kd_temp; ind++) indices[ind] = ind;
                    random_shuffle(indices.begin(), indices.end());
                    int del_conn = indices[rand()/(RAND_MAX/(indices.size()) + 1)];
                    kdvec_temp.erase(kdvec_temp.begin()+del_conn);
                    kdval_temp.erase(kdval_temp.begin()+del_conn);
                    n_kd_temp = n_kd_temp - 1; 
                 }         
                 break;
              }
       case 3:
              break;
    }
};*/
