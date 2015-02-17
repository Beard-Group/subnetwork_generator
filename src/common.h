#include <string>
#include <cfloat>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <iomanip>
#include <math.h>
#include <numeric>
#include <cmath>
#include <stdio.h>
#include <string.h>
#include <cuda.h>
//#include <boost/timer.hpp>
//#include "./lsoda/cuLsoda_kernel.cu"
// Thrust and Boost libraries
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/random.h>
#include <thrust/copy.h>
#include <thrust/transform.h>
#include <thrust/functional.h>
#include <thrust/fill.h>
#include <thrust/extrema.h>
#include <thrust/sort.h>
#ifndef COMM_H
#define COMM_H
const int N_gene = 1000;
const int N_time_points = 12;
const int probSize =  1600;
const int N_att = 10;
const int threadsPerBlock = 16;
const int blocksPerGrid = (probSize + threadsPerBlock -1)/threadsPerBlock;
#endif
//typedef float value_type;
typedef double value_type;
typedef thrust::device_vector < value_type > state_type;
//typedef thrust::host_vector < value_type > state_type;
typedef thrust::host_vector < value_type > host_type;
