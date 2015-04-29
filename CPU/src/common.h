#include <string>
#include <cfloat>
#include <ctime>
#include <iostream>
#include <functional>
#include <algorithm>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <iomanip>
#include <math.h>
#include <stdio.h>
#include <numeric>
#include <random>
#include <boost/timer.hpp>
#ifndef COMM_H
#define COMM_H
const int N_gene = 1000;
const int N_time_points = 12;
const int probSize = 3; // Number of acceptable profiles need to be generated in the simulation run
#endif
typedef double value_type;
typedef std::vector < value_type > state_type;
