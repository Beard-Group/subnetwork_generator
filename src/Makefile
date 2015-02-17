# Copyright 2009-2012 Karsten Ahnert
# Copyright 2009-2012 Mario Mulansky
#
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE_1_0.txt or
# copy at http://www.boost.org/LICENSE_1_0.txt)

CUDA_ROOT = /usr/local/cuda-5.0
ARCH = sm_35
#BOOST_ROOT = /usr/include/boost
THRUST_ROOT = /usr/local/cuda-5.0/include/thrust
SRC = .

CC = gcc-4.6.3
CXX = gcc-4.6.3
NVCC = $(CUDA_ROOT)/bin/nvcc

INCLUDES += -I$(THRUST_ROOT) -I$(CUDA_ROOT)/include -I../../../../.. 

NVCCFLAGS = -O3 $(INCLUDES) -arch $(ARCH) --compiler-bindir=/usr/local/gcc-4.6.3/bin/g++-4.6.3 #-Xcompiler -fopenmp -DTHRUST_HOST_BACKEND=THRUST_HOST_BACKEND_OMP
#NVCCFLAGS = -O3 $(INCLUDES) -arch $(ARCH) --compiler-bindir=/usr/local/gcc-4.6.3/bin/g++-4.6.3 -Xcompiler -DTHRUST_HOST_BACKEND=THRUST_HOST_BACKEND_CPP
#NVCCFLAGS = -O3 $(INCLUDES) --compiler-bindir=/home/raghu/gcc-4.6/bin/g++-4.6  -Xcompiler -fopenmp -DTHRUST_DEVICE_BACKEND=THRUST_DEVICE_BACKEND_OMP

#LDLIBS = -lcudart -lgomp -lm -lstdc++
LDLIBS = -lcudart -lm -lstdc++
LDFLAGS = -L$(CUDA_ROOT)/lib64

%.co : %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

all : run_N1000_gene645

run_N1000_gene645 :   $(SRC)/io/io.co $(SRC)/spline/spline.co $(SRC)/subnetsearch/search.co $(SRC)/opt/opt.co main.co
#lsoda_run :   $(SRC)/io/io.co $(SRC)/lsoda/cuLsoda.co main.co
#lsoda_run :   $(SRC)/io/io.co $(SRC)/spline/spline.co $(SRC)/model/model.co main.co
#run :   $(SRC)/io/io.co $(SRC)/spline/spline.co $(SRC)/model/model.co $(SRC)/MCsim/mc.co $(SRC)/neldermead/nm.co main.co
	$(CC) -o run_N1000_gene645 $(LDFLAGS) $(LDLIBS) $(SRC)/io/io.co $(SRC)/spline/spline.co $(SRC)/subnetsearch/search.co $(SRC)/opt/opt.co main.co
#	$(CC) -o lsoda_run $(LDFLAGS) $(LDLIBS) $(SRC)/lsoda/cuLsoda.co $(SRC)/io/io.co main.co
#	$(CC) -o run $(LDFLAGS) $(LDLIBS) $(SRC)/io/io.co $(SRC)/spline/spline.co $(SRC)/model/model.co main.co
#	$(CC) -o run $(LDFLAGS) $(LDLIBS) $(SRC)/io/io.co $(SRC)/spline/spline.co $(SRC)/model/model.co $(SRC)/MCsim/mc.co $(SRC)/neldermead/nm.co main.co
$(SRC)/io/io.co : $(SRC)/io/io.cu
#$(SRC)/lsoda/cuLsoda.co : $(SRC)/lsoda/cuLsoda.cu
$(SRC)/spline/spline.co : $(SRC)/spline/spline.cu
$(SRC)/subnetsearch/search.co : $(SRC)/subnetsearch/search.cu
$(SRC)/opt/opt.co : $(SRC)/opt/opt.cu
main.co : main.cu

clean :
#	rm *~ *.o *.co run
#       rm -f  run_lsoda *.o  *.co  *~ $(SRC)/*/*.o $(SRC)/*/*.co
	rm -f  *.o  *.co  *~ $(SRC)/*/*.o $(SRC)/*/*.co
